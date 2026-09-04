"""Read-only BioMA installation and runtime diagnostics.

The doctor intentionally does not import optional scientific packages in the
launcher process.  It probes them in their configured Python/R environments,
which makes the command useful on clusters where each module has a separate
Conda environment.  A project INI is optional; when supplied, the enabled
modules determine which checks are required.
"""

from __future__ import annotations

import configparser
import os
import platform
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .runtime import subprocess_environment


FALSE_VALUES = {"", "0", "false", "no", "off", "disabled", "none"}
KNOWN_MODULES = {"gf", "rona", "mar", "load", "niche", "wfmoment", "vulnerability"}

# Package-to-module hints are deliberately kept here rather than importing
# optional packages in BioMA itself.  This table controls whether a missing
# package is a blocking error for a project configuration.
R_PACKAGE_MODULES = {
    "gradientForest": {"gf"},
    "FNN": {"gf"},
    "data.table": {"gf", "rona", "mar", "load", "vulnerability"},
    "geosphere": {"gf"},
    "mar": {"mar"},
    "terra": {"rona"},
    "raster": {"gf", "rona", "niche", "load", "vulnerability"},
    "plotrix": {"rona"},
    "dismo": {"niche"},
    "sf": {"gf", "niche", "load", "vulnerability"},
    "dplyr": {"mar", "load", "niche", "vulnerability"},
    "tidyr": {"mar", "load", "niche", "vulnerability"},
    "ggplot2": {"rona", "mar", "load", "niche", "vulnerability", "wfmoment"},
    "readr": {"load", "niche"},
    "tibble": {"mar", "load", "niche"},
    "patchwork": {"load", "niche"},
    "cowplot": {"gf", "load", "niche"},
    "ragg": {"gf", "niche"},
    "future": {"niche", "load"},
    "mgcv": {"niche", "rona"},
    "ranger": {"load", "vulnerability"},
    "rsample": {"load"},
    "tidymodels": {"load"},
}

PYTHON_PACKAGE_MODULES = {
    "numpy": {"wfmoment"},
    "pandas": {"wfmoment"},
    "rasterio": {"wfmoment"},
    "wfmoments": {"wfmoment"},
}

# These are only opportunistic fallbacks mirroring module resolvers.  They are
# added when the path is present; a portable installation still uses PATH or
# an explicit INI value.
# Runtime discovery is intentionally portable.  Cluster installations should
# put their reviewed executables in the module INI (or PATH); no private host
# paths are embedded in the public package.
KNOWN_R_DEFAULTS = {}
KNOWN_PYTHON_DEFAULTS = {}


@dataclass
class DoctorCheck:
    """One diagnostic result.

    ``status`` is one of ``ok``, ``missing``, ``failed``, ``skipped`` or
    ``warning``.  ``required`` is kept in the output so users can distinguish
    an absent optional MaxEnt jar from a missing dependency for an enabled
    module.
    """

    name: str
    status: str
    required: bool
    detail: str = ""
    path: Optional[str] = None
    version: Optional[str] = None
    kind: str = "runtime"
    component: Optional[str] = None

    def as_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "ok": self.status == "ok",
            "required": self.required,
            "detail": self.detail,
            "path": self.path,
            "version": self.version,
            "kind": self.kind,
            "component": self.component,
        }


def _read_ini(path: Path) -> Optional[configparser.ConfigParser]:
    parser = configparser.ConfigParser(interpolation=None)
    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            parser.read_file(handle)
    except (OSError, UnicodeError, configparser.Error):
        return None
    return parser


def _resolve_value(value: str, base: Path) -> str:
    value = str(value or "").strip()
    if not value or value.lower() in FALSE_VALUES or value.lower() == "auto":
        return ""
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = base / candidate
    return str(candidate.resolve())


def _runtime_value(value: str, base: Path) -> str:
    """Resolve an executable path while preserving bare command names."""
    value = str(value or "").strip()
    if not value or value.lower() in FALSE_VALUES or value.lower() == "auto":
        return ""
    # ``Rscript``, ``python3`` and similar names should be looked up on PATH;
    # resolving them relative to an INI would turn them into a false path.
    if not any(token in value for token in ("/", "\\")) and not value.startswith((".", "~")):
        return value
    return _resolve_value(value, base) or value


def _enabled_modules(parser: Optional[configparser.ConfigParser], base: Path) -> Tuple[Tuple[str, ...], Dict[str, Path]]:
    """Return enabled modules and their INI paths without strict validation."""
    if parser is None or not parser.has_section("modules"):
        return tuple(), {}
    modules: List[str] = []
    paths: Dict[str, Path] = {}
    for module, raw in parser.items("modules"):
        value = (raw or "").strip()
        if not value or value.lower() in FALSE_VALUES:
            continue
        module_name = module.strip().lower()
        modules.append(module_name)
        if value.lower() not in {"true", "yes", "on", "1"}:
            candidate = Path(value).expanduser()
            if not candidate.is_absolute():
                candidate = base / candidate
            paths[module_name] = candidate.resolve()
        else:
            # Project files normally require a concrete module INI.  Retain a
            # sentinel path so the caller can report the malformed entry.
            paths[module_name] = (base / (".bioma-missing-{}-config.ini".format(module_name))).resolve()
    return tuple(dict.fromkeys(modules)), paths


def _config_context(config_path: Optional[Path]) -> Dict[str, Any]:
    """Collect non-fatal runtime hints from a project or standalone INI."""
    context: Dict[str, Any] = {
        "config": None,
        "modules": tuple(),
        "module_configs": {},
        "pythons": [],
        "rscripts": [],
        "rscript_modules": {},
        "python_modules": {},
        "gdal": [],
        "java": [],
        "maxent_jar": "",
        "config_error": "",
        "config_issues": [],
    }
    if config_path is None:
        return context
    path = Path(config_path).expanduser().resolve()
    context["config"] = str(path)
    parser = _read_ini(path)
    if parser is None:
        context["config_error"] = "configuration could not be read: {}".format(path)
        return context
    base = path.parent
    modules, module_paths = _enabled_modules(parser, base)
    unknown_modules = sorted(set(modules).difference(KNOWN_MODULES))
    if unknown_modules:
        context["config_issues"].append("unknown module(s): {}".format(", ".join(unknown_modules)))
    # A standalone module INI has no [modules] section.  Infer its module from
    # characteristic sections/name so that package checks remain useful.
    if not modules:
        stem = path.stem.lower()
        for candidate in ("gf", "rona", "mar", "load", "niche", "wfmoment", "vulnerability"):
            if candidate in stem:
                modules = (candidate,)
                break
        # Names such as ``project.ini`` are common for a standalone file.  In
        # that case infer the module from its required input keys/sections.
        if not modules and parser.has_section("inputs"):
            keys = {key.lower() for key in parser["inputs"].keys()}
            if {"variables", "tuning", "projection"}.issubset(set(parser.sections())):
                modules = ("niche",)
            elif "present_climate" in keys or "coordinates" in keys:
                modules = ("gf",)
            elif "vcf_dir" in keys:
                modules = ("load",)
            elif "alt_frequency" in keys or "unld_dir" in keys:
                modules = ("rona",)
            elif "current_raster" in keys:
                modules = ("wfmoment",)
            elif any(key in keys for key in ("gf_offsets", "rona_ensemble", "load_predictors", "niche_raster")):
                modules = ("vulnerability",)
            elif "vcf" in keys and parser.has_section("parameters"):
                modules = ("mar",)
        context["standalone_parser"] = parser
        context["standalone_base"] = base
    context["modules"] = modules
    context["module_configs"] = module_paths

    def add_unique(key: str, value: str) -> None:
        if value and value not in context[key]:
            context[key].append(value)

    def add_runtime(key: str, value: str, runtime_modules: Iterable[str] = ()) -> None:
        if not value:
            return
        add_unique(key, value)
        mapping_name = "rscript_modules" if key == "rscripts" else "python_modules"
        mapping = context[mapping_name]
        mapping.setdefault(value, set()).update(runtime_modules)

    shared_runtime_r = ""
    shared_runtime_py = ""
    if parser.has_section("shared"):
        shared = parser["shared"]
        shared_runtime_r = _runtime_value(shared.get("rscript", ""), base)
        shared_runtime_py = _runtime_value(shared.get("compute_python", ""), base)
        add_runtime("rscripts", shared_runtime_r, modules)
        add_runtime("pythons", shared_runtime_py, modules)
    if parser.has_section("inputs"):
        inputs = parser["inputs"]
        context["maxent_jar"] = _resolve_value(inputs.get("maxent_jar", ""), base)
    if parser.has_section("doctor"):
        doctor = parser["doctor"]
        for key in ("python", "compute_python"):
            add_runtime("pythons", _runtime_value(doctor.get(key, ""), base), modules)
        for key in ("rscript", "plot_rscript"):
            add_runtime("rscripts", _runtime_value(doctor.get(key, ""), base), modules)
        for key in ("gdalinfo", "gdallocationinfo", "ogrinfo"):
            add_unique("gdal", _runtime_value(doctor.get(key, ""), base))
        add_unique("java", _runtime_value(doctor.get("java", ""), base))
        jar = _resolve_value(doctor.get("maxent_jar", ""), base)
        if jar:
            context["maxent_jar"] = jar

    # Read each module INI, tolerating a missing/invalid path so doctor can
    # report all other checks in one invocation.
    for module, module_path in module_paths.items():
        module_parser = _read_ini(module_path)
        if module_parser is None:
            context["config_issues"].append("module {} configuration could not be read: {}".format(module, module_path))
            continue
        module_base = module_path.parent
        for section_name in ("parameters", "analysis"):
            if not module_parser.has_section(section_name):
                continue
            section = module_parser[section_name]
            for key in ("rscript", "plot_rscript"):
                raw = section.get(key, "").strip()
                if raw and not shared_runtime_r:
                    add_runtime("rscripts", _runtime_value(raw, module_base), (module,))
            for key in ("compute_python", "calc_python", "python"):
                raw = section.get(key, "").strip()
                if raw and not shared_runtime_py:
                    add_runtime("pythons", _runtime_value(raw, module_base), (module,))
            for key in ("gdalinfo", "gdallocationinfo", "ogrinfo"):
                raw = section.get(key, "").strip()
                if raw:
                    add_unique("gdal", _runtime_value(raw, module_base))
            raw_java = section.get("java", "").strip()
            if raw_java:
                add_unique("java", _runtime_value(raw_java, module_base))
        if module_parser.has_section("inputs"):
            raw_jar = module_parser["inputs"].get("maxent_jar", "").strip()
            if raw_jar and not context["maxent_jar"]:
                context["maxent_jar"] = _resolve_value(raw_jar, module_base)

    # For a standalone INI, process its own sections exactly as a module file.
    standalone = context.get("standalone_parser")
    if standalone is not None:
        module_base = context["standalone_base"]
        for section_name in ("parameters", "analysis"):
            if not standalone.has_section(section_name):
                continue
            section = standalone[section_name]
            for key in ("rscript", "plot_rscript"):
                raw = section.get(key, "").strip()
                if raw:
                    add_runtime("rscripts", _runtime_value(raw, module_base), modules)
            for key in ("compute_python", "calc_python", "python"):
                raw = section.get(key, "").strip()
                if raw:
                    add_runtime("pythons", _runtime_value(raw, module_base), modules)
        if standalone.has_section("inputs"):
            raw_jar = standalone["inputs"].get("maxent_jar", "").strip()
            if raw_jar:
                context["maxent_jar"] = _resolve_value(raw_jar, module_base)

    # Mirror the defaults used by the workflow loaders when those paths are
    # actually installed on the host.  This is especially useful on the
    # cluster, while remaining invisible on a normal desktop installation.
    assigned_r_modules = set().union(*(set(value) for value in context["rscript_modules"].values())) if context["rscript_modules"] else set()
    assigned_py_modules = set().union(*(set(value) for value in context["python_modules"].values())) if context["python_modules"] else set()
    for module in modules:
        if module not in assigned_r_modules:
            for candidate in KNOWN_R_DEFAULTS.get(module, ()):
                if Path(candidate).is_file():
                    add_runtime("rscripts", candidate, (module,))
        if module not in assigned_py_modules:
            for candidate in KNOWN_PYTHON_DEFAULTS.get(module, ()):
                if Path(candidate).is_file():
                    add_runtime("pythons", candidate, (module,))
    return context


def _which(value: str, fallback: str) -> Optional[str]:
    candidate = str(value or fallback or "").strip()
    if not candidate:
        return None
    # shutil.which handles both a bare command and an absolute executable.
    found = shutil.which(candidate)
    if found:
        return str(Path(found).resolve())
    path = Path(candidate).expanduser()
    if path.is_file():
        return str(path.resolve())
    return None


def _run(command: Sequence[str], timeout: int = 15) -> Tuple[Optional[int], str, str]:
    try:
        result = subprocess.run(
            list(command),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
            # Keep Conda libraries/data isolated when the user points to an
            # executable inside an environment, matching workflow runners.
            env=subprocess_environment(str(command[0])) if command else os.environ.copy(),
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return None, "", str(error)
    return result.returncode, (getattr(result, "stdout", "") or "").strip(), (getattr(result, "stderr", "") or "").strip()


def _version_from_output(stdout: str, stderr: str = "") -> str:
    text = (stdout or stderr or "").strip().replace("\r", "")
    if not text:
        return ""
    return text.splitlines()[0][:240]


def _command_check(name: str, value: str, fallback: str, required: bool, args: Sequence[str] = ("--version",)) -> DoctorCheck:
    resolved = _which(value, fallback)
    if not resolved:
        return DoctorCheck(name, "missing" if required else "skipped", required, "executable not found", str(value or fallback))
    code, stdout, stderr = _run([resolved] + list(args))
    if code is None:
        return DoctorCheck(name, "failed" if required else "warning", required, stderr or "could not execute", resolved)
    if code != 0:
        detail = stderr or stdout or "exit status {}".format(code)
        return DoctorCheck(name, "failed" if required else "warning", required, detail[:240], resolved)
    return DoctorCheck(name, "ok", required, "available", resolved, _version_from_output(stdout, stderr))


def _python_package_check(name: str, executable: str, package: str, required: bool) -> DoctorCheck:
    resolved = _which(executable, sys.executable)
    if not resolved:
        return DoctorCheck(name, "missing" if required else "skipped", required, "Python executable not found", str(executable) if executable else None, kind="python-package", component=package)
    expression = (
        "import importlib.util; s=importlib.util.find_spec({!r}); "
        "raise SystemExit(0 if s is not None else 1)"
    ).format(package)
    code, stdout, stderr = _run([resolved, "-c", expression])
    if code == 0:
        version_expr = (
            "import importlib.metadata as m; "
            "print(m.version({!r}))"
        ).format(package)
        vcode, vout, _ = _run([resolved, "-c", version_expr])
        version = _version_from_output(vout) if vcode == 0 else "installed"
        return DoctorCheck(name, "ok", required, "Python package available", resolved, version, "python-package", package)
    detail = stderr or stdout or "Python package '{}' is not installed".format(package)
    return DoctorCheck(name, "missing" if required else "warning", required, detail[:240], resolved, kind="python-package", component=package)


def _r_package_check(
    name: str,
    rscript: Optional[str],
    package: str,
    required: bool,
    symbol: Optional[str] = None,
) -> DoctorCheck:
    if not rscript:
        return DoctorCheck(name, "missing" if required else "skipped", required, "Rscript was not found", None, kind="r-package", component=package)
    resolved = _which(rscript, "Rscript")
    if not resolved:
        return DoctorCheck(name, "missing" if required else "skipped", required, "Rscript executable not found", str(rscript), kind="r-package", component=package)
    symbol_test = ""
    if symbol:
        symbol_test = "; if (!exists('{}', envir=asNamespace('{}'), inherits=FALSE)) quit(status=1)".format(symbol, package)
    expression = (
        "ok <- requireNamespace('{0}', quietly=TRUE); "
        "if (!ok) quit(status=1){1}; cat(as.character(packageVersion('{0}')))"
    ).format(package, symbol_test)
    code, stdout, stderr = _run([resolved, "--vanilla", "-e", expression])
    if code == 0:
        return DoctorCheck(name, "ok", required, "R package available", resolved, _version_from_output(stdout, stderr), "r-package", package)
    detail = stderr or stdout or "R package '{}' is not installed".format(package)
    return DoctorCheck(name, "missing" if required else "warning", required, detail[:240], resolved, kind="r-package", component=package)


def _maxent_check(path_value: str, required: bool) -> DoctorCheck:
    if not path_value:
        return DoctorCheck("maxent_jar", "missing" if required else "skipped", required, "maxent.jar is not configured", kind="file", component="maxent.jar")
    path = Path(path_value).expanduser()
    if not path.is_file():
        return DoctorCheck("maxent_jar", "missing" if required else "warning", required, "file does not exist", str(path), kind="file", component="maxent.jar")
    try:
        with zipfile.ZipFile(str(path), "r") as archive:
            bad = archive.testzip()
            if bad:
                return DoctorCheck("maxent_jar", "failed", required, "corrupt archive member: {}".format(bad), str(path), kind="file", component="maxent.jar")
            members = len(archive.namelist())
    except (OSError, zipfile.BadZipFile) as error:
        return DoctorCheck("maxent_jar", "failed", required, "invalid JAR/ZIP: {}".format(error), str(path), kind="file", component="maxent.jar")
    return DoctorCheck("maxent_jar", "ok", required, "valid JAR archive ({} entries)".format(members), str(path), kind="file", component="maxent.jar")


def _discover_maxent_jar(rscript: Optional[str]) -> str:
    """Ask dismo for its bundled jar when the user did not configure one."""
    resolved = _which(rscript or "", "Rscript")
    if not resolved:
        return ""
    expression = "cat(system.file('java', 'maxent.jar', package='dismo'))"
    code, stdout, _ = _run([resolved, "--vanilla", "-e", expression])
    if code != 0:
        return ""
    # R may print startup text around the path.  Select the first existing
    # path rather than trusting arbitrary output from a site-wide profile.
    for token in stdout.splitlines():
        candidate = token.strip()
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    return ""


def run_doctor(
    config_path: Optional[Path] = None,
    *,
    python: Optional[str] = None,
    rscript: Optional[str] = None,
    gdalinfo: Optional[str] = None,
    gdallocationinfo: Optional[str] = None,
    ogrinfo: Optional[str] = None,
    java: Optional[str] = None,
    maxent_jar: Optional[Path] = None,
    strict: bool = False,
) -> Dict[str, Any]:
    """Run read-only dependency checks and return a JSON-serialisable report.

    ``config_path`` may be a project INI or any standalone module INI.  CLI
    overrides take precedence over values discovered in the configuration.
    In a project, dependencies used by enabled modules are required.  Without
    a project all optional scientific components are reported but do not make
    the command fail unless ``strict=True``.
    """
    context = _config_context(config_path)
    modules = set(context.get("modules", tuple()))
    configured_rs = list(context.get("rscripts", []))
    configured_py = list(context.get("pythons", []))
    configured_gdal = list(context.get("gdal", []))
    configured_java = list(context.get("java", []))
    if rscript:
        configured_rs = [rscript]
        context["rscript_modules"] = {rscript: set(modules), str(Path(rscript).expanduser().resolve()): set(modules)}
    if python:
        configured_py = [python]
        context["python_modules"] = {python: set(modules), str(Path(python).expanduser().resolve()): set(modules)}
    if gdalinfo:
        configured_gdal.insert(0, gdalinfo)
    if java:
        configured_java.insert(0, java)
    configured_rs = list(dict.fromkeys(x for x in configured_rs if x))
    configured_py = list(dict.fromkeys(x for x in configured_py if x))
    configured_gdal = list(dict.fromkeys(x for x in configured_gdal if x))
    configured_java = list(dict.fromkeys(x for x in configured_java if x))

    # Modules with no explicit interpreter use the ambient command.  Keep it
    # as a separate candidate when another module has a dedicated environment,
    # so a missing package in one environment cannot hide a usable one in the
    # other.
    r_assignments = context.setdefault("rscript_modules", {})
    assigned_r_modules = set().union(*(set(value) for value in r_assignments.values())) if r_assignments else set()
    unassigned_r_modules = modules.difference(assigned_r_modules)
    if unassigned_r_modules and "Rscript" not in configured_rs:
        configured_rs.append("Rscript")
        r_assignments["Rscript"] = set(unassigned_r_modules)
    py_assignments = context.setdefault("python_modules", {})
    assigned_py_modules = set().union(*(set(value) for value in py_assignments.values())) if py_assignments else set()
    unassigned_py_modules = modules.difference(assigned_py_modules)
    if unassigned_py_modules and sys.executable not in configured_py:
        configured_py.append(sys.executable)
        py_assignments[sys.executable] = set(unassigned_py_modules)

    # Base interpreter is always required: otherwise the command could not run.
    checks: List[DoctorCheck] = []
    py_value = configured_py[0] if configured_py else sys.executable
    checks.append(_command_check("python", py_value, sys.executable, True, ("--version",)))

    r_required = bool(modules.intersection({"gf", "rona", "mar", "load", "niche", "vulnerability", "wfmoment"}))
    primary_r = configured_rs[0] if configured_rs else None
    r_check = _command_check("rscript", primary_r or "", "Rscript", r_required, ("--version",))
    checks.append(r_check)
    for index, candidate in enumerate(configured_rs[1:], start=2):
        assigned = r_assignments.get(candidate, set())
        checks.append(_command_check("rscript[{}]".format(index), candidate, "Rscript", bool(assigned), ("--version",)))
    # R is normally adjacent to Rscript; explicitly checking it catches broken
    # installations while remaining harmless for custom wrapper scripts.
    r_path = None
    if r_check.path:
        candidate = Path(r_check.path).with_name("R" + Path(r_check.path).suffix)
        r_path = str(candidate) if candidate.exists() else "R"
    checks.append(_command_check("r", r_path or "", "R", r_required, ("--version",)))

    gdal_required = bool(modules.intersection({"gf", "rona", "niche", "load", "vulnerability"}))
    gdal_names = (("gdalinfo", gdalinfo), ("gdallocationinfo", gdallocationinfo), ("ogrinfo", ogrinfo))
    for index, (name, override) in enumerate(gdal_names):
        configured = override or (configured_gdal[index] if index < len(configured_gdal) else "")
        checks.append(_command_check(name, configured, name, gdal_required, ("--version",)))

    java_required = "niche" in modules
    checks.append(_command_check("java", java or (configured_java[0] if configured_java else ""), "java", java_required, ("-version",)))

    jar_value = str(Path(maxent_jar).expanduser().resolve()) if maxent_jar else context.get("maxent_jar", "")
    if not jar_value and java_required:
        for candidate in (configured_rs or [primary_r]):
            jar_value = _discover_maxent_jar(candidate)
            if jar_value:
                break
    checks.append(_maxent_check(jar_value, java_required))

    # Probe each configured Rscript for the package checks.  A project can use
    # separate R installations for GF/RONA/MAR; checking each one avoids a
    # false green result from the first environment only.
    r_candidates = configured_rs or ([r_check.path] if r_check.path else [])
    if not r_candidates:
        r_candidates = [None]
    package_requirements = list(R_PACKAGE_MODULES.items())
    for label, required_modules in package_requirements:
        applicable: List[Tuple[int, Optional[str], set]] = []
        for index, candidate in enumerate(r_candidates):
            assigned = context.get("rscript_modules", {}).get(candidate, set())
            if not assigned and index < len(configured_rs):
                assigned = context.get("rscript_modules", {}).get(configured_rs[index], set())
            # A package is checked in an environment assigned to a module that
            # uses it.  With no project assignment, probe the first ambient
            # environment for every package as an optional inventory.
            if assigned and not assigned.intersection(required_modules):
                continue
            applicable.append((index, candidate, assigned))
        for position, (index, candidate, assigned) in enumerate(applicable):
            check_name = label if len(applicable) == 1 else "{}[{}]".format(label, position + 1)
            package_required = bool(assigned.intersection(required_modules)) if assigned else bool(modules.intersection(required_modules) and position == 0)
            symbol = "MARPIPELINE" if label == "mar" else None
            checks.append(_r_package_check(check_name, candidate, label, package_required, symbol=symbol))

    py_candidates = configured_py or [py_value]
    for package, required_modules in PYTHON_PACKAGE_MODULES.items():
        applicable = []
        for index, candidate in enumerate(py_candidates):
            assigned = context.get("python_modules", {}).get(candidate, set())
            if not assigned and index < len(configured_py):
                assigned = context.get("python_modules", {}).get(configured_py[index], set())
            if assigned and not assigned.intersection(required_modules):
                continue
            applicable.append((candidate, assigned))
        for position, (candidate, assigned) in enumerate(applicable):
            suffix = "" if len(applicable) == 1 else "[{}]".format(position + 1)
            check_name = "{}{}".format(package, suffix)
            required = bool(assigned.intersection(required_modules)) if assigned else bool(modules.intersection(required_modules) and position == 0)
            checks.append(_python_package_check(check_name, candidate, package, required))

    if context.get("config_error"):
        # A supplied configuration that cannot be read is actionable.  Keep it
        # required so `bioma doctor project.ini` exits non-zero rather than
        # silently checking an unrelated ambient environment.
        checks.append(DoctorCheck("configuration", "failed", True, context["config_error"], context.get("config"), kind="configuration", component="project.ini"))
    for issue in context.get("config_issues", []):
        checks.append(DoctorCheck("configuration", "failed", True, str(issue), context.get("config"), kind="configuration", component="module.ini"))

    failures = [
        item
        for item in checks
        if (
            item.status in {"missing", "failed"} and (item.required or strict)
        ) or (item.status in {"skipped", "warning"} and strict)
    ]
    warnings = [item for item in checks if item.status in {"warning", "missing", "failed"} and item not in failures]
    check_rows = [item.as_dict() for item in checks]
    check_map = {item.name: item.as_dict() for item in checks}
    runtime_versions = {
        row["name"]: row.get("version")
        for row in check_rows
        if row.get("kind") == "runtime" and row.get("version")
    }
    report: Dict[str, Any] = {
        "tool": "bioma-doctor",
        "bioma_version": _bioma_version(),
        "status": "fail" if failures else "pass",
        "result": "failed" if failures else "passed",
        "ok": not failures,
        "passed": not failures,
        "strict": bool(strict),
        "config": context.get("config"),
        "modules": sorted(modules),
        "module_configs": {name: str(path) for name, path in context.get("module_configs", {}).items()},
        "platform": platform.platform(),
        "checks": check_rows,
        "check_results": check_rows,
        "checks_by_name": check_map,
        "check_map": check_map,
        "components": {name: row.get("status") for name, row in check_map.items()},
        "runtime_versions": runtime_versions,
        "next_steps": (
            ["Install or configure the required components listed above."]
            if failures
            else (["Run with --strict or provide a project INI to require every module dependency."] if not config_path else [])
        ),
        "summary": {
            "total": len(checks),
            "ok": sum(item.status == "ok" for item in checks),
            "warnings": len(warnings),
            "failures": len(failures),
        },
    }
    return report


def _bioma_version() -> str:
    try:
        from . import __version__

        return str(__version__)
    except Exception:
        return "unknown"


def format_doctor_report(report: Mapping[str, Any]) -> str:
    """Render a compact human-readable report for terminal users."""
    status = "PASS" if report.get("status") == "pass" else "FAIL"
    lines = ["BioMA doctor: {}".format(status)]
    if report.get("config"):
        lines.append("Configuration: {}".format(report["config"]))
    modules = report.get("modules") or []
    if modules:
        lines.append("Enabled modules: {}".format(", ".join(modules)))
    symbols = {"ok": "OK", "missing": "MISSING", "failed": "FAIL", "warning": "WARN", "skipped": "SKIP"}
    for item in report.get("checks", []):
        marker = symbols.get(str(item.get("status")), str(item.get("status", "")).upper())
        required = " required" if item.get("required") else " optional"
        detail = item.get("detail") or ""
        version = item.get("version")
        if version:
            detail = "{} ({})".format(detail, version)
        lines.append("[{}] {:18s}{}{}".format(marker, str(item.get("name", "")), required, ": " + detail if detail else ""))
    summary = report.get("summary", {})
    lines.append("Summary: {ok} OK, {warnings} warning(s), {failures} failure(s)".format(**{
        "ok": summary.get("ok", 0),
        "warnings": summary.get("warnings", 0),
        "failures": summary.get("failures", 0),
    }))
    return "\n".join(lines)


check_environment = run_doctor


__all__ = ["DoctorCheck", "check_environment", "format_doctor_report", "run_doctor"]
