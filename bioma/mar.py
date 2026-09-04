"""Configuration-driven wrapper for the MAR population-genomics workflow."""

from __future__ import annotations

import configparser
import gzip
import hashlib
import json
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Optional, Sequence, Tuple

from .gf_frequency import InputError
from .runtime import subprocess_environment


@dataclass(frozen=True)
class MarConfig:
    config_path: Path
    vcf: Path
    lonlat: Path
    scenario_file: Optional[Path]
    output_dir: Path
    name: str
    geom_id: int
    scheme: str
    nrep: int
    xfrac: float
    quorum: bool
    randseed: int
    maxsnps: Optional[int]
    marsteps: Tuple[str, ...]
    rscript: str

    def payload(self) -> Dict[str, object]:
        return {
            "inputs": {"vcf": str(self.vcf), "lonlat": str(self.lonlat), "scenario_file": str(self.scenario_file) if self.scenario_file else None},
            "analysis": {"output_dir": str(self.output_dir), "name": self.name, "geom_id": self.geom_id},
            "parameters": {
                "scheme": self.scheme,
                "nrep": self.nrep,
                "xfrac": self.xfrac,
                "quorum": self.quorum,
                "randseed": self.randseed,
                "maxsnps": self.maxsnps,
                "marsteps": list(self.marsteps),
                "rscript": self.rscript,
            },
        }


def _path(value: str, base: Path, label: str) -> Path:
    value = value.strip()
    if not value:
        raise InputError("Missing [inputs] {}".format(label))
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _rscript(value: str) -> str:
    if value.strip():
        return value.strip()
    return "Rscript"


def load_mar_config(path: Path) -> MarConfig:
    path = path.expanduser().resolve()
    parser = configparser.ConfigParser(interpolation=None)
    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            parser.read_file(handle)
    except (OSError, configparser.Error, UnicodeError) as error:
        raise InputError("Cannot read MAR configuration: {}".format(error))
    if not parser.has_section("inputs") or not parser.has_section("analysis"):
        raise InputError("MAR configuration requires [inputs] and [analysis]")
    inputs, analysis = parser["inputs"], parser["analysis"]
    params = parser["parameters"] if parser.has_section("parameters") else {}
    base = path.parent
    vcf = _path(inputs.get("vcf", ""), base, "vcf")
    lonlat = _path(inputs.get("lonlat", ""), base, "lonlat")
    scenario_value = inputs.get("scenario_file", "").strip()
    scenario_file = _path(scenario_value, base, "scenario_file") if scenario_value else None
    output = _path(analysis.get("output_dir", ""), base, "output_dir")
    if not vcf.is_file() or not lonlat.is_file():
        raise InputError("MAR input file does not exist")
    scheme = params.get("scheme", "random").strip().lower()
    if scheme not in {"random", "eastwest", "westeast", "northsouth", "southnorth"}:
        raise InputError("scheme must be random, eastwest, westeast, northsouth, or southnorth")
    nrep = int(params.get("nrep", "10")); xfrac = float(params.get("xfrac", "0.01"))
    if nrep < 1 or not (0 < xfrac <= 1):
        raise InputError("nrep must be positive and xfrac must be in (0, 1]")
    steps = tuple(x.strip() for x in params.get("marsteps", "data,gm,sfs,mar,ext,plot").split(",") if x.strip())
    maxsnps_value = params.get("maxsnps", "auto").strip().lower()
    maxsnps = None if maxsnps_value in {"", "auto", "all"} else int(maxsnps_value)
    if maxsnps is not None and maxsnps < 1:
        raise InputError("maxsnps must be positive or auto")
    geom_id = int(analysis.get("geom_id", "7"))
    return MarConfig(path, vcf, lonlat, scenario_file, output, analysis.get("name", "bioma_mar").strip() or "bioma_mar", geom_id, scheme, nrep, xfrac, params.get("quorum", "true").strip().lower() not in {"0", "false", "no"}, int(params.get("randseed", "123")), maxsnps, steps, _rscript(params.get("rscript", "")))


def count_vcf_sites(path: Path) -> int:
    opener = gzip.open if path.suffix.lower() in {".gz", ".bgz", ".bgzip"} else open
    count = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line and not line.startswith("#"):
                count += 1
    if count < 1:
        raise InputError("VCF contains no variant records: {}".format(path))
    return count


def run_mar_workflow(config_path: Path, dry_run: bool = False, progress: Optional[Callable[[str], None]] = None) -> Dict[str, object]:
    config = load_mar_config(config_path)
    maxsnps = config.maxsnps if config.maxsnps is not None else count_vcf_sites(config.vcf)
    manifest: Dict[str, object] = {"module": "mar-workflow", "status": "planned" if dry_run else "running", "config": config.payload(), "resolved_maxsnps": maxsnps}
    config.output_dir.mkdir(parents=True, exist_ok=True)
    if dry_run:
        (config.output_dir / "mar_dry_run.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        return manifest
    script = Path(__file__).resolve().parent / "scripts" / "mar_compute.R"
    plot_script = Path(__file__).resolve().parent / "scripts" / "mar_plot.R"
    started = time.time()
    if progress: progress("MAR pipeline: scheme={}, maxsnps={}, nrep={}".format(config.scheme, maxsnps, config.nrep))
    args = [config.rscript, str(script), str(config.name), str(config.output_dir), str(config.vcf), str(config.lonlat), config.scheme, str(config.nrep), str(config.xfrac), "TRUE" if config.quorum else "FALSE", str(config.randseed), str(maxsnps), ",".join(config.marsteps)]
    result = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(config.rscript))
    (config.output_dir / "mar_compute.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise InputError("MAR computation failed; see {}".format(config.output_dir / "mar_compute.log"))
    if progress: progress("MAR plotting")
    plot_args = [config.rscript, str(plot_script), str(config.output_dir), str(config.name), str(config.scenario_file or ""), str(config.geom_id)]
    result = subprocess.run(plot_args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(config.rscript))
    (config.output_dir / "mar_plot.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise InputError("MAR plotting failed; see {}".format(config.output_dir / "mar_plot.log"))
    manifest.update({"status": "complete", "elapsed_seconds": round(time.time() - started, 3), "outputs": {"output_dir": str(config.output_dir)}})
    manifest["config_sha256"] = hashlib.sha256(config.config_path.read_bytes()).hexdigest()
    (config.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest
