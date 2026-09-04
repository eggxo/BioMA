"""Content fingerprints shared by standalone BioMA module manifests.

The project runner has a detailed provenance implementation of its own.  The
standalone commands need the same guarantees when they are run outside a
project, so this small dependency-free module exposes the stable parts of that
contract without importing the project orchestrator.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional

from .gf_frequency import InputError


HASH_ALGORITHM = "sha256"
DIRECTORY_HASH_DESCRIPTION = "sha256(canonical-json-lines)"


def sha256_file(path: Path) -> str:
    """Hash a file in bounded-size chunks."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as error:
        raise InputError("Cannot read input for provenance: {} ({})".format(path, error))
    return digest.hexdigest()


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _aggregate(entries: Iterable[Mapping[str, object]]) -> str:
    digest = hashlib.sha256()
    for entry in entries:
        digest.update(_canonical_bytes(entry))
        digest.update(b"\n")
    return digest.hexdigest()


def _shapefile_components(path: Path):
    """Return the primary shapefile and common same-stem sidecars."""
    if path.suffix.lower() != ".shp":
        return [path]
    stem = path.stem.casefold()
    suffixes = {
        ".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj", ".qix", ".fix",
        ".sbn", ".sbx", ".fbn", ".fbx", ".ain", ".aih", ".atx", ".ixs",
        ".mxs", ".xml",
    }
    try:
        candidates = [
            item for item in path.parent.iterdir()
            if item.is_file()
            and (
                item == path
                or item.name.casefold() == path.name.casefold()
                or (item.stem.casefold() == stem and item.suffix.casefold() in suffixes)
                or item.name.casefold() == stem + ".shp.xml"
            )
        ]
    except OSError as error:
        raise InputError("Cannot inspect input directory {}: {}".format(path.parent, error))
    if path not in candidates:
        candidates.append(path)
    return sorted({item.resolve() for item in candidates}, key=lambda item: item.name.casefold())


def _directory_entries(root: Path):
    """Yield deterministic records for every member of *root*."""
    root = root.resolve()

    def visit(directory: Path, relative: str):
        try:
            entries = sorted(directory.iterdir(), key=lambda item: item.name.casefold())
        except OSError as error:
            raise InputError("Cannot read input directory {}: {}".format(directory, error))
        for item in entries:
            rel = "{}/{}".format(relative, item.name) if relative else item.name
            try:
                if item.is_symlink():
                    yield {
                        "relative_path": rel.replace(os.sep, "/"),
                        "kind": "symlink",
                        "target": os.readlink(str(item)).replace(os.sep, "/"),
                    }
                elif item.is_dir():
                    children = iter(visit(item, rel))
                    try:
                        first = next(children)
                    except StopIteration:
                        yield {"relative_path": rel.replace(os.sep, "/"), "kind": "directory"}
                    else:
                        yield first
                        yield from children
                elif item.is_file():
                    stat = item.stat()
                    yield {
                        "relative_path": rel.replace(os.sep, "/"),
                        "kind": "file",
                        "size_bytes": int(stat.st_size),
                        "sha256": sha256_file(item),
                    }
                else:
                    yield {"relative_path": rel.replace(os.sep, "/"), "kind": "special"}
            except OSError as error:
                raise InputError("Cannot inspect input {}: {}".format(item, error))

    yield from visit(root, "")


def fingerprint_path(path: Path, *, required: bool = True) -> Dict[str, object]:
    """Return a stable content fingerprint for a file or directory.

    ``required=False`` is useful for optional configuration inputs: an omitted
    value is represented explicitly instead of silently disappearing from the
    manifest.  A configured-but-missing path remains an error because it would
    otherwise make a run non-reproducible.
    """
    if path is None:  # type: ignore[comparison-overlap]
        if required:
            raise InputError("Required provenance input is empty")
        return {"path": None, "kind": "absent", "sha256": None}
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        if required:
            raise InputError("Input does not exist: {}".format(resolved))
        return {"path": str(resolved), "kind": "missing", "sha256": None}
    if resolved.is_symlink():
        return {
            "path": str(resolved),
            "kind": "symlink",
            "target": os.readlink(str(resolved)).replace(os.sep, "/"),
            "sha256": _aggregate(
                [{"kind": "symlink", "target": os.readlink(str(resolved)).replace(os.sep, "/")}]
            ),
        }
    if resolved.is_dir():
        entries = list(_directory_entries(resolved))
        entries.sort(key=lambda item: (str(item.get("relative_path", "")).casefold(), str(item.get("relative_path", ""))))
        return {
            "path": str(resolved),
            "kind": "directory",
            "sha256": _aggregate(entries),
            "entry_count": len(entries),
            "total_size_bytes": sum(int(item.get("size_bytes", 0) or 0) for item in entries),
            "files": entries,
        }
    if not resolved.is_file():
        return {"path": str(resolved), "kind": "special", "sha256": None}
    components = []
    for component in _shapefile_components(resolved):
        try:
            stat = component.stat()
        except OSError as error:
            raise InputError("Cannot inspect input {}: {}".format(component, error))
        components.append({
            "relative_path": component.name,
            "kind": "file",
            "size_bytes": int(stat.st_size),
            "sha256": sha256_file(component),
        })
    components.sort(key=lambda item: str(item["relative_path"]).casefold())
    primary = next((row for row in components if str(row["relative_path"]).casefold() == resolved.name.casefold()), components[0])
    return {
        "path": str(resolved),
        "kind": "file",
        "size_bytes": int(primary["size_bytes"]),
        "sha256": _aggregate(components) if len(components) > 1 else str(primary["sha256"]),
        "primary_sha256": str(primary["sha256"]),
        "component_count": len(components),
        **({"components": components} if len(components) > 1 else {}),
    }


def fingerprint_inputs(inputs: Mapping[str, Optional[Path]]) -> Dict[str, Dict[str, object]]:
    """Fingerprint named inputs in stable key order."""
    return {
        str(label): fingerprint_path(path, required=path is not None)
        for label, path in sorted(inputs.items(), key=lambda item: str(item[0]))
    }


def attach_input_fingerprints(
    manifest: Dict[str, object], inputs: Mapping[str, Optional[Path]]
) -> Dict[str, object]:
    """Attach complete input fingerprints and an aggregate signature."""
    fingerprints = fingerprint_inputs(inputs)
    manifest["input_fingerprints"] = fingerprints
    # Keep the records discoverable beside the legacy semantic ``inputs``
    # payload as well as at the top level.  Existing consumers can continue to
    # read their old keys, while new consumers have one consistent location.
    if isinstance(manifest.get("inputs"), dict):
        manifest["inputs"]["fingerprints"] = fingerprints  # type: ignore[index]
    input_signature = _aggregate(
        {"label": label, "fingerprint": fingerprints[label]}
        for label in sorted(fingerprints)
    )
    manifest["input_signature_sha256"] = input_signature
    manifest["provenance"] = {
        "hash_algorithm": HASH_ALGORITHM,
        "directory_hash": DIRECTORY_HASH_DESCRIPTION,
        "input_signature_sha256": input_signature,
    }
    return manifest
