"""Subprocess environment helpers for reproducible Conda runtimes."""

from __future__ import annotations

import os
import platform
import shutil
from pathlib import Path
from typing import Dict


def subprocess_environment(executable: str) -> Dict[str, str]:
    """Use a Conda executable with libraries and geospatial data from its prefix."""
    environment = os.environ.copy()
    discovered = shutil.which(executable)
    path = Path(discovered or executable).expanduser()
    try:
        path = path.resolve()
    except OSError:
        return environment
    prefix = path.parent.parent
    if not (prefix / "conda-meta").is_dir():
        return environment
    environment["CONDA_PREFIX"] = str(prefix)
    environment["PATH"] = str(prefix / "bin") + os.pathsep + environment.get("PATH", "")
    library_dir = prefix / "lib"
    if library_dir.is_dir() and platform.system() == "Linux":
        environment["LD_LIBRARY_PATH"] = str(library_dir)
    proj_dir = prefix / "share" / "proj"
    if (proj_dir / "proj.db").is_file():
        environment["PROJ_DATA"] = str(proj_dir)
        environment["PROJ_LIB"] = str(proj_dir)
    gdal_dir = prefix / "share" / "gdal"
    if gdal_dir.is_dir():
        environment["GDAL_DATA"] = str(gdal_dir)
    return environment
