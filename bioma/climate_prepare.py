"""Validate and standardize present/future bioclimatic inputs for BioMA."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from . import __version__
from .gf_frequency import InputError
from .runtime import subprocess_environment


BIO_COUNT = 19
SCENARIO_RE = re.compile(
    r"^(?P<period>\d{4}-\d{4})-ssp(?P<ssp>\d+)-(?P<model>.+)$",
    re.IGNORECASE,
)
BIO_RE = re.compile(r"bio[_-]?(?P<number>\d{1,2})(?=\D|$)", re.IGNORECASE)


@dataclass(frozen=True)
class Scenario:
    name: str
    period: str
    ssp: str
    model: str
    directory: Path
    rasters: Tuple[Path, ...]


@dataclass(frozen=True)
class RasterGeometry:
    width: int
    height: int
    transform: Tuple[float, ...]
    projection: str

    @property
    def xmin(self) -> float:
        return self.transform[0]

    @property
    def ymax(self) -> float:
        return self.transform[3]

    @property
    def xmax(self) -> float:
        return self.xmin + self.width * self.transform[1]

    @property
    def ymin(self) -> float:
        return self.ymax + self.height * self.transform[5]


@dataclass(frozen=True)
class PopulationPoint:
    population_id: str
    group: str
    lon: float
    lat: float
    supplied_bio: Tuple[Optional[float], ...]


def _run(
    command: Sequence[str],
    stdin: Optional[str] = None,
    description: Optional[str] = None,
    isolated: bool = False,
) -> subprocess.CompletedProcess:
    try:
        result = subprocess.run(
            list(command),
            input=stdin,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=subprocess_environment(command[0]) if isolated else None,
        )
    except FileNotFoundError:
        raise InputError("Required program was not found: {}".format(command[0]))
    if result.returncode != 0:
        label = description or "Command"
        detail = result.stderr.strip() or result.stdout.strip() or "unknown error"
        raise InputError("{} failed: {}".format(label, detail))
    return result


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def discover_bioclim_rasters(directory: Path) -> Tuple[Path, ...]:
    """Return BIO1..BIO19 in numeric order and reject incomplete collections."""
    if not directory.is_dir():
        raise InputError("Climate directory does not exist: {}".format(directory))
    found: Dict[int, Path] = {}
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.suffix.lower() not in (".tif", ".tiff"):
            continue
        matches = list(BIO_RE.finditer(path.stem))
        if not matches:
            continue
        number = int(matches[-1].group("number"))
        if not 1 <= number <= BIO_COUNT:
            continue
        if number in found:
            raise InputError(
                "More than one BIO{} raster in {}: {}, {}".format(
                    number, directory, found[number].name, path.name
                )
            )
        found[number] = path.resolve()
    missing = [number for number in range(1, BIO_COUNT + 1) if number not in found]
    if missing:
        raise InputError(
            "Climate directory {} is missing BIO rasters: {}".format(
                directory, ", ".join(map(str, missing))
            )
        )
    return tuple(found[number] for number in range(1, BIO_COUNT + 1))


def discover_scenarios(
    future_root: Path,
    models: Sequence[str] = (),
    ssps: Sequence[str] = (),
    periods: Sequence[str] = (),
) -> List[Scenario]:
    if not future_root.is_dir():
        raise InputError("Future-climate root does not exist: {}".format(future_root))
    model_filter = set(models)
    ssp_filter = {str(value).lower().replace("ssp", "") for value in ssps}
    period_filter = set(periods)
    available: Dict[Tuple[str, str, str], Path] = {}
    selected: List[Scenario] = []
    for directory in sorted(path for path in future_root.iterdir() if path.is_dir()):
        match = SCENARIO_RE.match(directory.name)
        if not match:
            continue
        period = match.group("period")
        ssp = match.group("ssp")
        model = match.group("model")
        key = (model, ssp, period)
        if key in available:
            raise InputError("Duplicate future scenario: {}".format(key))
        available[key] = directory
        if model_filter and model not in model_filter:
            continue
        if ssp_filter and ssp not in ssp_filter:
            continue
        if period_filter and period not in period_filter:
            continue
        selected.append(
            Scenario(
                name=directory.name,
                period=period,
                ssp=ssp,
                model=model,
                directory=directory.resolve(),
                rasters=discover_bioclim_rasters(directory),
            )
        )
    if not available:
        raise InputError("No future scenario directories were found under {}".format(future_root))
    if not selected:
        raise InputError("No future scenarios match the requested model/SSP/period filters")

    available_models = {key[0] for key in available}
    available_ssps = {key[1] for key in available}
    available_periods = {key[2] for key in available}
    missing_models = model_filter - available_models
    missing_ssps = ssp_filter - available_ssps
    missing_periods = period_filter - available_periods
    if missing_models or missing_ssps or missing_periods:
        parts = []
        if missing_models:
            parts.append("models={}".format(",".join(sorted(missing_models))))
        if missing_ssps:
            parts.append("SSPs={}".format(",".join(sorted(missing_ssps))))
        if missing_periods:
            parts.append("periods={}".format(",".join(sorted(missing_periods))))
        raise InputError("Requested future data are absent: {}".format("; ".join(parts)))

    if model_filter and ssp_filter and period_filter:
        selected_keys = {(item.model, item.ssp, item.period) for item in selected}
        expected_keys = {
            (model, ssp, period)
            for model in model_filter
            for ssp in ssp_filter
            for period in period_filter
        }
        absent = sorted(expected_keys - selected_keys)
        if absent:
            preview = ", ".join("{}/{}/{}".format(*key) for key in absent[:10])
            raise InputError("Requested scenario combinations are absent: {}".format(preview))
    return selected


def read_population_points(path: Path) -> List[PopulationPoint]:
    if not path.is_file():
        raise InputError("Coordinate file does not exist: {}".format(path))
    lines = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    if len(lines) < 2:
        raise InputError("Coordinate file contains no population rows")
    header = re.split(r"\s+", lines[0])
    normalized = {name.lower(): index for index, name in enumerate(header)}
    id_index = normalized.get("population_id", normalized.get("id"))
    group_index = normalized.get("group", normalized.get("pop"))
    lon_index = normalized.get("lon", normalized.get("longitude"))
    lat_index = normalized.get("lat", normalized.get("latitude"))
    if id_index is None or lon_index is None or lat_index is None:
        raise InputError(
            "Coordinate file must contain ID/population_id, lon, and lat columns"
        )
    bio_indexes = [normalized.get("bio{}".format(number)) for number in range(1, BIO_COUNT + 1)]
    points: List[PopulationPoint] = []
    seen = set()
    for line_number, line in enumerate(lines[1:], start=2):
        fields = re.split(r"\s+", line)
        if len(fields) != len(header):
            raise InputError(
                "Coordinate line {} has {} fields; expected {}".format(
                    line_number, len(fields), len(header)
                )
            )
        population_id = fields[id_index]
        if population_id in seen:
            raise InputError("Duplicate population ID in coordinates: {}".format(population_id))
        seen.add(population_id)
        try:
            lon = float(fields[lon_index])
            lat = float(fields[lat_index])
        except ValueError:
            raise InputError("Non-numeric longitude/latitude at line {}".format(line_number))
        if not (-180 <= lon <= 180 and -90 <= lat <= 90):
            raise InputError("Invalid longitude/latitude at line {}".format(line_number))
        supplied: List[Optional[float]] = []
        for index in bio_indexes:
            if index is None or fields[index].upper() in ("NA", "NAN", "."):
                supplied.append(None)
            else:
                try:
                    supplied.append(float(fields[index]))
                except ValueError:
                    raise InputError("Non-numeric BIO value at line {}".format(line_number))
        group = fields[group_index] if group_index is not None else ""
        points.append(PopulationPoint(population_id, group, lon, lat, tuple(supplied)))
    return points


def _gdal_geometry(gdalinfo: str, raster: Path) -> RasterGeometry:
    # GDAL 1.x, which is still common on HPC servers, has no -json option.
    try:
        json_result = subprocess.run(
            [gdalinfo, "-json", str(raster)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        raise InputError("Required program was not found: {}".format(gdalinfo))
    if json_result.returncode == 0:
        try:
            info = json.loads(json_result.stdout)
            width, height = info["size"]
            transform = tuple(float(value) for value in info["geoTransform"])
            projection = info.get("coordinateSystem", {}).get("wkt", "")
            band_count = len(info["bands"])
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
            raise InputError("Cannot parse raster metadata for {}: {}".format(raster, error))
    else:
        text_result = _run([gdalinfo, str(raster)], description="gdalinfo")
        size_match = re.search(r"^Size is\s+(\d+)\s*,\s*(\d+)\s*$", text_result.stdout, re.MULTILINE)
        origin_match = re.search(
            r"^Origin\s*=\s*\(([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\)\s*$",
            text_result.stdout,
            re.MULTILINE,
        )
        pixel_match = re.search(
            r"^Pixel Size\s*=\s*\(([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\)\s*$",
            text_result.stdout,
            re.MULTILINE,
        )
        projection_match = re.search(
            r"Coordinate System is:\s*(.*?)\nOrigin\s*=",
            text_result.stdout,
            re.DOTALL,
        )
        if not size_match or not origin_match or not pixel_match or not projection_match:
            raise InputError("Cannot parse legacy gdalinfo output for {}".format(raster))
        width, height = int(size_match.group(1)), int(size_match.group(2))
        origin_x, origin_y = float(origin_match.group(1)), float(origin_match.group(2))
        pixel_x, pixel_y = float(pixel_match.group(1)), float(pixel_match.group(2))
        transform = (origin_x, pixel_x, 0.0, origin_y, 0.0, pixel_y)
        projection = re.sub(r"\s+", "", projection_match.group(1))
        band_count = len(re.findall(r"^Band\s+\d+\b", text_result.stdout, re.MULTILINE))
    if band_count != 1:
        raise InputError("BIO raster must have one band: {}".format(raster))
    if len(transform) != 6 or transform[2] != 0 or transform[4] != 0:
        raise InputError("Rotated or invalid raster grid is unsupported: {}".format(raster))
    if not projection:
        raise InputError("Raster has no coordinate reference system: {}".format(raster))
    return RasterGeometry(int(width), int(height), transform, projection)


def _close(left: float, right: float, tolerance: float = 1e-9) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def _same_grid(left: RasterGeometry, right: RasterGeometry) -> bool:
    return (
        left.width == right.width
        and left.height == right.height
        and all(_close(a, b) for a, b in zip(left.transform, right.transform))
        and left.projection == right.projection
    )


def _aligned_to(reference: RasterGeometry, candidate: RasterGeometry) -> bool:
    if reference.projection != candidate.projection:
        return False
    if not _close(reference.transform[1], candidate.transform[1]):
        return False
    if not _close(reference.transform[5], candidate.transform[5]):
        return False
    x_steps = (candidate.xmin - reference.xmin) / reference.transform[1]
    y_steps = (candidate.ymax - reference.ymax) / reference.transform[5]
    return _close(x_steps, round(x_steps), 1e-7) and _close(y_steps, round(y_steps), 1e-7)


def _validate_collection(
    rasters: Sequence[Path], gdalinfo: str
) -> Tuple[RasterGeometry, List[RasterGeometry]]:
    geometries = [_gdal_geometry(gdalinfo, path) for path in rasters]
    first = geometries[0]
    for raster, geometry in zip(rasters[1:], geometries[1:]):
        if not _same_grid(first, geometry):
            raise InputError("BIO rasters do not share one grid: {}".format(raster))
    return first, geometries


def _extract_population_environment(
    points: Sequence[PopulationPoint],
    rasters: Sequence[Path],
    gdallocationinfo: str,
) -> List[List[float]]:
    coordinates = "".join("{:.12g} {:.12g}\n".format(point.lon, point.lat) for point in points)
    columns: List[List[float]] = []
    for raster in rasters:
        result = _run(
            [gdallocationinfo, "-valonly", "-geoloc", str(raster)],
            stdin=coordinates,
            description="Population climate extraction",
        )
        values = result.stdout.splitlines()
        if len(values) != len(points):
            raise InputError(
                "Raster extraction returned {} values for {} population points: {}".format(
                    len(values), len(points), raster
                )
            )
        try:
            columns.append([float(value.strip()) for value in values])
        except ValueError:
            raise InputError("A population point falls outside or on NoData in {}".format(raster))
    return [[columns[bio][row] for bio in range(BIO_COUNT)] for row in range(len(points))]


def _resolve_rscript(requested: Optional[str]) -> str:
    candidates = []
    if requested:
        candidates.append(requested)
    discovered = shutil.which("Rscript")
    if discovered:
        candidates.append(discovered)
    checked = set()
    for candidate in candidates:
        if candidate in checked:
            continue
        checked.add(candidate)
        try:
            result = subprocess.run(
                [candidate, "--vanilla", "-e", "suppressPackageStartupMessages(library(raster))"],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=subprocess_environment(candidate),
            )
        except FileNotFoundError:
            continue
        if result.returncode == 0:
            return candidate
    raise InputError(
        "No Rscript with the 'raster' package was found; provide one with --rscript"
    )


def _validate_mask(mask: Path, ogrinfo: str) -> None:
    if not mask.is_file() or mask.suffix.lower() != ".shp":
        raise InputError("Mask must be an existing ESRI Shapefile (.shp): {}".format(mask))
    for suffix in (".dbf", ".shx"):
        companion = mask.with_suffix(suffix)
        if not companion.is_file():
            raise InputError("Mask companion file is missing: {}".format(companion))
    _run(
        [ogrinfo, "-ro", "-so", str(mask), mask.stem],
        description="Mask validation",
    )


def _extract_masked_background(
    rscript: str,
    helper_script: Path,
    mask: Path,
    rasters: Sequence[Path],
    output: Path,
) -> None:
    command = [
        rscript,
        "--vanilla",
        str(helper_script),
        "mask",
        str(mask),
        str(output),
    ] + [str(path) for path in rasters]
    _run(command, description="Masked climate extraction", isolated=True)


def _extract_at_coordinates(
    rscript: str,
    helper_script: Path,
    coordinates: Path,
    rasters: Sequence[Path],
    output: Path,
) -> None:
    command = [
        rscript,
        "--vanilla",
        str(helper_script),
        "extract",
        str(coordinates),
        str(output),
    ] + [str(path) for path in rasters]
    _run(command, description="Climate extraction at canonical grid cells", isolated=True)


def _background_profile(path: Path) -> Tuple[int, str, int]:
    digest = hashlib.sha256()
    count = 0
    missing_rows = 0
    with gzip.open(str(path), "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        expected = ["lon", "lat"] + ["bio{}".format(i) for i in range(1, BIO_COUNT + 1)]
        if header != expected:
            raise InputError("Unexpected background table columns in {}".format(path))
        for row in reader:
            if len(row) != len(expected):
                raise InputError("Malformed background row in {}".format(path))
            digest.update((row[0] + "\t" + row[1] + "\n").encode("utf-8"))
            if any(value.strip().upper() in ("NA", "NAN", "") for value in row[2:]):
                missing_rows += 1
            count += 1
    if count == 0:
        raise InputError("Mask produced no climate cells: {}".format(path))
    return count, digest.hexdigest(), missing_rows


def prepare_climate_inputs(
    present_dir: Path,
    future_root: Path,
    coordinates_path: Path,
    current_mask: Path,
    output_dir: Path,
    future_mask: Optional[Path] = None,
    models: Sequence[str] = (),
    ssps: Sequence[str] = (),
    periods: Sequence[str] = (),
    validate_only: bool = False,
    gdalinfo: str = "gdalinfo",
    gdallocationinfo: str = "gdallocationinfo",
    ogrinfo: str = "ogrinfo",
    rscript: Optional[str] = None,
    supplied_bio_tolerance: float = 1e-6,
) -> Dict[str, object]:
    """Create standardized climate tables and manifests without modifying sources."""
    present_dir = present_dir.expanduser().resolve()
    future_root = future_root.expanduser().resolve()
    coordinates_path = coordinates_path.expanduser().resolve()
    current_mask = current_mask.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    future_mask = (future_mask or current_mask).expanduser().resolve()
    if output_dir.exists():
        raise InputError("Output directory already exists; choose a new directory: {}".format(output_dir))
    if supplied_bio_tolerance < 0:
        raise InputError("supplied_bio_tolerance cannot be negative")

    present_rasters = discover_bioclim_rasters(present_dir)
    scenarios = discover_scenarios(future_root, models, ssps, periods)
    points = read_population_points(coordinates_path)
    _validate_mask(current_mask, ogrinfo)
    if future_mask != current_mask:
        _validate_mask(future_mask, ogrinfo)

    present_geometry, present_geometries = _validate_collection(present_rasters, gdalinfo)
    for point in points:
        if not (
            present_geometry.xmin <= point.lon <= present_geometry.xmax
            and present_geometry.ymin <= point.lat <= present_geometry.ymax
        ):
            raise InputError("Population {} falls outside the present climate grid".format(point.population_id))

    scenario_geometries: Dict[str, Tuple[RasterGeometry, List[RasterGeometry]]] = {}
    for scenario in scenarios:
        geometry, all_geometries = _validate_collection(scenario.rasters, gdalinfo)
        if not _aligned_to(present_geometry, geometry):
            raise InputError("Future grid is not aligned to the present grid: {}".format(scenario.name))
        scenario_geometries[scenario.name] = (geometry, all_geometries)

    population_environment = _extract_population_environment(
        points, present_rasters, gdallocationinfo
    )
    warnings: List[Tuple[str, str, str]] = []
    max_difference = 0.0
    compared_values = 0
    for point, extracted in zip(points, population_environment):
        for supplied, actual in zip(point.supplied_bio, extracted):
            if supplied is not None:
                max_difference = max(max_difference, abs(supplied - actual))
                compared_values += 1
    if compared_values and max_difference > supplied_bio_tolerance:
        warnings.append(
            (
                "SUPPLIED_BIO_MISMATCH",
                "coordinates",
                "Largest difference between supplied and re-extracted BIO values is {:.12g}".format(max_difference),
            )
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(tempfile.mkdtemp(prefix=".{}-tmp-".format(output_dir.name), dir=str(output_dir.parent)))
    try:
        _write_tsv(
            temporary_dir / "population_environment.tsv",
            ["population_id", "group", "lon", "lat"]
            + ["bio{}".format(i) for i in range(1, BIO_COUNT + 1)],
            (
                [point.population_id, point.group, "{:.12g}".format(point.lon), "{:.12g}".format(point.lat)]
                + ["{:.12g}".format(value) for value in values]
                for point, values in zip(points, population_environment)
            ),
        )

        _write_tsv(
            temporary_dir / "future_manifest.tsv",
            ["scenario", "model", "ssp", "period", "directory"],
            ([s.name, s.model, "ssp{}".format(s.ssp), s.period, str(s.directory)] for s in scenarios),
        )

        raster_rows = []
        for bio_index, (path, geometry) in enumerate(zip(present_rasters, present_geometries), start=1):
            stat = path.stat()
            raster_rows.append(
                ["present", "present", "", "", "", "bio{}".format(bio_index), str(path), geometry.width, geometry.height,
                 "{:.12g}".format(geometry.xmin), "{:.12g}".format(geometry.ymin), "{:.12g}".format(geometry.xmax),
                 "{:.12g}".format(geometry.ymax), "{:.12g}".format(geometry.transform[1]),
                 "{:.12g}".format(abs(geometry.transform[5])), stat.st_size, int(stat.st_mtime)]
            )
        for scenario in scenarios:
            _, geometries = scenario_geometries[scenario.name]
            for bio_index, (path, geometry) in enumerate(zip(scenario.rasters, geometries), start=1):
                stat = path.stat()
                raster_rows.append(
                    ["future", scenario.name, scenario.model, "ssp{}".format(scenario.ssp), scenario.period,
                     "bio{}".format(bio_index), str(path), geometry.width, geometry.height,
                     "{:.12g}".format(geometry.xmin), "{:.12g}".format(geometry.ymin), "{:.12g}".format(geometry.xmax),
                     "{:.12g}".format(geometry.ymax), "{:.12g}".format(geometry.transform[1]),
                     "{:.12g}".format(abs(geometry.transform[5])), stat.st_size, int(stat.st_mtime)]
                )
        _write_tsv(
            temporary_dir / "raster_manifest.tsv",
            ["climate", "scenario", "model", "ssp", "period", "variable", "path", "width", "height",
             "xmin", "ymin", "xmax", "ymax", "resolution_x", "resolution_y", "size_bytes", "mtime_epoch"],
            raster_rows,
        )

        background_counts: Dict[str, int] = {}
        background_nodata_rows: Dict[str, int] = {}
        if not validate_only:
            resolved_rscript = _resolve_rscript(rscript)
            helper_script = Path(__file__).resolve().parent / "scripts" / "extract_masked_background.R"
            if not helper_script.is_file():
                raise InputError("Bundled background-extraction helper is missing: {}".format(helper_script))
            current_output = temporary_dir / "current_background.tsv.gz"
            _extract_masked_background(resolved_rscript, helper_script, current_mask, present_rasters, current_output)
            current_count, current_coordinate_hash, current_missing = _background_profile(current_output)
            background_counts["current"] = current_count
            background_nodata_rows["current"] = current_missing
            local_dir = temporary_dir / "future_local"
            search_dir = temporary_dir / "future_search"
            local_dir.mkdir()
            search_dir.mkdir()
            same_masks = current_mask == future_mask
            if same_masks:
                search_reference = current_output
            else:
                search_reference = temporary_dir / ".future_search_reference.tsv.gz"
                _extract_masked_background(
                    resolved_rscript,
                    helper_script,
                    future_mask,
                    present_rasters,
                    search_reference,
                )
                _background_profile(search_reference)
            for scenario in scenarios:
                local_output = local_dir / (scenario.name + ".tsv.gz")
                _extract_at_coordinates(
                    resolved_rscript,
                    helper_script,
                    current_output,
                    scenario.rasters,
                    local_output,
                )
                local_count, local_coordinate_hash, local_missing = _background_profile(local_output)
                if local_count != current_count or local_coordinate_hash != current_coordinate_hash:
                    raise InputError(
                        "Future local grid does not match current background cells: {}".format(scenario.name)
                    )
                background_counts["local/{}".format(scenario.name)] = local_count
                background_nodata_rows["local/{}".format(scenario.name)] = local_missing
                if local_missing:
                    warnings.append(
                        (
                            "FUTURE_LOCAL_NODATA_CELLS",
                            scenario.name,
                            "{} of {} current-mask cells contain future NoData; local offset will be NA there and incomplete search candidates must be excluded".format(
                                local_missing, local_count
                            ),
                        )
                    )
                search_output = search_dir / (scenario.name + ".tsv.gz")
                if same_masks:
                    os.link(str(local_output), str(search_output))
                    background_counts["search/{}".format(scenario.name)] = local_count
                    background_nodata_rows["search/{}".format(scenario.name)] = local_missing
                else:
                    _extract_at_coordinates(
                        resolved_rscript,
                        helper_script,
                        search_reference,
                        scenario.rasters,
                        search_output,
                    )
                    search_count, _, search_missing = _background_profile(search_output)
                    background_counts["search/{}".format(scenario.name)] = search_count
                    background_nodata_rows["search/{}".format(scenario.name)] = search_missing
                    if search_missing:
                        warnings.append(
                            (
                                "FUTURE_SEARCH_NODATA_CELLS",
                                scenario.name,
                                "{} of {} future-search cells contain NoData and must be excluded as offset destinations".format(
                                    search_missing, search_count
                                ),
                            )
                        )
            if not same_masks:
                search_reference.unlink()
        else:
            resolved_rscript = None

        _write_tsv(temporary_dir / "warnings.tsv", ["warning_code", "scope", "message"], warnings)
        manifest: Dict[str, object] = {
            "module": "climate-prepare",
            "bioma_version": __version__,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "engineering_status": "completed",
            "scientific_status": "completed_with_warning" if warnings else "completed_validated",
            "inputs": {
                "present_dir": str(present_dir),
                "future_root": str(future_root),
                "coordinates": str(coordinates_path),
                "coordinates_sha256": _sha256(coordinates_path),
                "current_mask": str(current_mask),
                "current_mask_sha256": _sha256(current_mask),
                "future_mask": str(future_mask),
                "future_mask_sha256": _sha256(future_mask),
            },
            "parameters": {
                "bio_variables": ["bio{}".format(i) for i in range(1, BIO_COUNT + 1)],
                "models": list(models),
                "ssps": list(ssps),
                "periods": list(periods),
                "validate_only": validate_only,
                "supplied_bio_tolerance": supplied_bio_tolerance,
            },
            "counts": {
                "populations": len(points),
                "future_scenarios": len(scenarios),
                "raster_files": len(raster_rows),
                "warnings": len(warnings),
                "background_cells": background_counts,
                "background_nodata_rows": background_nodata_rows,
            },
            "checks": {
                "grids_aligned": True,
                "population_bio_values_compared": compared_values,
                "population_bio_max_absolute_difference": max_difference,
            },
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "rscript": resolved_rscript,
            },
        }
        with (temporary_dir / "run_manifest.json").open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        os.replace(str(temporary_dir), str(output_dir))
        return manifest
    except Exception:
        shutil.rmtree(str(temporary_dir), ignore_errors=True)
        raise
