#!/usr/bin/env python3
"""Compute 2-D deme WFmoments projections.

The command is intentionally independent of the BioMA launcher environment:
the caller can point it at the Python environment that contains ``wfmoments``
and rasterio.  Inputs are ordinary binary habitat masks and small metadata
tables, so every intermediate decision can be audited from the output TSVs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import traceback
from collections import deque
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import rasterio
import wfmoments as wfm


SCENARIO_RE = re.compile(r"(\d{4}-\d{4}).*?(ssp\d+)", re.IGNORECASE)


def read_table_auto(path: Path) -> pd.DataFrame:
    """Read comma- or tab-delimited tables without relying on the suffix."""
    errors: List[Exception] = []
    for sep in (",", "\t"):
        try:
            frame = pd.read_csv(path, sep=sep)
            if frame.shape[1] > 1:
                return frame
        except Exception as exc:  # pragma: no cover - fallback diagnostics
            errors.append(exc)
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception as exc:
        errors.append(exc)
        raise ValueError("Cannot read table {}: {}".format(path, errors[-1]))


def finite_float(value: object, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        raise ValueError("{} must be numeric, got {!r}".format(label, value))
    if not np.isfinite(result):
        raise ValueError("{} must be finite".format(label))
    return result


def scalarize(value: object) -> float:
    arr = np.asarray(value, dtype=float)
    if arr.size == 0:
        return float("nan")
    if arr.size == 1:
        return float(arr.reshape(-1)[0])
    return float(np.nanmean(arr))


def read_species_value(path: Path, species: str, field: str) -> float:
    frame = read_table_auto(path)
    if "species" not in frame.columns:
        raise ValueError("{} is missing the 'species' column".format(path))
    if field not in frame.columns:
        raise ValueError("{} is missing the '{}' column".format(path, field))
    hit = frame.loc[frame["species"].astype(str) == str(species), field]
    if hit.empty:
        raise ValueError("Species '{}' is not present in {}".format(species, path))
    return finite_float(hit.iloc[0], "{} for {}".format(field, species))


def read_binary_raster(path: Path) -> Tuple[np.ndarray, np.ndarray, rasterio.Affine, Dict[str, object]]:
    with rasterio.open(path) as src:
        raw = src.read(1, masked=True)
        values = np.asarray(raw.filled(0), dtype=float)
        valid = ~np.asarray(raw.mask, dtype=bool)
        valid &= np.isfinite(values)
        occupied = (values > 0) & valid
        transform = src.transform
        info = {
            "width": int(src.width),
            "height": int(src.height),
            "crs": str(src.crs) if src.crs else None,
            "transform": tuple(float(x) for x in transform),
            "nodata": None if src.nodata is None else float(src.nodata),
        }
    return occupied, valid, transform, info


def split_edges(length: int, parts: int) -> np.ndarray:
    if parts < 1 or parts > length:
        raise ValueError("Grid dimension {} must be between 1 and raster size {}".format(parts, length))
    edges = np.rint(np.linspace(0, length, parts + 1)).astype(int)
    edges[0] = 0
    edges[-1] = length
    if np.any(np.diff(edges) <= 0):
        raise ValueError("Grid dimension creates empty deme blocks")
    return edges


def build_deme_table(
    occupied: np.ndarray,
    valid: np.ndarray,
    transform: rasterio.Affine,
    nx: int,
    ny: int,
    threshold: float,
    min_valid: int,
) -> pd.DataFrame:
    nrows, ncols = occupied.shape
    row_edges = split_edges(nrows, ny)
    col_edges = split_edges(ncols, nx)
    records: List[Dict[str, object]] = []
    for row in range(ny):
        for col in range(nx):
            r0, r1 = int(row_edges[row]), int(row_edges[row + 1])
            c0, c1 = int(col_edges[col]), int(col_edges[col + 1])
            block_valid = valid[r0:r1, c0:c1]
            block_occupied = occupied[r0:r1, c0:c1]
            n_valid = int(block_valid.sum())
            n_occupied = int(block_occupied.sum())
            fraction = (n_occupied / n_valid) if n_valid else float("nan")
            is_occupied = int(
                n_valid >= min_valid
                and np.isfinite(fraction)
                and fraction >= threshold
            )
            rr = (r0 + r1 - 1) / 2.0
            cc = (c0 + c1 - 1) / 2.0
            x_center, y_center = rasterio.transform.xy(transform, rr, cc)
            records.append(
                {
                    "row": row,
                    "col": col,
                    "deme_index": row * nx + col,
                    "deme_id": "deme_{:02d}_{:02d}".format(row, col),
                    "pixel_r0": r0,
                    "pixel_r1": r1,
                    "pixel_c0": c0,
                    "pixel_c1": c1,
                    "n_valid": n_valid,
                    "n_occupied": n_occupied,
                    "frac_current": fraction,
                    "occupied_current": is_occupied,
                    "x_center": x_center,
                    "y_center": y_center,
                }
            )
    return pd.DataFrame.from_records(records)


def normalize_direction(direction: str) -> str:
    value = direction.strip().lower().replace("-", "_")
    aliases = {
        "eastwest": "east_to_west",
        "e2w": "east_to_west",
        "westeast": "west_to_east",
        "w2e": "west_to_east",
        "northsouth": "north_to_south",
        "n2s": "north_to_south",
        "southnorth": "south_to_north",
        "s2n": "south_to_north",
    }
    value = aliases.get(value, value)
    allowed = {"east_to_west", "west_to_east", "north_to_south", "south_to_north"}
    if value not in allowed:
        raise ValueError("direction must be one of {}".format(", ".join(sorted(allowed))))
    return value


def order_demes(deme_table: pd.DataFrame, direction: str) -> pd.DataFrame:
    direction = normalize_direction(direction)
    occupied = deme_table.loc[deme_table["occupied_current"] == 1].copy()
    if direction == "east_to_west":
        keys, ascending = ["col", "row", "deme_index"], [False, True, True]
    elif direction == "west_to_east":
        keys, ascending = ["col", "row", "deme_index"], [True, True, True]
    elif direction == "north_to_south":
        keys, ascending = ["row", "col", "deme_index"], [True, True, True]
    else:
        keys, ascending = ["row", "col", "deme_index"], [False, True, True]
    ordered = occupied.sort_values(keys, ascending=ascending, kind="mergesort").reset_index(drop=True)
    ordered["loss_rank"] = np.arange(1, len(ordered) + 1, dtype=int)
    return ordered


def active_mask(active: Iterable[int], nx: int, ny: int) -> np.ndarray:
    mask = np.ones((ny, nx), dtype=bool)
    for index in active:
        index = int(index)
        row, col = divmod(index, nx)
        if row < 0 or row >= ny or col < 0 or col >= nx:
            raise ValueError("Deme index {} is outside {}x{} grid".format(index, nx, ny))
        mask[row, col] = False
    return mask


def count_components(active: Iterable[int], nx: int, ny: int) -> int:
    remaining = set(int(x) for x in active)
    components = 0
    while remaining:
        components += 1
        start = remaining.pop()
        queue = deque([start])
        while queue:
            index = queue.popleft()
            row, col = divmod(index, nx)
            neighbours = []
            if row:
                neighbours.append((row - 1) * nx + col)
            if row + 1 < ny:
                neighbours.append((row + 1) * nx + col)
            if col:
                neighbours.append(row * nx + col - 1)
            if col + 1 < nx:
                neighbours.append(row * nx + col + 1)
            for neighbour in neighbours:
                if neighbour in remaining:
                    remaining.remove(neighbour)
                    queue.append(neighbour)
    return components


def extant_indices(extinct: np.ndarray, nx: int, ny: int) -> List[int]:
    index_map, _ = wfm.build_2d_index(nx, ny, extinct_demes=extinct)
    return list(range(len(index_map)))


def build_model(theta: float, migration: float, nx: int, ny: int, extinct: np.ndarray):
    return wfm.build_2d_spatial(
        theta=theta,
        migration_rate=migration,
        xlen=nx,
        ylen=ny,
        extinct_demes=extinct,
    )


def split_halves(extinct: np.ndarray, nx: int, ny: int) -> Tuple[List[int], List[int]]:
    index_map, _ = wfm.build_2d_index(nx, ny, extinct_demes=extinct)
    coords = list(index_map)
    if len(coords) < 2:
        raise ValueError("At least two extant demes are required for FST calibration")
    columns = np.asarray([int(coord[1]) for coord in coords])
    cutoff = float(np.median(columns))
    west = [i for i, coord in enumerate(coords) if int(coord[1]) <= cutoff]
    east = [i for i, coord in enumerate(coords) if int(coord[1]) > cutoff]
    if not west or not east:
        half = max(1, len(coords) // 2)
        west, east = list(range(half)), list(range(half, len(coords)))
    return west, east


def evaluate_migration(
    target_fst: float,
    nx: int,
    ny: int,
    extinct: np.ndarray,
    migration: float,
    metric: str,
    theta_probe: float,
) -> Dict[str, float]:
    try:
        moment_mat, const_vec = build_model(theta_probe, migration, nx, ny, extinct)
        equilibrium = wfm.compute_equilibrium(moment_mat, const_vec)
        west, east = split_halves(extinct, nx, ny)
        if metric == "nei":
            fst = scalarize(wfm.compute_fst_nei(equilibrium, west, east))
        else:
            fst = scalarize(wfm.compute_fst_hudson(equilibrium, west, east))
        difference = abs(fst - target_fst) if np.isfinite(fst) else float("inf")
        return {
            "migration_rate": float(migration),
            "fst_model": float(fst),
            "fst_abs_diff": float(difference),
            "status": "ok" if np.isfinite(difference) else "nonfinite",
        }
    except Exception as exc:
        return {
            "migration_rate": float(migration),
            "fst_model": float("nan"),
            "fst_abs_diff": float("inf"),
            "status": "failed: {}".format(str(exc)[:180]),
        }


def calibrate_migration(
    target_fst: float,
    nx: int,
    ny: int,
    extinct: np.ndarray,
    metric: str,
    theta_probe: float,
    grid: Sequence[float],
) -> Tuple[float, pd.DataFrame]:
    rows = [evaluate_migration(target_fst, nx, ny, extinct, value, metric, theta_probe) for value in grid]
    scan = pd.DataFrame(rows)
    usable = scan[np.isfinite(scan["fst_abs_diff"].to_numpy())]
    if usable.empty:
        raise RuntimeError("All migration calibration candidates failed")
    best = usable.sort_values(["fst_abs_diff", "migration_rate"], kind="mergesort").iloc[0]
    return float(best["migration_rate"]), scan.sort_values(["fst_abs_diff", "migration_rate"], kind="mergesort")


def calibrate_theta(
    pi_obs: float,
    nx: int,
    ny: int,
    extinct: np.ndarray,
    migration: float,
    theta_probe: float,
) -> Tuple[float, float, object]:
    probe_m, probe_v = build_model(theta_probe, migration, nx, ny, extinct)
    probe_eq = wfm.compute_equilibrium(probe_m, probe_v)
    demes = extant_indices(extinct, nx, ny)
    pi_probe = scalarize(wfm.compute_pi(probe_eq, demes))
    if not np.isfinite(pi_probe) or pi_probe <= 0:
        raise RuntimeError("The theta calibration model returned invalid pi ({})".format(pi_probe))
    theta = theta_probe * pi_obs / pi_probe
    model_m, model_v = build_model(theta, migration, nx, ny, extinct)
    equilibrium = wfm.compute_equilibrium(model_m, model_v)
    pi_current = scalarize(wfm.compute_pi(equilibrium, demes))
    if not np.isfinite(pi_current) or pi_current <= 0:
        raise RuntimeError("The calibrated model returned invalid current pi ({})".format(pi_current))
    return float(theta), float(pi_current), equilibrium


def safe_pi(moment: object, demes: Sequence[int]) -> float:
    if moment is None:
        return float("nan")
    try:
        result = scalarize(wfm.compute_pi(moment, demes))
        return float(result) if np.isfinite(result) else float("nan")
    except Exception:
        return float("nan")


def safe_evolve(moment_mat, const_vec, moments, time: float):
    try:
        return wfm.evolve_forward(moment_mat, const_vec, moments, time)
    except Exception:
        return None


def parse_float_list(value: str) -> List[float]:
    result = []
    for item in str(value).split(","):
        item = item.strip()
        if item:
            result.append(finite_float(item, "migration_grid"))
    if not result:
        raise ValueError("migration_grid must contain at least one value")
    return result


def parse_auto_float(value: str, label: str) -> Optional[float]:
    text = str(value).strip().lower()
    if text in {"", "auto", "none", "na", "nan"}:
        return None
    return finite_float(text, label)


def scenario_name(value: object) -> str:
    return str(value).strip().replace("-", "_")


def scenario_time(scenario: str, time3: float, time5: float) -> float:
    match = SCENARIO_RE.search(scenario.replace("_", "-"))
    if match and match.group(1).startswith("2061"):
        return time3
    if match and match.group(1).startswith("2081"):
        return time5
    return time3


def read_area_table(path: Path) -> pd.DataFrame:
    frame = read_table_auto(path)
    if "Scenario" not in frame.columns or "Area_km2" not in frame.columns:
        raise ValueError("area_file must contain Scenario and Area_km2 columns")
    frame = frame[["Scenario", "Area_km2"]].copy()
    frame["Scenario"] = frame["Scenario"].map(scenario_name)
    frame["Area_km2"] = pd.to_numeric(frame["Area_km2"], errors="coerce")
    frame = frame[np.isfinite(frame["Area_km2"]) & (frame["Area_km2"] > 0)]
    current = frame.loc[frame["Scenario"].str.lower() == "current", "Area_km2"]
    if current.empty:
        raise ValueError("area_file has no Scenario=current row")
    current_area = float(current.iloc[0])
    future = frame[frame["Scenario"].str.lower() != "current"].copy()
    future["A_remaining"] = future["Area_km2"] / current_area
    future["A_remaining"] = future["A_remaining"].clip(lower=0, upper=1)
    future["habitat_loss_pct"] = (1.0 - future["A_remaining"]) * 100.0
    return future.reset_index(drop=True)


def area_from_future_masks(current_path: Path, mapping_path: Path) -> pd.DataFrame:
    """Derive area ratios from a JSON object mapping scenario to binary TIFF."""
    with mapping_path.open("r", encoding="utf-8") as handle:
        mapping = json.load(handle)
    _, current_valid, _, _ = read_binary_raster(current_path)
    current_count = int(current_valid.sum())
    # Area is based on suitable cells, not all cells in the rectangular raster.
    current_occ, _, _, _ = read_binary_raster(current_path)
    current_count = int(current_occ.sum())
    if current_count <= 0:
        raise ValueError("Current raster contains no occupied cells")
    rows = []
    for raw_name, raw_path in mapping.items():
        future_occ, _, _, _ = read_binary_raster(Path(raw_path))
        if future_occ.shape != current_occ.shape:
            raise ValueError("Future raster shape differs from current raster: {}".format(raw_path))
        remaining = float(future_occ.sum()) / current_count
        rows.append(
            {
                "Scenario": scenario_name(raw_name),
                "A_remaining": min(1.0, max(0.0, remaining)),
                "habitat_loss_pct": (1.0 - remaining) * 100.0,
            }
        )
    return pd.DataFrame(rows)


def normalize_curve(frame: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, float]]:
    frame = frame.sort_values("habitat_loss_pct").reset_index(drop=True).copy()
    zero = frame[np.isclose(frame["habitat_loss_pct"].to_numpy(dtype=float), 0.0)]
    if zero.empty:
        raise ValueError("Curve has no zero-loss baseline")
    reference = float(zero.iloc[0]["pi_remaining_immediate"])
    if not np.isfinite(reference) or reference <= 0:
        raise ValueError("Invalid zero-loss immediate pi: {}".format(reference))
    columns = [
        "pi_remaining_immediate",
        "pi_remaining_3gen",
        "pi_remaining_5gen",
        "pi_remaining_equilibrium",
    ]
    baseline: Dict[str, float] = {"baseline_raw_immediate": reference}
    for column in columns:
        if column not in frame.columns:
            continue
        raw_name = column + "_raw_modelscale"
        frame[raw_name] = frame[column]
        baseline["baseline_raw_{}".format(column)] = float(zero.iloc[0][column])
        frame[column] = pd.to_numeric(frame[column], errors="coerce") / reference * 100.0
    return frame, baseline


def summarize_replicates(replicates: pd.DataFrame) -> pd.DataFrame:
    value_columns = [
        "gdar_remaining",
        "pi_remaining_immediate",
        "pi_remaining_3gen",
        "pi_remaining_5gen",
        "pi_remaining_equilibrium",
    ]
    value_columns = [column for column in value_columns if column in replicates.columns]
    rows = []
    group_columns = ["habitat_loss_pct", "A_remaining", "occupied_demes"]
    for keys, group in replicates.groupby(group_columns, sort=True, dropna=False):
        for column in value_columns:
            values = pd.to_numeric(group[column], errors="coerce")
            finite = values[np.isfinite(values)]
            rows.append(
                {
                    "curve_type": {
                        "gdar_remaining": "GDAR_short_term",
                        "pi_remaining_immediate": "WF_Immediate",
                        "pi_remaining_3gen": "WF_3_gen",
                        "pi_remaining_5gen": "WF_5_gen",
                        "pi_remaining_equilibrium": "WF_equilibrium",
                    }[column],
                    "habitat_loss_pct": float(keys[0]),
                    "A_remaining": float(keys[1]),
                    "occupied_demes": int(keys[2]),
                    "mean_pi_remaining": float(finite.mean()) if len(finite) else float("nan"),
                    # WF columns are already normalised to percent; GDAR is
                    # intentionally retained as a fraction and carries its
                    # own curve type so plotting cannot confuse it with WF
                    # immediate diversity.
                    "value_scale": "fraction" if column == "gdar_remaining" else "percent",
                    "sd_pi_remaining": float(finite.std(ddof=1)) if len(finite) > 1 else 0.0 if len(finite) else float("nan"),
                    "n_nonmissing": int(len(finite)),
                    "n_total": int(len(values)),
                    "loss_mode": str(group["loss_mode"].iloc[0]) if "loss_mode" in group.columns else "",
                }
            )
    return pd.DataFrame(rows).sort_values(["curve_type", "habitat_loss_pct"]).reset_index(drop=True)


def interpolate_value(curve: pd.DataFrame, column: str, loss: float) -> float:
    values = curve[["habitat_loss_pct", column]].copy()
    values["habitat_loss_pct"] = pd.to_numeric(values["habitat_loss_pct"], errors="coerce")
    values[column] = pd.to_numeric(values[column], errors="coerce")
    values = values[np.isfinite(values["habitat_loss_pct"]) & np.isfinite(values[column])]
    if values.empty:
        return float("nan")
    values = values.sort_values("habitat_loss_pct").drop_duplicates("habitat_loss_pct")
    return float(np.interp(float(loss), values["habitat_loss_pct"], values[column]))


def make_scenario_points(curve: pd.DataFrame, area: pd.DataFrame, time3: float, time5: float) -> pd.DataFrame:
    rows = []
    for _, item in area.iterrows():
        scenario = scenario_name(item["Scenario"])
        loss = float(item["habitat_loss_pct"])
        rows.append(
            {
                "scenario": scenario,
                "scenario_raw": str(item["Scenario"]),
                "A_remaining": float(item["A_remaining"]),
                "habitat_loss_pct": loss,
                "time_forward": scenario_time(scenario, time3, time5),
                "gdar_remaining": interpolate_value(curve, "gdar_remaining", loss)
                if "gdar_remaining" in curve.columns
                else float("nan"),
                "GDAR_short_term": interpolate_value(curve, "gdar_remaining", loss)
                if "gdar_remaining" in curve.columns
                else float("nan"),
                "pi_remaining_immediate": interpolate_value(curve, "pi_remaining_immediate", loss),
                "pi_remaining_3gen": interpolate_value(curve, "pi_remaining_3gen", loss),
                "pi_remaining_5gen": interpolate_value(curve, "pi_remaining_5gen", loss),
                "pi_remaining_equilibrium": interpolate_value(curve, "pi_remaining_equilibrium", loss)
                if "pi_remaining_equilibrium" in curve.columns
                else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compute_state(
    remaining: Sequence[int],
    current_extinct: np.ndarray,
    eq_current: object,
    theta: float,
    migration: float,
    nx: int,
    ny: int,
    pi_obs: float,
    z_gdar: float,
    a_remaining: float,
    time3: float,
    time5: float,
    include_equilibrium: bool,
) -> Dict[str, object]:
    remaining = sorted(int(x) for x in remaining)
    if not remaining:
        return {
            "occupied_demes": 0,
            "n_components": 0,
            "pi_remaining_immediate": 0.0,
            "pi_remaining_3gen": 0.0,
            "pi_remaining_5gen": 0.0,
            "pi_remaining_equilibrium": 0.0 if include_equilibrium else float("nan"),
            "gdar_remaining": 0.0 if a_remaining <= 0 else float(a_remaining) ** z_gdar,
            "step_status": "extinct_zero_deme",
        }
    extinct_now = active_mask(remaining, nx, ny)
    demes_now = extant_indices(extinct_now, nx, ny)
    try:
        moments_now = wfm.get_moments_2d(eq_current, current_extinct, extinct_now)
        model_m, model_v = build_model(theta, migration, nx, ny, extinct_now)
        pi0 = safe_pi(moments_now, demes_now)
        mom3 = safe_evolve(model_m, model_v, moments_now, time3)
        mom5 = safe_evolve(model_m, model_v, moments_now, time5)
        pi3 = safe_pi(mom3, demes_now)
        pi5 = safe_pi(mom5, demes_now)
        pieq = float("nan")
        status = "ok"
        if include_equilibrium:
            try:
                equilibrium = wfm.compute_equilibrium(model_m, model_v)
                pieq = safe_pi(equilibrium, demes_now)
            except Exception:
                status = "equilibrium_failed"
        if not np.isfinite(pi0):
            status = "failed_immediate"
        elif not np.isfinite(pi3) or not np.isfinite(pi5):
            status = "failed_forward"
        return {
            "occupied_demes": len(remaining),
            "n_components": count_components(remaining, nx, ny),
            "pi_remaining_immediate": pi0 / pi_obs if np.isfinite(pi0) else float("nan"),
            "pi_remaining_3gen": pi3 / pi_obs if np.isfinite(pi3) else float("nan"),
            "pi_remaining_5gen": pi5 / pi_obs if np.isfinite(pi5) else float("nan"),
            "pi_remaining_equilibrium": pieq / pi_obs if np.isfinite(pieq) else float("nan"),
            "gdar_remaining": 0.0 if a_remaining <= 0 else float(a_remaining) ** z_gdar,
            "step_status": status,
        }
    except Exception as exc:
        return {
            "occupied_demes": len(remaining),
            "n_components": count_components(remaining, nx, ny),
            "pi_remaining_immediate": float("nan"),
            "pi_remaining_3gen": float("nan"),
            "pi_remaining_5gen": float("nan"),
            "pi_remaining_equilibrium": float("nan"),
            "gdar_remaining": 0.0 if a_remaining <= 0 else float(a_remaining) ** z_gdar,
            "step_status": "failed_exception: {}".format(str(exc)[:180]),
        }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="2-D deme WFmoments projection")
    parser.add_argument("--current-raster", required=True, type=Path)
    parser.add_argument("--pi-file", required=True, type=Path)
    parser.add_argument("--structure-file", type=Path)
    parser.add_argument("--param-file", type=Path)
    parser.add_argument("--z-gdar", type=float)
    parser.add_argument("--area-file", type=Path)
    parser.add_argument("--future-masks-json", type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--nx", type=int, default=20)
    parser.add_argument("--ny", type=int, default=20)
    parser.add_argument("--threshold", type=float, default=0.25)
    parser.add_argument("--min-valid", type=int, default=0)
    parser.add_argument("--loss-mode", choices=["edge", "random"], default="edge")
    parser.add_argument("--direction", default="east_to_west")
    parser.add_argument("--migration", default="25")
    parser.add_argument("--migration-grid", default="0.1,0.3,1,3,10,25,50")
    parser.add_argument("--fst-metric", choices=["hudson", "nei"], default="hudson")
    parser.add_argument("--theta", default="auto")
    parser.add_argument("--theta-probe", type=float, default=1e-4)
    parser.add_argument("--time3", type=float, default=3.0)
    parser.add_argument("--time5", type=float, default=5.0)
    parser.add_argument("--mu", type=float, default=3.75e-8)
    parser.add_argument("--midterm-generations", default="auto")
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--include-equilibrium", action="store_true")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args(argv)

    for path in [args.current_raster, args.pi_file]:
        if not path.is_file():
            raise ValueError("Input does not exist: {}".format(path))
    if args.structure_file is not None and not args.structure_file.is_file():
        raise ValueError("structure_file does not exist: {}".format(args.structure_file))
    if args.param_file is not None and not args.param_file.is_file():
        raise ValueError("param_file does not exist: {}".format(args.param_file))
    if args.area_file is not None and not args.area_file.is_file():
        raise ValueError("area_file does not exist: {}".format(args.area_file))
    if args.future_masks_json is not None and not args.future_masks_json.is_file():
        raise ValueError("future_masks_json does not exist: {}".format(args.future_masks_json))
    if args.nx < 1 or args.ny < 1 or args.replicates < 1:
        raise ValueError("nx, ny and replicates must be positive")
    if not 0 <= args.threshold <= 1:
        raise ValueError("threshold must be between 0 and 1")
    if args.min_valid < 0:
        raise ValueError("min_valid must not be negative")
    direction = normalize_direction(args.direction)
    pi_obs = read_species_value(args.pi_file, args.species, "pi_obs")
    z_gdar = (
        finite_float(args.z_gdar, "z_gdar")
        if args.z_gdar is not None
        else read_species_value(args.param_file, args.species, "z_gdar")
        if args.param_file is not None
        else None
    )
    if z_gdar is None:
        raise ValueError("Provide --param-file or --z-gdar")
    if args.mu <= 0 or args.theta_probe <= 0:
        raise ValueError("mu and theta_probe must be positive")

    args.outdir.mkdir(parents=True, exist_ok=True)
    occupied, valid, transform, raster_info = read_binary_raster(args.current_raster)
    deme_table = build_deme_table(occupied, valid, transform, args.nx, args.ny, args.threshold, args.min_valid)
    current = sorted(int(x) for x in deme_table.loc[deme_table["occupied_current"] == 1, "deme_index"])
    if not current:
        raise ValueError("No occupied demes under threshold/min_valid")
    current_set = set(current)
    current_extinct = active_mask(current, args.nx, args.ny)

    migration_value = parse_auto_float(args.migration, "migration")
    migration_scan = None
    target_fst = None
    if args.structure_file is not None:
        target_fst = read_species_value(args.structure_file, args.species, "fst_global_est")
    if migration_value is None:
        if target_fst is None:
            raise ValueError("migration=auto requires --structure-file")
        migration_value, migration_scan = calibrate_migration(
            target_fst,
            args.nx,
            args.ny,
            current_extinct,
            args.fst_metric,
            args.theta_probe,
            parse_float_list(args.migration_grid),
        )
    elif migration_value <= 0:
        raise ValueError("migration must be positive")

    theta_value = parse_auto_float(args.theta, "theta")
    if theta_value is None:
        theta_value, pi_current, eq_current = calibrate_theta(
            pi_obs, args.nx, args.ny, current_extinct, migration_value, args.theta_probe
        )
    else:
        if theta_value <= 0:
            raise ValueError("theta must be positive")
        model_m, model_v = build_model(theta_value, migration_value, args.nx, args.ny, current_extinct)
        eq_current = wfm.compute_equilibrium(model_m, model_v)
        pi_current = safe_pi(eq_current, extant_indices(current_extinct, args.nx, args.ny))
        if not np.isfinite(pi_current) or pi_current <= 0:
            raise ValueError("Explicit theta produced invalid current pi")

    midterm = parse_auto_float(args.midterm_generations, "midterm_generations")
    if str(args.midterm_generations).strip().lower() == "auto":
        midterm = pi_obs / (4.0 * args.mu) / 2.0

    area = None
    if args.area_file is not None:
        area = read_area_table(args.area_file)
    elif args.future_masks_json is not None:
        area = area_from_future_masks(args.current_raster, args.future_masks_json)

    replicate_rows: List[Dict[str, object]] = []
    for replicate in range(1, args.replicates + 1):
        if args.loss_mode == "random":
            rng = np.random.default_rng(args.seed + replicate - 1)
            order = rng.permutation(current).tolist()
        else:
            order = order_demes(deme_table, direction)["deme_index"].astype(int).tolist()
        for removed_count in range(0, len(current) + 1):
            remaining = current_set - set(order[:removed_count])
            loss = 100.0 * removed_count / len(current)
            state = compute_state(
                sorted(remaining),
                current_extinct,
                eq_current,
                theta_value,
                migration_value,
                args.nx,
                args.ny,
                pi_obs,
                z_gdar,
                1.0 - removed_count / len(current),
                args.time3,
                args.time5,
                args.include_equilibrium,
            )
            row = {
                "replicate": replicate,
                "seed": args.seed + replicate - 1,
                "habitat_loss_pct": loss,
                "A_remaining": 1.0 - removed_count / len(current),
                "direction": direction,
                "loss_mode": args.loss_mode,
                "z_gdar": z_gdar,
            }
            row.update(state)
            replicate_rows.append(row)

    replicates = pd.DataFrame(replicate_rows)
    # The model is calibrated to pi_obs, but the immediate zero-loss state is
    # the reproducible normalization reference used in the paper workflow.
    normalized_parts = []
    for replicate, group in replicates.groupby("replicate", sort=True):
        norm, baseline = normalize_curve(group)
        norm["baseline_raw_immediate"] = baseline["baseline_raw_immediate"]
        normalized_parts.append(norm)
    replicates = pd.concat(normalized_parts, ignore_index=True)
    summary = summarize_replicates(replicates)

    if area is not None:
        # Convert the long summary to a wide interpolation frame.
        curve_wide = summary.pivot_table(
            index=["habitat_loss_pct", "A_remaining", "occupied_demes"],
            columns="curve_type",
            values="mean_pi_remaining",
            aggfunc="first",
        ).reset_index()
        curve_wide = curve_wide.rename(
            columns={
                "GDAR_short_term": "gdar_remaining",
                "WF_Immediate": "pi_remaining_immediate",
                "WF_3_gen": "pi_remaining_3gen",
                "WF_5_gen": "pi_remaining_5gen",
                "WF_equilibrium": "pi_remaining_equilibrium",
            }
        )
        points = make_scenario_points(curve_wide, area, args.time3, args.time5)
    else:
        points = pd.DataFrame()

    deme_table.to_csv(args.outdir / "deme_table.tsv", sep="\t", index=False)
    order_table = order_demes(deme_table, direction)
    order_table.to_csv(args.outdir / "deme_order.tsv", sep="\t", index=False)
    replicates.to_csv(args.outdir / "curve_replicates.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "curve_summary.tsv", sep="\t", index=False)
    if not points.empty:
        points.to_csv(args.outdir / "scenario_points.tsv", sep="\t", index=False)
    if migration_scan is not None:
        migration_scan.to_csv(args.outdir / "migration_scan.tsv", sep="\t", index=False)

    meta = {
        "module": "wfmoment-2d-deme",
        "species": args.species,
        "pi_obs": pi_obs,
        "z_gdar": z_gdar,
        "gdar_definition": "A_remaining^z_gdar",
        "mu": args.mu,
        "Ne_est_from_pi": pi_obs / (4.0 * args.mu),
        "midterm_generations": midterm,
        "nx": args.nx,
        "ny": args.ny,
        "threshold": args.threshold,
        "min_valid": args.min_valid,
        "occupied_demes_current": len(current),
        "loss_mode": args.loss_mode,
        "direction": direction,
        "replicates": args.replicates,
        "seed": args.seed,
        "migration_rate": migration_value,
        "migration_input": args.migration,
        "target_fst": target_fst,
        "fst_metric": args.fst_metric,
        "theta": theta_value,
        "pi_current_model": pi_current,
        "theta_calibration_relative_error": (pi_current - pi_obs) / pi_obs,
        "time3": args.time3,
        "time5": args.time5,
        "include_equilibrium": bool(args.include_equilibrium),
        "raster": raster_info,
        "inputs": {
            "current_raster": str(args.current_raster),
            "pi_file": str(args.pi_file),
            "structure_file": str(args.structure_file) if args.structure_file else None,
            "param_file": str(args.param_file) if args.param_file else None,
            "area_file": str(args.area_file) if args.area_file else None,
            "future_masks_json": str(args.future_masks_json) if args.future_masks_json else None,
        },
        "input_sha256": {
            str(path): sha256(path)
            for path in [args.current_raster, args.pi_file, args.structure_file, args.param_file, args.area_file]
            if path is not None and path.is_file()
        },
    }
    (args.outdir / "metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    (args.outdir / "metadata.tsv").write_text(
        "\t".join(meta.keys()) + "\n" + "\t".join(str(meta[key]) for key in meta.keys()) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"outdir": str(args.outdir), "meta": meta, "scenarios": len(points)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print("wfmoment_compute error: {}".format(exc), file=sys.stderr)
        traceback.print_exc()
        raise
