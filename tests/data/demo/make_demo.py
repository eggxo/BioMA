#!/usr/bin/env python3
"""Generate a small, server-independent BioMA project fixture.

The generated files are intentionally deterministic.  Raster placeholders are
used unless ``--valid-rasters`` is supplied; this keeps the basic demo usable
without importing rasterio while preserving a path-complete project layout.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


BIO_NAMES = ["bio{}".format(index) for index in range(1, 20)]
MODELS = ("demoA", "demoB")
SCENARIO = "2061-2080-ssp245-{}"


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _vcf(path: Path, samples) -> None:
    rows = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    genotypes = [
        ["0/0", "0/1", "1/1", "0/1", "0/0", "1/1"],
        ["0/1", "0/1", "0/0", "1/1", "0/1", "0/0"],
        ["1/1", "0/0", "0/1", "0/0", "1/1", "0/1"],
        ["0/0", "1/1", "0/1", "1/1", "0/0", "0/1"],
    ]
    for index, values in enumerate(genotypes, start=1):
        rows.append("1\t{}\tadaptive_{}\tA\tG\t.\tPASS\t.\tGT\t{}".format(index * 100, index, "\t".join(values)))
    _write(path, "\n".join(rows) + "\n")


def _raster(path: Path, value: float, valid: bool) -> None:
    if valid:
        try:
            import numpy as np
            import rasterio
            from rasterio.transform import from_origin
        except ImportError as error:
            raise SystemExit("--valid-rasters requires rasterio and numpy: {}".format(error))
        path.parent.mkdir(parents=True, exist_ok=True)
        data = np.full((5, 5), value, dtype="float32")
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            height=5,
            width=5,
            count=1,
            dtype="float32",
            crs="EPSG:4326",
            transform=from_origin(99.5, 25.5, 1.0, 1.0),
            nodata=-9999,
        ) as dst:
            dst.write(data, 1)
    else:
        _write(path, "BioMA demo raster placeholder {}\n".format(value))


def _module_configs(root: Path) -> None:
    configs = root / "configs"
    data = root / "data"
    common = {
        "adaptive": data / "adaptive_sites.vcf",
        "samples": data / "population_samples.tsv",
        "environment": data / "population_environment.tsv",
        "climate_current": data / "climate" / "current",
        "climate_future": data / "climate" / "future",
        "mask": data / "species_mask.geojson",
    }
    _write(
        configs / "gf.ini",
        """[inputs]\nvcf = {adaptive}\nsamples = {samples}\ncoordinates = {environment}\npresent_climate = {current}\nfuture_climate = {future}\ncurrent_mask = {mask}\n\n[analysis]\noutput_dir = {out}\nmodels = demoA,demoB\nssps = 245\nperiods = 2061-2080\n\n[parameters]\n""".format(adaptive=common["adaptive"], samples=common["samples"], environment=common["environment"], current=common["climate_current"], future=common["climate_future"], mask=common["mask"], out=root / "results" / "01_gf"),
    )
    _write(
        configs / "rona.ini",
        """[inputs]\nalt_frequency = {frequency}\nunld_dir = {unld}\nenvironment = {environment}\nfuture_climate = {future}\nmask = {mask}\n\n[analysis]\noutput_dir = {out}\nmodels = demoA,demoB\nssps = 245\nperiods = 2061-2080\n\n[parameters]\nrscript = Rscript\n""".format(frequency=data / "adaptive_frequency.tsv", unld=data / "rona_ld", environment=common["environment"], future=common["climate_future"], mask=common["mask"], out=root / "results" / "02_rona"),
    )
    _write(
        configs / "mar.ini",
        """[inputs]\nvcf = {vcf}\nlonlat = {lonlat}\n\n[analysis]\noutput_dir = {out}\nname = demo_mar\n\n[parameters]\nmaxsnps = auto\nscheme = random\n""".format(vcf=data / "whole_genome.mis0.9.maf0.00001.vcf", lonlat=data / "mar_lonlat.tsv", out=root / "results" / "03_mar"),
    )
    _write(
        configs / "load.ini",
        """[inputs]\nvcf_dir = {vcfs}\npopulation_dir = {popdir}\npredictors = {predictors}\nfuture_dir = {future}\nmask = {mask}\n\n[analysis]\noutput_dir = {out}\n\n[parameters]\nexpected_future_files = 1\n""".format(vcfs=data / "load_3vcf", popdir=data / "population_lists", predictors=data / "load_predictors.csv", future=data / "load_future", mask=common["mask"], out=root / "results" / "04_load"),
    )
    _write(
        configs / "niche.ini",
        """[inputs]\nwork_dir = {root}\noccurrence_csv = {occurrence}\ncurrent_env_dir = {current}\nfuture_root = {future}\nmask_shp = {mask}\n\n[variables]\ncorrelation_threshold = 0.8\ncandidate_subset_sizes = 1\n\n[tuning]\nmax_models = 1\nbackground_n = 10\ncv_folds = 2\nfeature_classes = L\nbeta_multipliers = 1\n\n[projection]\nperiods = 2061-2080\nscenarios = ssp245\ngcms = demoA,demoB\nensemble_method = mean\nbinary_output = false\n\n[analysis]\noutput_dir = {out}\nrscript = Rscript\nplot_rscript = Rscript\n""".format(root=root, occurrence=data / "occurrence.csv", current=common["climate_current"], future=common["climate_future"], mask=common["mask"], out=root / "results" / "05_niche"),
    )
    _write(
        configs / "wfmoment.ini",
        """[inputs]\ncurrent_raster = {raster}\npi_file = {pi}\narea_file = {area}\n\n[analysis]\nspecies = demo_species\noutput_dir = {out}\n\n[parameters]\ncompute_python = python3\nrscript = Rscript\nnx = 2\nny = 2\nreplicates = 1\n""".format(raster=data / "habitat_current.tif", pi=data / "pi_observed.csv", area=data / "area_summary.csv", out=root / "results" / "06_wfmoment"),
    )
    _write(
        configs / "vulnerability.ini",
        """[inputs]\ngf_offsets = {gf}\nrona_ensemble = {rona}\nload_predictors = {load}\nniche_raster = {niche}\n\n[analysis]\noutput_dir = {out}\nscenario_label = BioMA demo\n\n[parameters]\ngf_field = local_offset_mean\nrona_variables = auto\nrona_summary = mean\nloadM_field = mean_loadM_relax\nloadD_field = mean_loadD_relax\ngroup_column = cluster\ngroup_order =\nexclude_groups =\n""".format(gf=data / "gf_offsets.tsv", rona=data / "rona_ensemble.tsv", load=data / "load_predictors.csv", niche=data / "niche_maladaptation.tif", out=root / "results" / "07_vulnerability"),
    )
    modules = "\n".join("{} = {}".format(name, configs / "{}.ini".format(name)) for name in ("gf", "rona", "mar", "load", "niche", "wfmoment", "vulnerability"))
    _write(
        root / "project.ini",
        """[project]\nname = bioma-demo\noutput_dir = {out}\nresume = true\nstop_on_error = true\nreport = true\n\n[modules]\n{modules}\n\n[shared]\nspecies = demo_species\nmodels = demoA,demoB\nssps = 245\nperiods = 2061-2080\n\n[integration]\nperiod = 2061-2080\nssp = 245\nensemble_method = mean\nrona_variables = auto\nrona_summary = mean\ngroup_order =\nexclude_groups =\n""".format(out=root / "results", modules=modules),
    )


def generate(root: Path, valid_rasters: bool = False) -> Path:
    if root.resolve() == Path(__file__).resolve().parent:
        raise ValueError("Refusing to overwrite the checked-in demo directory; choose --outdir elsewhere")
    if root.exists():
        marker = root / "demo_manifest.json"
        if any(root.iterdir()) and not marker.is_file():
            raise ValueError("Refusing to remove an existing directory without demo_manifest.json: {}".format(root))
        if any(root.iterdir()):
            shutil.rmtree(root)
    data = root / "data"
    samples = ["S1", "S2", "S3", "S4", "S5", "S6"]
    _vcf(data / "adaptive_sites.vcf", samples)
    _vcf(data / "whole_genome.mis0.9.maf0.00001.vcf", samples)
    _write(data / "population_samples.tsv", "sample_id\tpopulation_id\n" + "\n".join("{}\tP{}".format(sample, (index // 2) + 1) for index, sample in enumerate(samples)) + "\n")
    env_header = ["ID", "pop", "lon", "lat", "cluster"] + BIO_NAMES
    env_rows = []
    for index in range(1, 4):
        env_rows.append(["P{}".format(index), "G{}".format(index), str(100 + index), str(20 + index), "G{}".format(index)] + [str(index + bio / 100.0) for bio in range(1, 20)])
    _write(data / "population_environment.tsv", "\t".join(env_header) + "\n" + "\n".join("\t".join(row) for row in env_rows) + "\n")
    _write(data / "adaptive_frequency.tsv", "Pop\tadaptive_1\tadaptive_2\tadaptive_3\tadaptive_4\n" + "P1\t0.1\t0.2\t0.3\t0.4\nP2\t0.2\t0.3\t0.4\t0.5\nP3\t0.3\t0.4\t0.5\t0.6\n")
    for bio in range(1, 20):
        _write(data / "rona_ld" / "LD_BIO{}.prune.in".format(bio), "adaptive_1\nadaptive_2\n")
    future_directories = [data / "climate" / "future" / SCENARIO.format(model) for model in MODELS]
    for directory in (data / "climate" / "current", *future_directories):
        for bio in range(1, 20):
            name = "bio{}.tif".format(bio) if directory.name == "current" else "{}bio{}.cut.tif".format(directory.name, bio)
            _raster(directory / name, float(bio), valid_rasters)
    _raster(data / "habitat_current.tif", 1.0, valid_rasters)
    _raster(data / "niche_maladaptation.tif", 0.2, valid_rasters)
    _write(data / "species_mask.geojson", json.dumps({"type": "FeatureCollection", "features": []}) + "\n")
    _write(data / "mar_lonlat.tsv", "ID\tLONGITUDE\tLATITUDE\nS1\t101\t21\nS2\t101\t21\nS3\t102\t22\nS4\t102\t22\nS5\t103\t23\nS6\t103\t23\n")
    _write(data / "load_predictors.csv", "pop,lat,lon,cluster," + ",".join(BIO_NAMES) + ",mean_loadM_relax,mean_loadD_relax\n" + "\n".join(",".join(["P{}".format(i), str(20 + i), str(100 + i), "G{}".format(i)] + [str(i + bio / 100.0) for bio in range(1, 20)] + [str(i / 10), str(i / 20)]) for i in range(1, 4)) + "\n")
    for name in ("strict.synonymous.vcf", "strict.nonsynonymous.vcf", "strict.deleterious_plain.vcf", "relaxed.deleterious_with_warning.vcf"):
        _vcf(data / "load_3vcf" / name, samples)
    for index in range(1, 4):
        _write(data / "population_lists" / "P{}".format(index), "S{}\nS{}\n".format(index * 2 - 1, index * 2))
    _write(data / "load_future" / "future_env_ssp245_mean.csv", "pop,bio1\nP1,1\nP2,2\nP3,3\n")
    _write(data / "occurrence.csv", "species,lon,lat\ndemo_species,101,21\ndemo_species,102,22\ndemo_species,103,23\n")
    _write(data / "pi_observed.csv", "species,pi_obs\ndemo_species,0.01\n")
    _write(data / "area_summary.csv", "Scenario,Area_km2\ncurrent,100\n2061-2080_ssp245,80\n")
    _write(data / "gf_offsets.tsv", "lon\tlat\tlocal_offset_mean\n101\t21\t0.1\n102\t22\t0.2\n103\t23\t0.3\n")
    _write(data / "rona_ensemble.tsv", "ID\tBIO1_RONA\tBIO2_RONA\nP1\t0.1\t0.2\nP2\t0.2\t0.3\nP3\t0.3\t0.4\n")
    _module_configs(root)
    _write(root / "demo_manifest.json", json.dumps({"schema": 1, "scenarios": [SCENARIO.format(model) for model in MODELS], "server_independent": True, "valid_rasters": bool(valid_rasters)}, indent=2) + "\n")
    return root


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--valid-rasters", action="store_true")
    args = parser.parse_args()
    root = generate(args.outdir.expanduser().resolve(), args.valid_rasters)
    print(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
