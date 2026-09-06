import configparser
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GUIDES = (ROOT / "docs" / "USER_GUIDE.en.md", ROOT / "docs" / "USER_GUIDE.zh-CN.md")

SECTION_KEYS = {
    3: {
        "name", "output_dir", "resume", "stop_on_error", "report",
        "gf", "rona", "mar", "load", "niche", "wfmoment", "vulnerability",
        "adaptive_vcf", "adaptive_frequency", "whole_genome_vcf", "load_vcf_dir",
        "population_samples", "population_environment", "population_dir",
        "rona_unld_dir", "climate_current", "climate_future", "species_mask",
        "load_predictors", "load_future_dir", "mar_lonlat", "occurrence_csv",
        "maxent_jar", "wfmoment_current_raster", "pi_file", "structure_file",
        "species_parameters", "area_file", "future_masks_json", "species", "models",
        "ssps", "periods", "seed", "rscript", "compute_python", "period", "ssp",
        "ensemble_method", "rona_variables", "rona_summary", "rona_weights",
        "group_order", "exclude_groups", "python", "gdalinfo", "gdallocationinfo",
        "ogrinfo", "java",
    },
    4: {
        "vcf", "samples", "sample_groups_dir", "coordinates", "present_climate",
        "future_climate", "current_mask", "future_mask", "output_dir", "models",
        "ssps", "periods", "expected_sites", "predictors", "min_population_samples",
        "warn_population_samples", "ntree", "nbin", "corr_threshold", "max_level",
        "seed", "forward_radii", "initial_knn_k", "batch_size", "verify_sample",
        "tie_tolerance", "map_radius", "minimum_models", "supplied_bio_tolerance",
        "rscript", "gdalinfo", "gdallocationinfo", "ogrinfo",
    },
    5: {
        "alt_frequency", "unld_dir", "environment", "future_climate", "mask",
        "output_dir", "models", "ssps", "periods", "expected_populations", "rscript",
        "interpolate", "grid_step", "grid_k",
    },
    6: {
        "vcf", "lonlat", "scenario_file", "output_dir", "name", "geom_id",
        "maxsnps", "scheme", "nrep", "xfrac", "quorum", "randseed", "marsteps",
        "rscript",
    },
    7: {
        "vcf_dir", "population_dir", "predictors", "future_dir", "mask", "output_dir",
        "n_splits", "cv_v", "cv_repeats", "grid_size", "workers", "scheme", "seed",
        "rscript", "calc_python", "calc_script", "rf_script", "predict_script",
        "future_pattern", "expected_future_files",
    },
    8: {
        "work_dir", "occurrence_csv", "current_env_dir", "future_root", "mask_shp",
        "maxent_jar", "correlation_threshold", "correlation_method", "correlation_scope",
        "candidate_subset_sizes", "max_models", "background_n", "cv_folds",
        "feature_classes", "beta_multipliers", "selection_metric", "seed", "periods",
        "scenarios", "gcms", "ensemble_method", "binary_output", "output_dir", "rscript",
        "plot_rscript",
    },
    9: {
        "current_raster", "pi_file", "structure_file", "param_file", "area_file",
        "future_masks_json", "species", "output_dir", "plot_title", "compute_python",
        "rscript", "nx", "ny", "threshold", "min_valid", "loss_mode", "direction",
        "migration", "migration_grid", "fst_metric", "theta", "z_gdar", "theta_probe",
        "time3", "time5", "mu", "midterm_generations", "replicates", "seed",
        "include_equilibrium", "plot_equilibrium", "plot_raw",
    },
    10: {
        "gf_offsets", "rona_ensemble", "load_predictors", "niche_raster", "output_dir",
        "scenario_label", "gf_field", "rona_variables", "rona_summary", "rona_weights",
        "rona_bio3_field", "rona_bio15_field", "loadM_field", "loadD_field",
        "group_column", "group_order", "exclude_groups", "nearest_tolerance", "rscript",
    },
}


def section(text: str, number: int) -> str:
    start = text.index("## {}.".format(number))
    next_marker = "## {}.".format(number + 1)
    end = text.find(next_marker, start + 1)
    return text[start:] if end < 0 else text[start:end]


def documents_field(text: str, key: str) -> bool:
    pattern = r"`[^`\n]*\b{}\b[^`\n]*`".format(re.escape(key))
    return re.search(pattern, text, flags=re.IGNORECASE) is not None


class DocumentationContractTest(unittest.TestCase):
    def test_bilingual_guides_cover_every_supported_field(self):
        for guide in GUIDES:
            text = guide.read_text(encoding="utf-8")
            for number, keys in SECTION_KEYS.items():
                body = section(text, number)
                missing = sorted(key for key in keys if not documents_field(body, key))
                self.assertEqual(missing, [], "{} section {}".format(guide.name, number))

    def test_every_example_ini_field_is_documented_in_both_languages(self):
        examples = sorted(ROOT.glob("workflow*.example.ini"))
        for example in examples:
            parser = configparser.ConfigParser(interpolation=None)
            parser.read(example, encoding="utf-8")
            keys = {key for section_name in parser.sections() for key in parser[section_name]}
            for guide in GUIDES:
                text = guide.read_text(encoding="utf-8")
                missing = sorted(key for key in keys if not documents_field(text, key))
                self.assertEqual(missing, [], "{} missing fields from {}".format(guide.name, example.name))

    def test_readme_and_source_manifest_publish_both_guides(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
        for guide in GUIDES:
            relative = guide.relative_to(ROOT).as_posix()
            self.assertIn(relative, readme)
        self.assertIn("recursive-include docs *.md", manifest)

    def test_method_definitions_and_core_references_are_explicit(self):
        required = (
            "local offset", "forward offset", "reverse offset",
            "abs(alpha + beta*E_future - p_obs)", "S = c*A^z",
            "loadM = Pn / (Pn + Ps)", "maladaptation",
            "A_remaining^z_gdar", "10.1073/pnas.2514371123",
        )
        for guide in GUIDES:
            text = guide.read_text(encoding="utf-8")
            for phrase in required:
                self.assertIn(phrase, text, guide.name)


if __name__ == "__main__":
    unittest.main()
