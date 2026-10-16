"""Static contracts for edge cases in the MaxEnt projection script."""

import unittest
import tempfile
from pathlib import Path

from bioma.niche import load_config


class NicheProjectionContractTest(unittest.TestCase):
    def test_single_gcm_is_used_directly_as_ensemble(self):
        script = (
            Path(__file__).resolve().parents[1]
            / "bioma"
            / "scripts"
            / "niche_pipeline.R"
        ).read_text(encoding="utf-8")
        self.assertIn("if (length(keys) == 1L)", script)
        self.assertIn("ensemble <- future_preds[[keys[[1L]]]]", script)

    def test_range_change_configuration_is_optional_and_recorded(self):
        with tempfile.TemporaryDirectory(prefix="bioma-niche-range-contract-") as tmp:
            root = Path(tmp)
            config = root / "niche.ini"
            config.write_text(
                """[inputs]\nwork_dir = .\noccurrence_csv = occurrence.csv\ncurrent_env_dir = current\nfuture_root = future\nmask_shp = mask.shp\n\n[variables]\ncorrelation_threshold = 0.8\n\n[tuning]\nmax_models = 1\nbackground_n = 10\ncv_folds = 2\nfeature_classes = L\nbeta_multipliers = 1\n\n[projection]\nperiods = 2061-2080\nscenarios = ssp245\ngcms = A,B,C\nensemble_method = mean\nbinary_output = false\n\n[range_change]\nenabled = true\nthreshold = auto\nconsensus = majority\nminimum_gcm_agreement = auto\n\n[analysis]\noutput_dir = output\nrscript = Rscript\nplot_rscript = Rscript\n""",
                encoding="utf-8",
            )
            for path in (root / "occurrence.csv", root / "mask.shp"):
                path.write_text("placeholder\n", encoding="utf-8")
            (root / "current").mkdir()
            (root / "future").mkdir()
            cfg = load_config(config)
            self.assertTrue(cfg.range_change_enabled)
            self.assertEqual(cfg.range_threshold, "auto")
            self.assertEqual(cfg.range_consensus, "majority")
            self.assertEqual(cfg.range_minimum_gcm_agreement, "auto")

    def test_range_change_script_defines_transition_classes_and_area_summary(self):
        script = (
            Path(__file__).resolve().parents[1]
            / "bioma"
            / "scripts"
            / "niche_range_change.R"
        ).read_text(encoding="utf-8")
        self.assertIn("transition_class_", script)
        self.assertIn("lost_pct_current", script)
        self.assertIn("future_binary", script)


if __name__ == "__main__":
    unittest.main()
