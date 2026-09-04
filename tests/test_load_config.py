import tempfile
import unittest
from pathlib import Path

from bioma.gf_frequency import InputError
from bioma.load import load_config


class LoadConfigTest(unittest.TestCase):
    def test_bundled_scripts_are_defaults(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for directory in ("vcfs", "populations", "future"):
                (root / directory).mkdir()
            (root / "predictors.csv").write_text("pop,longitude,latitude\n", encoding="utf-8")
            config = root / "load.ini"
            config.write_text(
                "[inputs]\nvcf_dir = vcfs\npopulation_dir = populations\n"
                "predictors = predictors.csv\nfuture_dir = future\n\n"
                "[analysis]\noutput_dir = results\n\n"
                "[parameters]\nexpected_future_files = 0\n",
                encoding="utf-8",
            )

            loaded = load_config(config)

            self.assertEqual(loaded.calc_script.name, "calc_loadM_loadD_from_3vcf.py")
            self.assertEqual(loaded.rf_script.name, "rf_one_target_svd.R")
            self.assertEqual(loaded.predict_script.name, "predict_future_one_target.R")
            self.assertTrue(loaded.calc_script.is_file())
            self.assertTrue(loaded.rf_script.is_file())
            self.assertTrue(loaded.predict_script.is_file())
            self.assertEqual(loaded.expected_future_files, 0)

    def test_unsupported_scheme_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for directory in ("vcfs", "populations", "future"):
                (root / directory).mkdir()
            (root / "predictors.csv").write_text("pop,longitude,latitude\n", encoding="utf-8")
            config = root / "load.ini"
            config.write_text(
                "[inputs]\nvcf_dir = vcfs\npopulation_dir = populations\n"
                "predictors = predictors.csv\nfuture_dir = future\n\n"
                "[analysis]\noutput_dir = results\n\n"
                "[parameters]\nscheme = unsupported\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(InputError, "scheme must be svd"):
                load_config(config)

    def test_rf_script_uses_configured_grid_and_seed(self):
        script = Path(__file__).resolve().parents[1] / "bioma" / "scripts" / "rf_one_target_svd.R"
        script_text = script.read_text(encoding="utf-8")
        self.assertIn("size = grid_size", script_text)
        self.assertIn("seed_base  = seed_base", script_text)


if __name__ == "__main__":
    unittest.main()
