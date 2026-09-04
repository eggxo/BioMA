"""Static contracts for edge cases in the MaxEnt projection script."""

import unittest
from pathlib import Path


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


if __name__ == "__main__":
    unittest.main()
