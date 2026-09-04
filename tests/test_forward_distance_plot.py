from __future__ import annotations

import unittest
from pathlib import Path

from bioma.forward_distance_plot import _normalize_ssp, plot_forward_distance
from bioma.gf_frequency import InputError


class ForwardDistancePlotInputTest(unittest.TestCase):
    def test_plot_uses_ordered_categories_without_axis_break_slashes(self):
        helper = Path(__file__).resolve().parents[1] / "bioma" / "scripts" / "plot_forward_by_distance.R"
        source = helper.read_text(encoding="utf-8")
        self.assertIn("x_pos = seq_along(finite_values)", source)
        self.assertNotIn('label = "//"', source)

    def test_numeric_ssp_is_normalized(self):
        self.assertEqual(_normalize_ssp("245"), "ssp245")

    def test_prefixed_ssp_is_normalized(self):
        self.assertEqual(_normalize_ssp("SSP585"), "ssp585")

    def test_invalid_ssp_is_blocked(self):
        with self.assertRaises(InputError):
            _normalize_ssp("RCP85")

    def test_single_model_is_blocked(self):
        with self.assertRaises(InputError):
            plot_forward_distance(
                inputs=[Path("unused.tsv.gz")],
                output_dir=Path("unused-output"),
                models=["BCC-CSM2-MR"],
            )

    def test_minimum_models_cannot_enable_single_model_plot(self):
        with self.assertRaises(InputError):
            plot_forward_distance(
                inputs=[Path("unused.tsv.gz")],
                output_dir=Path("unused-output"),
                minimum_models=1,
            )


if __name__ == "__main__":
    unittest.main()
