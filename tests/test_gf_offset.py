from __future__ import annotations

import math
import unittest
from pathlib import Path

from bioma.gf_frequency import InputError
from bioma.gf_offset import parse_forward_radii


class GFOffsetInputTest(unittest.TestCase):
    def test_runtime_probe_checks_compressed_table_dependency(self):
        source = (Path(__file__).resolve().parents[1] / "bioma" / "gf_offset.py").read_text(encoding="utf-8")
        self.assertIn("library(R.utils)", source)

    def test_default_style_radii(self):
        values = parse_forward_radii("100,250,500,1000,inf")
        self.assertEqual(values[:4], [100.0, 250.0, 500.0, 1000.0])
        self.assertTrue(math.isinf(values[4]))

    def test_duplicate_radii_are_removed(self):
        self.assertEqual(parse_forward_radii("100,100,250"), [100.0, 250.0])

    def test_invalid_radius_is_blocked(self):
        with self.assertRaises(InputError):
            parse_forward_radii("100,0")


if __name__ == "__main__":
    unittest.main()
