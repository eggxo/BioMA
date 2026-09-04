from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bioma.climate_prepare import discover_bioclim_rasters, discover_scenarios
from bioma.gf_frequency import InputError


class ClimateDiscoveryTest(unittest.TestCase):
    def test_bio_files_are_sorted_numerically(self):
        with tempfile.TemporaryDirectory(prefix="bioma-climate-") as root:
            directory = Path(root)
            for number in range(19, 0, -1):
                (directory / "scenario-bio{}.cut.tif".format(number)).touch()
            rasters = discover_bioclim_rasters(directory)
            self.assertEqual(rasters[0].name, "scenario-bio1.cut.tif")
            self.assertEqual(rasters[18].name, "scenario-bio19.cut.tif")

    def test_missing_bio_is_blocked(self):
        with tempfile.TemporaryDirectory(prefix="bioma-climate-") as root:
            directory = Path(root)
            for number in range(1, 19):
                (directory / "bio_{}.tif".format(number)).touch()
            with self.assertRaises(InputError):
                discover_bioclim_rasters(directory)

    def test_requested_absent_model_is_blocked(self):
        with tempfile.TemporaryDirectory(prefix="bioma-climate-") as root:
            scenario = Path(root) / "2061-2080-ssp245-BCC-CSM2-MR"
            scenario.mkdir()
            for number in range(1, 20):
                (scenario / "bio{}.tif".format(number)).touch()
            with self.assertRaises(InputError):
                discover_scenarios(Path(root), models=["MPI-ESM1-2-HR"])


if __name__ == "__main__":
    unittest.main()
