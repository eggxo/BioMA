from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from bioma.gf_frequency import InputError
from bioma.offset_plot import _r_subprocess_env, _resolve_input, _validate_header, plot_offsets


VALID_HEADER = [
    "scenario",
    "period",
    "ssp",
    "model",
    "radius_km",
    "lon",
    "lat",
    "local_offset",
    "forward_offset",
    "reverse_offset",
]


class OffsetPlotInputTest(unittest.TestCase):
    def test_single_model_plot_is_blocked(self):
        with self.assertRaises(InputError):
            plot_offsets(
                inputs=[Path("unused.tsv.gz")],
                output_dir=Path("unused-output"),
                models=["MPI-ESM1-2-HR"],
            )

    def test_duplicate_model_names_do_not_satisfy_ensemble_requirement(self):
        with self.assertRaises(InputError):
            plot_offsets(
                inputs=[Path("unused.tsv.gz")],
                output_dir=Path("unused-output"),
                models=["MPI-ESM1-2-HR", "MPI-ESM1-2-HR"],
            )

    def test_directory_resolves_to_final_offset_table(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table = root / "all_offsets.tsv.gz"
            with gzip.open(str(table), "wt", encoding="utf-8") as handle:
                handle.write("\t".join(VALID_HEADER) + "\n")
            self.assertEqual(_resolve_input(root), table.resolve())

    def test_gzip_header_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory) / "offsets.tsv.gz"
            with gzip.open(str(table), "wt", encoding="utf-8") as handle:
                handle.write("\t".join(VALID_HEADER) + "\n")
            _validate_header(table)

    def test_missing_final_column_is_blocked(self):
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory) / "offsets.tsv"
            table.write_text("\t".join(VALID_HEADER[:-1]) + "\n", encoding="utf-8")
            with self.assertRaises(InputError):
                _validate_header(table)

    def test_conda_r_runtime_is_isolated_from_login_libraries(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory)
            rscript = prefix / "bin" / "Rscript"
            rscript.parent.mkdir()
            rscript.touch()
            (prefix / "conda-meta").mkdir()
            (prefix / "lib").mkdir()
            proj = prefix / "share" / "proj"
            proj.mkdir(parents=True)
            (proj / "proj.db").touch()
            gdal = prefix / "share" / "gdal"
            gdal.mkdir()
            with mock.patch("bioma.runtime.platform.system", return_value="Linux"):
                environment = _r_subprocess_env(str(rscript))
            self.assertEqual(environment["LD_LIBRARY_PATH"], str(prefix / "lib"))
            self.assertEqual(environment["PROJ_DATA"], str(proj))
            self.assertEqual(environment["GDAL_DATA"], str(gdal))


if __name__ == "__main__":
    unittest.main()
