import tempfile
import unittest
from pathlib import Path

from bioma.gf_frequency import InputError
from bioma.mar import load_mar_config


class MarConfigurationTest(unittest.TestCase):
    def _config(self, root: Path, scheme: str = "random", steps: str = "data,gm,ext") -> Path:
        vcf = root / "whole_genome.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
            "1\t1\tv1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )
        lonlat = root / "lonlat.tsv"
        lonlat.write_text("ID\tLONGITUDE\tLATITUDE\ns1\t100\t30\n", encoding="utf-8")
        config = root / "mar.ini"
        config.write_text(
            "[inputs]\n"
            "vcf = whole_genome.vcf\n"
            "lonlat = lonlat.tsv\n"
            "scenario_file =\n\n"
            "[analysis]\n"
            "output_dir = results\n\n"
            "[parameters]\n"
            "scheme = {}\n"
            "marsteps = {}\n".format(scheme, steps),
            encoding="utf-8",
        )
        return config

    def test_schemes_match_pinned_mar_package(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for scheme in ("random", "inwards", "outwards", "northsouth", "southnorth"):
                self.assertEqual(load_mar_config(self._config(root, scheme)).scheme, scheme)

    def test_unsupported_direction_alias_is_rejected_during_preflight(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(InputError, "scheme must be"):
                load_mar_config(self._config(Path(temporary), "eastwest"))

    def test_full_workflow_requires_extinction_objects(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(InputError, "require marsteps"):
                load_mar_config(self._config(Path(temporary), steps="data,gm,sfs,mar,plot"))

    def test_optional_scenario_file_is_validated(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self._config(root)
            text = config.read_text(encoding="utf-8").replace(
                "scenario_file =", "scenario_file = missing.tsv"
            )
            config.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(InputError, "scenario_file does not exist"):
                load_mar_config(config)


if __name__ == "__main__":
    unittest.main()
