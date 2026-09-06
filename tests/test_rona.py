import json
import tempfile
import unittest
from pathlib import Path

from bioma.rona import _ld_prune_records, run_rona_workflow


class RonaInputSmokeTest(unittest.TestCase):
    @staticmethod
    def _write_valid_inputs(root, *, expected="", ld_overlap=True):
        unld = root / "unld"
        future = root / "future" / "2061-2080-ssp245-demo"
        unld.mkdir(parents=True)
        future.mkdir(parents=True)
        loci = ["chr1:{}".format(i) for i in range(1, 20)]
        for bio, locus in enumerate(loci, start=1):
            listed = locus if ld_overlap else "absent:{}".format(bio)
            (unld / "LD_BIO{}.prune.in".format(bio)).write_text(
                listed + "\n", encoding="utf-8"
            )
            (future / "2061-2080-ssp245-demobio{}.cut.tif".format(bio)).write_text(
                "demo", encoding="utf-8"
            )
        header = "Pop\t" + "\t".join(loci) + "\n"
        rows = "".join(
            pop + "\t" + "\t".join(["0.1"] * 19) + "\n"
            for pop in ("P1", "P2", "P3")
        )
        (root / "frequency.tsv").write_text(header + rows, encoding="utf-8")
        env_header = "ID\tpop\tlon\tlat\t" + "\t".join(
            "bio{}".format(i) for i in range(1, 20)
        ) + "\n"
        env_rows = "".join(
            "{}\t{}\t{}\t{}\t{}\n".format(
                pop, group, index, index, "\t".join([str(index)] * 19)
            )
            for index, (pop, group) in enumerate(
                (("P1", "A"), ("P2", "B"), ("P3", "C")), start=1
            )
        )
        (root / "environment.tsv").write_text(env_header + env_rows, encoding="utf-8")
        config = root / "rona.ini"
        config.write_text(
            "[inputs]\n"
            "alt_frequency = frequency.tsv\n"
            "unld_dir = unld\n"
            "environment = environment.tsv\n"
            "future_climate = future\n\n"
            "[analysis]\noutput_dir = result\nmodels = demo\nssps = 245\n"
            "periods = 2061-2080\nexpected_populations = {}\n\n".format(expected)
            + "[parameters]\nrscript = Rscript\n",
            encoding="utf-8",
        )
        return config

    def test_ld_lists_are_required_and_fingerprinted(self):
        with tempfile.TemporaryDirectory(prefix="bioma-rona-") as temporary:
            root = Path(temporary)
            unld = root / "unld"
            unld.mkdir()
            for bio in range(1, 20):
                (unld / "LD_BIO{}.prune.in".format(bio)).write_text(
                    "chr1:1\nchr1:2\n", encoding="utf-8"
                )
            records = _ld_prune_records(unld)
            self.assertEqual(len(records), 19)
            self.assertTrue(all(row["sha256"] for row in records))
            self.assertEqual(records[0]["loci"], 2)
            unique_records, unique_count = _ld_prune_records(unld, include_unique=True)
            self.assertEqual(len(unique_records), 19)
            self.assertEqual(unique_count, 2)
            (unld / "LD_BIO19.prune.in").unlink()
            with self.assertRaisesRegex(Exception, "missing LD-pruning"):
                _ld_prune_records(unld)

    def test_dry_run_records_ld_inputs_without_server_paths(self):
        with tempfile.TemporaryDirectory(prefix="bioma-rona-") as temporary:
            root = Path(temporary)
            config = self._write_valid_inputs(root)
            manifest = run_rona_workflow(config, dry_run=True)
            self.assertEqual(manifest["status"], "planned")
            self.assertEqual(len(manifest["inputs"]["ld_pruning"]["files"]), 19)
            self.assertEqual(manifest["counts"]["overlapping_populations"], 3)
            self.assertEqual(
                manifest["inputs"]["ld_pruning"]["files"][4]["frequency_overlap_loci"],
                1,
            )
            saved = json.loads((root / "result" / "rona_dry_run.json").read_text(encoding="utf-8"))
            self.assertEqual(saved["inputs"]["ld_pruning"]["files"][4]["bio"], 5)

    def test_expected_population_count_is_enforced(self):
        with tempfile.TemporaryDirectory(prefix="bioma-rona-") as temporary:
            config = self._write_valid_inputs(Path(temporary), expected="4")
            with self.assertRaisesRegex(Exception, "found 3.*expected 4"):
                run_rona_workflow(config, dry_run=True)

    def test_ld_lists_must_overlap_frequency_columns(self):
        with tempfile.TemporaryDirectory(prefix="bioma-rona-") as temporary:
            config = self._write_valid_inputs(Path(temporary), ld_overlap=False)
            with self.assertRaisesRegex(Exception, "no frequency-column overlap"):
                run_rona_workflow(config, dry_run=True)


if __name__ == "__main__":
    unittest.main()
