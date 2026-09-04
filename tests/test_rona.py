import json
import tempfile
import unittest
from pathlib import Path

from bioma.rona import _ld_prune_records, run_rona_workflow


class RonaInputSmokeTest(unittest.TestCase):
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
            unld = root / "unld"
            future = root / "future" / "2061-2080-ssp245-demo"
            unld.mkdir(parents=True)
            future.mkdir(parents=True)
            for bio in range(1, 20):
                (unld / "LD_BIO{}.prune.in".format(bio)).write_text("v{}\n".format(bio), encoding="utf-8")
                (future / "2061-2080-ssp245-demobio{}.cut.tif".format(bio)).write_text("demo", encoding="utf-8")
            (root / "frequency.tsv").write_text("Pop\tv1\nP1\t0.1\n", encoding="utf-8")
            (root / "environment.tsv").write_text("ID\tpop\tlon\tlat\tbio1\nP1\tA\t1\t1\t1\n", encoding="utf-8")
            config = root / "rona.ini"
            config.write_text(
                "[inputs]\n"
                "alt_frequency = frequency.tsv\n"
                "unld_dir = unld\n"
                "environment = environment.tsv\n"
                "future_climate = future\n\n"
                "[analysis]\noutput_dir = result\nmodels = demo\nssps = 245\nperiods = 2061-2080\n\n"
                "[parameters]\nrscript = Rscript\n",
                encoding="utf-8",
            )
            manifest = run_rona_workflow(config, dry_run=True)
            self.assertEqual(manifest["status"], "planned")
            self.assertEqual(len(manifest["inputs"]["ld_pruning"]["files"]), 19)
            saved = json.loads((root / "result" / "rona_dry_run.json").read_text(encoding="utf-8"))
            self.assertEqual(saved["inputs"]["ld_pruning"]["files"][4]["bio"], 5)


if __name__ == "__main__":
    unittest.main()
