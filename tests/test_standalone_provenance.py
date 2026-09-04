"""Regression checks for standalone module input provenance."""

import tempfile
import unittest
import sys
from pathlib import Path

from bioma.gf_workflow import run_gf_workflow
from bioma.load import run_load_workflow
from bioma.mar import run_mar_workflow
from bioma.niche import run_niche_workflow
from bioma.rona import run_rona_workflow
from bioma.vulnerability import run_vulnerability_workflow
from bioma.wfmoment import run_wfmoment_workflow

sys.path.insert(0, str(Path(__file__).resolve().parent / "data" / "demo"))
from make_demo import generate  # noqa: E402


class StandaloneProvenanceTest(unittest.TestCase):
    def test_each_module_records_content_fingerprints_in_dry_run(self):
        runners = {
            "gf": run_gf_workflow,
            "rona": run_rona_workflow,
            "mar": run_mar_workflow,
            "load": run_load_workflow,
            "niche": run_niche_workflow,
            "wfmoment": run_wfmoment_workflow,
            "vulnerability": run_vulnerability_workflow,
        }
        with tempfile.TemporaryDirectory(prefix="bioma-provenance-") as temporary:
            root = generate(Path(temporary) / "demo")
            for module, runner in runners.items():
                with self.subTest(module=module):
                    manifest = runner(root / "configs" / (module + ".ini"), dry_run=True)
                    fingerprints = manifest.get("input_fingerprints")
                    self.assertIsInstance(fingerprints, dict)
                    self.assertTrue(manifest.get("input_signature_sha256"))
                    self.assertEqual(
                        manifest["provenance"]["input_signature_sha256"],
                        manifest["input_signature_sha256"],
                    )
                    self.assertTrue(all("kind" in value for value in fingerprints.values()))
                    self.assertTrue(
                        all(
                            value.get("kind") in {"file", "directory", "symlink"}
                            and value.get("sha256")
                            for value in fingerprints.values()
                            if value.get("kind") != "absent"
                        )
                    )

            # Optional inputs remain explicit and do not get mistaken for a
            # missing required file.
            self.assertEqual(
                run_mar_workflow(root / "configs" / "mar.ini", dry_run=True)[
                    "input_fingerprints"
                ]["scenario_file"]["kind"],
                "absent",
            )
            self.assertEqual(
                run_wfmoment_workflow(root / "configs" / "wfmoment.ini", dry_run=True)[
                    "input_fingerprints"
                ]["structure_file"]["kind"],
                "absent",
            )

    def test_directory_content_change_changes_signature(self):
        with tempfile.TemporaryDirectory(prefix="bioma-provenance-") as temporary:
            root = generate(Path(temporary) / "demo")
            config = root / "configs" / "rona.ini"
            first = run_rona_workflow(config, dry_run=True)
            (root / "data" / "rona_ld" / "LD_BIO1.prune.in").write_text(
                "adaptive_1\nadaptive_new\n", encoding="utf-8"
            )
            second = run_rona_workflow(config, dry_run=True)
            self.assertNotEqual(
                first["input_fingerprints"]["unld_dir"]["sha256"],
                second["input_fingerprints"]["unld_dir"]["sha256"],
            )
            self.assertNotEqual(
                first["input_signature_sha256"], second["input_signature_sha256"]
            )


if __name__ == "__main__":
    unittest.main()
