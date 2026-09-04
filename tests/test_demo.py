import tempfile
import unittest
from pathlib import Path

import sys

from bioma.gf_workflow import run_gf_workflow
from bioma.load import run_load_workflow
from bioma.mar import run_mar_workflow
from bioma.niche import run_niche_workflow
from bioma.rona import run_rona_workflow
from bioma.vulnerability import run_vulnerability_workflow
from bioma.wfmoment import run_wfmoment_workflow
sys.path.insert(0, str(Path(__file__).resolve().parent / "data" / "demo"))
from make_demo import generate


class PublicDemoSmokeTest(unittest.TestCase):
    def test_generator_has_no_server_private_paths_and_all_modules_validate(self):
        with tempfile.TemporaryDirectory(prefix="bioma-demo-") as temporary:
            root = generate(Path(temporary))
            project_text = (root / "project.ini").read_text(encoding="utf-8")
            self.assertNotIn("/usr_storage", project_text)
            self.assertTrue((root / "data" / "adaptive_sites.vcf").is_file())
            self.assertEqual(len(list((root / "data" / "rona_ld").glob("LD_BIO*.prune.in"))), 19)
            self.assertEqual(len(list((root / "data" / "climate" / "current").glob("*.tif"))), 19)

            runners = {
                "gf": run_gf_workflow,
                "rona": run_rona_workflow,
                "mar": run_mar_workflow,
                "load": run_load_workflow,
                "niche": run_niche_workflow,
                "wfmoment": run_wfmoment_workflow,
                "vulnerability": run_vulnerability_workflow,
            }
            for module, runner in runners.items():
                result = runner(root / "configs" / (module + ".ini"), dry_run=True)
                self.assertEqual(result.get("status") or result.get("engineering_status"), "planned", module)


if __name__ == "__main__":
    unittest.main()
