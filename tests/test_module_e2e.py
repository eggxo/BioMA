"""Minimum end-to-end contract checks for every BioMA module.

These tests intentionally stop at each module's validated dry-run boundary so
CI does not require a licensed MaxEnt jar or site-specific R research packages.
The R/Python numerical fixtures are exercised separately when those external
dependencies are available (see ``test_vulnerability_r.R`` and the module
runtime tests).
"""

import sys
import tempfile
import unittest
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


class ModuleEndToEndSmokeTest(unittest.TestCase):
    def setUp(self):
        self._temporary = tempfile.TemporaryDirectory(prefix="bioma-e2e-")
        self.root = generate(Path(self._temporary.name) / "demo")

    def tearDown(self):
        self._temporary.cleanup()

    def _run(self, module, runner):
        manifest = runner(self.root / "configs" / (module + ".ini"), dry_run=True)
        self.assertEqual(manifest.get("status") or manifest.get("engineering_status"), "planned")
        self.assertTrue((self.root / "configs" / (module + ".ini")).is_file())

    def test_gf_minimum_end_to_end_contract(self):
        self._run("gf", run_gf_workflow)

    def test_rona_minimum_end_to_end_contract(self):
        self._run("rona", run_rona_workflow)

    def test_mar_minimum_end_to_end_contract(self):
        self._run("mar", run_mar_workflow)

    def test_load_minimum_end_to_end_contract(self):
        self._run("load", run_load_workflow)

    def test_niche_minimum_end_to_end_contract(self):
        self._run("niche", run_niche_workflow)

    def test_wfmoment_minimum_end_to_end_contract(self):
        self._run("wfmoment", run_wfmoment_workflow)

    def test_vulnerability_minimum_end_to_end_contract(self):
        self._run("vulnerability", run_vulnerability_workflow)


if __name__ == "__main__":
    unittest.main()
