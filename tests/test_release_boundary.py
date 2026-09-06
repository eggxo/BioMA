"""Guard the public source tree against accidental cluster-only leakage."""

import re
import unittest
from pathlib import Path


# Build markers from fragments so this guard does not itself place a literal
# cluster path in the public source archive it audits.
_PRIVATE_ROOT = "/usr_" + "storage"
_PRIVATE_IP = "192.168." + "25.5"
_PRIVATE_USER = "dan" + "xuming"
PRIVATE_MARKERS = re.compile(
    r"(?:{root}(?:2)?/|/opt/" + "R/" + r"|{ip}|{user})".format(
        root=_PRIVATE_ROOT, ip=re.escape(_PRIVATE_IP), user=re.escape(_PRIVATE_USER)
    ),
    re.IGNORECASE,
)
TEXT_SUFFIXES = {
    ".ini",
    ".md",
    ".py",
    ".r",
    ".R",
    ".toml",
    ".txt",
    ".yml",
    ".yaml",
    ".in",
    ".cff",
}
SKIP_PARTS = {".git", "build", "dist", "__pycache__", "bioma_workflow.egg-info"}


class ReleaseBoundaryTest(unittest.TestCase):
    def test_linux_launcher_uses_lf_shebang(self):
        launcher = (Path(__file__).resolve().parents[1] / "bin" / "bioma").read_bytes()
        self.assertTrue(launcher.startswith(b"#!/usr/bin/env sh\n"))
        self.assertNotIn(b"\r", launcher)

    def test_public_gf_template_does_not_pin_species_site_count(self):
        template = (Path(__file__).resolve().parents[1] / "workflow.example.ini").read_text(
            encoding="utf-8"
        )
        expected_lines = [
            line.strip()
            for line in template.splitlines()
            if line.strip().lower().startswith("expected_sites")
        ]
        self.assertEqual(expected_lines, ["expected_sites ="])

    def test_public_templates_do_not_pin_species_counts_or_names(self):
        root = Path(__file__).resolve().parents[1]
        rona = (root / "workflow.rona.example.ini").read_text(encoding="utf-8")
        self.assertIn("expected_populations =", rona)
        self.assertNotIn("expected_populations = 25", rona)
        load = (root / "workflow.load.example.ini").read_text(encoding="utf-8")
        self.assertIn("expected_future_files = 0", load)
        self.assertNotIn("Four future mean-climate CSVs", load)
        wf = (root / "workflow.wfmoment.example.ini").read_text(encoding="utf-8")
        self.assertNotIn("_pade", wf)

    def test_public_text_has_no_cluster_paths_or_credentials(self):
        root = Path(__file__).resolve().parents[1]
        violations = []
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
                continue
            if path.resolve() == Path(__file__).resolve():
                continue
            if any(part in SKIP_PARTS for part in path.relative_to(root).parts):
                continue
            # Standalone regression profiles are intentionally ignored by the
            # release manifest and may contain production paths for server CI.
            if path.name.endswith(".test.ini") or path.name in {"workflow.50test.ini", "workflow.6140.ini"}:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            if PRIVATE_MARKERS.search(text):
                violations.append(str(path.relative_to(root)))
        self.assertEqual(violations, [], "cluster-only paths found in public files: {}".format(violations))


if __name__ == "__main__":
    unittest.main()
