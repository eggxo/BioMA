from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "data" / "gf_frequency_50"


def read_frequency(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    header = rows[0]
    values = {
        row[0]: [float(value) for value in row[1:]]
        for row in rows[1:]
    }
    return header, values


class GFFrequencyIntegrationTest(unittest.TestCase):
    def setUp(self):
        self.temp_root = Path(tempfile.mkdtemp(prefix="bioma-test-"))

    def tearDown(self):
        shutil.rmtree(str(self.temp_root), ignore_errors=True)

    def run_module(self, vcf: Path, name: str):
        outdir = self.temp_root / name
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "bioma.cli",
                "gf-frequency",
                "--vcf",
                str(vcf),
                "--samples",
                str(FIXTURE / "samples.tsv"),
                "--outdir",
                str(outdir),
                "--expected-sites",
                "50",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        return outdir

    def assert_matches_reference(self, outdir: Path):
        actual_header, actual = read_frequency(
            outdir / "population_alt_frequency.tsv"
        )
        expected_header, expected = read_frequency(
            FIXTURE / "expected_population_alt_frequency.tsv"
        )
        self.assertEqual(actual_header, expected_header)
        self.assertEqual(set(actual), set(expected))
        for population_id in expected:
            self.assertEqual(len(actual[population_id]), 50)
            for observed, reference in zip(actual[population_id], expected[population_id]):
                # The historical vcftools matrix stores six decimal places;
                # BioMA retains the more precise genotype-derived frequency.
                self.assertAlmostEqual(observed, reference, delta=5e-7)

        manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["counts"]["variants"], 50)
        self.assertEqual(manifest["counts"]["populations"], 25)
        self.assertEqual(manifest["counts"]["assigned_samples"], 175)

    def test_50_site_regression(self):
        outdir = self.run_module(FIXTURE / "adaptive_50.vcf", "plain")
        self.assert_matches_reference(outdir)

    def test_file_named_gz_is_detected_by_content(self):
        misleading = self.temp_root / "adaptive_50.vcf.gz"
        shutil.copyfile(str(FIXTURE / "adaptive_50.vcf"), str(misleading))
        outdir = self.run_module(misleading, "misleading-extension")
        self.assert_matches_reference(outdir)

    def test_duplicate_sample_is_blocked(self):
        duplicate = self.temp_root / "duplicate_samples.tsv"
        text = (FIXTURE / "samples.tsv").read_text(encoding="utf-8")
        first_data_line = text.splitlines()[1]
        duplicate.write_text(text + first_data_line + "\n", encoding="utf-8")
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "bioma.cli",
                "gf-frequency",
                "--vcf",
                str(FIXTURE / "adaptive_50.vcf"),
                "--samples",
                str(duplicate),
                "--outdir",
                str(self.temp_root / "duplicate-output"),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("occurs more than once", result.stderr)


if __name__ == "__main__":
    unittest.main()
