from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bioma.gf_frequency import InputError
from bioma.gf_train import read_training_inputs


class GFTrainInputTest(unittest.TestCase):
    def _write(self, path: Path, text: str) -> Path:
        path.write_text(text, encoding="utf-8")
        return path

    def test_population_rows_are_joined_by_id(self):
        with tempfile.TemporaryDirectory(prefix="bioma-gf-train-") as root:
            root_path = Path(root)
            frequency = self._write(
                root_path / "frequency.tsv",
                "population_id\tv1\np2\t0.8\np1\t0.2\np3\t0.4\np4\t0.6\np5\t0.1\n",
            )
            environment = self._write(
                root_path / "environment.tsv",
                "population_id\tbio1\np1\t1\np2\t2\np3\t3\np4\t4\np5\t5\n",
            )
            inputs = read_training_inputs(frequency, environment, predictors=["bio1"])
            self.assertEqual(inputs.population_ids, ["p1", "p2", "p3", "p4", "p5"])
            self.assertEqual(inputs.frequencies["p1"], [0.2])
            self.assertEqual(inputs.frequency_row_indexes["p1"], 2)

    def test_population_mismatch_is_blocked(self):
        with tempfile.TemporaryDirectory(prefix="bioma-gf-train-") as root:
            root_path = Path(root)
            frequency = self._write(
                root_path / "frequency.tsv",
                "population_id\tv1\np1\t0.2\np2\t0.3\np3\t0.4\np4\t0.5\np5\t0.6\n",
            )
            environment = self._write(
                root_path / "environment.tsv",
                "population_id\tbio1\np1\t1\np2\t2\np3\t3\np4\t4\np6\t5\n",
            )
            with self.assertRaises(InputError):
                read_training_inputs(frequency, environment, predictors=["bio1"])


if __name__ == "__main__":
    unittest.main()
