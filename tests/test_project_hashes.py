import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bioma.project import (
    _fingerprint_path,
    _fingerprint_paths,
    run_project_workflow,
)


class ProjectInputFingerprintTest(unittest.TestCase):
    def test_directory_manifest_is_stable_and_content_addressed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first"
            second = root / "second"
            for directory in (first, second):
                (directory / "nested").mkdir(parents=True)
                (directory / "nested" / "b.bin").write_bytes(b"second")
                (directory / "a.txt").write_bytes(b"first")

            left = _fingerprint_path(first)
            right = _fingerprint_path(second)
            self.assertEqual(left["kind"], "directory")
            self.assertEqual(left["sha256"], right["sha256"])
            self.assertEqual(left["entry_count"], 2)
            self.assertEqual(left["total_size_bytes"], len(b"first") + len(b"second"))
            self.assertEqual(
                [entry["relative_path"] for entry in left["files"]],
                ["a.txt", "nested/b.bin"],
            )

            # A byte-level change invalidates the aggregate, independently of
            # mtime and absolute directory location.
            (second / "nested" / "b.bin").write_bytes(b"changed")
            self.assertNotEqual(left["sha256"], _fingerprint_path(second)["sha256"])

    def test_shapefile_fingerprint_includes_sidecars(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            shp = root / "range.shp"
            shp.write_bytes(b"geometry")
            (root / "range.shx").write_bytes(b"index")
            (root / "range.dbf").write_bytes(b"attributes")
            (root / "range.prj").write_bytes(b"projection")
            # Similar-looking files with another stem are not companions.
            (root / "ranger.txt").write_bytes(b"unrelated")

            record = _fingerprint_path(shp)
            names = [item["relative_path"] for item in record["components"]]
            self.assertEqual(names, ["range.dbf", "range.prj", "range.shp", "range.shx"])
            before = record["sha256"]
            (root / "range.prj").write_bytes(b"changed projection")
            changed = _fingerprint_path(shp)
            self.assertNotEqual(before, changed["sha256"])
            self.assertEqual(changed["primary_sha256"], record["primary_sha256"])

    def test_fingerprint_cache_reuses_shared_path(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.txt"
            path.write_text("fixture", encoding="utf-8")
            records = _fingerprint_paths({"one": path, "two": path})
            self.assertEqual(records["one"]["sha256"], records["two"]["sha256"])

    @staticmethod
    def _fake_runner(config_path, dry_run=False, progress=None):
        return {"module": "fake", "status": "planned" if dry_run else "complete"}

    def test_content_change_is_detected_before_metadata_mutation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "data.tsv"
            input_path.write_text("one\n", encoding="utf-8")
            module = root / "gf.ini"
            module.write_text(
                "[inputs]\ncustom_input = data.tsv\n\n"
                "[analysis]\noutput_dir = ignored\n\n[parameters]\n",
                encoding="utf-8",
            )
            project = root / "project.ini"
            output = root / "results"
            project.write_text(
                "[project]\nname = hash-test\noutput_dir = {}\nresume = true\nreport = false\n\n"
                "[modules]\ngf = {}\n".format(output, module)
                + "\n[shared]\n",
                encoding="utf-8",
            )

            with patch("bioma.project._runner", return_value=self._fake_runner):
                run_project_workflow(project, dry_run=True)

            project_dir = output / "00_project"
            old_resolved = (project_dir / "resolved_project.json").read_bytes()
            old_contract = (project_dir / "input_contract.tsv").read_bytes()
            old_manifest = (project_dir / "input_manifest.json").read_bytes()
            input_path.write_text("two\n", encoding="utf-8")

            with patch("bioma.project._runner", return_value=self._fake_runner), patch(
                "bioma.project._materialize_project_inputs",
                side_effect=AssertionError("materialization must not run after a signature mismatch"),
            ):
                with self.assertRaisesRegex(Exception, "Project configuration changed"):
                    run_project_workflow(project, dry_run=True)

            self.assertEqual(old_resolved, (project_dir / "resolved_project.json").read_bytes())
            self.assertEqual(old_contract, (project_dir / "input_contract.tsv").read_bytes())
            self.assertEqual(old_manifest, (project_dir / "input_manifest.json").read_bytes())

            resolved = json.loads(old_resolved.decode("utf-8"))
            self.assertIn("input_fingerprints", resolved)
            self.assertTrue(resolved["input_fingerprints"])

    def test_manifest_tampering_is_detected_before_reuse(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module = root / "gf.ini"
            module.write_text("[analysis]\noutput_dir = ignored\n", encoding="utf-8")
            project = root / "project.ini"
            output = root / "results"
            project.write_text(
                "[project]\nname = manifest-test\noutput_dir = {}\nreport = false\n\n"
                "[modules]\ngf = {}\n".format(output, module)
                + "\n[shared]\n",
                encoding="utf-8",
            )
            with patch("bioma.project._runner", return_value=self._fake_runner):
                run_project_workflow(project, dry_run=True)
                manifest_path = output / "00_project" / "input_manifest.json"
                manifest_path.write_text(manifest_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
                with self.assertRaisesRegex(Exception, "input manifest changed"):
                    run_project_workflow(project, dry_run=True)


if __name__ == "__main__":
    unittest.main()
