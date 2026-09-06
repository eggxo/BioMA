import io
import json
import tempfile
import unittest
import zipfile
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

from bioma.cli import main
from bioma.doctor import DoctorCheck, R_PACKAGE_MODULES, _config_context, _maxent_check, run_doctor


class DoctorConfigTest(unittest.TestCase):
    def test_gf_requires_compressed_table_runtime(self):
        self.assertEqual(R_PACKAGE_MODULES["R.utils"], {"gf"})

    def test_standalone_gf_config_is_inferred(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = Path(temporary) / "config.ini"
            config.write_text(
                "[inputs]\nvcf = adaptive.vcf\ncoordinates = env.tsv\n"
                "present_climate = current\n\n[analysis]\noutput_dir = out\n",
                encoding="utf-8",
            )
            context = _config_context(config)
            self.assertEqual(context["modules"], ("gf",))

    def test_project_runtime_hints_and_module_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            niche = root / "niche.ini"
            niche.write_text(
                "[inputs]\nmaxent_jar = jars/maxent.jar\n\n"
                "[analysis]\nrscript = Rscript\n\n",
                encoding="utf-8",
            )
            project = root / "project.ini"
            project.write_text(
                "[project]\noutput_dir = results\n\n"
                "[modules]\nniche = niche.ini\nmar = false\n\n"
                "[shared]\ncompute_python = envs/wf/bin/python\n\n"
                "[doctor]\ngdalinfo = gdalinfo\njava = java\n",
                encoding="utf-8",
            )
            context = _config_context(project)
            self.assertEqual(context["modules"], ("niche",))
            self.assertIn(str((root / "envs/wf/bin/python").resolve()), context["pythons"])
            self.assertIn("Rscript", context["rscripts"])
            self.assertEqual(context["maxent_jar"], str((root / "jars/maxent.jar").resolve()))
            self.assertIn("gdalinfo", context["gdal"])
            self.assertIn("java", context["java"])

    def test_shared_runtime_overrides_module_runtime_hints(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            module = root / "gf.ini"
            module.write_text("[inputs]\n\n[analysis]\n\n[parameters]\nrscript = module-r\n", encoding="utf-8")
            project = root / "project.ini"
            project.write_text(
                "[project]\noutput_dir = results\n\n[modules]\ngf = gf.ini\n\n"
                "[shared]\nrscript = shared-r\n",
                encoding="utf-8",
            )
            context = _config_context(project)
            self.assertEqual(context["rscripts"], ["shared-r"])

    def test_sparse_gdal_hints_keep_their_command_names(self):
        """A single configured GDAL utility must not shift the others."""
        with tempfile.TemporaryDirectory() as temporary:
            config = Path(temporary) / "mar.ini"
            config.write_text(
                "[inputs]\nvcf = input.vcf\n\n"
                "[parameters]\ncompute_python = custom-python\n\n"
                "[doctor]\n"
                "ogrinfo = custom-ogrinfo\n",
                encoding="utf-8",
            )
            context = _config_context(config)
            self.assertEqual(context["gdal_commands"], {"ogrinfo": "custom-ogrinfo"})
            self.assertIn("custom-python", context["pythons"])

    def test_maxent_archive_check_is_read_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            jar = Path(temporary) / "maxent.jar"
            with zipfile.ZipFile(str(jar), "w") as archive:
                archive.writestr("manifest.mf", "Manifest-Version: 1.0\n")
            before = jar.read_bytes()
            result = _maxent_check(str(jar), True)
            self.assertEqual(result.status, "ok")
            self.assertEqual(jar.read_bytes(), before)


class DoctorRuntimeTest(unittest.TestCase):
    def test_cli_runtime_override_is_forwarded(self):
        calls = []

        def probe(name, value, fallback, required, args=("--version",)):
            calls.append((name, value))
            return DoctorCheck(name, "ok", required, "mocked", value, "test")

        with patch("bioma.doctor._command_check", side_effect=probe), patch(
            "bioma.doctor._r_package_check", side_effect=lambda *a, **k: DoctorCheck("pkg", "ok", False)
        ), patch("bioma.doctor._python_package_check", side_effect=lambda *a, **k: DoctorCheck("pkg", "ok", False)), patch(
            "bioma.doctor._maxent_check", return_value=DoctorCheck("maxent_jar", "skipped", False)
        ):
            run_doctor(python="/tmp/custom-python", rscript="/tmp/custom-rscript")
        self.assertEqual(calls[0][0], "python")
        self.assertEqual(calls[0][1], "/tmp/custom-python")
        self.assertEqual(calls[1][0], "rscript")
        self.assertEqual(calls[1][1], "/tmp/custom-rscript")

    def test_unreadable_config_is_reported_as_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            missing = Path(temporary) / "missing.ini"
            report = run_doctor(str(missing))
        self.assertFalse(report["ok"])
        self.assertTrue(any(row["name"] == "configuration" for row in report["checks"]))

    def test_strict_turns_optional_failures_into_failures(self):
        def missing(*args, **kwargs):
            return DoctorCheck("probe", "missing", False, "not installed")

        with patch("bioma.doctor._command_check", side_effect=missing), patch(
            "bioma.doctor._r_package_check", side_effect=missing
        ), patch("bioma.doctor._python_package_check", side_effect=missing), patch(
            "bioma.doctor._maxent_check", side_effect=missing
        ):
            relaxed = run_doctor()
            strict = run_doctor(strict=True)
        self.assertTrue(relaxed["ok"])
        self.assertFalse(strict["ok"])
        self.assertGreater(strict["summary"]["failures"], 0)

    def test_strict_turns_optional_warnings_into_failures(self):
        def warning(*args, **kwargs):
            return DoctorCheck("probe", "warning", False, "not installed")

        with patch("bioma.doctor._command_check", side_effect=warning), patch(
            "bioma.doctor._r_package_check", side_effect=warning
        ), patch("bioma.doctor._python_package_check", side_effect=warning), patch(
            "bioma.doctor._maxent_check", side_effect=warning
        ):
            report = run_doctor(strict=True)
        self.assertFalse(report["ok"])


class DoctorCliTest(unittest.TestCase):
    def test_json_output_and_exit_status(self):
        report = {
            "status": "pass",
            "ok": True,
            "checks": [],
            "summary": {"ok": 0, "warnings": 0, "failures": 0},
        }
        output = io.StringIO()
        with patch("bioma.cli.run_doctor", return_value=report) as mocked:
            with redirect_stdout(output):
                code = main(["doctor", "--json"])
        self.assertEqual(code, 0)
        self.assertEqual(json.loads(output.getvalue())["status"], "pass")
        mocked.assert_called_once()

    def test_doctor_run_alias(self):
        report = {"status": "pass", "ok": True, "checks": [], "summary": {}}
        with patch("bioma.cli.run_doctor", return_value=report):
            with redirect_stdout(io.StringIO()):
                self.assertEqual(main(["doctor-run"]), 0)


if __name__ == "__main__":
    unittest.main()
