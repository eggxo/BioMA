import configparser
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bioma.gf_frequency import InputError
from bioma.project import run_project_workflow


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ProjectWorkflowTest(unittest.TestCase):
    def _module_config(self, root, name):
        path = root / "{}.ini".format(name)
        path.write_text(
            "[inputs]\nrelative_input = data/input.tsv\n\n"
            "[analysis]\noutput_dir = ignored\nmodels = old\nssps = 126\nperiods = 2000-2020\n\n"
            "[parameters]\nseed = 1\n",
            encoding="utf-8",
        )
        return path

    def _project_config(self, root, modules):
        output = root / "results"
        module_lines = "\n".join("{} = {}".format(key, value) for key, value in modules.items())
        path = root / "project.ini"
        path.write_text(
            "[project]\nname = test-project\noutput_dir = {}\nresume = true\nreport = true\n\n"
            "[modules]\n{}\n\n"
            "[shared]\nspecies = test_species\nmodels = model-a,model-b\n"
            "ssps = ssp245\nperiods = 2061-2080\nseed = 42\n\n"
            "[integration]\nperiod = 2061-2080\nssp = 245\nensemble_method = mean\n".format(output, module_lines),
            encoding="utf-8",
        )
        return path, output

    @staticmethod
    def _fake_runner(config_path, dry_run=False, progress=None):
        parser = configparser.ConfigParser(interpolation=None)
        parser.read(config_path, encoding="utf-8")
        output = Path(parser.get("analysis", "output_dir"))
        if dry_run:
            output.mkdir(parents=True, exist_ok=True)
            return {"module": "fake", "status": "planned"}
        output.mkdir(parents=True, exist_ok=True)
        manifest = {
            "module": "fake",
            "status": "complete",
            "config_sha256": _sha256(Path(config_path)),
        }
        (output / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return manifest

    def test_standard_output_shared_inheritance_and_resume(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gf = self._module_config(root, "gf")
            project, output = self._project_config(root, {"gf": gf})
            with patch("bioma.project._runner", return_value=self._fake_runner):
                first = run_project_workflow(project)
                second = run_project_workflow(project)

            self.assertEqual(first["status"], "complete")
            self.assertEqual(first["modules"][0]["status"], "completed")
            self.assertEqual(second["modules"][0]["status"], "reused")
            effective = configparser.ConfigParser(interpolation=None)
            effective.read(output / "00_project" / "effective_configs" / "gf.ini", encoding="utf-8")
            self.assertEqual(Path(effective.get("analysis", "output_dir")), output / "01_gf")
            self.assertEqual(effective.get("analysis", "models"), "model-a,model-b")
            self.assertEqual(effective.get("analysis", "ssps"), "245")
            self.assertEqual(effective.get("parameters", "seed"), "42")
            self.assertTrue((output / "report" / "index.html").is_file())
            self.assertTrue((output / "00_project" / "project_summary.tsv").is_file())

    def test_integration_vulnerability_controls_are_forwarded(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            vulnerability = self._module_config(root, "vulnerability")
            project, output = self._project_config(root, {"vulnerability": vulnerability})
            parser = configparser.ConfigParser(interpolation=None)
            parser.read(project, encoding="utf-8")
            parser.set("integration", "rona_variables", "BIO1_RONA,BIO7_RONA")
            parser.set("integration", "rona_summary", "max")
            parser.set("integration", "rona_weights", "BIO1_RONA:0.2,BIO7_RONA:0.8")
            parser.set("integration", "group_order", "group-a,group-b")
            parser.set("integration", "exclude_groups", "excluded")
            with project.open("w", encoding="utf-8") as handle:
                parser.write(handle)

            with patch("bioma.project._runner", return_value=self._fake_runner):
                run_project_workflow(project, dry_run=True)

            effective = configparser.ConfigParser(interpolation=None)
            effective.read(output / "00_project" / "effective_configs" / "vulnerability.ini", encoding="utf-8")
            self.assertEqual(effective.get("parameters", "rona_variables"), "BIO1_RONA,BIO7_RONA")
            self.assertEqual(effective.get("parameters", "rona_summary"), "max")
            self.assertEqual(effective.get("parameters", "rona_weights"), "BIO1_RONA:0.2,BIO7_RONA:0.8")
            self.assertEqual(effective.get("parameters", "group_order"), "group-a,group-b")
            self.assertEqual(effective.get("parameters", "exclude_groups"), "excluded")

    def test_configuration_change_requires_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gf = self._module_config(root, "gf")
            project, output = self._project_config(root, {"gf": gf})
            with patch("bioma.project._runner", return_value=self._fake_runner):
                run_project_workflow(project, dry_run=True)
                effective_path = output / "00_project" / "effective_configs" / "gf.ini"
                resolved_path = output / "00_project" / "resolved_project.json"
                report_text = (output / "report" / "index.html").read_text(encoding="utf-8")
                self.assertIn("validated", report_text)
                self.assertIn("01_gf", report_text)
                old_effective = effective_path.read_bytes()
                old_resolved = resolved_path.read_bytes()
                gf.write_text(gf.read_text(encoding="utf-8") + "\n# changed\n", encoding="utf-8")
                with self.assertRaises(InputError):
                    run_project_workflow(project, dry_run=True)
                self.assertEqual(effective_path.read_bytes(), old_effective)
                self.assertEqual(resolved_path.read_bytes(), old_resolved)

    def test_overwrite_archives_results_and_project_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gf = self._module_config(root, "gf")
            project, output = self._project_config(root, {"gf": gf})
            with patch("bioma.project._runner", return_value=self._fake_runner):
                run_project_workflow(project)
                gf.write_text(gf.read_text(encoding="utf-8") + "\n# changed\n", encoding="utf-8")
                rerun = run_project_workflow(project, overwrite=True)

            self.assertEqual(rerun["modules"][0]["status"], "completed")
            backups = output / "_bioma_backups"
            self.assertTrue(any(path.name == "00_project" for path in backups.rglob("00_project")))
            self.assertTrue(any(path.name == "01_gf" for path in backups.rglob("01_gf")))

    def test_vulnerability_auto_inputs_are_wired(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            configs = {key: self._module_config(root, key) for key in ("gf", "rona", "load", "niche")}
            vulnerability = root / "vulnerability.ini"
            vulnerability.write_text(
                "[inputs]\ngf_offsets = auto\nrona_ensemble = auto\n"
                "load_predictors = auto\nniche_raster = auto\n\n"
                "[analysis]\noutput_dir = ignored\nscenario_label =\n\n[parameters]\n",
                encoding="utf-8",
            )
            configs["vulnerability"] = vulnerability
            project, output = self._project_config(root, configs)
            with patch("bioma.project._runner", return_value=self._fake_runner):
                plan = run_project_workflow(project, dry_run=True)

            self.assertEqual(plan["preflight"]["vulnerability"]["status"], "deferred")
            effective = configparser.ConfigParser(interpolation=None)
            effective.read(output / "00_project" / "effective_configs" / "vulnerability.ini", encoding="utf-8")
            self.assertEqual(
                Path(effective.get("inputs", "gf_offsets")),
                output / "01_gf" / "05_offset_plots" / "2061-2080" / "ssp245" / "ensemble_mean_offsets.tsv.gz",
            )
            self.assertEqual(
                Path(effective.get("inputs", "niche_raster")),
                output / "05_niche" / "rasters" / "maladaptation_2061-2080_ssp245_ensemble_mean.tif",
            )

    def test_failed_upstream_blocks_auto_vulnerability(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gf = self._module_config(root, "gf")
            vulnerability = root / "vulnerability.ini"
            vulnerability.write_text(
                "[inputs]\ngf_offsets = auto\nrona_ensemble = existing/rona.tsv\n"
                "load_predictors = existing/load.csv\nniche_raster = existing/niche.tif\n\n"
                "[analysis]\noutput_dir = ignored\nscenario_label =\n\n[parameters]\n",
                encoding="utf-8",
            )
            for relative in ("existing/rona.tsv", "existing/load.csv", "existing/niche.tif"):
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("fixture", encoding="utf-8")
            project, _ = self._project_config(root, {"gf": gf, "vulnerability": vulnerability})
            parser = configparser.ConfigParser(interpolation=None)
            parser.read(project, encoding="utf-8")
            parser.set("project", "stop_on_error", "false")
            with project.open("w", encoding="utf-8") as handle:
                parser.write(handle)

            def runner(module):
                if module == "gf":
                    def fail(config_path, dry_run=False, progress=None):
                        if dry_run:
                            return {"status": "planned"}
                        raise InputError("expected test failure")
                    return fail
                return self._fake_runner

            with patch("bioma.project._runner", side_effect=runner):
                result = run_project_workflow(project)

            self.assertEqual(result["status"], "complete_with_failures")
            self.assertEqual([row["status"] for row in result["modules"]], ["failed", "blocked"])
            self.assertEqual(result["counts"]["blocked"], 1)

    def test_non_directory_module_output_is_reported_as_invalid(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gf = self._module_config(root, "gf")
            project, output = self._project_config(root, {"gf": gf})
            output.mkdir(parents=True)
            (output / "01_gf").write_text("not a directory", encoding="utf-8")
            with patch("bioma.project._runner", return_value=self._fake_runner):
                with self.assertRaisesRegex(InputError, "invalid"):
                    run_project_workflow(project)

    def test_auto_vulnerability_requires_selected_or_existing_upstream(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            configs = {key: self._module_config(root, key) for key in ("gf", "rona", "load", "niche")}
            vulnerability = root / "vulnerability.ini"
            vulnerability.write_text(
                "[inputs]\ngf_offsets = auto\nrona_ensemble = auto\n"
                "load_predictors = auto\nniche_raster = auto\n\n"
                "[analysis]\noutput_dir = ignored\nscenario_label =\n\n[parameters]\n",
                encoding="utf-8",
            )
            configs["vulnerability"] = vulnerability
            project, _ = self._project_config(root, configs)
            with patch("bioma.project._runner", return_value=self._fake_runner):
                with self.assertRaisesRegex(InputError, "requires selected module"):
                    run_project_workflow(project, dry_run=True, only_modules=("vulnerability",))

    def test_unified_inputs_generate_compatibility_files_and_wire_modules(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            adaptive = root / "adaptive.vcf.gz"
            adaptive.write_text("fixture", encoding="utf-8")
            whole_genome = root / "whole_genome_mis0.9_maf0.00001.vcf.gz"
            whole_genome.write_text("fixture", encoding="utf-8")
            load_vcfs = root / "load_3vcf"
            load_vcfs.mkdir()
            samples = root / "population_samples.tsv"
            samples.write_text(
                "sample_id\tpopulation_id\nS1\tP1\nS2\tP1\nS3\tP2\n",
                encoding="utf-8",
            )
            environment = root / "population_environment.tsv"
            environment.write_text(
                "ID\tpop\tlon\tlat\t" + "\t".join("bio{}".format(i) for i in range(1, 20)) + "\n"
                + "P1\tNorthwest\t105\t32\t" + "\t".join(["1"] * 19) + "\n"
                + "P2\tEast\t119\t30\t" + "\t".join(["2"] * 19) + "\n",
                encoding="utf-8",
            )
            current = root / "climate_current"; current.mkdir()
            future = root / "climate_future"; future.mkdir()
            mask = root / "mask.shp"; mask.write_text("fixture", encoding="utf-8")
            unld = root / "unld"; unld.mkdir()
            project_modules = {key: self._module_config(root, key) for key in ("gf", "rona", "mar", "load", "niche")}
            project, output = self._project_config(root, project_modules)
            parser = configparser.ConfigParser(interpolation=None)
            parser.read(project, encoding="utf-8")
            parser.add_section("inputs")
            values = {
                "adaptive_vcf": adaptive,
                "whole_genome_vcf": whole_genome,
                "load_vcf_dir": load_vcfs,
                "population_samples": samples,
                "population_environment": environment,
                "rona_unld_dir": unld,
                "climate_current": current,
                "climate_future": future,
                "species_mask": mask,
            }
            for key, value in values.items():
                parser.set("inputs", key, str(value))
            with project.open("w", encoding="utf-8") as handle:
                parser.write(handle)

            with patch("bioma.project._runner", return_value=self._fake_runner):
                run_project_workflow(project, dry_run=True)

            effective = configparser.ConfigParser(interpolation=None)
            effective.read(output / "00_project" / "effective_configs" / "gf.ini", encoding="utf-8")
            self.assertEqual(Path(effective.get("inputs", "vcf")), adaptive)
            self.assertTrue((output / "00_project" / "generated_inputs" / "gf_samples.tsv").is_file())
            self.assertTrue((output / "00_project" / "generated_inputs" / "load_population_lists" / "P1").is_file())

            rona = configparser.ConfigParser(interpolation=None)
            rona.read(output / "00_project" / "effective_configs" / "rona.ini", encoding="utf-8")
            self.assertEqual(
                Path(rona.get("inputs", "alt_frequency")),
                output / "01_gf" / "01_frequency" / "population_alt_frequency.tsv",
            )
            mar = configparser.ConfigParser(interpolation=None)
            mar.read(output / "00_project" / "effective_configs" / "mar.ini", encoding="utf-8")
            self.assertTrue(Path(mar.get("inputs", "lonlat")).is_file())
            niche = configparser.ConfigParser(interpolation=None)
            niche.read(output / "00_project" / "effective_configs" / "niche.ini", encoding="utf-8")
            self.assertTrue(Path(niche.get("inputs", "occurrence_csv")).is_file())


if __name__ == "__main__":
    unittest.main()
