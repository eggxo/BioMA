from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bioma.gf_workflow import (
    _materialize_sample_design,
    build_workflow_plan,
    load_workflow_config,
    run_gf_workflow,
)


class GFWorkflowTest(unittest.TestCase):
    def test_sample_group_directory_becomes_two_column_design(self):
        with tempfile.TemporaryDirectory(prefix="bioma-workflow-") as root:
            root_path = Path(root)
            groups = root_path / "groups"
            groups.mkdir()
            (groups / "population-B").write_text("sample-3\n", encoding="utf-8")
            (groups / "population-A").write_text("sample-1\nsample-2\n", encoding="utf-8")
            (groups / ".ignored.swp").write_text("junk\n", encoding="utf-8")
            output = _materialize_sample_design(groups, root_path / "samples.tsv")
            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "sample_id\tpopulation_id\n"
                "sample-1\tpopulation-A\n"
                "sample-2\tpopulation-A\n"
                "sample-3\tpopulation-B\n",
            )

    def _fixture(self, root: Path) -> Path:
        for name in ("adaptive.vcf", "samples.tsv", "coordinates.tsv", "mask.shp"):
            (root / name).write_text("fixture\n", encoding="utf-8")
        present = root / "present"
        future = root / "future"
        present.mkdir()
        future.mkdir()
        for number in range(1, 20):
            (present / "bio{}.tif".format(number)).touch()
        for model in ("Model-A", "Model-B"):
            scenario = future / "2061-2080-ssp245-{}".format(model)
            scenario.mkdir()
            for number in range(1, 20):
                (scenario / "bio{}.tif".format(number)).touch()
        config = root / "workflow.ini"
        config.write_text(
            """[inputs]
vcf = adaptive.vcf
samples = samples.tsv
coordinates = coordinates.tsv
present_climate = present
future_climate = future
current_mask = mask.shp

[analysis]
output_dir = output
models = Model-A, Model-B
ssps = 245
periods = 2061-2080
expected_sites = 50
""",
            encoding="utf-8",
        )
        return config

    def test_config_builds_complete_scenario_plan(self):
        with tempfile.TemporaryDirectory(prefix="bioma-workflow-") as root:
            config = load_workflow_config(self._fixture(Path(root)))
            plan = build_workflow_plan(config)
            self.assertEqual(plan.models, ("Model-A", "Model-B"))
            self.assertEqual(plan.ssps, ("245",))
            self.assertEqual(len(plan.scenarios), 2)

    def test_runner_executes_all_stages_and_then_resumes(self):
        with tempfile.TemporaryDirectory(prefix="bioma-workflow-") as root:
            config_path = self._fixture(Path(root))
            calls = []

            def fake(module, outputs):
                def run(**kwargs):
                    output_dir = kwargs["output_dir"]
                    output_dir.mkdir(parents=True)
                    for name in outputs:
                        path = output_dir / name
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_text("fixture\n", encoding="utf-8")
                    manifest = {
                        "module": module,
                        "engineering_status": "completed",
                        "scientific_status": "completed_validated",
                    }
                    (output_dir / "run_manifest.json").write_text(
                        json.dumps(manifest), encoding="utf-8"
                    )
                    calls.append(module)
                    return manifest

                return run

            with patch(
                "bioma.gf_workflow.build_population_frequency",
                side_effect=fake("gf-frequency", ("population_alt_frequency.tsv", "variant_manifest.tsv")),
            ), patch(
                "bioma.gf_workflow.prepare_climate_inputs",
                side_effect=fake(
                    "climate-prepare",
                    ("population_environment.tsv", "current_background.tsv.gz", "future_manifest.tsv"),
                ),
            ), patch(
                "bioma.gf_workflow.train_gradient_forest",
                side_effect=fake("gf-train", ("all_gfmod.data", "training_alignment.tsv")),
            ), patch(
                "bioma.gf_workflow.compute_gf_offsets",
                side_effect=fake("gf-offset", ("all_offsets.tsv.gz", "scenario_summary.tsv")),
            ), patch(
                "bioma.gf_workflow.plot_offsets",
                side_effect=fake(
                    "offset-plot",
                    ("ensemble_mean_offsets.tsv.gz", "all_offsets_maps.png", "offset_relationships.png"),
                ),
            ), patch(
                "bioma.gf_workflow.plot_forward_distance",
                side_effect=fake(
                    "forward-distance-plot",
                    ("forward_offset_by_distance.tsv", "forward_offset_by_distance.png"),
                ),
            ):
                first = run_gf_workflow(config_path)
                first_call_count = len(calls)
                second = run_gf_workflow(config_path)

            self.assertEqual(first["engineering_status"], "completed")
            self.assertEqual(first["counts"]["stages"], 7)
            self.assertEqual(first["counts"]["reused_stages"], 0)
            self.assertEqual(first_call_count, 7)
            self.assertEqual(len(calls), first_call_count)
            self.assertEqual(second["counts"]["reused_stages"], 7)


if __name__ == "__main__":
    unittest.main()
