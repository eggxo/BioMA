import io
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

from bioma.cli import main


class ProjectCliTest(unittest.TestCase):
    def test_project_with_recorded_failures_returns_nonzero(self):
        with patch(
            "bioma.cli.run_project_workflow",
            return_value={"status": "complete_with_failures"},
        ):
            with redirect_stdout(io.StringIO()):
                exit_code = main(["project", "project.ini"])
        self.assertEqual(exit_code, 1)


if __name__ == "__main__":
    unittest.main()
