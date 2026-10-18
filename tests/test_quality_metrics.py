import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tests.quality_metrics import (
    DIRECT_PATH_TOLERANCE_SAMPLES,
    FIXTURE_ATOL,
    evaluate_fixtures,
    repeat_is_identical,
)
from tests.scenarios import omni_rir


class QualityMetricTests(unittest.TestCase):
    def setUp(self):
        self._previous = os.getcwd()
        os.chdir(ROOT)

    def tearDown(self):
        os.chdir(self._previous)

    def test_fixtures_within_atol(self):
        for row in evaluate_fixtures():
            self.assertLessEqual(
                row["max_abs"],
                FIXTURE_ATOL,
                "%s max abs %s exceeds %s" % (row["key"], row["max_abs"], FIXTURE_ATOL),
            )
            self.assertTrue(row["within_atol"])

    def test_repeat_call_is_identical(self):
        self.assertTrue(repeat_is_identical())

    def test_omni_direct_path_within_one_sample(self):
        rows = evaluate_fixtures()
        omni = next(row for row in rows if row["key"] == "omni")
        report = omni["direct_paths"][0]
        self.assertLessEqual(abs(report["error_samples"]), DIRECT_PATH_TOLERANCE_SAMPLES, report)
        self.assertTrue(report["within_tolerance"])
        # The scenario helper and the metric agree on the same impulse.
        self.assertEqual(omni["samples"], omni_rir().shape[0])


if __name__ == "__main__":
    unittest.main()
