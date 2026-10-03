import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from src.acoustics import rt60_to_absorption
from tests.scenarios import FIXTURES, directional_rir, facade_rir, omni_rir


class RirExactTests(unittest.TestCase):
    def setUp(self):
        self._previous = os.getcwd()
        os.chdir(ROOT)

    def tearDown(self):
        os.chdir(self._previous)

    def test_omni_rir_matches_reference(self):
        actual = omni_rir()
        expected = np.load(FIXTURES / "rir_omni.npy")
        np.testing.assert_array_equal(actual, expected)

    def test_directional_rir_matches_reference(self):
        actual = directional_rir()
        expected = np.load(FIXTURES / "rir_directional.npy")
        np.testing.assert_array_equal(actual, expected)

    def test_facade_rir_matches_reference(self):
        actual = facade_rir()
        expected = np.load(FIXTURES / "rir_facade.npy")
        np.testing.assert_array_equal(actual, expected)

    def test_repeated_call_is_identical(self):
        first = omni_rir()
        second = omni_rir()
        np.testing.assert_array_equal(first, second)

    def test_rt60_to_absorption_matches_reference(self):
        actual = np.asarray(rt60_to_absorption([3.5, 4.0, 2.5], 0.3))
        expected = np.load(FIXTURES / "absorption.npy")
        np.testing.assert_array_equal(actual, expected)


if __name__ == "__main__":
    unittest.main()
