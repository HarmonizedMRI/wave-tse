"""Regression tests for an unidentifiable low-amplitude Wave constant phase."""

import unittest
from unittest.mock import patch

import numpy as np

from utils import wave_psf_coefficient_processing as processing


class LowAmplitudeConstantPhaseTests(unittest.TestCase):
    """Require explicit opt-in, accurate constrained recovery and amplitude rejection."""

    def evaluate(self, amplitude, enabled):
        """Fit synthetic independent projections after simulating a c frequency failure."""
        x = np.arange(512, dtype=float)
        w = 2 * np.pi * 4 / len(x)
        a, b = np.sin(w * x), .7 * np.cos(w * x + .2)
        c = .03 + amplitude * np.sin(w * x + .4)
        original = processing._fit_one_coefficient

        def fail_unconstrained_c(name, *args):
            """Represent a c term with insufficient independent frequency information."""
            if name == "c":
                raise ValueError("c frequency is unidentifiable")
            return original(name, *args)

        with patch.object(processing, "_fit_one_coefficient", side_effect=fail_unconstrained_c):
            return processing.process_wave_yz_coefficients(
                a, b, c / 2, c / 2, 128, 384, expected_cycles=4,
                allow_low_amplitude_c_frequency_locked_fallback=enabled), c

    def test_default_rejects_unidentifiable_c(self):
        """Preserve default validation behavior when fallback has not been requested."""
        with self.assertRaisesRegex(ValueError, "unidentifiable"):
            self.evaluate(.02, False)

    def test_opt_in_recovers_small_constant_phase(self):
        """Recover the small programmed-frequency curve and identify its frequency provenance."""
        result, truth = self.evaluate(.02, True)
        self.assertLess(np.max(abs(result["c_processed"].numpy() - truth)), .002)
        report = result["report"]
        self.assertTrue(report["validation_passed"])
        self.assertTrue(report["coefficients"]["c"]["frequency_locked"])
        self.assertIsNone(report["frequency_relative_errors"]["c"])

    def test_opt_in_still_rejects_large_constant_phase(self):
        """Prevent the weak-signal branch from admitting a strong misidentified c term."""
        with self.assertRaisesRegex(ValueError, "amplitude_at_most"):
            self.evaluate(.2, True)


if __name__ == "__main__":
    unittest.main()
