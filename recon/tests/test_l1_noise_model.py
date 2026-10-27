"""Tests of the noise-model primitives in utils/l1_noise_model.py on synthetic data."""

import unittest

import numpy as np
import torch

from utils.l1_noise_model import (
    colored_noise, far_margin_rows, margin_covariance, margin_leakage, margin_whitener, normal_diagonal,
    outer_lines, region_rms, whitened_hybrid_power,
)
from utils.mb2_cg_sense_wave import mb2_wave_forward
from utils.wave_cg_sense_precondition import fftc_dim, ifftc_dim


def random_covariance(ncoil, seed):
    """Return a random Hermitian positive-definite coil covariance with correlated coils."""
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((ncoil, ncoil)) + 1j * rng.standard_normal((ncoil, ncoil))
    return a @ a.conj().T / ncoil + np.eye(ncoil)


def receiver_profile(nx_os):
    """Return a plateau readout profile that rolls off linearly to 0.3 at both edges."""
    u = abs(np.arange(nx_os) - nx_os / 2) / (nx_os / 2)
    return np.where(u < .6, 1.3, 1.3 - (u - .6) / .4)


class SelectionTests(unittest.TestCase):
    """Margin and outer-line selections."""

    def test_far_margins_and_outer_lines(self):
        """Margins are the outer 20% at both ends; outer lines keep |line - Ny/2| > 0.35 Ny."""
        np.testing.assert_array_equal(far_margin_rows(10, .2), [0, 1, 8, 9])
        self.assertEqual(len(far_margin_rows(2760, .2)), 2 * 552)
        for fraction in (.01, .5):
            with self.assertRaises(ValueError):
                far_margin_rows(10, fraction)
        lines = np.arange(0, 920, 4)
        selected = outer_lines(lines, 920, .3)
        self.assertTrue(np.all(abs(selected - 460) > 322))
        self.assertEqual(len(selected), np.sum(abs(lines - 460) > 322))


class NoiseModelTests(unittest.TestCase):
    """Whitening, colored noise and the margin covariance."""

    def test_whitener_keeps_the_mean_level(self):
        """W cov W^H is the mean coil variance times the identity."""
        cov = random_covariance(6, 1)
        w = margin_whitener(cov)
        np.testing.assert_allclose(w @ cov @ w.conj().T, np.real(np.diag(cov)).mean() * np.eye(6), atol=1e-10)

    def test_colored_noise_follows_the_model(self):
        """Hybrid noise has covariance q(x) * cov, the margin covariance recovers it, and seeds repeat."""
        cov = random_covariance(4, 2)
        profile = receiver_profile(60)
        noise = colored_noise(cov, profile, (4, 60, 4000), seed=7)
        self.assertTrue(torch.equal(noise, colored_noise(cov, profile, (4, 60, 4000), seed=7)))
        hybrid = ifftc_dim(noise, dim=1).to(torch.complex128).numpy()
        for region in (slice(0, 6), slice(25, 35)):
            samples = hybrid[:, region].reshape(4, -1)
            measured = samples @ samples.conj().T / samples.shape[1]
            np.testing.assert_allclose(measured, profile[region].mean() * cov, atol=.03 * abs(cov).max())
        margin = far_margin_rows(60, .2)
        estimate = margin_covariance(noise, np.arange(4000), margin)
        np.testing.assert_allclose(estimate, profile[margin].mean() * cov, atol=.03 * abs(cov).max())
        whitener = margin_whitener(estimate)
        power = whitened_hybrid_power(noise, whitener, np.arange(4000))
        level = np.real(np.diag(estimate)).mean()
        self.assertAlmostEqual(power[margin].mean() / level, 1, delta=.01)
        self.assertAlmostEqual(power[25:35].mean() / level, 1.3 / profile[margin].mean(), delta=.03)


class OperatorScaleTests(unittest.TestCase):
    """Normal-operator diagonal, region RMS and margin leakage."""

    def test_normal_diagonal_matches_the_dense_operator(self):
        """line_fraction * sum_c |S|^2 equals the diagonal of A^H A of a Wave MB2 operator."""
        rng = np.random.default_rng(3)
        nx_os, ny, rows = 24, 16, slice(8, 16)
        sens = np.zeros((2, 3, nx_os, ny), np.complex64)
        sens[:, :, rows] = rng.standard_normal((2, 3, 8, ny)) + 1j * rng.standard_normal((2, 3, 8, ny))
        sens, phase, psf = (torch.as_tensor(v) for v in (
            sens, np.exp(1j * rng.uniform(0, 2 * np.pi, (2, ny))).astype(np.complex64),
            np.exp(1j * rng.uniform(-2, 2, (2, nx_os, ny))).astype(np.complex64)))
        mask = torch.zeros((1, 1, ny))
        mask[..., ::3] = 1
        diagonal = np.zeros((2, nx_os, ny))
        for b in range(2):
            for x in range(nx_os):
                for y in range(ny):
                    e = torch.zeros((2, nx_os, ny), dtype=torch.complex64)
                    e[b, x, y] = 1
                    diagonal[b, x, y] = float(torch.sum(abs(mb2_wave_forward(e, sens, phase, psf, mask))**2))
        regions = rng.uniform(size=(2, 8, ny)) > .3
        expected = [diagonal[b, rows][regions[b]].mean() for b in range(2)]
        np.testing.assert_allclose(normal_diagonal(sens, rows, regions, float(mask.mean())), expected, rtol=1e-5)

    def test_region_rms_and_margin_leakage(self):
        """RMS uses only the region; leakage is small for native-row signal and large otherwise."""
        images = np.zeros((2, 4, 4), np.complex64)
        images[0, :2] = 3j
        images[1] = 2
        regions = np.zeros((2, 4, 4), bool)
        regions[0, :2] = regions[1, 0] = True
        self.assertEqual(region_rms(images, regions), [3, 2])
        rng = np.random.default_rng(4)
        nx_os, ny = 60, 32
        noise = colored_noise(np.eye(4), np.ones(nx_os), (4, nx_os, ny), seed=5)
        whitener, lines, margin = np.eye(4), np.arange(ny), far_margin_rows(nx_os, .2)
        inside = np.zeros((4, nx_os, ny), np.complex64)
        inside[:, 20:40] = 50 * rng.standard_normal((4, 20, ny))
        everywhere = 50 * rng.standard_normal((4, nx_os, ny)).astype(np.complex64)
        for image, small in ((inside, True), (everywhere, False)):
            prediction = fftc_dim(torch.as_tensor(image), dim=1)
            leakage = margin_leakage(noise + prediction, prediction, whitener, lines, margin)
            self.assertEqual(leakage < 1e-6, small)


if __name__ == "__main__":
    unittest.main()
