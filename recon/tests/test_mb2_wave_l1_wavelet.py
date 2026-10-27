"""Test the L1-wavelet FISTA solver of utils/mb2_wave_l1_wavelet.py.

The cases are copied unchanged from the research test
test_terrax_240um_recon_side_tests.py of pulseq_basic (commit 6c01ffa); its
tests of dataset-specific side-test helpers are not part of this repository.
"""

import unittest

import numpy as np
import torch

from utils.mb2_cg_sense_wave import cg_sense_mb2_wave_2d, mb2_wave_adjoint, mb2_wave_forward
from utils.mb2_wave_l1_wavelet import (
    complex_soft_threshold, fista_l1_wavelet_mb2_wave, normal_operator_norm, wavelet_analysis,
    wavelet_l1, wavelet_prox, wavelet_synthesis,
)


def small_problem(nx_os=48, ny=32, rows=slice(16, 32), ncoil=4, seed=0):
    """Return a small MB2 Wave problem (sens, phase, psf, mask, truth) with sensitivity only on rows."""
    rng = np.random.default_rng(seed)
    x, y = np.mgrid[:nx_os, :ny]
    sens = np.zeros((2, ncoil, nx_os, ny), np.complex64)
    for b in range(2):
        for c in range(ncoil):
            cx, cy = rng.uniform(0, nx_os), rng.uniform(0, ny)
            sens[b, c] = np.exp(-((x - cx)**2 + (y - cy)**2) / 400) * np.exp(1j * rng.uniform(0, 2 * np.pi))
    sens[:, :, :rows.start] = 0
    sens[:, :, rows.stop:] = 0
    phase = np.exp(1j * rng.uniform(0, 2 * np.pi, (2, ny))).astype(np.complex64)
    psf = np.exp(1j * rng.uniform(-1, 1, (2, nx_os, ny))).astype(np.complex64)
    mask = np.zeros((1, 1, ny), np.float32)
    mask[..., ::2] = 1
    mask[..., ny // 2 - 2:ny // 2 + 2] = 1
    truth = np.zeros((2, nx_os, ny), np.complex64)
    truth[:, rows] = rng.standard_normal((2, rows.stop - rows.start, ny)) + 1j * rng.standard_normal((2, rows.stop - rows.start, ny))
    return tuple(torch.as_tensor(v) for v in (sens, phase, psf, mask, truth))


class WaveletTests(unittest.TestCase):
    """Verify orthonormal band-separate transforms and complex shrinkage."""

    def test_perfect_reconstruction_norm_and_band_separation(self):
        """Zero threshold reproduces rows exactly; a large band-0 threshold leaves band 1 unchanged."""
        rng = np.random.default_rng(1)
        image = (rng.standard_normal((64, 48)) + 1j * rng.standard_normal((64, 48))).astype(np.complex64)
        array, slices = wavelet_analysis(image)
        self.assertAlmostEqual(np.linalg.norm(array) / np.linalg.norm(image), 1, places=5)
        np.testing.assert_allclose(wavelet_synthesis(array, slices), image, atol=1e-5)
        x = torch.zeros((2, 96, 48), dtype=torch.complex64)
        x[:, 16:80] = torch.as_tensor(np.stack([image, 2 * image]))
        x[:, :16] = 1
        same = wavelet_prox(x, [0, 0], slice(16, 80))
        np.testing.assert_allclose(same[:, 16:80].numpy(), x[:, 16:80].numpy(), atol=1e-5)
        self.assertTrue(torch.all(same[:, :16] == 0))
        shrunk = wavelet_prox(x, [1e9, 0], slice(16, 80))
        np.testing.assert_allclose(shrunk[1, 16:80].numpy(), x[1, 16:80].numpy(), atol=1e-5)
        detail = wavelet_analysis(shrunk[0, 16:80].numpy())[0]
        detail[slices[0]] = 0
        self.assertLess(np.abs(detail).max(), 1e-5)
        self.assertGreater(wavelet_l1(x, slice(16, 80))[1], wavelet_l1(x, slice(16, 80))[0])

    def test_complex_soft_threshold(self):
        """Magnitudes shrink by tau with phase kept; protected coefficients and invalid tau are handled."""
        array = np.array([[3 + 4j, .5, 0], [1j, -2, 6]], np.complex64)
        result = complex_soft_threshold(array, 1., (slice(1, 2), slice(2, 3)))
        np.testing.assert_allclose(result, [[(3 + 4j) * .8, 0, 0], [0, -1, 6]], atol=1e-6)
        for tau in (-1, np.nan):
            with self.assertRaises(ValueError):
                complex_soft_threshold(array, tau, (slice(0, 0), slice(0, 0)))


class FistaTests(unittest.TestCase):
    """Verify FISTA against CG, objective descent and the power-iteration norm."""

    def test_lambda0_matches_cg(self):
        """With zero penalty FISTA converges to the CG least-squares solution on the sensitivity rows."""
        sens, phase, psf, mask, truth = small_problem()
        y = mb2_wave_forward(truth, sens, phase, psf, mask)
        cg = cg_sense_mb2_wave_2d(y, sens, phase, psf, mask, n_iter=400, tol=1e-8, verbose=False)
        fista, info = fista_l1_wavelet_mb2_wave(y, sens, phase, psf, mask, [0, 0], slice(16, 32), n_iter=4000, tol=1e-9)
        error = torch.linalg.vector_norm(fista - cg) / torch.linalg.vector_norm(cg)
        self.assertLess(float(error), 2e-3)
        self.assertTrue(torch.all(fista[:, :16] == 0) and torch.all(fista[:, 32:] == 0))
        self.assertGreater(info["lipschitz"], 0)

    def test_penalty_lowers_regularized_objective(self):
        """The FISTA objective is below that of the least-squares start on noisy data."""
        sens, phase, psf, mask, truth = small_problem(seed=3)
        generator = torch.Generator().manual_seed(4)
        noise = torch.complex(torch.randn((4, 48, 32), generator=generator), torch.randn((4, 48, 32), generator=generator))
        y = (mb2_wave_forward(truth, sens, phase, psf, mask) + .3 * noise.to(torch.complex64)) * mask
        start = cg_sense_mb2_wave_2d(y, sens, phase, psf, mask, n_iter=400, tol=1e-8, verbose=False)
        lam = [.05, .05]
        x, info = fista_l1_wavelet_mb2_wave(y, sens, phase, psf, mask, lam, slice(16, 32), x0=start, n_iter=500, tol=1e-7)

        def objective(v):
            """Return the regularized objective for one band pair."""
            data = .5 * float(torch.sum(abs(mb2_wave_forward(v, sens, phase, psf, mask) - y)**2))
            return data + sum(l * n for l, n in zip(lam, wavelet_l1(v, slice(16, 32))))

        self.assertLess(objective(x), objective(start))
        self.assertAlmostEqual(info["data_term"] + sum(info["penalty_terms"]), objective(x), delta=1e-3 * objective(x))
        with self.assertRaises(ValueError):
            fista_l1_wavelet_mb2_wave(y, sens, phase, psf, mask, [-1, 0], slice(16, 32))

    def test_power_iteration_bounds_dense_normal_matrix(self):
        """The Rayleigh estimate approaches the largest eigenvalue of A^H A from below.

        FISTA relies on 1.05 x the default 30-iteration estimate exceeding the true
        largest eigenvalue, so that the gradient step 1/L is safe.
        """
        sens, phase, psf, mask, _ = small_problem(nx_os=24, ny=16, rows=slice(8, 16), ncoil=3, seed=5)
        basis = torch.eye(2 * 24 * 16, dtype=torch.complex64).reshape(-1, 2, 24, 16)
        columns = [mb2_wave_adjoint(mb2_wave_forward(e, sens, phase, psf, mask), sens, phase, psf, mask).reshape(-1)
                   for e in basis]
        dense = torch.stack(columns, dim=1).numpy().astype(np.complex128)
        largest = np.linalg.eigvalsh((dense + dense.conj().T) / 2).max()
        converged = normal_operator_norm(sens, phase, psf, mask, n_iter=300)
        default = normal_operator_norm(sens, phase, psf, mask)
        self.assertLessEqual(converged, largest * (1 + 1e-5))
        self.assertGreater(converged / largest, .99)
        self.assertGreaterEqual(1.05 * default, largest)


if __name__ == "__main__":
    unittest.main()
