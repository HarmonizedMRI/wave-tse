"""L1-wavelet regularized MB2 Wave reconstruction with FISTA.

For one SMS group the solver minimizes

``0.5 * ||mb2_wave_forward(x) - y||^2 + sum_b lam[b] * ||W(x[b, rows])||_1``

with the unchanged operator of :mod:`utils.mb2_cg_sense_wave`. ``W`` is an
orthonormal periodized 2D DWT applied to each band image separately on the
non-oversampled readout rows; it never mixes the two physical slices, and the
coarsest approximation band is not penalized. Readout rows outside ``rows``
carry no sensitivity and stay exactly zero.
"""

import numpy as np
import pywt
import torch

from utils.mb2_cg_sense_wave import mb2_wave_forward, mb2_wave_adjoint


def wavelet_analysis(image, wavelet="db4", levels=3):
    """Return the packed orthonormal periodized 2D DWT of one complex image.

    Parameters
    ----------
    image : numpy.ndarray
        Complex 2D image ``(Nx, Ny)``; both sizes divisible by ``2**levels``.
    wavelet : str
        PyWavelets orthogonal wavelet name.
    levels : int
        Decomposition depth.

    Returns
    -------
    tuple
        Packed coefficient array with the image shape, and the PyWavelets
        slice list required by :func:`wavelet_synthesis`. ``slices[0]``
        indexes the unpenalized coarsest approximation band.

    Raises
    ------
    ValueError
        If the image is not 2D or its sizes are not divisible by ``2**levels``.
    """
    if image.ndim != 2 or any(n % 2**levels for n in image.shape):
        raise ValueError(f"Image shape {image.shape} is incompatible with {levels} periodized levels")
    coefficients = pywt.wavedec2(image, wavelet, mode="periodization", level=levels)
    return pywt.coeffs_to_array(coefficients)


def wavelet_synthesis(array, slices, wavelet="db4"):
    """Invert :func:`wavelet_analysis`.

    Parameters
    ----------
    array : numpy.ndarray
        Packed complex coefficient array.
    slices : list
        Slice list returned by :func:`wavelet_analysis`.
    wavelet : str
        Wavelet used for the analysis.

    Returns
    -------
    numpy.ndarray
        Complex 2D image with the packed array shape.
    """
    coefficients = pywt.array_to_coeffs(array, slices, output_format="wavedec2")
    return pywt.waverec2(coefficients, wavelet, mode="periodization")


def complex_soft_threshold(array, tau, protected):
    """Shrink complex magnitudes by ``tau`` and keep one protected block unchanged.

    Parameters
    ----------
    array : numpy.ndarray
        Complex coefficients.
    tau : float
        Nonnegative threshold in coefficient units.
    protected : tuple
        Index (slice tuple) of coefficients returned unchanged.

    Returns
    -------
    numpy.ndarray
        New array ``max(|a| - tau, 0) * a / |a|`` outside ``protected``.

    Raises
    ------
    ValueError
        If ``tau`` is negative or nonfinite.
    """
    if not np.isfinite(tau) or tau < 0:
        raise ValueError("Threshold must be finite and nonnegative")
    scale = np.maximum(1 - tau / np.maximum(np.abs(array), 1e-30), 0)
    result = array * scale.astype(array.real.dtype)
    result[protected] = array[protected]
    return result


def wavelet_prox(x, thresholds, rows, wavelet="db4", levels=3):
    """Apply per-band wavelet soft thresholding on the native readout rows.

    Parameters
    ----------
    x : torch.Tensor
        Band images ``(2, Nx_os, Ny)``.
    thresholds : sequence of float
        One coefficient-domain threshold per band.
    rows : slice
        Native (non-oversampled) readout rows that carry sensitivity.

    Returns
    -------
    torch.Tensor
        Thresholded images on ``x.device``; rows outside ``rows`` are zero.
    """
    if len(thresholds) != x.shape[0]:
        raise ValueError("Need one threshold per band")
    result = torch.zeros_like(x)
    block = x[:, rows].detach().cpu().numpy()
    for band in range(x.shape[0]):
        array, slices = wavelet_analysis(block[band], wavelet, levels)
        array = complex_soft_threshold(array, float(thresholds[band]), slices[0])
        image = wavelet_synthesis(array, slices, wavelet).astype(np.complex64)
        result[band, rows] = torch.as_tensor(image, device=x.device)
    return result


def wavelet_l1(x, rows, wavelet="db4", levels=3):
    """Return the penalized (non-approximation) wavelet L1 norm of each band, as floats."""
    norms = []
    block = x[:, rows].detach().cpu().numpy()
    for band in range(x.shape[0]):
        array, slices = wavelet_analysis(block[band], wavelet, levels)
        penalized = np.ones(array.shape, bool)
        penalized[slices[0]] = False
        norms.append(float(np.abs(array[penalized]).astype(np.float64).sum()))
    return norms


def normal_operator_norm(sens, caipi_phase, wave_psf, mask, n_iter=30, seed=0):
    """Estimate the largest eigenvalue of the MB2 Wave normal operator by power iteration.

    Parameters
    ----------
    sens, caipi_phase, wave_psf, mask : torch.Tensor
        Encoding terms of :func:`utils.mb2_cg_sense_wave.mb2_wave_forward`.
    n_iter : int
        Power iterations.
    seed : int
        CPU generator seed for the deterministic start vector.

    Returns
    -------
    float
        Final Rayleigh quotient, a lower bound that converges to the norm.
    """
    generator = torch.Generator(device="cpu").manual_seed(seed)
    shape = (2, *sens.shape[2:])
    v = torch.complex(torch.randn(shape, generator=generator), torch.randn(shape, generator=generator)).to(sens.device)
    v = v / torch.linalg.vector_norm(v)
    value = 0.
    for _ in range(int(n_iter)):
        w = mb2_wave_adjoint(mb2_wave_forward(v, sens, caipi_phase, wave_psf, mask), sens, caipi_phase, wave_psf, mask)
        value = float(torch.vdot(v.reshape(-1), w.reshape(-1)).real)
        v = w / torch.linalg.vector_norm(w)
    return value


def fista_l1_wavelet_mb2_wave(y, sens, caipi_phase, wave_psf, mask, lam, rows, x0=None,
                              n_iter=300, tol=1e-4, lipschitz=None, wavelet="db4", levels=3):
    """Solve the L1-wavelet regularized MB2 Wave problem with restarted FISTA.

    Parameters
    ----------
    y : torch.Tensor
        Masked measured k-space ``(Ncoil, Nx_os, Ny)``; for fractional masks
        pass ``mask * y`` so that the data term is weighted least squares.
    sens, caipi_phase, wave_psf, mask : torch.Tensor
        Encoding terms of :func:`utils.mb2_cg_sense_wave.mb2_wave_forward`.
    lam : sequence of float
        Nonnegative L1 weight per band.
    rows : slice
        Native readout rows on which the wavelet penalty acts.
    x0 : torch.Tensor, optional
        Start image ``(2, Nx_os, Ny)``; zeros when omitted.
    n_iter : int
        Maximum iterations.
    tol : float
        Stop when ``||x_k+1 - x_k|| / ||x_k+1||`` falls below this value.
    lipschitz : float, optional
        Gradient Lipschitz constant; 1.05 x power-iteration estimate when omitted.

    Returns
    -------
    tuple
        ``(x, info)``: the solution and a dictionary with iterations,
        convergence flag, Lipschitz constant, restart count, relative-change
        history and final data/penalty objective terms.

    Raises
    ------
    ValueError
        If ``lam`` does not hold two finite nonnegative values.
    """
    lam = [float(v) for v in lam]
    if len(lam) != 2 or not all(np.isfinite(v) and v >= 0 for v in lam):
        raise ValueError("lam must hold two finite nonnegative band weights")
    if lipschitz is None:
        lipschitz = 1.05 * normal_operator_norm(sens, caipi_phase, wave_psf, mask)
    step = 1 / lipschitz
    thresholds = [v * step for v in lam]
    x = torch.zeros((2, *sens.shape[2:]), dtype=torch.complex64, device=sens.device) if x0 is None else x0.clone()
    z, t = x.clone(), 1.
    history, restarts, converged = [], 0, False
    for iteration in range(int(n_iter)):
        gradient = mb2_wave_adjoint(mb2_wave_forward(z, sens, caipi_phase, wave_psf, mask) - y,
                                    sens, caipi_phase, wave_psf, mask)
        x_new = wavelet_prox(z - step * gradient, thresholds, rows, wavelet, levels)
        # Gradient-based adaptive restart removes momentum oscillations.
        if float(torch.vdot((z - x_new).reshape(-1), (x_new - x).reshape(-1)).real) > 0:
            t, restarts = 1., restarts + 1
        t_new = (1 + np.sqrt(1 + 4 * t * t)) / 2
        z = x_new + ((t - 1) / t_new) * (x_new - x)
        change = float(torch.linalg.vector_norm(x_new - x) / torch.clamp(torch.linalg.vector_norm(x_new), min=1e-30))
        x, t = x_new, t_new
        history.append(change)
        if change < tol:
            converged = True
            break
    data_term = 0.5 * float(torch.sum(abs(mb2_wave_forward(x, sens, caipi_phase, wave_psf, mask) - y).double()**2))
    penalty = [v * n for v, n in zip(lam, wavelet_l1(x, rows, wavelet, levels))]
    info = {"iterations": iteration + 1, "converged": converged, "lipschitz": float(lipschitz),
            "restarts": restarts, "relative_change_history": history,
            "data_term": data_term, "penalty_terms": penalty, "lambda": lam}
    return x, info
