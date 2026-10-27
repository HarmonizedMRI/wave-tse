"""Measured TSE noise model that scales the L1-wavelet penalty.

The L1-wavelet weight of band ``b`` is ``lambda_b = c * sigma_b * d_b``, with
``sigma_b`` the RMS of the lambda-0 reconstruction of noise alone and ``d_b``
the mean diagonal of the normal operator, both inside the same region. The
noise model follows the TerraX research reconstruction of pulseq_basic
(``experiment_recon/tse_2d_recon_update``, commit 6c01ffa):

* the virtual-coil covariance ``C`` is measured on the far readout margins of
  the TSE data, which receive no signal;
* in the hybrid (readout position ``x``, PE line) domain the noise of sample
  ``x`` has the covariance ``q(x) * C``. The receiver profile ``q`` is about 1
  on average over the far margins and larger inside the FOV, because the
  receiver filter rolls off towards the edges of the sampled band.

The research reconstruction took ``q`` from a fully sampled reference scan.
The reconstruction command estimates ``q`` from the scan itself: the lambda-0
residual of the outer PE lines contains noise only, but inside the FOV the fit
absorbs part of that noise, so the residual power is divided by the residual
power that the same fit leaves of synthetic noise with a flat profile. The
functions here are the shared numerical steps; :func:`margin_whitener`,
:func:`colored_noise`, :func:`normal_diagonal` and :func:`region_rms` keep the
arithmetic of the research helpers.
"""

import numpy as np
import torch

from utils.wave_cg_sense_precondition import fftc_dim, ifftc_dim


def far_margin_rows(nx_os, fraction):
    """Return the readout samples of both far margins of the oversampled readout.

    Parameters
    ----------
    nx_os : int
        Oversampled readout size.
    fraction : float
        Width of each margin as a fraction of ``nx_os``.

    Returns
    -------
    numpy.ndarray
        Indices of the first and of the last ``round(fraction * nx_os)``
        readout samples.

    Raises
    ------
    ValueError
        If a margin would be empty or the two margins would meet.
    """
    width = int(round(fraction * nx_os))
    if width < 1 or 2 * width >= nx_os:
        raise ValueError(f"Margin fraction {fraction} does not fit a {nx_os}-sample readout")
    return np.r_[0:width, nx_os - width:nx_os]


def outer_lines(lines, ny, fraction):
    """Return the sampled PE lines that lie in the outer part of the PE grid.

    Parameters
    ----------
    lines : numpy.ndarray
        Sampled zero-based PE lines.
    ny : int
        PE matrix size; line ``ny // 2`` is the k-space center.
    fraction : float
        Outer fraction of the PE grid; 0.3 keeps the lines with
        ``|line - ny // 2| > 0.35 * ny``.

    Returns
    -------
    numpy.ndarray
        The selected subset of ``lines``.
    """
    lines = np.asarray(lines)
    return lines[abs(lines - ny // 2) > (1 - fraction) * ny / 2]


def margin_covariance(data, lines, margin):
    """Return the virtual-coil covariance on the far readout margins.

    Parameters
    ----------
    data : torch.Tensor
        Measured k-space ``(coil, Nx_os, Ny)``.
    lines : numpy.ndarray
        Sampled PE lines.
    margin : numpy.ndarray
        Readout samples from :func:`far_margin_rows`.

    Returns
    -------
    numpy.ndarray
        complex128 ``(coil, coil)``: the mean of ``h h^H`` over the margin
        samples of the sampled lines, where ``h`` is the readout-IFFT
        (hybrid) data.
    """
    selected = data[:, :, torch.as_tensor(lines, device=data.device)]
    hybrid = ifftc_dim(selected, dim=1)[:, torch.as_tensor(margin, device=data.device)]
    samples = hybrid.reshape(hybrid.shape[0], -1).to(torch.complex128)
    return (samples @ samples.conj().T / samples.shape[1]).cpu().numpy()


def margin_whitener(cov):
    """Return ``W`` with ``W cov W^H = mean(diag(cov)) I``.

    ``W`` is the scaled inverse Cholesky factor. The mean noise level is kept,
    so that solver tolerances and residual scales stay comparable with
    unwhitened data.

    Parameters
    ----------
    cov : numpy.ndarray
        Hermitian positive-definite covariance ``(coil, coil)``.

    Returns
    -------
    numpy.ndarray
        complex128 whitener ``(coil, coil)``, applied as ``W @ data``.
    """
    return np.linalg.inv(np.linalg.cholesky(cov)) * np.sqrt(np.real(np.diag(cov)).mean())


def colored_noise(cov, relative_variance, shape, seed):
    """Return reproducible k-space noise that follows the measured TSE noise model.

    In the hybrid domain, the noise of readout sample ``x`` has the covariance
    ``relative_variance[x] * cov``; it is white along PE.

    Parameters
    ----------
    cov : numpy.ndarray
        Virtual-coil covariance ``(coil, coil)``.
    relative_variance : numpy.ndarray
        Nonnegative readout profile ``(Nx_os,)``.
    shape : tuple of int
        Output shape ``(coil, Nx_os, Ny)``.
    seed : int
        Seed of the CPU random generator.

    Returns
    -------
    torch.Tensor
        complex64 noise ``(coil, Nx_os, Ny)`` on the CPU.
    """
    generator = torch.Generator(device="cpu").manual_seed(int(seed))
    white = torch.complex(torch.randn(shape, generator=generator),
                          torch.randn(shape, generator=generator)) / np.sqrt(2)
    chol = torch.as_tensor(np.linalg.cholesky(cov).astype(np.complex64))
    scale = torch.as_tensor(np.sqrt(relative_variance).astype(np.float32))
    shaped = fftc_dim(white * scale[None, :, None], dim=1)
    return torch.einsum("cd,dkl->ckl", chol, shaped)


def whitened_hybrid_power(kspace, whitener, lines):
    """Return the mean whitened hybrid power of every readout sample.

    Parameters
    ----------
    kspace : torch.Tensor
        k-space ``(coil, Nx_os, Ny)``, for example data or a residual.
    whitener : numpy.ndarray
        Coil whitener ``(coil, coil)`` from :func:`margin_whitener`.
    lines : numpy.ndarray
        PE lines to average over.

    Returns
    -------
    numpy.ndarray
        float64 ``(Nx_os,)``: the mean over coils and ``lines`` of
        ``|IFFT_x(W k)|^2``.
    """
    w = torch.as_tensor(whitener.astype(np.complex64), device=kspace.device)
    selected = kspace[:, :, torch.as_tensor(lines, device=kspace.device)]
    hybrid = ifftc_dim(torch.einsum("cd,dkl->ckl", w, selected), dim=1)
    return (abs(hybrid)**2).double().mean(dim=(0, 2)).cpu().numpy()


def margin_leakage(data, prediction, whitener, lines, margin):
    """Return the predicted signal power in the far readout margins relative to the data power there.

    The lambda-0 prediction reaches the margins only where the Wave PSF
    spreads the FOV signal into them; a small value confirms that the margin
    covariance measures noise.

    Parameters
    ----------
    data : torch.Tensor
        Measured k-space ``(coil, Nx_os, Ny)``.
    prediction : torch.Tensor
        k-space predicted by the lambda-0 solution, same shape.
    whitener : numpy.ndarray
        Coil whitener ``(coil, coil)``.
    lines : numpy.ndarray
        Sampled PE lines.
    margin : numpy.ndarray
        Readout samples of the far margins.

    Returns
    -------
    float
        Mean predicted power divided by mean data power over the margins.
    """
    predicted = whitened_hybrid_power(prediction, whitener, lines)
    measured = whitened_hybrid_power(data, whitener, lines)
    return float(predicted[margin].mean() / measured[margin].mean())


def normal_diagonal(sens, rows, regions, line_fraction):
    """Return the mean diagonal of the normal operator of each band inside its region.

    With unit-magnitude CAIPI phase and Wave PSF and orthonormal FFTs, the
    diagonal of ``A^H A`` is ``line_fraction * sum_c |S_c|^2``.

    Parameters
    ----------
    sens : torch.Tensor
        Sensitivities ``(band, coil, Nx_os, Ny)``.
    rows : slice
        Native readout rows.
    regions : numpy.ndarray
        Boolean regions ``(band, Nx, Ny)``.
    line_fraction : float
        Fraction of the PE lines that were sampled.

    Returns
    -------
    list of float
        One mean diagonal per band.
    """
    power = (abs(sens[:, :, rows])**2).sum(dim=1).cpu().numpy()
    return [float(line_fraction * power[b][regions[b]].mean()) for b in range(power.shape[0])]


def region_rms(images, regions):
    """Return the complex RMS of each band image inside its region.

    Parameters
    ----------
    images : numpy.ndarray
        Complex band images ``(band, Nx, Ny)``.
    regions : numpy.ndarray
        Boolean regions of the same shape.

    Returns
    -------
    list of float
        One RMS per band.
    """
    return [float(np.sqrt(np.mean(abs(images[b][regions[b]])**2))) for b in range(len(images))]
