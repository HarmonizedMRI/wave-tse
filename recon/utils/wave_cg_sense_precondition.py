"""
wave_cg_sense_precondition.py

Centered orthonormal 1D FFT helpers shared by the calibration import and the
MB2 Wave operators. The module keeps the name of its Wave-MPRAGE and
Wave-GRE counterparts; their CG-SENSE classes are not used here.
"""

import torch


def fftc_dim(x, dim):
    """Centered orthonormal 1D FFT along one dimension."""
    return torch.fft.fftshift(
        torch.fft.fft(torch.fft.ifftshift(x, dim=(dim,)), dim=dim, norm='ortho'),
        dim=(dim,))


def ifftc_dim(x, dim):
    """Centered orthonormal 1D IFFT along one dimension."""
    return torch.fft.fftshift(
        torch.fft.ifft(torch.fft.ifftshift(x, dim=(dim,)), dim=dim, norm='ortho'),
        dim=(dim,))
