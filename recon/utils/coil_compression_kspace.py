"""
MRI k-space helpers for slice-wise coil compression and readout cropping.

The reconstruction keeps the coil dimension first for SigPy:

    kspace_sigpy.shape == (Ncoil, Nx, Ny)

Functions
---------
crop_readout_oversampling_kspace(kspace, readout_os_factor, axis=0)
    Remove readout oversampling by cropping the centered image-domain FOV.

apply_cc_coilfirst_np_2d(kspace_cxy, W)
    Apply one coil-compression matrix W, shaped (Ncoil, Ncc), to coil-first
    2D k-space shaped (Ncoil, Nx, Ny), returning (Ncc, Nx, Ny).

Use the same compression matrix for every calibration slice and every TSE
group so that the sensitivities and the TSE data share one coil space.
"""

import numpy as np
import torch


def crop_readout_oversampling_kspace(kspace, readout_os_factor, axis=0):
    """Remove readout oversampling by cropping in the image domain.

    Parameters
    ----------
    kspace : np.ndarray or torch.Tensor
        Centered k-space with the oversampled readout dimension on ``axis``.
    readout_os_factor : int
        Integer ratio between the acquired and retained readout FOV.
    axis : int
        Readout axis.

    Returns
    -------
    np.ndarray or torch.Tensor
        Centered k-space with the same type and dtype as the input and a
        readout length reduced by ``readout_os_factor``.

    Raises
    ------
    ValueError
        If the factor is invalid or does not divide the readout length.

    Notes
    -----
    Striding through k-space folds the enlarged readout FOV onto the retained
    FOV. Cropping the centered readout image first removes out-of-FOV support
    before returning to the encoded k-space grid.
    """
    if isinstance(readout_os_factor, bool):
        raise ValueError("readout_os_factor must be a positive integer.")
    factor = int(readout_os_factor)
    if factor < 1 or factor != readout_os_factor:
        raise ValueError("readout_os_factor must be a positive integer.")

    ndim = kspace.ndim
    axis = int(axis) % ndim
    readout_size = int(kspace.shape[axis])
    if readout_size % factor != 0:
        raise ValueError(
            f"Readout length {readout_size} is not divisible by factor {factor}."
        )
    if factor == 1:
        return kspace

    retained_size = readout_size // factor
    crop_start = (readout_size - retained_size) // 2
    crop_stop = crop_start + retained_size
    crop_slices = [slice(None)] * ndim
    crop_slices[axis] = slice(crop_start, crop_stop)

    if isinstance(kspace, torch.Tensor):
        image = torch.fft.fftshift(
            torch.fft.ifft(
                torch.fft.ifftshift(kspace, dim=(axis,)),
                dim=axis,
                norm="ortho",
            ),
            dim=(axis,),
        )
        image = image[tuple(crop_slices)]
        return torch.fft.fftshift(
            torch.fft.fft(
                torch.fft.ifftshift(image, dim=(axis,)),
                dim=axis,
                norm="ortho",
            ),
            dim=(axis,),
        )

    array = np.asarray(kspace)
    image = np.fft.fftshift(
        np.fft.ifft(
            np.fft.ifftshift(array, axes=(axis,)),
            axis=axis,
            norm="ortho",
        ),
        axes=(axis,),
    )
    image = image[tuple(crop_slices)]
    cropped = np.fft.fftshift(
        np.fft.fft(
            np.fft.ifftshift(image, axes=(axis,)),
            axis=axis,
            norm="ortho",
        ),
        axes=(axis,),
    )
    return cropped.astype(array.dtype, copy=False)


# -----------------------------------------------------------------------------
# 2D multislice helpers
# -----------------------------------------------------------------------------
# These functions treat the slice dimension as an external batch index.  They are
# intended for non-MB 2D acquisitions stored as (Nx, Ny, Nslices, Ncoil), where
# every physical slice has a 2D k-space matrix (kx, ky) and there is no FFT along
# the slice dimension.


def apply_cc_coilfirst_np_2d(kspace_cxy, W):
    """
    Apply coil compression to one 2D coil-first k-space array.

    Parameters
    ----------
    kspace_cxy : np.ndarray
        Shape (Ncoil, Nx, Ny).
    W : np.ndarray
        Shape (Ncoil, Ncc).

    Returns
    -------
    np.ndarray
        Shape (Ncc, Nx, Ny), complex64.
    """
    arr = np.asarray(kspace_cxy)
    if arr.ndim != 3:
        raise ValueError(f"kspace_cxy must have shape (Ncoil, Nx, Ny), got {arr.shape}")
    return np.einsum("cxy,cn->nxy", arr, W, optimize=True).astype(np.complex64, copy=False)
