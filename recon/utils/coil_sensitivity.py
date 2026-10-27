"""Slice-wise 2D ESPIRiT coil-sensitivity calibration."""

import gc

import numpy as np
import sigpy as sp
import sigpy.mri as mr
import torch

from utils.coil_compression_kspace import (
    apply_cc_coilfirst_np_2d,
    crop_readout_oversampling_kspace,
)


def estimate_espirit_csm_2d_slicewise(
    kspace_acs,
    Wcc,
    encoded_matrix,
    acs_start_y,
    readout_os_factor,
    sigpy_device,
    calib_width=24,
    crop=0.98,
    return_diagnostics=False,
):
    """Estimate matched 2D ESPIRiT maps independently for every slice.

    Parameters
    ----------
    kspace_acs : torch.Tensor
        Compact central ACS shaped ``(Nx_os, Nky, Nslices, Ncoil)``.
    Wcc : np.ndarray
        One global compression basis shaped ``(Ncoil, Ncc)``.
    return_diagnostics : bool
        Also return largest eigenvalues and unmasked calibration RSS images.
        This opt-in does not change calibration or the default array return.

    Returns
    -------
    np.ndarray or tuple[np.ndarray, dict]
        Unit-RSS maps shaped ``(Ncc, Nx, Ny, Nslices)``. The compact ACS is
        embedded at its physical ky indices before calibration, and no
        operation mixes the physical-slice dimension.
        When requested, the second return contains ``eigenvalues`` and
        ``calibration_rss``, each shaped ``(Nx, Ny, Nslices)``. RSS is derived
        from the acquired ACS, not from the cropped CSM support.
    """
    if not isinstance(kspace_acs, torch.Tensor) or kspace_acs.ndim != 4:
        raise ValueError("kspace_acs must be a 4D coil-last torch tensor.")

    encoded_nx, encoded_ny = (int(value) for value in encoded_matrix)
    acs_start_y = int(acs_start_y)
    nx_os, nky, nslices, ncoil = kspace_acs.shape
    Wcc = np.asarray(Wcc, dtype=np.complex64)
    if nx_os != encoded_nx * readout_os_factor:
        raise ValueError(
            f"Nx_os={nx_os} does not equal encoded_nx*readout_os_factor="
            f"{encoded_nx * readout_os_factor}."
        )
    if Wcc.ndim != 2 or Wcc.shape[0] != ncoil:
        raise ValueError(f"Wcc must have shape ({ncoil}, Ncc), got {Wcc.shape}.")
    if not 0 <= acs_start_y <= encoded_ny - nky:
        raise ValueError("ACS ky block lies outside the encoded phase matrix.")

    ncc = Wcc.shape[1]
    csm = np.empty((ncc, encoded_nx, encoded_ny, nslices), dtype=np.complex64)
    if return_diagnostics:
        diagnostics = {key: np.empty((encoded_nx, encoded_ny, nslices), dtype=np.float32)
                       for key in ("eigenvalues", "calibration_rss")}
    for slice_index in range(nslices):
        print(f"Estimating matched 2D ESPIRiT CSM {slice_index + 1}/{nslices}")
        kspace_coilfirst = (
            crop_readout_oversampling_kspace(
                kspace_acs[:, :, slice_index, :].permute(2, 0, 1).contiguous(),
                readout_os_factor=readout_os_factor,
                axis=1,
            )
            .cpu()
            .numpy()
            .astype(np.complex64, copy=False)
        )
        if kspace_coilfirst.shape != (ncoil, encoded_nx, nky):
            raise RuntimeError(f"Unexpected decimated ACS shape {kspace_coilfirst.shape}.")

        kspace_cc = apply_cc_coilfirst_np_2d(kspace_coilfirst, Wcc)
        kspace_full = np.zeros((ncc, encoded_nx, encoded_ny), dtype=np.complex64)
        kspace_full[:, :, acs_start_y : acs_start_y + nky] = kspace_cc

        kspace_device = sp.to_device(kspace_full, sigpy_device)
        result = mr.app.EspiritCalib(
            kspace_device,
            calib_width=min(calib_width, nky),
            device=sigpy_device,
            crop=crop,
            output_eigenvalue=return_diagnostics,
            show_pbar=True,
        ).run()
        if return_diagnostics:
            csm_device, eigenvalues_device = result
            diagnostics["eigenvalues"][:, :, slice_index] = np.asarray(
                sp.to_device(eigenvalues_device, sp.Device(-1))).reshape(encoded_nx, encoded_ny)
            diagnostics["calibration_rss"][:, :, slice_index] = np.sqrt(
                np.sum(abs(sp.ifft(kspace_full, axes=(1, 2)))**2, axis=0))
            del eigenvalues_device
        else:
            csm_device = result
        csm[:, :, :, slice_index] = sp.to_device(
            csm_device,
            sp.Device(-1),
        ).astype(np.complex64, copy=False)
        _normalize_csm_inplace(csm[:, :, :, slice_index])

        del kspace_coilfirst, kspace_cc, kspace_full, kspace_device, csm_device, result
        # Slice temporaries are acyclic and release through reference counts;
        # retain periodic collection for any backend-created reference cycles.
        if (slice_index + 1) % 10 == 0:
            gc.collect()

    return (csm, diagnostics) if return_diagnostics else csm


def _normalize_csm_inplace(csm, eps=1e-8):
    """Normalize finite CSM support and zero numerically empty voxels.

    Parameters
    ----------
    csm : np.ndarray
        Coil-first complex sensitivity maps modified in place.
    eps : float
        RSS values below this threshold are treated as empty support.

    Returns
    -------
    None
        The input array is normalized in place.
    """

    rss = np.sqrt(np.sum(np.abs(csm) ** 2, axis=0, keepdims=True))
    support = rss >= eps
    csm /= np.where(support, rss, 1.0)
    csm *= support
