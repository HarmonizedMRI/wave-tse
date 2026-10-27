"""ESPIRiT CSM calibration and geometry-aware sampling for 2D slices."""

import gc

import numpy as np
import sigpy as sp
import sigpy.mri as mr
import torch
from scipy.ndimage import map_coordinates

from utils.coil_compression_kspace import (
    apply_cc_coilfirst_np,
    apply_cc_coilfirst_np_2d,
    crop_readout_oversampling_kspace,
    estimate_cc_matrix_coillast,
)


def estimate_espirit_csm_3d(
    kspace_acs,
    encoded_matrix,
    acs_start_yz,
    readout_os_factor,
    sigpy_device,
    ncc=24,
    calib_width=24,
    crop=0.98,
):
    """Estimate a global compression basis and 3D ESPIRiT maps.

    Parameters
    ----------
    kspace_acs : torch.Tensor
        Compact central ACS in coil-last layout ``(Nx_os, Nky, Nkz, Ncoil)``.

    Returns
    -------
    Wcc : np.ndarray
        Global coil-compression matrix shaped ``(Ncoil, Ncc)``.
    cc_energy : np.ndarray
        Cumulative retained coil-compression energy.
    csm_3d : np.ndarray
        Compressed-coil ESPIRiT maps on the full encoded matrix.
    """
    if not isinstance(kspace_acs, torch.Tensor) or kspace_acs.ndim != 4:
        raise ValueError("kspace_acs must be a 4D coil-last torch tensor.")

    encoded_nx, encoded_ny, encoded_nz = (int(value) for value in encoded_matrix)
    acs_start_y, acs_start_z = (int(value) for value in acs_start_yz)
    nx_os, nky, nkz, ncoil = kspace_acs.shape
    if nx_os != encoded_nx * readout_os_factor:
        raise ValueError(
            f"Nx_os={nx_os} does not equal encoded_nx*readout_os_factor="
            f"{encoded_nx * readout_os_factor}."
        )
    if ncc > ncoil:
        raise ValueError(f"ncc={ncc} exceeds the physical coil count {ncoil}.")

    Wcc, _cc_svals, cc_energy = estimate_cc_matrix_coillast(
        kspace_acs,
        ncc=ncc,
        acs=min(nky, nkz),
        x_step=readout_os_factor,
    )

    # Remove enlarged readout FOV support before returning to encoded k-space.
    kspace_coilfirst = (
        crop_readout_oversampling_kspace(
            kspace_acs.permute(3, 0, 1, 2).contiguous(),
            readout_os_factor=readout_os_factor,
            axis=1,
        )
        .cpu()
        .numpy()
        .astype(np.complex64, copy=False)
    )
    if kspace_coilfirst.shape != (ncoil, encoded_nx, nky, nkz):
        raise RuntimeError(f"Unexpected decimated ACS shape {kspace_coilfirst.shape}.")

    kspace_cc_compact = apply_cc_coilfirst_np(kspace_coilfirst, Wcc)
    del kspace_coilfirst
    gc.collect()

    if not 0 <= acs_start_y <= encoded_ny - nky:
        raise ValueError("ACS ky block lies outside the encoded phase matrix.")
    if not 0 <= acs_start_z <= encoded_nz - nkz:
        raise ValueError("ACS kz block lies outside the encoded partition matrix.")

    # REFSCAN/SET contains a compact central ky-kz block. Embed it at its
    # physical encoded indices before ESPIRiT so the maps retain the full FOV.
    # Treating the compact 32x32 block as a complete grid would shrink its FOV
    # and then incorrectly stretch the resulting sensitivities.
    kspace_cc_full = np.zeros(
        (ncc, encoded_nx, encoded_ny, encoded_nz),
        dtype=np.complex64,
    )
    kspace_cc_full[
        :,
        :,
        acs_start_y : acs_start_y + nky,
        acs_start_z : acs_start_z + nkz,
    ] = kspace_cc_compact
    del kspace_cc_compact
    gc.collect()

    kspace_cc_device = sp.to_device(kspace_cc_full, sigpy_device)
    del kspace_cc_full
    gc.collect()
    csm_device = mr.app.EspiritCalib(
        kspace_cc_device,
        calib_width=min(calib_width, nky, nkz),
        device=sigpy_device,
        crop=crop,
        show_pbar=True,
    ).run()
    csm_3d = sp.to_device(csm_device, sp.Device(-1)).astype(np.complex64, copy=False)

    _normalize_csm_inplace(csm_3d)
    return Wcc.astype(np.complex64, copy=False), cc_energy, csm_3d


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


def resample_csm_to_2d_slices(
    csm_3d,
    source_fov_mm,
    target_matrix_xy,
    target_fov_xy_mm,
    target_slice_offsets_mm,
    coordinate_tolerance_mm=1e-5,
):
    """Sample a centered 3D CSM grid onto physical 2D slice locations.

    Source and target orientation and FOV center must already be validated.
    Interpolation is skipped when x, y, and all requested z coordinates match
    the source grid exactly.

    Returns
    -------
    csm_target : np.ndarray
        Array shaped ``(Ncc, target_nx, target_ny, Nslice)``.
    interpolation_applied : bool
        Whether coordinate interpolation was required.
    """
    csm_3d = np.asarray(csm_3d, dtype=np.complex64)
    if csm_3d.ndim != 4:
        raise ValueError(f"csm_3d must be (Ncc, Nx, Ny, Nz), got {csm_3d.shape}.")

    source_fov_mm = np.asarray(source_fov_mm, dtype=np.float64)
    target_fov_xy_mm = np.asarray(target_fov_xy_mm, dtype=np.float64)
    target_slice_offsets_mm = np.asarray(target_slice_offsets_mm, dtype=np.float64).reshape(-1)
    target_nx, target_ny = (int(target_matrix_xy[0]), int(target_matrix_xy[1]))
    if source_fov_mm.shape != (3,) or target_fov_xy_mm.shape != (2,):
        raise ValueError("source_fov_mm and target_fov_xy_mm must have lengths 3 and 2.")

    _, source_nx, source_ny, source_nz = csm_3d.shape
    source_x = _centered_coordinates(source_nx, source_fov_mm[0])
    source_y = _centered_coordinates(source_ny, source_fov_mm[1])
    source_z = _centered_coordinates(source_nz, source_fov_mm[2])
    target_x = _centered_coordinates(target_nx, target_fov_xy_mm[0])
    target_y = _centered_coordinates(target_ny, target_fov_xy_mm[1])

    z_indices = _matching_coordinate_indices(
        source_z,
        target_slice_offsets_mm,
        atol=coordinate_tolerance_mm,
    )
    exact_grid = (
        source_nx == target_nx
        and source_ny == target_ny
        and np.allclose(source_x, target_x, rtol=0.0, atol=coordinate_tolerance_mm)
        and np.allclose(source_y, target_y, rtol=0.0, atol=coordinate_tolerance_mm)
        and z_indices is not None
    )
    if exact_grid:
        csm_target = csm_3d[:, :, :, z_indices].copy()
        _normalize_csm_inplace(csm_target)
        return csm_target, False

    target_bounds = (
        ("readout", -0.5 * target_fov_xy_mm[0], 0.5 * target_fov_xy_mm[0]),
        ("phase", -0.5 * target_fov_xy_mm[1], 0.5 * target_fov_xy_mm[1]),
        ("slice", target_slice_offsets_mm.min(), target_slice_offsets_mm.max()),
    )
    for axis, (name, target_min, target_max) in enumerate(target_bounds):
        source_min = -0.5 * source_fov_mm[axis]
        source_max = 0.5 * source_fov_mm[axis]
        if (
            target_min < source_min - coordinate_tolerance_mm
            or target_max > source_max + coordinate_tolerance_mm
        ):
            raise ValueError(f"Target {name} coordinates exceed the source CSM coverage.")

    x_index = _coordinates_to_indices(target_x, source_nx, source_fov_mm[0])
    y_index = _coordinates_to_indices(target_y, source_ny, source_fov_mm[1])
    z_index = _coordinates_to_indices(target_slice_offsets_mm, source_nz, source_fov_mm[2])
    coordinate_grid = np.stack(
        np.broadcast_arrays(
            x_index[:, None, None],
            y_index[None, :, None],
            z_index[None, None, :],
        ),
        axis=0,
    )

    csm_target = np.empty(
        (csm_3d.shape[0], target_nx, target_ny, target_slice_offsets_mm.size),
        dtype=np.complex64,
    )
    for coil in range(csm_3d.shape[0]):
        real = map_coordinates(
            csm_3d[coil].real,
            coordinate_grid,
            order=1,
            mode="nearest",
            prefilter=False,
        )
        imag = map_coordinates(
            csm_3d[coil].imag,
            coordinate_grid,
            order=1,
            mode="nearest",
            prefilter=False,
        )
        csm_target[coil] = real + 1j * imag

    _normalize_csm_inplace(csm_target)
    return csm_target, True


def _centered_coordinates(size, fov_mm):
    voxel_size = float(fov_mm) / size
    return (np.arange(size, dtype=np.float64) + 0.5 - 0.5 * size) * voxel_size


def _coordinates_to_indices(coordinates_mm, source_size, source_fov_mm):
    return (
        np.asarray(coordinates_mm, dtype=np.float64) / (source_fov_mm / source_size)
        + 0.5 * source_size
        - 0.5
    )


def _matching_coordinate_indices(source, target, atol):
    indices = []
    for coordinate in target:
        index = int(np.argmin(np.abs(source - coordinate)))
        if abs(source[index] - coordinate) > atol:
            return None
        indices.append(index)
    return np.asarray(indices, dtype=np.int64)


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
