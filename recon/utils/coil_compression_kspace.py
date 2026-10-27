"""
MRI k-space coil-compression helpers.

This module provides utilities for preparing Siemens Twix/mapVBVD k-space data
for low-resolution ESPIRiT coil sensitivity calibration and reconstruction.

Assumed data conventions
------------------------
1. Raw k-space tensors are assumed to use coil-last layout:

       kspace.shape == (Nx, Ny, Nz, Ncoil)

   where:
       Nx    : readout dimension, possibly oversampled
       Ny    : first phase-encode dimension
       Nz    : second phase-encode / partition dimension
       Ncoil : physical receive channels, e.g. 32

2. SigPy MRI routines generally expect coil-first layout:

       kspace_sigpy.shape == (Ncoil, Nx, Ny, Nz)

3. Complex data should be stored as complex64 whenever possible to reduce
   memory use.

4. For accelerated acquisitions with separate ACS/reference scans, the image
   data may contain only the undersampled Ry x Rz grid, while the fully sampled
   ACS block may live in the Twix REF/REFSCAN stream. The ACS region should be
   merged into the image k-space before ESPIRiT calibration.

Functions
---------

estimate_cc_matrix_coillast(kspace, ncc=12, acs=24, x_step=4, eps=1e-8)
    Estimate a coil-compression matrix W from a small central ACS calibration
    block using the coil covariance matrix. The input is coil-last k-space
    shaped (Nx, Ny, Nz, Ncoil). The returned W has shape (Ncoil, Ncc).

apply_cc_coillast_torch(kspace, W, x_chunk=8)
    Apply a previously estimated coil-compression matrix to a full coil-last
    PyTorch tensor in readout chunks to avoid excessive memory use. Returns
    compressed k-space shaped (Nx, Ny, Nz, Ncc).

apply_cc_coilfirst_np(kspace_cf, W)
    Apply the same compression matrix to a coil-first NumPy array shaped
    (Ncoil, Nx, Ny, Nz), returning a compressed array shaped
    (Ncc, Nx, Ny, Nz). This is useful before passing low-resolution k-space
    to SigPy ESPIRiT.

Typical workflow
----------------
1. Load and merge image + ACS/reference k-space:

       kspace_echo, mask_2d = load_merged_img_ref(...)

2. Estimate coil compression from the merged full k-space:

       Wcc, svals, energy = estimate_cc_matrix_coillast(
           kspace_echo,
           ncc=12,
           acs=24,
           x_step=4,
       )

3. Build a low-resolution coil-first k-space for ESPIRiT:

       kspace_np = (
           kspace_echo
           .permute(3, 0, 1, 2)
           .contiguous()
           .numpy()
           .astype(np.complex64, copy=False)
       )

       kspace_low_np = sp.resize(kspace_np, (Ncoil, 64, 64, 48))
       kspace_low_cc_np = apply_cc_coilfirst_np(kspace_low_np, Wcc)

4. Run ESPIRiT on the compressed low-resolution k-space:

       kspace_low_cc_sp = sp.to_device(kspace_low_cc_np, sp.Device(0))

       csm_low_cc = mr.app.EspiritCalib(
           kspace_low_cc_sp,
           calib_width=24,
           device=sp.Device(0),
           crop=0,
           show_pbar=True,
       ).run()

5. Compress full-resolution k-space, if needed:

       kspace_cc = apply_cc_coillast_torch(kspace_echo, Wcc, x_chunk=8)

Notes
-----
- Coil compression is estimated from the central ACS region but can be applied
  to the full k-space.
- Use the same compression matrix W for wave and no-wave datasets if they were
  acquired with the same coil configuration and ordering.
- ESPIRiT memory scales roughly as:

       Nx * Ny * Nz * Ncoil^2

  so reducing 32 coils to 12 coils greatly reduces GPU memory usage.
- For a 32-channel complex64 dataset shaped (32, 1024, 256, 192), the full
  coil-first array is approximately 12.9 GB. Avoid unnecessary full copies.
"""

import numpy as np
import torch
import scipy.linalg as la


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

def estimate_cc_matrix_coillast(kspace, ncc=12, acs=24, x_step=4, eps=1e-8):
    """
    Estimate SVD coil compression matrix from central ACS.

    kspace: torch or numpy, shape (Nx, Ny, Nz, ncoil)
    returns W: numpy complex64, shape (ncoil, ncc)
    """
    Nx, Ny, Nz, ncoil = kspace.shape

    cy, cz = Ny // 2, Nz // 2
    y0, y1 = cy - acs // 2, cy + acs // 2
    z0, z1 = cz - acs // 2, cz + acs // 2

    calib = crop_readout_oversampling_kspace(
        kspace[:, y0:y1, z0:z1, :],
        readout_os_factor=x_step,
        axis=0,
    )

    if isinstance(calib, torch.Tensor):
        calib = calib.detach().cpu().numpy()

    calib = np.asarray(calib, dtype=np.complex64)
    X = calib.reshape(-1, ncoil)

    # Remove zero-filled samples
    power = np.sum(np.abs(X) ** 2, axis=1)
    keep = power > eps * power.max()
    X = X[keep]

    # 32 x 32 covariance, much faster than full SVD
    C = X.conj().T @ X
    C = 0.5 * (C + C.conj().T)

    evals, evecs = la.eigh(C)
    idx = np.argsort(evals)[::-1]

    evals = evals[idx]
    evecs = evecs[:, idx]

    W = evecs[:, :ncc].astype(np.complex64, copy=False)
    svals = np.sqrt(np.maximum(evals, 0))

    energy = np.cumsum(evals) / np.sum(evals)

    return W, svals, energy


def apply_cc_coillast_torch(kspace, W, x_chunk=8):
    """
    Apply coil compression to full k-space in chunks.

    kspace: torch tensor, shape (Nx, Ny, Nz, ncoil)
    W: numpy array, shape (ncoil, ncc)
    returns: torch tensor, shape (Nx, Ny, Nz, ncc)
    """
    assert isinstance(kspace, torch.Tensor)

    Nx, Ny, Nz, ncoil = kspace.shape
    ncc = W.shape[1]

    W_t = torch.as_tensor(W, dtype=kspace.dtype, device=kspace.device)
    out = torch.empty((Nx, Ny, Nz, ncc), dtype=kspace.dtype, device=kspace.device)

    with torch.no_grad():
        for x0 in range(0, Nx, x_chunk):
            x1 = min(x0 + x_chunk, Nx)

            block = kspace[x0:x1]  # (chunk, Ny, Nz, ncoil)
            X = block.reshape(-1, ncoil)
            Xcc = X @ W_t

            out[x0:x1] = Xcc.reshape(x1 - x0, Ny, Nz, ncc)

    return out


def apply_cc_coilfirst_np(kspace_cf, W):
    """
    kspace_cf: numpy, shape (ncoil, x, y, z)
    W: shape (ncoil, ncc)
    returns: shape (ncc, x, y, z)
    """
    return np.einsum(
        "cxyz,cn->nxyz",
        kspace_cf,
        W,
        optimize=True,
    ).astype(np.complex64, copy=False)


# -----------------------------------------------------------------------------
# 2D multislice helpers
# -----------------------------------------------------------------------------
# These functions treat the slice dimension as an external batch index.  They are
# intended for non-MB 2D acquisitions stored as (Nx, Ny, Nslices, Ncoil), where
# every physical slice has a 2D k-space matrix (kx, ky) and there is no FFT along
# the slice dimension.


def estimate_cc_matrix_coillast_2d(kspace_xyc, ncc=12, acs=24, x_step=1, eps=1e-8):
    """
    Estimate an SVD coil-compression matrix from a central 2D ACS block.

    Parameters
    ----------
    kspace_xyc : torch.Tensor or np.ndarray
        Coil-last 2D k-space, shape (Nx, Ny, Ncoil). Nx may be readout
        oversampled.
    ncc : int
        Number of compressed virtual coils.
    acs : int
        Number of central ky lines to use.
    x_step : int
        Readout oversampling factor removed by an image-domain center crop.
    eps : float
        Relative power threshold for discarding zero-filled samples.

    Returns
    -------
    W : np.ndarray, complex64, shape (Ncoil, Ncc)
    svals : np.ndarray
    energy : np.ndarray
        Cumulative retained energy for each virtual coil.
    """
    if isinstance(kspace_xyc, torch.Tensor):
        arr = kspace_xyc.detach().cpu().numpy()
    else:
        arr = np.asarray(kspace_xyc)

    if arr.ndim != 3:
        raise ValueError(f"kspace_xyc must have shape (Nx, Ny, Ncoil), got {arr.shape}")

    Nx, Ny, ncoil = arr.shape
    if ncc > ncoil:
        raise ValueError(f"ncc={ncc} cannot exceed ncoil={ncoil}")

    acs_eff = min(int(acs), Ny)
    cy = Ny // 2
    y0 = max(0, cy - acs_eff // 2)
    y1 = min(Ny, y0 + acs_eff)
    y0 = max(0, y1 - acs_eff)

    x_step = max(1, int(x_step))
    calib = crop_readout_oversampling_kspace(
        arr[:, y0:y1, :],
        readout_os_factor=x_step,
        axis=0,
    )
    calib = np.asarray(calib, dtype=np.complex64)
    X = calib.reshape(-1, ncoil)

    power = np.sum(np.abs(X) ** 2, axis=1)
    if power.size == 0 or power.max() <= 0:
        raise ValueError("No nonzero samples available for 2D coil-compression calibration.")

    keep = power > eps * power.max()
    X = X[keep]
    if X.shape[0] < ncoil:
        raise ValueError(
            f"Too few nonzero calibration samples ({X.shape[0]}) for {ncoil} coils. "
            "Check ACS selection or data loading."
        )

    C = X.conj().T @ X
    C = 0.5 * (C + C.conj().T)

    evals, evecs = la.eigh(C)
    idx = np.argsort(evals)[::-1]
    evals = np.maximum(evals[idx], 0)
    evecs = evecs[:, idx]

    W = evecs[:, :ncc].astype(np.complex64, copy=False)
    svals = np.sqrt(evals)
    denom = np.sum(evals)
    energy = np.cumsum(evals) / denom if denom > 0 else np.zeros_like(evals)
    return W, svals, energy


def estimate_cc_matrix_coillast_multislice_2d(
    kspace_xyzc,
    ncc=12,
    acs=24,
    x_step=1,
    eps=1e-8,
):
    """Estimate one global compression basis from all physical 2D slices.

    The physical-slice dimension contributes calibration samples but is never
    Fourier transformed or interpolated. The returned basis can therefore be
    applied unchanged to every calibration slice and every SMS group.
    """
    if isinstance(kspace_xyzc, torch.Tensor):
        arr = kspace_xyzc.detach().cpu().numpy()
    else:
        arr = np.asarray(kspace_xyzc)
    if arr.ndim != 4:
        raise ValueError(
            f"kspace_xyzc must have shape (Nx, Ny, Nslices, Ncoil), got {arr.shape}"
        )

    x_step = max(1, int(x_step))
    cropped = crop_readout_oversampling_kspace(
        arr,
        readout_os_factor=x_step,
        axis=0,
    )
    nx, ny, nslices, ncoil = cropped.shape
    combined = (
        np.asarray(cropped, dtype=np.complex64)
        .transpose(0, 2, 1, 3)
        .reshape((nx * nslices, ny, ncoil))
    )
    return estimate_cc_matrix_coillast_2d(
        combined,
        ncc=ncc,
        acs=acs,
        x_step=1,
        eps=eps,
    )


def apply_cc_coillast_torch_2d(kspace_xyc, W, x_chunk=8):
    """
    Apply coil compression to one 2D coil-last slice.

    Parameters
    ----------
    kspace_xyc : torch.Tensor
        Shape (Nx, Ny, Ncoil).
    W : np.ndarray or torch.Tensor
        Shape (Ncoil, Ncc).

    Returns
    -------
    torch.Tensor
        Shape (Nx, Ny, Ncc).
    """
    if not isinstance(kspace_xyc, torch.Tensor):
        kspace_xyc = torch.as_tensor(kspace_xyc)

    if kspace_xyc.ndim != 3:
        raise ValueError(f"kspace_xyc must have shape (Nx, Ny, Ncoil), got {tuple(kspace_xyc.shape)}")

    Nx, Ny, ncoil = kspace_xyc.shape
    W_t = torch.as_tensor(W, dtype=kspace_xyc.dtype, device=kspace_xyc.device)
    if W_t.shape[0] != ncoil:
        raise ValueError(f"W first dim {W_t.shape[0]} must match ncoil={ncoil}")

    ncc = W_t.shape[1]
    out = torch.empty((Nx, Ny, ncc), dtype=kspace_xyc.dtype, device=kspace_xyc.device)

    x_chunk = max(1, int(x_chunk))
    with torch.no_grad():
        for x0 in range(0, Nx, x_chunk):
            x1 = min(x0 + x_chunk, Nx)
            X = kspace_xyc[x0:x1].reshape(-1, ncoil)
            Xcc = X @ W_t
            out[x0:x1] = Xcc.reshape(x1 - x0, Ny, ncc)
    return out


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


def estimate_cc_matrices_coillast_slicewise(kspace_xyzc, ncc=12, acs=24, x_step=1, eps=1e-8):
    """
    Estimate one 2D coil-compression matrix per acquired slice.

    kspace_xyzc: shape (Nx, Ny, Nslices, Ncoil)
    returns:
        W_slices: shape (Nslices, Ncoil, Ncc)
        svals_slices: shape (Nslices, Ncoil)
        energy_slices: shape (Nslices, Ncoil)
    """
    if isinstance(kspace_xyzc, torch.Tensor):
        shape = tuple(kspace_xyzc.shape)
    else:
        shape = np.asarray(kspace_xyzc).shape
    if len(shape) != 4:
        raise ValueError(f"kspace_xyzc must have shape (Nx, Ny, Nslices, Ncoil), got {shape}")

    Nz = shape[2]
    W_list, svals_list, energy_list = [], [], []
    for iz in range(Nz):
        W, svals, energy = estimate_cc_matrix_coillast_2d(
            kspace_xyzc[:, :, iz, :],
            ncc=ncc,
            acs=acs,
            x_step=x_step,
            eps=eps,
        )
        W_list.append(W)
        svals_list.append(svals)
        energy_list.append(energy)
    return (
        np.stack(W_list, axis=0).astype(np.complex64, copy=False),
        np.stack(svals_list, axis=0),
        np.stack(energy_list, axis=0),
    )


def apply_cc_coillast_torch_slicewise(kspace_xyzc, W_slices, x_chunk=8):
    """
    Apply slice-wise coil compression to a 2D multislice acquisition.

    Parameters
    ----------
    kspace_xyzc : torch.Tensor or np.ndarray
        Shape (Nx, Ny, Nslices, Ncoil).
    W_slices : np.ndarray or torch.Tensor
        Shape (Nslices, Ncoil, Ncc), or shape (Ncoil, Ncc) to reuse the same
        matrix for every slice.

    Returns
    -------
    torch.Tensor
        Shape (Nx, Ny, Nslices, Ncc).
    """
    if not isinstance(kspace_xyzc, torch.Tensor):
        kspace_xyzc = torch.as_tensor(kspace_xyzc)
    if kspace_xyzc.ndim != 4:
        raise ValueError(f"kspace_xyzc must have shape (Nx, Ny, Nslices, Ncoil), got {tuple(kspace_xyzc.shape)}")

    Nx, Ny, Nz, ncoil = kspace_xyzc.shape
    W_arr = torch.as_tensor(W_slices, dtype=kspace_xyzc.dtype, device=kspace_xyzc.device)
    if W_arr.ndim == 2:
        W_arr = W_arr.unsqueeze(0).expand(Nz, -1, -1)
    if W_arr.shape[0] != Nz or W_arr.shape[1] != ncoil:
        raise ValueError(
            f"W_slices must be (Nslices,Ncoil,Ncc) or (Ncoil,Ncc). "
            f"Got {tuple(W_arr.shape)}, expected Nslices={Nz}, Ncoil={ncoil}."
        )

    ncc = W_arr.shape[2]
    out = torch.empty((Nx, Ny, Nz, ncc), dtype=kspace_xyzc.dtype, device=kspace_xyzc.device)
    for iz in range(Nz):
        out[:, :, iz, :] = apply_cc_coillast_torch_2d(kspace_xyzc[:, :, iz, :], W_arr[iz], x_chunk=x_chunk)
    return out
