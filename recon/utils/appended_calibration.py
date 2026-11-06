"""Import and coil-map estimation for the 3D FLASH calibration appended to a TSE.

The generators ``seq/tse_2d_wave_mb2v_v3.m`` and ``seq/tse_2d_wave_mb1.m``
append the calibration of ``seq/utils/designFlashCalibration.m`` after the
TSE in the same sequence file (``appended_ref_lin_v1``). Every calibration
ADC is a REF readout with a compact zero-based LIN in acquisition order, so
the calibration arrives in the TWIX ``refscan`` stream and the TSE in the
``image`` stream. The TSE sidecar (``*_PEorder.mat``) holds the calibration
tables in its ``calibration`` struct.

The z convention is the one the mixed calibration measured on the scanner:
the data of an encoding moment ``k`` [1/m] are ``s(k) = sum_z m(z)
exp(+2i*pi*k*z)``. :func:`embed_kz` turns the ACS kz lines into a centered
k-space grid whose SigPy/centered inverse FFT places the image planes at the
physical slice centers.
"""

import numpy as np
import sigpy as sp
import sigpy.mri as mr

from utils.compact_calibration_import import read_compact_records


STORAGE_LAYOUT = "appended_ref_lin_v1"


def load_appended_metadata(pe):
    """Return the validated calibration metadata of a TSE sidecar.

    Args:
        pe: TSE sidecar loaded with ``scipy.io.loadmat(simplify_cells=True)``.

    Returns:
        dict: The sidecar ``calibration`` struct plus the TSE geometry used by
        the reconstruction (``Nslices``, ``slicePositions`` ascending,
        ``sliceThickness``, ``smsGroups`` for MB2) and ``kind`` =
        ``"appended"``. Index maps are integer arrays.

    Raises:
        ValueError: If the sidecar has no appended calibration or its tables
            are not a consistent bijection over the calibration ADCs.
    """
    if not bool(pe.get("isAppendCalibration", 0)) or "calibration" not in pe:
        raise ValueError("The TSE sidecar has no appended calibration; pass --calibration-twix/--calibration-seq")
    meta = dict(pe["calibration"])
    if meta.get("storageLayoutVersion") != STORAGE_LAYOUT:
        raise ValueError(f"Appended calibration layout {meta.get('storageLayoutVersion')} is not {STORAGE_LAYOUT}")
    count = int(meta["expectedAdcCount"])
    records = meta["acqTable"]
    if isinstance(records, dict):
        records = [records]
    if len(records) != count or not np.array_equal([r["rawLin0"] for r in records], np.arange(count)):
        raise ValueError("acqTable is not a contiguous compact LIN record table")
    if any(int(r["isRef"]) != 1 for r in records):
        raise ValueError("Every appended calibration readout must be a REF readout")
    psf = np.atleast_2d(np.asarray(meta["psfAdcIndices"], dtype=np.int64))
    valid = np.asarray(meta["psfValidMask"], dtype=bool).reshape(psf.shape)
    acs = np.atleast_2d(np.asarray(meta["acsAdcIndices"], dtype=np.int64))
    lengths = np.atleast_1d(np.asarray(meta["psfValidLengths"], dtype=np.int64))
    if psf.shape[0] != 4 or tuple(acs.shape) != tuple(np.atleast_1d(meta["Nacs"]).astype(int)):
        raise ValueError("Unexpected PSF/ACS map dimensions")
    if np.any(psf[~valid] != 0) or not np.array_equal(valid.sum(axis=1), lengths):
        raise ValueError("PSF padding/mask/valid lengths disagree")
    if not np.array_equal(np.sort(np.r_[psf[valid], acs.ravel()]), np.arange(1, count + 1)):
        raise ValueError("PSF and ACS maps are not a bijection over all calibration records")
    for part in range(4):
        if any(records[i - 1]["setID"] != part for i in psf[part, valid[part]]):
            raise ValueError("PSF map and scientific part IDs differ")
    ky = np.atleast_1d(np.asarray(meta["acsKyIndices"], dtype=np.int64))
    kz = np.atleast_1d(np.asarray(meta["acsKzIndices"], dtype=np.int64))
    for (row, column), index in np.ndenumerate(acs):
        record = records[index - 1]
        if record["setID"] != 4 or record["kyPhysicalIndex"] != ky[row] or record["kzPhysicalIndex"] != kz[column]:
            raise ValueError("ACS map, part ID and ky/kz indices differ")
    if np.any(np.diff(ky) != 1) or np.any(np.diff(kz) != 1):
        raise ValueError("ACS ky/kz lines are not contiguous")
    grid = int(meta["acsGridNz"])
    areas = np.atleast_1d(np.asarray(meta["acsKzAreas"], float))
    if not np.allclose(areas, (kz - 1 - grid // 2) / float(meta["acsFovZ"]), rtol=0, atol=1e-9):
        raise ValueError("ACS kz moments differ from the ACS kz grid")

    meta["kind"] = "appended"
    meta["acqTable"] = records
    meta["psfAdcIndices"], meta["psfValidMask"], meta["psfValidLengths"] = psf, valid, lengths
    meta["acsAdcIndices"], meta["acsKyIndices"], meta["acsKzIndices"] = acs, ky, kz
    meta["acsKzAreas"] = areas
    meta["slicePositions"] = np.sort(np.atleast_1d(np.asarray(meta["slicePositions"], float)))
    meta["Nslices"] = len(meta["slicePositions"])
    if int(pe.get("SMSFactor", 1)) == 2:
        meta["smsGroups"] = np.asarray(pe["smsGroups"], int)
    for key in ("Nx", "Ny", "Nx_os", "Ncycles", "tseAdcCount", "expectedAdcCount", "acsGridNz"):
        meta[key] = int(meta[key])
    meta["fov"] = np.atleast_1d(np.asarray(meta["fov"], float))[:2]
    return meta


def index_appended_stream(scan, meta):
    """Return the compact-LIN record index of the appended calibration.

    Args:
        scan: Header-only mapVBVD measurement of the TSE with its calibration.
        meta: :func:`load_appended_metadata` result.

    Returns:
        tuple: ``(names, local)`` arrays over the calibration LIN values, the
        stream name (``"refscan"``) and the record position in that stream,
        as used by ``read_compact_records``.

    Raises:
        ValueError: If the refscan stream is missing, does not hold every
            compact LIN exactly once, has nonzero storage counters or another
            readout length.
    """
    count = int(meta["expectedAdcCount"])
    if "refscan" not in scan:
        raise ValueError("The TSE TWIX file has no refscan stream; the appended calibration is missing")
    stream = scan["refscan"]
    raw_lin = np.asarray(stream.Lin, dtype=np.int64)
    if not np.array_equal(np.sort(raw_lin), np.arange(count)):
        raise ValueError("refscan LIN coverage differs from the appended calibration")
    if any(np.any(np.asarray(getattr(stream, key)) != 0) for key in ("Par", "Sli", "Set", "Eco", "Rep", "Ave", "Seg")):
        raise ValueError("refscan has unexpected nonzero storage counters")
    if int(stream.NCol) != int(meta["Nx_os"]):
        raise ValueError("Calibration ADC sample count differs from the sidecar")
    names = np.full(count, "refscan", dtype="U7")
    local = np.full(count, -1, dtype=np.int64)
    local[raw_lin] = np.arange(count)
    return names, local


def acs_z_grid(meta):
    """Return the z grid on which the ACS image planes fall on the slice centers.

    The image grid has ``acsGridNz`` planes spaced ``acsGridSpacingZ`` with
    centered coordinates shifted by ``z0`` (half a plane for an even slice
    count), so that every slice center is a grid plane.

    Args:
        meta: :func:`load_appended_metadata` result.

    Returns:
        dict: ``n`` (grid size), ``spacing`` [m], ``z0`` [m], ``positions``
        (n,) physical plane coordinates [m], ``slice_planes`` (Nslices,)
        plane index of every slice and ``line_index`` (Nacs_z,) grid index of
        every ACS kz line.

    Raises:
        ValueError: If a slice center or a kz moment is not on the grid.
    """
    n, spacing = int(meta["acsGridNz"]), float(meta["acsGridSpacingZ"])
    slices = np.asarray(meta["slicePositions"], float)
    z0 = slices[0] - np.round(slices[0] / spacing) * spacing
    steps = (slices - z0) / spacing
    planes = np.round(steps).astype(int) + n // 2
    if np.max(abs(steps - np.round(steps))) > 1e-6 or planes.min() < 0 or planes.max() >= n:
        raise ValueError("Slice centers are not planes of the ACS z grid")
    moments = np.asarray(meta["acsKzAreas"], float) * n * spacing
    if np.max(abs(moments - np.round(moments))) > 1e-6:
        raise ValueError("ACS kz moments are not multiples of the ACS kz step")
    # Reversing kz turns exp(+2i*pi*k*z) data into the exp(-2i*pi*k*z) kernel
    # of a centered inverse FFT. Only a complete grid reaches the Nyquist line
    # k = -n/2, whose reversal +n/2 is the same grid line modulo n.
    line_index = np.mod(n // 2 - np.round(moments).astype(int), n)
    if len(np.unique(line_index)) != len(line_index):
        raise ValueError("ACS kz lines do not fit on the ACS grid")
    return {"n": n, "spacing": spacing, "z0": float(z0), "positions": z0 + (np.arange(n) - n // 2) * spacing,
            "slice_planes": planes, "line_index": line_index}


def embed_kz(lines, meta, grid=None):
    """Place ACS kz lines on the centered grid of :func:`acs_z_grid`.

    Args:
        lines: Complex array whose last axis holds the ACS kz lines in sidecar
            order.
        meta: :func:`load_appended_metadata` result.
        grid: Optional precomputed :func:`acs_z_grid` result.

    Returns:
        numpy.ndarray: complex64 array with the last axis replaced by the
        ``n`` grid positions; a centered inverse FFT along it gives the image
        at ``grid["positions"]``.
    """
    grid = acs_z_grid(meta) if grid is None else grid
    lines = np.asarray(lines)
    shift = np.exp(-2j * np.pi * np.asarray(meta["acsKzAreas"], float) * grid["z0"]).astype(np.complex64)
    out = np.zeros(lines.shape[:-1] + (grid["n"],), np.complex64)
    out[..., grid["line_index"]] = lines * shift
    return out


def embed_ky(lines, meta, ny, axis):
    """Place ACS ky lines at their physical indices of the TSE PE grid.

    Args:
        lines: Complex array with the ACS ky lines on ``axis``.
        meta: :func:`load_appended_metadata` result.
        ny: TSE PE matrix size.
        axis: Axis of ``lines`` that holds the ky lines.

    Returns:
        numpy.ndarray: complex64 array with that axis of length ``ny``.
    """
    lines = np.moveaxis(np.asarray(lines), axis, -1)
    out = np.zeros(lines.shape[:-1] + (ny,), np.complex64)
    out[..., np.asarray(meta["acsKyIndices"]) - 1] = lines
    return np.moveaxis(out, -1, axis)


def read_acs(scan, record_index, meta, basis):
    """Return the compressed ACS hybrid space on the native readout rows.

    Args:
        scan: Header-only measurement with the calibration.
        record_index: :func:`index_appended_stream` result.
        meta: :func:`load_appended_metadata` result.
        basis: Coil-compression basis (coil, virtual coil).

    Returns:
        numpy.ndarray: complex64 (Nx, virtual coil, Nacs_y, Nacs_z): readout in
        image space (centered inverse FFT, central Nx samples), PE in k-space.
    """
    nx, nx_os = int(meta["Nx"]), int(meta["Nx_os"])
    nacs_y, nacs_z = np.asarray(meta["acsAdcIndices"]).shape
    rows = []
    for row in np.asarray(meta["acsAdcIndices"]):
        data = read_compact_records(scan, record_index, row) @ basis  # (Nx_os, Nacs_z, virtual coil)
        rows.append(data)
    data = np.stack(rows, axis=1)  # (Nx_os, Nacs_y, Nacs_z, virtual coil)
    start = (nx_os - nx) // 2
    hybrid = sp.ifft(data, axes=(0,))[start:start + nx]
    return np.ascontiguousarray(hybrid.transpose(0, 3, 1, 2), dtype=np.complex64)


def espirit_planes(hybrid_planes, meta, ny, calib_width, crop, device):
    """Estimate unit-RSS ESPIRiT maps of readout planes at the slice centers.

    Each readout position is calibrated independently in (ky, kz), as the
    slice2d backend of Wave-GRE does.

    Args:
        hybrid_planes: complex64 (plane, virtual coil, Nacs_y, Nacs_z).
        meta: :func:`load_appended_metadata` result.
        ny: TSE PE matrix size.
        calib_width, crop: ESPIRiT settings.
        device: SigPy device.

    Returns:
        tuple: ``(maps, eigenvalues, rss)`` shaped (plane, virtual coil, Ny,
        Nslices), (plane, Ny, Nslices) and (plane, Ny, Nslices); ``rss`` is the
        coil root-sum-of-squares of the zero-filled ACS image.
    """
    grid = acs_z_grid(meta)
    planes = grid["slice_planes"]
    nplane, ncoil = hybrid_planes.shape[:2]
    maps = np.zeros((nplane, ncoil, ny, len(planes)), np.complex64)
    eigen = np.zeros((nplane, ny, len(planes)), np.float32)
    rss = np.zeros((nplane, ny, len(planes)), np.float32)
    width = min(calib_width, *hybrid_planes.shape[2:])
    for index, plane in enumerate(hybrid_planes):
        kspace = embed_ky(embed_kz(plane, meta, grid), meta, ny, axis=1)  # (virtual coil, Ny, n)
        image = sp.ifft(kspace, axes=(1, 2))
        rss[index] = np.sqrt(np.sum(abs(image) ** 2, axis=0))[:, planes]
        if not np.any(kspace):
            continue
        result, value = mr.app.EspiritCalib(sp.to_device(kspace, device), calib_width=width, crop=crop,
                                            device=device, output_eigenvalue=True, show_pbar=False).run()
        result = np.asarray(sp.to_device(result, sp.Device(-1)), np.complex64)[:, :, planes]
        power = np.sqrt(np.sum(abs(result) ** 2, axis=0))
        result /= np.where(power >= 1e-8, power, 1)
        result *= power >= 1e-8
        maps[index] = result
        eigen[index] = np.asarray(sp.to_device(value, sp.Device(-1))).reshape(ny, grid["n"])[:, planes]
    return maps, eigen, rss


def compare_acs_with_z_reference(scan, record_index, meta):
    """Check the z orientation of the ACS planes against the z projection.

    The no-Wave z projection (part 2) is the established reference: its
    centered inverse FFT has physical z at minus its coordinate, as the mixed
    calibration showed with direct slices. It is low-pass filtered to the ACS
    kz extent and compared, at the center kx sample and ky = 0, with the ACS
    planes computed by :func:`embed_kz` at the stated (+1) and reversed (-1)
    slice positions.

    Args:
        scan: Header-only measurement with the calibration.
        record_index: :func:`index_appended_stream` result.
        meta: :func:`load_appended_metadata` result.

    Returns:
        dict: Per hypothesis the complex receive-vector correlation of every
        slice and its median and minimum over supported slices, the
        supported slices and a method note. No files are written.

    Raises:
        ValueError: Without a z projection, a ky = 0 ACS line or enough signal.
    """
    from utils.wave_cg_sense_precondition import ifftc_dim
    import torch

    n_ref = int(meta["psfValidLengths"][2])
    if n_ref == 0:
        raise ValueError("The appended calibration has no z projection reference")
    kx_center = int(meta["Nx_os"]) // 2
    ky = np.asarray(meta["acsKyIndices"])
    center = np.flatnonzero(ky == int(meta["Ny"]) // 2 + 1)
    if len(center) != 1:
        raise ValueError("The ACS has no unique ky = 0 line")
    grid = acs_z_grid(meta)
    acs = read_compact_records(scan, record_index, np.asarray(meta["acsAdcIndices"])[center[0]])[kx_center]
    acs_profile = sp.ifft(embed_kz(acs.T, meta, grid), axes=(-1,))[:, grid["slice_planes"]].T  # (slice, coil)

    indices = np.asarray(meta["psfAdcIndices"])[2, :n_ref]
    reference = read_compact_records(scan, record_index, indices)[kx_center]  # (kz line, coil)
    kz = np.array([meta["acqTable"][i - 1]["kzArea"] for i in indices])
    reference = np.where((abs(kz) <= np.max(abs(meta["acsKzAreas"])) + 1e-9)[:, None], reference, 0)
    profile = ifftc_dim(torch.from_numpy(np.ascontiguousarray(reference)), dim=0).numpy()
    coordinate = np.asarray(meta["psfZCoordinates"], float)
    positions = np.asarray(meta["slicePositions"], float)
    power = np.sum(abs(acs_profile) ** 2, axis=-1)
    valid = power > .05 * power.max()
    if valid.sum() < 10:
        raise ValueError("Insufficient ACS signal to check the z orientation")
    scores = {}
    for sign in (1, -1):
        # Physical z is minus the z-projection coordinate.
        sample = np.stack([np.interp(-sign * positions, coordinate, profile[:, c].real)
                           + 1j * np.interp(-sign * positions, coordinate, profile[:, c].imag)
                           for c in range(profile.shape[1])], axis=-1)
        correlation = abs(np.sum(acs_profile.conj() * sample, axis=-1)) / np.maximum(
            np.linalg.norm(acs_profile, axis=-1) * np.linalg.norm(sample, axis=-1), 1e-30)
        scores[str(sign)] = {"correlation_by_slice": correlation.tolist(),
                             "median_supported_correlation": float(np.median(correlation[valid])),
                             "minimum_supported_correlation": float(np.min(correlation[valid]))}
    return {"hypotheses": scores, "support_slices_0based": np.flatnonzero(valid).tolist(),
            "method": "ACS ky=0 planes versus the low-passed no-Wave z projection, complex receive vectors "
                      "at the center kx sample; +1 is the stated slice orientation"}
