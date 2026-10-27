"""Bounded-memory import for compact_adc_lin_v1 mixed FLASH calibration.

Python mapVBVD counters are zero-based, unlike MATLAB mapVBVD counters.
Sidecar ADC maps are one-based and refer to global records, not REF array
positions. Never interpret scanner SLC as physical slice identity here.
"""

from pathlib import Path

import mapvbvd
import numpy as np
from scipy.io import loadmat


def open_twix(path):
    """Return the last measurement header/streams without loading sample arrays."""
    measurements = mapvbvd.mapVBVD(str(Path(path).resolve(strict=True)), quiet=True)
    scan = measurements[-1] if isinstance(measurements, list) else measurements
    for name in ("image", "refscan"):
        if name in scan:
            scan[name].flagRemoveOS = False
            scan[name].squeeze = False
    return scan


def read_records(stream, indices, batch_size=32):
    """Read selected zero-based acquisition indices as (sample, record, coil).

    The returned complex64 array contains only requested records. Internal
    reads are bounded by batch_size; per-record reflection/RO corrections use
    the original acquisition indices. No sorted Cartesian grid is allocated.
    Raises ValueError for invalid indices and IOError for broken/nonfinite data.
    """
    indices = np.asarray(indices, dtype=np.int64).reshape(-1)
    if not len(indices) or np.any(indices < 0) or np.any(indices >= stream.NAcq):
        raise ValueError("Requested acquisition indices are empty or out of range.")
    if batch_size < 1:
        raise ValueError("batch_size must be positive.")
    stream.flagRemoveOS = False
    nx, nc = int(stream.NCol), int(stream.NCha)
    output = np.empty((nx, indices.size, nc), dtype=np.complex64)
    for start in range(0, indices.size, batch_size):
        selection = indices[start:start + batch_size]
        shape = np.array([nx, nc, len(selection)], dtype=np.int64)
        # Explicit cIxToRaw prevents unsorted(subset)'s local-index correction
        # bug: reflection/RO offsets must belong to these exact MDH records.
        block = stream.readData(
            np.asarray(stream.memPos)[selection],
            cIxToTarg=np.arange(len(selection)), cIxToRaw=selection,
            selRange=[np.arange(nx), np.arange(nc)],
            selRangeSz=shape, outSize=shape,
        ).reshape(tuple(shape))
        if getattr(stream, "isBrokenFile", False) or not np.isfinite(block).all():
            raise IOError("Broken or nonfinite Twix sample data.")
        output[:, start:start + len(selection), :] = block.transpose(0, 2, 1)
    return output


def load_compact_metadata(path):
    """Load and validate sidecar bijections, physical slices and ky coordinates.

    Returns the original MATLAB metadata dictionary (simplified scipy format).
    Legacy files must use their existing importer; malformed maps raise
    ValueError before any raw sample allocation.
    """
    meta = loadmat(Path(path).resolve(strict=True), simplify_cells=True)
    if meta.get("storageLayoutVersion") != "compact_adc_lin_v1":
        raise ValueError("Not compact_adc_lin_v1; use the legacy importer.")
    count = int(meta["expectedAdcCount"])
    records = meta["acqTable"]
    if len(records) != count or not np.array_equal(
        [r["rawLin0"] for r in records], np.arange(count)
    ):
        raise ValueError("acqTable is not a contiguous global LIN record table.")
    psf = np.asarray(meta["psfAdcIndices"], dtype=np.int64)
    valid = np.asarray(meta["psfValidMask"], dtype=bool)
    csm = np.asarray(meta["csmAdcIndices"], dtype=np.int64)
    expected_shape = (int(meta["Nslices"]), int(meta["Nacs_y"]))
    if psf.shape != valid.shape or csm.shape != expected_shape:
        raise ValueError("Unexpected PSF/CSM map dimensions.")
    if np.any(psf[~valid] != 0) or not np.array_equal(
        valid.sum(axis=1), meta["psfValidLengths"]
    ):
        raise ValueError("PSF padding/mask/valid lengths disagree.")
    if not np.array_equal(np.sort(np.r_[psf[valid], csm.ravel()]),
                          np.arange(1, count + 1)):
        raise ValueError("PSF and CSM maps are not a bijection over all records.")
    for part in range(4):
        for index in psf[part, valid[part]]:
            r = records[index - 1]
            if r["setID"] != part or r["isRef"] != 0:
                raise ValueError("PSF map and scientific part/REF routing differ.")
    for physical, row in enumerate(csm):
        selected = [records[index - 1] for index in row]
        if any(r["setID"] != 4 or r["isRef"] != 1 or
               r["slcLabel"] != physical for r in selected):
            raise ValueError("CSM map and physical slice/REF routing differ.")
        if not np.allclose([r["zPosition"] for r in selected],
                           meta["slicePositions"][physical], atol=1e-10):
            raise ValueError("CSM physical positions differ from slicePositions.")
        if not np.array_equal([r["linLocal"] for r in selected],
                              np.arange(csm.shape[1])):
            raise ValueError("CSM local ky is not ordered.")
        ky = np.array([r["kyPhysicalIndex"] for r in selected])
        if not np.all(np.diff(ky) == 1):
            raise ValueError("CSM ky coordinates are not contiguous.")
    return meta


def index_compact_streams(scan, meta):
    """Return global-LIN to (stream, local-record) maps, validating raw coverage.

    Input scan is a header-only mapVBVD measurement; meta is validated sidecar
    metadata. Scanner label extents, routing and unique complete coverage must
    match. This function does not load any sample data.
    """
    count = int(meta["expectedAdcCount"])
    names = np.full(count, "", dtype="U7")
    local = np.full(count, -1, dtype=np.int64)
    ncoil = None
    for name, map_name in (("image", "imageAdcIndices"),
                           ("refscan", "refAdcIndices")):
        expected = np.atleast_1d(meta[map_name]).astype(np.int64) - 1
        if expected.size == 0:
            if name in scan and scan[name].NAcq:
                raise ValueError(f"Unexpected {name} records.")
            continue
        if name not in scan:
            raise ValueError(f"Missing {name} stream.")
        stream = scan[name]
        raw_lin = np.asarray(stream.Lin, dtype=np.int64)
        if not np.array_equal(np.sort(raw_lin), np.sort(expected)):
            raise ValueError(f"{name} global LIN coverage differs from sidecar.")
        if any(np.any(np.asarray(getattr(stream, key)) != 0)
               for key in ("Par", "Sli", "Set", "Eco", "Rep", "Ave", "Seg")):
            raise ValueError(f"{name} has unexpected nonzero storage counters.")
        if int(stream.NCol) != int(meta["Nx_os"]):
            raise ValueError("ADC sample count differs from calibration metadata.")
        if ncoil is not None and ncoil != int(stream.NCha):
            raise ValueError("Calibration streams have different receive counts.")
        ncoil = int(stream.NCha)
        names[raw_lin] = name
        local[raw_lin] = np.arange(len(raw_lin))
    if np.any(local < 0):
        raise ValueError("Global compact calibration has missing ADC records.")
    return names, local


def read_compact_records(scan, record_index, adc_indices_1based, batch_size=32):
    """Return (sample, requested record, coil) data from a one-based sidecar map.

    Requests must belong to one stream, as each physical CSM slice and each
    PSF projection does. The order supplied by the sidecar is preserved.
    """
    names, local = record_index
    indices = np.asarray(adc_indices_1based, dtype=np.int64).reshape(-1) - 1
    if np.any(indices < 0) or np.any(indices >= len(local)):
        raise ValueError("Sidecar ADC index outside global record range.")
    selected = np.unique(names[indices])
    if len(selected) != 1:
        raise ValueError("Request must contain one nonempty scanner stream.")
    return read_records(scan[selected[0]], local[indices], batch_size)


def compare_direct_slice_fourier_z(scan, record_index, meta):
    """Return a calibration-only direction audit using complex receive profiles.

    Inputs are a compact Twix measurement, validated global record index and
    matching sidecar. Compares center-kx/ky direct-slice coil vectors against
    the center-kx Fourier-z reference at positive/negative sidecar positions.
    Only the diagnostic 3D profile is interpolated; no CSM slice interpolation
    is performed. No output files are written by this helper.
    """
    import torch
    from utils.wave_cg_sense_precondition import ifftc_dim

    kx_center = int(meta["Nx_os"]) // 2
    first = meta["csmAdcIndices"][0]
    ky = np.array([meta["acqTable"][index - 1]["kyPhysicalIndex"] for index in first])
    center = np.flatnonzero(ky == int(meta["Ny"]) // 2 + 1)
    if len(center) != 1:
        raise ValueError("Direct-slice CSM table has no unique ky=0 record.")
    direct = []
    for row in meta["csmAdcIndices"]:
        data = read_compact_records(scan, record_index, row[center])
        direct.append(data[kx_center, 0])
    direct = np.asarray(direct)
    indices = meta["psfAdcIndices"][2, :int(meta["psfValidLengths"][2])]
    data = read_compact_records(scan, record_index, indices)
    profile = ifftc_dim(torch.from_numpy(data), dim=1)[kx_center].numpy()
    coordinate = np.asarray(meta["psfZCoordinates"])
    positions = np.asarray(meta["slicePositions"])
    power = np.sum(abs(direct)**2, axis=-1)
    valid = power > .05 * power.max()
    if valid.sum() < 10:
        raise ValueError("Insufficient direct-slice signal to establish Fourier-z direction.")
    scores = {}
    for sign in (1, -1):
        sample = np.stack([np.interp(sign * positions, coordinate, profile[:, c].real)
                           + 1j * np.interp(sign * positions, coordinate, profile[:, c].imag)
                           for c in range(profile.shape[1])], axis=-1)
        correlation = abs(np.sum(direct.conj() * sample, axis=-1)) / np.maximum(
            np.linalg.norm(direct, axis=-1) * np.linalg.norm(sample, axis=-1), 1e-30)
        scores[str(sign)] = {"correlation_by_slice": correlation.tolist(),
                             "median_supported_correlation": float(np.median(correlation[valid])),
                             "minimum_supported_correlation": float(np.min(correlation[valid]))}
    return {"hypotheses": scores, "support_slices_0based": np.flatnonzero(valid).tolist(),
            "method": "complex receive-vector similarity; calibration only, no TSE residual fitting"}
