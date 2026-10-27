"""Noise covariance estimation and receive-coil prewhitening utilities."""

import mapvbvd
import numpy as np


def load_twix_noise_samples_coillast(filename):
    """Load the Twix noise-adjustment scan as ``(Nsamples, Ncoil)``.

    The noise scan is normally the first measurement in a multi-RAID Siemens
    Twix file. All non-coil dimensions are flattened into independent noise
    samples; no Fourier transform is applied.
    """
    twix_obj = mapvbvd.mapVBVD(filename, quiet=True)
    scans = twix_obj if isinstance(twix_obj, list) else [twix_obj]
    noise_obj = next(
        (scan.noise for scan in scans if getattr(scan, "noise", None) is not None),
        None,
    )
    if noise_obj is None:
        raise ValueError(f"Twix file has no noise-adjustment scan: {filename}")

    noise_obj.flagRemoveOS = False
    data = np.asarray(noise_obj[""])
    data_dims = list(noise_obj.dataDims)
    full_size = tuple(int(value) for value in noise_obj.fullSize)
    retained_dims = [name for name, size in zip(data_dims, full_size) if size > 1]
    data = np.squeeze(data)
    if "Cha" not in retained_dims:
        raise ValueError("Twix noise scan has no non-singleton channel dimension.")
    if data.ndim != len(retained_dims):
        raise RuntimeError(
            f"Unexpected noise array shape {data.shape} for retained dimensions {retained_dims}."
        )

    channel_axis = retained_dims.index("Cha")
    data = np.moveaxis(data, channel_axis, -1)
    samples = np.ascontiguousarray(data.reshape(-1, data.shape[-1]), dtype=np.complex64)
    if samples.shape[0] < samples.shape[1] or not np.isfinite(samples).all():
        raise ValueError(f"Invalid Twix noise samples with shape {samples.shape}.")
    return samples


def estimate_noise_covariance(samples_coillast, normalize_trace=True):
    """Estimate a Hermitian receive-noise covariance matrix.

    Parameters
    ----------
    samples_coillast : np.ndarray
        Noise samples shaped ``(..., Ncoil)``.
    normalize_trace : bool
        Normalize the mean diagonal variance to one. This removes an arbitrary
        global noise scale without changing relative channel weighting.
    """
    samples = np.asarray(samples_coillast)
    if samples.ndim < 2:
        raise ValueError("Noise samples must have a final receive-coil dimension.")
    samples = samples.reshape(-1, samples.shape[-1]).astype(np.complex128, copy=False)
    if samples.shape[0] < samples.shape[1] or not np.isfinite(samples).all():
        raise ValueError(f"Invalid noise sample matrix with shape {samples.shape}.")

    samples = samples - samples.mean(axis=0, keepdims=True)
    covariance = samples.conj().T @ samples / samples.shape[0]
    covariance = 0.5 * (covariance + covariance.conj().T)
    if normalize_trace:
        mean_variance = float(np.trace(covariance).real / covariance.shape[0])
        if not np.isfinite(mean_variance) or mean_variance <= 0:
            raise ValueError("Noise covariance has a nonpositive mean variance.")
        covariance /= mean_variance
    return covariance.astype(np.complex64)


def make_noise_whitener(covariance, eigenvalue_floor_relative=1e-8):
    """Return ``A`` for row-coil data such that ``A.H @ Psi @ A = I``."""
    covariance = np.asarray(covariance, dtype=np.complex128)
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError(f"Noise covariance must be square, got {covariance.shape}.")
    if not np.isfinite(covariance).all():
        raise ValueError("Noise covariance contains nonfinite values.")

    covariance = 0.5 * (covariance + covariance.conj().T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    largest = float(eigenvalues[-1])
    floor = float(eigenvalue_floor_relative) * largest
    if largest <= 0 or eigenvalues[0] <= floor:
        raise ValueError(
            "Noise covariance is singular or below the requested eigenvalue floor: "
            f"min={eigenvalues[0]:.6e}, floor={floor:.6e}."
        )

    whitener = (
        eigenvectors * (1.0 / np.sqrt(eigenvalues))[None, :]
    ) @ eigenvectors.conj().T
    return whitener.astype(np.complex64), eigenvalues.astype(np.float64)


def validate_noise_whitener(covariance, whitener, atol=2e-5):
    """Validate the row-coil whitening convention and return maximum error."""
    covariance = np.asarray(covariance, dtype=np.complex128)
    whitener = np.asarray(whitener, dtype=np.complex128)
    if whitener.shape != covariance.shape:
        raise ValueError(
            f"Whitener shape {whitener.shape} does not match covariance {covariance.shape}."
        )
    whitened = whitener.conj().T @ covariance @ whitener
    error = float(np.max(np.abs(whitened - np.eye(covariance.shape[0]))))
    if not np.isfinite(error) or error > atol:
        raise ValueError(f"Noise whitener validation failed with max error {error:.6e}.")
    return error
