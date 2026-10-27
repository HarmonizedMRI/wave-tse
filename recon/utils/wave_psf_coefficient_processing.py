"""Process noisy 2D Wave PSF coefficients with central sine-line fits."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.ndimage import median_filter
from scipy.optimize import least_squares
from scipy.signal import lombscargle
import torch


TWO_PI = 2.0 * np.pi
NONLINEAR_FIT_MAX_NFEV = 500
AUTO_RANGE_ALGORITHM = "psf-sine-line-auto-range"
AUTO_RANGE_VERSION = 9
AUTO_FIT_PREFILTER_WINDOW = 9
CENTRAL_CORE_SAMPLES = 80
MINIMUM_CORE_SUPPORT_FRACTION = 0.75
QUALITY_DENSITY_WINDOW_SAMPLES = 21
MINIMUM_LOCAL_SUPPORT_FRACTION = 0.75
REGIONAL_MINIMUM_READOUT_FRACTION = 0.20
SUSTAINED_CORRUPTION_READOUT_FRACTION = 0.05
GLOBAL_MINIMUM_SUPPORT_FRACTION = 0.30
REGIONAL_MINIMUM_SUPPORT_FRACTION = 0.60
CENTER_SLOPE_LIMIT_MULTIPLIER = 2.5
CENTER_VARIANCE_LIMIT_MULTIPLIER = 6.0
CENTER_REFERENCE_READOUT_FRACTION = 0.20


def nan_moving_average(values, window):
    """Return a centered moving average that ignores non-finite samples.

    Parameters
    ----------
    values : array-like
        One-dimensional coefficient samples.
    window : int
        Positive odd moving-average width.

    Returns
    -------
    np.ndarray
        Smoothed values, with NaN where a window has no finite support.

    Raises
    ------
    ValueError
        If the input is not one-dimensional or the window is invalid.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("Coefficient values must be one-dimensional.")
    if window < 1 or window % 2 == 0:
        raise ValueError("Moving-average window must be a positive odd integer.")
    finite = np.isfinite(values)
    kernel = np.ones(window, dtype=np.float64)
    numerator = np.convolve(np.where(finite, values, 0.0), kernel, mode="same")
    denominator = np.convolve(finite.astype(np.float64), kernel, mode="same")
    result = np.full(values.shape, np.nan, dtype=np.float64)
    np.divide(numerator, denominator, out=result, where=denominator > 0)
    return result


def sine_line_model(t, amplitude, angular_frequency, phase, slope, intercept):
    """Evaluate a sine plus linear trend at readout coordinates.

    Parameters
    ----------
    t : array-like
        Readout sample coordinates.
    amplitude, angular_frequency, phase, slope, intercept : float
        Parameters of ``A*sin(w*t+phi)+C1*t+C2``.

    Returns
    -------
    np.ndarray
        Model values at the supplied coordinates.
    """
    t = np.asarray(t, dtype=np.float64)
    return amplitude * np.sin(angular_frequency * t + phase) + slope * t + intercept


def fit_sine_plus_line(t, values):
    """Fit a sine plus line and report optimizer and stability diagnostics.

    This implementation follows the Wave-MPRAGE upstream fitter: a
    Lomb-Scargle frequency initialization is refined with bounded nonlinear
    least squares in centered readout coordinates.

    Parameters
    ----------
    t, values : array-like
        One-dimensional finite readout coordinates and coefficient samples.

    Returns
    -------
    dict[str, Any]
        Fitted parameters and JSON-compatible numerical diagnostics.

    Raises
    ------
    ValueError
        If there is insufficient support or the optimizer fails.
    """
    t = np.asarray(t, dtype=np.float64).ravel()
    values = np.asarray(values, dtype=np.float64).ravel()
    valid = np.isfinite(t) & np.isfinite(values)
    t = t[valid]
    values = values[valid]
    if t.size < 6:
        raise ValueError("At least 6 finite samples are required for sine-line fitting.")
    order = np.argsort(t)
    t = t[order]
    values = values[order]
    unique_t = np.unique(t)
    if unique_t.size < 2 or np.ptp(unique_t) == 0:
        raise ValueError("Fit coordinates must contain distinct samples.")

    t_ref = float(np.mean(t))
    x = t - t_ref
    span = float(np.ptp(x))
    median_dt = float(np.median(np.diff(unique_t)))
    if not np.isfinite(median_dt) or median_dt <= 0:
        raise ValueError("Fit coordinates must have positive finite spacing.")
    w_min = TWO_PI / span
    w_max = np.pi / median_dt
    if not w_min < w_max:
        raise ValueError("Fit interval is too narrow to identify a sine frequency.")

    slope_initial, intercept_initial = np.polyfit(x, values, 1)
    detrended = values - (slope_initial * x + intercept_initial)
    detrended -= np.mean(detrended)
    frequency_grid = np.linspace(w_min, w_max, 10000)
    power = lombscargle(x, detrended, frequency_grid, normalize=True)
    w_initial = float(frequency_grid[int(np.argmax(power))])
    design = np.column_stack(
        [np.sin(w_initial * x), np.cos(w_initial * x), x, np.ones_like(x)]
    )
    sine_coef, cosine_coef, slope_initial, intercept_initial = np.linalg.lstsq(
        design, values, rcond=None
    )[0]
    initial = np.array(
        [
            max(float(np.hypot(sine_coef, cosine_coef)), np.finfo(float).eps),
            w_initial,
            float(np.arctan2(cosine_coef, sine_coef)),
            slope_initial,
            intercept_initial,
        ]
    )

    def residuals(parameters):
        """Evaluate residuals in centered coordinates for least squares.

        Parameters
        ----------
        parameters : array-like
            Candidate amplitude, frequency, phase, slope, and intercept.

        Returns
        -------
        np.ndarray
            Candidate model minus observed coefficient values.
        """
        amplitude, angular_frequency, phase_ref, slope, intercept_ref = parameters
        return (
            amplitude * np.sin(angular_frequency * x + phase_ref)
            + slope * x
            + intercept_ref
            - values
        )

    result = least_squares(
        residuals,
        initial,
        bounds=(
            [0.0, w_min, -np.inf, -np.inf, -np.inf],
            [np.inf, w_max, np.inf, np.inf, np.inf],
        ),
        method="trf",
        x_scale="jac",
        loss="linear",
        max_nfev=NONLINEAR_FIT_MAX_NFEV,
    )
    amplitude, angular_frequency, phase_ref, slope, intercept_ref = result.x
    phase = (phase_ref - angular_frequency * t_ref + np.pi) % TWO_PI - np.pi
    intercept = intercept_ref - slope * t_ref
    residual = np.asarray(result.fun, dtype=np.float64)
    jacobian = np.asarray(result.jac, dtype=np.float64)
    column_norm = np.linalg.norm(jacobian, axis=0)
    if np.any(~np.isfinite(column_norm)) or np.any(column_norm <= np.finfo(float).eps):
        condition_number = float("inf")
    else:
        condition_number = float(np.linalg.cond(jacobian / column_norm[None, :]))
    frequency_span = w_max - w_min
    boundary_fraction = float(
        min(angular_frequency - w_min, w_max - angular_frequency) / frequency_span
    )
    diagnostics: dict[str, Any] = {
        "A": float(amplitude),
        "w": float(angular_frequency),
        "phi": float(phase),
        "C1": float(slope),
        "C2": float(intercept),
        "period_samples": float(TWO_PI / angular_frequency),
        "cycles_in_fit_interval": float(angular_frequency * span / TWO_PI),
        "frequency_search_bounds": [float(w_min), float(w_max)],
        "frequency_boundary_fraction": boundary_fraction,
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "n_samples": int(t.size),
        "n_function_evaluations": int(result.nfev),
        "maximum_function_evaluations": NONLINEAR_FIT_MAX_NFEV,
        "cost": float(result.cost),
        "optimality": float(result.optimality),
        "residual_rmse": float(np.sqrt(np.mean(residual**2))),
        "residual_rmse_relative_to_range": float(
            np.sqrt(np.mean(residual**2))
            / max(float(np.ptp(values)), np.finfo(float).eps)
        ),
        "residual_median_absolute": float(np.median(np.abs(residual))),
        "residual_max_absolute": float(np.max(np.abs(residual))),
        "standardized_jacobian_condition_number": condition_number,
    }
    if not result.success or not np.all(np.isfinite(result.x)):
        raise ValueError(f"Sine-line optimizer failed: {result.message}")
    return diagnostics


def _smoothed_shape_metrics(values):
    """Measure local slope and variance after coefficient smoothing.

    Parameters
    ----------
    values : array-like
        One-dimensional raw coefficient vector.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Per-sample maximum adjacent slope and local variance.
    """
    smoothed = nan_moving_average(values, AUTO_FIT_PREFILTER_WINDOW)
    adjacent_slope = np.abs(np.diff(smoothed))
    slope = np.full(smoothed.shape, np.nan, dtype=np.float64)
    if adjacent_slope.size:
        slope[0] = adjacent_slope[0]
        slope[-1] = adjacent_slope[-1]
        if slope.size > 2:
            slope[1:-1] = np.maximum(adjacent_slope[:-1], adjacent_slope[1:])
    local_mean = nan_moving_average(smoothed, QUALITY_DENSITY_WINDOW_SAMPLES)
    local_variance = nan_moving_average(
        (smoothed - local_mean) ** 2,
        QUALITY_DENSITY_WINDOW_SAMPLES,
    )
    return slope, local_variance


def _quality_vector(quality, name, length, dtype):
    """Return one validated per-readout projection-quality vector.

    Parameters
    ----------
    quality : dict
        Projection-quality mapping.
    name : str
        Required vector key.
    length : int
        Expected vector length.
    dtype : numpy dtype
        Output conversion dtype.

    Returns
    -------
    np.ndarray
        Validated one-dimensional quality vector.
    """
    if name not in quality:
        raise ValueError(f"Projection quality is missing {name!r}.")
    vector = np.asarray(quality[name], dtype=dtype).reshape(-1)
    if vector.size != length:
        raise ValueError(
            f"Projection quality {name!r} has length {vector.size}; "
            f"expected {length}."
        )
    return vector


def _robust_upper_limit(values, absolute_limit):
    """Return a robust upper limit capped by a scientific absolute limit.

    Parameters
    ----------
    values : array-like
        Finite quality values.
    absolute_limit : float
        Maximum permitted threshold.

    Returns
    -------
    float
        Median/MAD-derived upper threshold.
    """
    values = np.asarray(values, dtype=np.float64)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    robust = max(median + 4.0 * 1.4826 * mad, 1.5 * median + 0.02)
    return float(min(absolute_limit, robust))


def _contiguous_true_runs(mask):
    """Return half-open intervals for contiguous true samples.

    Parameters
    ----------
    mask : array-like
        One-dimensional boolean mask.

    Returns
    -------
    list[tuple[int, int]]
        Half-open intervals in ascending order.
    """
    padded = np.pad(np.asarray(mask, dtype=bool), (1, 1), constant_values=False)
    changes = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1)
    return [
        (int(start), int(stop))
        for start, stop in zip(starts, stops, strict=True)
    ]


def _coefficient_stability_mask(values):
    """Reject transient and sustained raw-coefficient instabilities.

    Parameters
    ----------
    values : array-like
        One-dimensional raw coefficient vector.

    Returns
    -------
    tuple[np.ndarray, dict]
        Stability mask and JSON-compatible threshold diagnostics.
    """
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = np.isfinite(values)
    stable = finite.copy()
    if np.count_nonzero(finite) < 9:
        return stable, {"rejected_total_samples": 0}

    indices = np.arange(values.size, dtype=np.float64)
    filled = np.interp(indices, indices[finite], values[finite])
    local_median = median_filter(filled, size=9, mode="reflect")
    deviation = np.abs(filled - local_median)
    finite_deviation = deviation[finite]
    deviation_median = float(np.median(finite_deviation))
    deviation_mad = float(
        np.median(np.abs(finite_deviation - deviation_median))
    )
    deviation_p90 = float(np.percentile(finite_deviation, 90.0))
    finite_pairs = finite[:-1] & finite[1:]
    finite_steps = np.abs(np.diff(values)[finite_pairs])
    median_step = float(np.median(finite_steps)) if finite_steps.size else 0.0
    step_p90 = (
        float(np.percentile(finite_steps, 90.0)) if finite_steps.size else 0.0
    )
    deviation_limit = max(
        deviation_median + 8.0 * 1.4826 * deviation_mad,
        3.0 * deviation_p90,
        8.0 * median_step,
        4.0 * step_p90,
        100.0 * np.finfo(float).eps,
    )
    transient = finite & (deviation > deviation_limit)
    stable &= ~transient

    center_width = min(
        values.size,
        max(
            CENTRAL_CORE_SAMPLES,
            int(np.ceil(CENTER_REFERENCE_READOUT_FRACTION * values.size)),
        ),
    )
    center_start = (values.size - center_width) // 2
    center_values = values[center_start : center_start + center_width]
    center_values = center_values[np.isfinite(center_values)]
    if not center_values.size:
        center_values = values[finite]
    center_median = float(np.median(center_values))
    center_mad = float(np.median(np.abs(center_values - center_median)))
    amplitude_limit = max(
        5.0 * 1.4826 * center_mad,
        20.0 * median_step,
        0.25,
    )
    gross_amplitude = finite & (np.abs(values - center_median) > amplitude_limit)
    stable &= ~gross_amplitude

    smoothed_slope, smoothed_variance = _smoothed_shape_metrics(values)
    center_slice = slice(center_start, center_start + center_width)
    central_slope = smoothed_slope[center_slice]
    central_slope = central_slope[np.isfinite(central_slope)]
    central_variance = smoothed_variance[center_slice]
    central_variance = central_variance[np.isfinite(central_variance)]
    central_max_slope = float(np.max(central_slope)) if central_slope.size else 0.0
    central_max_variance = (
        float(np.max(central_variance)) if central_variance.size else 0.0
    )
    slope_limit = max(
        CENTER_SLOPE_LIMIT_MULTIPLIER * central_max_slope,
        100.0 * np.finfo(float).eps,
    )
    variance_limit = max(
        CENTER_VARIANCE_LIMIT_MULTIPLIER * central_max_variance,
        100.0 * np.finfo(float).eps,
    )
    excessive_shape = finite & (
        (smoothed_slope > slope_limit) | (smoothed_variance > variance_limit)
    )
    stable &= ~excessive_shape
    return stable, {
        "window_samples": 9,
        "deviation_limit": float(deviation_limit),
        "central_reference_interval": [center_start, center_start + center_width],
        "central_reference_median": center_median,
        "gross_amplitude_deviation_limit": float(amplitude_limit),
        "central_smoothed_max_absolute_slope": central_max_slope,
        "smoothed_slope_limit": float(slope_limit),
        "central_smoothed_max_local_variance": central_max_variance,
        "smoothed_local_variance_limit": float(variance_limit),
        "rejected_gross_amplitude_samples": int(np.count_nonzero(gross_amplitude)),
        "rejected_center_slope_or_variance_samples": int(
            np.count_nonzero(excessive_shape)
        ),
        "rejected_transient_samples": int(np.count_nonzero(transient)),
        "rejected_total_samples": int(np.count_nonzero(finite & ~stable)),
    }


def select_automatic_kx_range(raw_coefficients, projection_quality):
    """Select a reliable center-containing interval from a and b stability.

    This follows Wave-MPRAGE auto-range version 9. It combines smoothed
    a/b stability with the y/z wrapped-line fit quality at every kx. The c
    stability is reported but does not determine the range; c is validated
    against the independently fitted a/b frequency after selection.

    Parameters
    ----------
    raw_coefficients : sequence
        Equal-length raw coefficient vectors ``a``, ``b``, and ``c``.
    projection_quality : dict
        ``y`` and ``z`` mappings containing ``skipped``, ``valid_pixels``,
        ``masked_ratio``, and ``wrapped_rms`` vectors.

    Returns
    -------
    tuple[tuple[int, int], dict]
        Selected half-open range and JSON-compatible diagnostics. The private
        ``_fit_sample_mask`` entry identifies accepted samples in that range.

    Raises
    ------
    ValueError
        If inputs are malformed or no adequate center-containing range exists.
    """
    coefficients = [
        np.asarray(value, dtype=np.float64).reshape(-1)
        for value in raw_coefficients
    ]
    if len(coefficients) != 3 or not coefficients[0].size:
        raise ValueError("Auto-range requires non-empty a, b, and c vectors.")
    nx = int(coefficients[0].size)
    if any(value.size != nx for value in coefficients):
        raise ValueError("Auto-range coefficient vectors must have equal length.")
    if nx < CENTRAL_CORE_SAMPLES:
        raise ValueError(
            f"Auto-range requires at least {CENTRAL_CORE_SAMPLES} readout samples."
        )

    coefficient_valid = np.ones(nx, dtype=bool)
    coefficient_stability = {}
    for name, values in zip(("a", "b", "c"), coefficients, strict=True):
        stable, diagnostics = _coefficient_stability_mask(values)
        if name in ("a", "b"):
            coefficient_valid &= stable
        coefficient_stability[name] = diagnostics

    common_valid = coefficient_valid.copy()
    projection_diagnostics = {}
    for projection in ("y", "z"):
        if projection not in projection_quality:
            raise ValueError(f"Auto-range requires {projection!r} fit quality.")
        quality = projection_quality[projection]
        skipped = _quality_vector(quality, "skipped", nx, bool)
        valid_pixels = _quality_vector(quality, "valid_pixels", nx, float)
        masked_ratio = _quality_vector(quality, "masked_ratio", nx, float)
        wrapped_rms = _quality_vector(quality, "wrapped_rms", nx, float)
        support_fraction = 1.0 - masked_ratio
        provisional = (
            ~skipped
            & np.isfinite(valid_pixels)
            & np.isfinite(support_fraction)
            & np.isfinite(wrapped_rms)
            & (valid_pixels >= 10)
            & (support_fraction > 0.0)
        )
        if not np.any(provisional):
            raise ValueError(f"The {projection} fit has no reliable samples.")
        median_support = float(np.median(support_fraction[provisional]))
        support_floor = float(max(0.15, 0.5 * median_support))
        rms_ceiling = _robust_upper_limit(
            wrapped_rms[provisional], absolute_limit=0.7
        )
        accepted = (
            provisional
            & (support_fraction >= support_floor)
            & (wrapped_rms <= rms_ceiling)
        )
        common_valid &= accepted
        projection_diagnostics[projection] = {
            "provisional_samples": int(np.count_nonzero(provisional)),
            "accepted_samples": int(np.count_nonzero(accepted)),
            "minimum_valid_pixels": 10,
            "support_fraction_floor": support_floor,
            "wrapped_rms_ceiling_rad": rms_ceiling,
            "median_support_fraction": median_support,
            "median_wrapped_rms_rad": float(np.median(wrapped_rms[provisional])),
        }

    edge_guard = max(5, int(np.ceil(0.02 * nx)))
    common_valid[:edge_guard] = False
    common_valid[nx - edge_guard :] = False
    regional_minimum_width = max(
        CENTRAL_CORE_SAMPLES,
        24,
        int(np.ceil(REGIONAL_MINIMUM_READOUT_FRACTION * nx)),
    )
    core_start = (nx - CENTRAL_CORE_SAMPLES) // 2
    core_stop = core_start + CENTRAL_CORE_SAMPLES
    core_valid_samples = int(np.count_nonzero(common_valid[core_start:core_stop]))
    core_support_fraction = core_valid_samples / CENTRAL_CORE_SAMPLES
    if core_support_fraction < MINIMUM_CORE_SUPPORT_FRACTION:
        raise ValueError(
            "Auto-range has insufficient central support: "
            f"{core_valid_samples}/{CENTRAL_CORE_SAMPLES}."
        )

    density_kernel = np.ones(QUALITY_DENSITY_WINDOW_SAMPLES, dtype=np.float64)
    local_total = np.convolve(
        np.ones(nx, dtype=np.float64), density_kernel, mode="same"
    )
    local_valid = np.convolve(
        common_valid.astype(np.float64), density_kernel, mode="same"
    )
    dense_support = local_valid / local_total >= MINIMUM_LOCAL_SUPPORT_FRACTION
    dense_support[:edge_guard] = False
    dense_support[nx - edge_guard :] = False

    coefficient_local_valid = np.convolve(
        coefficient_valid.astype(np.float64), density_kernel, mode="same"
    )
    coefficient_dense = (
        coefficient_local_valid / local_total >= MINIMUM_LOCAL_SUPPORT_FRACTION
    )
    coefficient_dense[:edge_guard] = False
    coefficient_dense[nx - edge_guard :] = False
    low_density_runs = _contiguous_true_runs(~coefficient_dense)
    corruption_minimum_width = max(
        16, int(np.ceil(SUSTAINED_CORRUPTION_READOUT_FRACTION * nx))
    )
    interior_start = max(edge_guard, int(np.floor(0.15 * nx)))
    interior_stop = min(nx - edge_guard, int(np.ceil(0.85 * nx)))
    sustained_corruption = []
    for run in low_density_runs:
        overlap = max(
            0,
            min(run[1], interior_stop) - max(run[0], interior_start),
        )
        if overlap >= corruption_minimum_width:
            sustained_corruption.append(run)

    near_global = (edge_guard, nx - edge_guard)
    if sustained_corruption:
        candidates = [
            run
            for run in _contiguous_true_runs(coefficient_dense)
            if run[0] <= core_start
            and run[1] >= core_stop
            and run[1] - run[0] >= regional_minimum_width
        ]
        if not candidates:
            raise ValueError(
                "Auto-range detected sustained coefficient corruption but no "
                "adequate center-containing stable interval."
            )
        selected_start, selected_stop = max(
            candidates, key=lambda run: (run[1] - run[0], -run[0])
        )
        while selected_start > edge_guard and coefficient_valid[selected_start - 1]:
            selected_start -= 1
        while selected_stop < nx - edge_guard and coefficient_valid[selected_stop]:
            selected_stop += 1
        selected = (int(selected_start), int(selected_stop))
        strategy = "centered-region-after-sustained-coefficient-corruption"
        minimum_support = REGIONAL_MINIMUM_SUPPORT_FRACTION
    else:
        selected = near_global
        strategy = "near-global-no-sustained-coefficient-corruption"
        minimum_support = GLOBAL_MINIMUM_SUPPORT_FRACTION

    fit_sample_mask = np.zeros(nx, dtype=bool)
    fit_sample_mask[selected[0] : selected[1]] = common_valid[
        selected[0] : selected[1]
    ]
    fit_count = int(np.count_nonzero(fit_sample_mask))
    fit_support = fit_count / (selected[1] - selected[0])
    minimum_fit_samples = max(80, int(np.ceil(0.10 * nx)))
    if fit_count < minimum_fit_samples or fit_support < minimum_support:
        raise ValueError(
            f"Auto-range retained {fit_count} samples ({fit_support:.1%}) in "
            f"[{selected[0]}, {selected[1]})."
        )

    diagnostics = {
        "name": AUTO_RANGE_ALGORITHM,
        "version": AUTO_RANGE_VERSION,
        "coordinate": "oversampled_readout_sample_index",
        "range_convention": "half-open [min, max)",
        "readout_length": nx,
        "readout_center_index": nx // 2,
        "edge_guard_samples": edge_guard,
        "near_global_interval": list(near_global),
        "regional_minimum_interval_samples": regional_minimum_width,
        "required_central_core_interval": [core_start, core_stop],
        "central_core_valid_samples": core_valid_samples,
        "central_core_support_fraction": core_support_fraction,
        "quality_density_window_samples": QUALITY_DENSITY_WINDOW_SAMPLES,
        "minimum_local_support_fraction": MINIMUM_LOCAL_SUPPORT_FRACTION,
        "selection_strategy": strategy,
        "sustained_corruption_minimum_samples": corruption_minimum_width,
        "sustained_coefficient_corruption_intervals": [
            list(run) for run in sustained_corruption
        ],
        "selected_interval": list(selected),
        "fit_sample_count": fit_count,
        "fit_support_fraction": fit_support,
        "minimum_fit_samples": minimum_fit_samples,
        "minimum_interval_support_fraction": minimum_support,
        "common_valid_samples": int(np.count_nonzero(common_valid)),
        "quality_density_intervals": [
            list(run) for run in _contiguous_true_runs(dense_support)
        ],
        "coefficient_quality_density_intervals": [
            list(run) for run in _contiguous_true_runs(coefficient_dense)
        ],
        "coefficient_stability": coefficient_stability,
        "projection_quality": projection_diagnostics,
        "excluded_sample_indices_within_interval": (
            np.flatnonzero(~fit_sample_mask[selected[0] : selected[1]])
            + selected[0]
        ).astype(int).tolist(),
        "_fit_sample_mask": fit_sample_mask.tolist(),
    }
    return selected, diagnostics


def align_constant_phase_branch(values, reference_slice):
    """Align one constant-phase coefficient by an integer multiple of 2 pi.

    Parameters
    ----------
    values : array-like
        One-dimensional wrapped-branch coefficient vector.
    reference_slice : slice
        Central interval used to estimate the representative phase branch.

    Returns
    -------
    tuple[np.ndarray, int]
        Gauge-aligned vector and the number of added 2-pi turns.

    Raises
    ------
    ValueError
        If no finite reference samples are present.
    """
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    reference = values[reference_slice]
    reference = reference[np.isfinite(reference)]
    if not reference.size:
        raise ValueError("Constant-phase branch alignment has no finite reference samples.")
    turns = int(np.rint(-float(np.median(reference)) / TWO_PI))
    return values + turns * TWO_PI, turns


def _fit_one_coefficient(name, raw, fit_input, fit_indices, kx_all, fit_range):
    """Fit and validate one coefficient over a reviewed central interval.

    Parameters
    ----------
    name : str
        Coefficient name used in validation errors.
    raw, fit_input : np.ndarray
        Gauge-aligned raw values and prefiltered values.
    fit_indices, kx_all : np.ndarray
        Accepted fit coordinates and complete readout coordinates.
    fit_range : tuple[int, int]
        Reviewed half-open readout interval.

    Returns
    -------
    tuple[np.ndarray, dict]
        Complete fitted curve and fit/validation diagnostics.
    """
    params = fit_sine_plus_line(fit_indices, fit_input[fit_indices])
    fitted = sine_line_model(
        kx_all, params["A"], params["w"], params["phi"], params["C1"], params["C2"]
    )
    raw_prediction = sine_line_model(
        fit_indices,
        params["A"],
        params["w"],
        params["phi"],
        params["C1"],
        params["C2"],
    )
    raw_values = raw[fit_indices]
    raw_residual = raw_prediction - raw_values
    params["raw_observation_residual_rmse"] = float(
        np.sqrt(np.mean(raw_residual**2))
    )
    params["raw_observation_residual_rmse_relative_to_range"] = float(
        np.sqrt(np.mean(raw_residual**2))
        / max(float(np.ptp(raw_values)), np.finfo(float).eps)
    )

    interval_width = fit_range[1] - fit_range[0]
    trim = max(2, int(np.ceil(0.05 * interval_width)))
    trimmed_indices = fit_indices[
        (fit_indices >= fit_range[0] + trim) & (fit_indices < fit_range[1] - trim)
    ]
    stability: dict[str, Any] = {
        "trim_samples_per_side": trim,
        "refit_success": False,
        "refit_standardized_jacobian_condition_number": None,
        "full_readout_relative_l2_difference": None,
        "full_readout_max_difference_relative_to_fit_range": None,
    }
    if trimmed_indices.size >= 6:
        try:
            trimmed = fit_sine_plus_line(trimmed_indices, fit_input[trimmed_indices])
            trimmed_fitted = sine_line_model(
                kx_all,
                trimmed["A"],
                trimmed["w"],
                trimmed["phi"],
                trimmed["C1"],
                trimmed["C2"],
            )
            difference = fitted - trimmed_fitted
            stability.update(
                {
                    "refit_success": bool(
                        trimmed["success"]
                        and np.isfinite(trimmed["standardized_jacobian_condition_number"])
                        and trimmed["standardized_jacobian_condition_number"] <= 1.0e12
                    ),
                    "refit_standardized_jacobian_condition_number": trimmed[
                        "standardized_jacobian_condition_number"
                    ],
                    "full_readout_relative_l2_difference": float(
                        np.linalg.norm(difference)
                        / max(np.linalg.norm(fitted), np.finfo(float).eps)
                    ),
                    "full_readout_max_difference_relative_to_fit_range": float(
                        np.max(np.abs(difference))
                        / max(float(np.ptp(fit_input[fit_indices])), np.finfo(float).eps)
                    ),
                }
            )
        except ValueError as error:
            stability["error"] = str(error)

    gates = {
        "optimizer_converged": bool(params["success"]),
        "condition_number_at_most_1e12": bool(
            params["standardized_jacobian_condition_number"] <= 1.0e12
        ),
        "residual_rmse_relative_to_range_at_most_0p5": bool(
            params["residual_rmse_relative_to_range"] <= 0.5
        ),
        "period_at_least_4_samples": bool(params["period_samples"] >= 4.0),
        "frequency_not_on_search_boundary": bool(
            params["frequency_boundary_fraction"] >= 0.005
        ),
        "endpoint_trim_refit_succeeded": bool(stability["refit_success"]),
        "endpoint_trim_full_readout_relative_l2_at_most_1": bool(
            stability["full_readout_relative_l2_difference"] is not None
            and stability["full_readout_relative_l2_difference"] <= 1.0
        ),
        "full_fitted_curve_within_plus_minus_2pi": bool(
            np.max(np.abs(fitted)) <= TWO_PI
        ),
    }
    params["endpoint_trim_stability"] = stability
    params["validation_gates"] = gates
    params["validation_passed"] = bool(all(gates.values()))
    if not params["validation_passed"]:
        failed = [gate for gate, passed in gates.items() if not passed]
        raise ValueError(f"{name} sine-line validation failed: {failed}.")
    return fitted, params


def _fit_low_amplitude_frequency_locked_sine_line(
    raw,
    fit_input,
    fit_indices,
    kx_all,
    fit_range,
    reference_angular_frequency,
):
    """Fit a low-amplitude sine-plus-line at the programmed Wave frequency.

    The oscillatory frequency is fixed because a low-amplitude constant-phase
    term cannot identify it reliably. The sine/cosine coefficients, phase,
    slope, and intercept remain data-driven.

    Parameters
    ----------
    raw, fit_input : np.ndarray
        Gauge-aligned raw values and prefiltered values.
    fit_indices, kx_all : np.ndarray
        Accepted fit coordinates and complete readout coordinates.
    fit_range : tuple[int, int]
        Reviewed half-open readout interval.
    reference_angular_frequency : float
        Programmed Wave frequency in radians per readout sample.

    Returns
    -------
    tuple[np.ndarray, dict]
        Complete fitted sine-plus-line curve and validation diagnostics.

    Raises
    ------
    ValueError
        If the constrained fit is unstable, inaccurate, too large for the
        low-amplitude branch, or outside the accepted phase range.
    """
    fit_indices = np.asarray(fit_indices, dtype=np.float64)
    values = np.asarray(fit_input, dtype=np.float64)[fit_indices.astype(int)]
    if not np.isfinite(reference_angular_frequency) or reference_angular_frequency <= 0:
        raise ValueError("The locked c frequency must be finite and positive.")

    reference_design = np.column_stack(
        [
            np.sin(reference_angular_frequency * fit_indices),
            np.cos(reference_angular_frequency * fit_indices),
            fit_indices,
            np.ones(fit_indices.size),
        ]
    )
    reference_coefficients = np.linalg.lstsq(reference_design, values, rcond=None)[0]
    reference_amplitude = float(np.hypot(*reference_coefficients[:2]))
    reference_phase = float(
        np.arctan2(reference_coefficients[1], reference_coefficients[0])
    )
    full_design = np.column_stack(
        [
            np.sin(reference_angular_frequency * kx_all),
            np.cos(reference_angular_frequency * kx_all),
            kx_all,
            np.ones(kx_all.size),
        ]
    )
    fitted = full_design @ reference_coefficients
    prediction = reference_design @ reference_coefficients
    residual = prediction - values
    raw_values = np.asarray(raw, dtype=np.float64)[fit_indices.astype(int)]
    raw_residual = prediction - raw_values
    value_range = max(float(np.ptp(values)), np.finfo(float).eps)
    raw_range = max(float(np.ptp(raw_values)), np.finfo(float).eps)

    interval_width = fit_range[1] - fit_range[0]
    trim = max(2, int(np.ceil(0.05 * interval_width)))
    trimmed = fit_indices[
        (fit_indices >= fit_range[0] + trim) & (fit_indices < fit_range[1] - trim)
    ]
    trimmed_difference = None
    if trimmed.size >= 2:
        trimmed_values = np.asarray(fit_input, dtype=np.float64)[trimmed.astype(int)]
        trimmed_design = np.column_stack(
            [
                np.sin(reference_angular_frequency * trimmed),
                np.cos(reference_angular_frequency * trimmed),
                trimmed,
                np.ones(trimmed.size),
            ]
        )
        trimmed_coefficients = np.linalg.lstsq(
            trimmed_design, trimmed_values, rcond=None
        )[0]
        trimmed_fitted = full_design @ trimmed_coefficients
        trimmed_difference = float(
            np.linalg.norm(fitted - trimmed_fitted)
            / max(np.linalg.norm(fitted), np.finfo(float).eps)
        )

    diagnostics = {
        "model": "frequency-locked-low-amplitude-sine-line-fallback",
        "frequency_locked": True,
        "frequency_source": "sequence-expected-wave-cycles",
        "A": reference_amplitude,
        "w": float(reference_angular_frequency),
        "phi": reference_phase,
        "C1": float(reference_coefficients[2]),
        "C2": float(reference_coefficients[3]),
        "period_samples": float(TWO_PI / reference_angular_frequency),
        "full_readout_cycles": float(
            reference_angular_frequency * len(kx_all) / TWO_PI
        ),
        "reference_angular_frequency": float(reference_angular_frequency),
        "reference_frequency_amplitude_rad": reference_amplitude,
        "n_samples": int(fit_indices.size),
        "residual_rmse": float(np.sqrt(np.mean(residual**2))),
        "residual_rmse_relative_to_range": float(
            np.sqrt(np.mean(residual**2)) / value_range
        ),
        "raw_observation_residual_rmse": float(
            np.sqrt(np.mean(raw_residual**2))
        ),
        "raw_observation_residual_rmse_relative_to_range": float(
            np.sqrt(np.mean(raw_residual**2)) / raw_range
        ),
        "endpoint_trim_samples_per_side": trim,
        "endpoint_trim_full_readout_relative_l2_difference": trimmed_difference,
    }
    gates = {
        "reference_frequency_amplitude_at_most_0p05_rad": bool(
            reference_amplitude <= 0.05
        ),
        "residual_rmse_relative_to_range_at_most_0p5": bool(
            diagnostics["residual_rmse_relative_to_range"] <= 0.5
        ),
        "raw_residual_rmse_relative_to_range_at_most_0p5": bool(
            diagnostics["raw_observation_residual_rmse_relative_to_range"] <= 0.5
        ),
        "endpoint_trim_full_readout_relative_l2_at_most_1": bool(
            trimmed_difference is not None and trimmed_difference <= 1.0
        ),
        "full_fitted_curve_within_plus_minus_2pi": bool(
            np.max(np.abs(fitted)) <= TWO_PI
        ),
    }
    diagnostics["validation_gates"] = gates
    diagnostics["validation_passed"] = bool(all(gates.values()))
    if not diagnostics["validation_passed"]:
        failed = [gate for gate, passed in gates.items() if not passed]
        raise ValueError(
            "c low-amplitude frequency-locked sine-line fallback validation "
            "failed: " + str(failed) + "."
        )
    return fitted, diagnostics


def process_wave_yz_coefficients(
    a_raw,
    b_raw,
    c_y_raw,
    c_z_raw,
    fit_kx_min,
    fit_kx_max,
    prefilter_window=9,
    expected_cycles=10,
    projection_quality=None,
    allow_low_amplitude_c_frequency_locked_fallback=False,
    minimum_automatic_support_fraction=0.60,
):
    """Process two-axis Wave coefficients with manual or automatic support.

    The two projection pairs independently provide ``a, c_y`` and ``b, c_z``.
    Following the upstream Wave-MPRAGE convention, the fitted constant is
    ``c = c_y + c_z``. Only the reviewed central readout interval contributes
    to the NaN-aware prefilter or sine-plus-line fits.

    Parameters
    ----------
    a_raw, b_raw, c_y_raw, c_z_raw : array-like
        Equal-length raw coefficient vectors in radians.
    fit_kx_min, fit_kx_max : int or None
        Manual half-open fit interval. Set both to ``None`` for automatic
        Wave-MPRAGE-style range selection.
    prefilter_window : int
        Positive odd NaN-aware moving-average width used before fitting.
    expected_cycles : int
        Sequence-defined Wave cycles over the complete readout.
    projection_quality : dict or None
        Full-readout y/z wrapped-line quality evidence required by automatic
        range selection.
    allow_low_amplitude_c_frequency_locked_fallback : bool
        If true, permit a validated low-amplitude ``c`` sine-plus-line model
        whose frequency is fixed to the sequence-defined Wave cycles when the
        unconstrained ``c`` fit fails stability or its frequency disagrees with
        the reliable ``a/b`` fit. All low-amplitude fallback gates still apply.
    minimum_automatic_support_fraction : float
        Required valid-sample fraction inside an automatically selected fit
        interval. The historical default is 0.60; dataset-specific callers
        may lower it only while retaining all subsequent frequency gates.

    Returns
    -------
    dict
        Processed coefficients, masked fit inputs, gauge information, and
        serializable fit diagnostics.

    Raises
    ------
    ValueError
        If shapes, support, coefficient bounds, periods, or fit gates fail.
    """
    raw_vectors = {}
    for name, values in (
        ("a", a_raw),
        ("b", b_raw),
        ("c_y", c_y_raw),
        ("c_z", c_z_raw),
    ):
        vector = np.asarray(
            torch.as_tensor(values).detach().cpu(), dtype=np.float64
        ).squeeze()
        if vector.ndim != 1:
            raise ValueError(f"Raw {name} coefficients must be one-dimensional.")
        raw_vectors[name] = vector
    shapes = {values.shape for values in raw_vectors.values()}
    if len(shapes) != 1:
        raise ValueError("Raw two-axis coefficient vectors must have equal length.")

    nx = int(raw_vectors["a"].size)
    if expected_cycles <= 0:
        raise ValueError("Expected Wave cycles must be positive.")
    minimum_automatic_support_fraction = float(
        minimum_automatic_support_fraction
    )
    if not 0.0 < minimum_automatic_support_fraction <= 1.0:
        raise ValueError(
            "minimum_automatic_support_fraction must lie in (0, 1]."
        )

    c_raw = raw_vectors["c_y"] + raw_vectors["c_z"]
    if (fit_kx_min is None) != (fit_kx_max is None):
        raise ValueError("Specify both manual fit bounds or neither.")
    if fit_kx_min is None:
        if projection_quality is None:
            raise ValueError(
                "Automatic two-axis fitting requires projection quality evidence."
            )
        fit_range, range_diagnostics = select_automatic_kx_range(
            (raw_vectors["a"], raw_vectors["b"], c_raw),
            projection_quality,
        )
        fit_kx_min, fit_kx_max = fit_range
        fit_mask = np.asarray(
            range_diagnostics.pop("_fit_sample_mask"), dtype=bool
        )
        fit_range_selection = "automatic"
    else:
        fit_kx_min = int(fit_kx_min)
        fit_kx_max = int(fit_kx_max)
        if not (0 <= fit_kx_min < fit_kx_max <= nx):
            raise ValueError(
                "PSF fit range must satisfy 0 <= min < max <= readout length; "
                f"got [{fit_kx_min}, {fit_kx_max}) for {nx} samples."
            )
        fit_range = (fit_kx_min, fit_kx_max)
        fit_mask = np.ones(nx, dtype=bool)
        range_diagnostics = {
            "name": "manual-half-open-range",
            "version": 1,
            "selected_interval": list(fit_range),
        }
        fit_range_selection = "manual"

    reference_slice = slice(*fit_range)
    c_aligned, c_gauge_turns = align_constant_phase_branch(c_raw, reference_slice)
    coefficients = {
        "a": raw_vectors["a"],
        "b": raw_vectors["b"],
        "c": c_aligned,
    }
    common_valid = np.ones(nx, dtype=bool)
    for values in coefficients.values():
        common_valid &= np.isfinite(values) & (np.abs(values) <= TWO_PI)
    interval_mask = np.zeros(nx, dtype=bool)
    interval_mask[fit_kx_min:fit_kx_max] = True
    fit_mask &= common_valid & interval_mask
    fit_indices = np.flatnonzero(fit_mask)
    support_fraction = fit_indices.size / (fit_kx_max - fit_kx_min)
    minimum_support_fraction = (
        0.90
        if fit_range_selection == "manual"
        else minimum_automatic_support_fraction
    )
    if fit_indices.size < 80 or support_fraction < minimum_support_fraction:
        raise ValueError(
            f"Two-axis PSF interval retained {fit_indices.size} samples "
            f"({support_fraction:.1%}); at least 80 and "
            f"{minimum_support_fraction:.0%} are required."
        )

    fit_inputs = {
        name: nan_moving_average(np.where(fit_mask, values, np.nan), prefilter_window)
        for name, values in coefficients.items()
    }
    for values in fit_inputs.values():
        fit_mask &= np.isfinite(values)
    fit_indices = np.flatnonzero(fit_mask)
    kx_all = np.arange(nx, dtype=np.float64)
    processed = {}
    diagnostics = {}
    expected_period = nx / float(expected_cycles)
    expected_frequency = TWO_PI / expected_period
    c_fit_failure = None
    for name in ("a", "b", "c"):
        try:
            processed[name], diagnostics[name] = _fit_one_coefficient(
                name, coefficients[name], fit_inputs[name], fit_indices, kx_all, fit_range
            )
        except ValueError as error:
            if name != "c" or not allow_low_amplitude_c_frequency_locked_fallback:
                raise
            # A near-zero c term may not identify an independent frequency.
            # This opt-in branch retains the existing amplitude, residual,
            # trim-stability and full-range gates of the constrained model.
            processed[name], diagnostics[name] = _fit_low_amplitude_frequency_locked_sine_line(
                coefficients[name], fit_inputs[name], fit_indices, kx_all,
                fit_range, expected_frequency,
            )
            c_fit_failure = str(error)
            diagnostics[name]["independent_sine_line_failure"] = c_fit_failure

    period_relative_errors = {
        name: abs(values["period_samples"] - expected_period) / expected_period
        for name, values in diagnostics.items()
    }
    a_frequency = diagnostics["a"]["w"]
    b_frequency = diagnostics["b"]["w"]
    c_frequency = diagnostics["c"]["w"]
    ab_reference_frequency = 0.5 * (a_frequency + b_frequency)
    frequency_relative_errors = {
        name: abs(values["w"] - expected_frequency) / expected_frequency
        for name, values in diagnostics.items()
    }
    ab_frequency_relative_difference = (
        abs(a_frequency - b_frequency) / expected_frequency
    )
    c_frequency_relative_difference_to_ab = (
        abs(c_frequency - ab_reference_frequency) / ab_reference_frequency
    )
    period_gates = {
        "a_frequency_within_10_percent_of_sequence": (
            frequency_relative_errors["a"] <= 0.10
        ),
        "b_frequency_within_10_percent_of_sequence": (
            frequency_relative_errors["b"] <= 0.10
        ),
        "a_b_frequency_difference_within_10_percent": (
            ab_frequency_relative_difference <= 0.10
        ),
        "c_frequency_within_10_percent_of_ab_reference": (
            c_frequency_relative_difference_to_ab <= 0.10
        ),
    }
    c_fallback_reason = c_fit_failure
    if c_fit_failure is not None:
        period_relative_errors["c"] = None
        frequency_relative_errors["c"] = None
        c_frequency_relative_difference_to_ab = None
        period_gates.pop("c_frequency_within_10_percent_of_ab_reference")
        period_gates["c_low_amplitude_frequency_locked_sine_line_fallback_validated"] = True
    if not all(period_gates.values()):
        failed = [gate for gate, passed in period_gates.items() if not passed]
        c_only_failure = failed == [
            "c_frequency_within_10_percent_of_ab_reference"
        ]
        if (
            not allow_low_amplitude_c_frequency_locked_fallback
            or not c_only_failure
        ):
            raise ValueError(
                f"Two-axis sine-line frequency validation failed: {failed}."
            )
        c_fallback_reason = (
            "Independent c sine-line cross-frequency validation failed: "
            f"{failed}."
        )
        processed["c"], diagnostics["c"] = (
            _fit_low_amplitude_frequency_locked_sine_line(
                coefficients["c"],
                fit_inputs["c"],
                fit_indices,
                kx_all,
                fit_range,
                expected_frequency,
            )
        )
        diagnostics["c"]["independent_sine_line_failure"] = c_fallback_reason
        period_relative_errors["c"] = None
        frequency_relative_errors["c"] = None
        c_frequency_relative_difference_to_ab = None
        period_gates.pop("c_frequency_within_10_percent_of_ab_reference")
        period_gates[
            "c_low_amplitude_frequency_locked_sine_line_fallback_validated"
        ] = True

    report = {
        "coefficient_processing": "central-sine-line",
        "upstream_reference": "wave-mprage/recon: a=a_y, b=b_z, c=c_y+c_z",
        "model": "A*sin(w*kx+phi)+C1*kx+C2",
        "fit_range_selection": fit_range_selection,
        "fit_range_selection_diagnostics": range_diagnostics,
        "kx_range": [fit_kx_min, fit_kx_max],
        "kx_range_convention": "half-open [min, max)",
        "fit_input_processing": {
            "name": "range-masked-nan-aware-moving-average",
            "window_samples": int(prefilter_window),
        },
        "accepted_samples": int(fit_indices.size),
        "fit_support_fraction": float(
            fit_indices.size / (fit_kx_max - fit_kx_min)
        ),
        "minimum_required_fit_support_fraction": float(
            minimum_support_fraction
        ),
        "coefficient_acceptance_range_rad": [-TWO_PI, TWO_PI],
        "constant_combination": "c_y + c_z",
        "c_gauge_added_turns_2pi": int(c_gauge_turns),
        "c_gauge_added_radians": float(c_gauge_turns * TWO_PI),
        "expected_wave_cycles": int(expected_cycles),
        "expected_period_samples": float(expected_period),
        "expected_angular_frequency_rad_per_sample": float(expected_frequency),
        "period_relative_errors": period_relative_errors,
        "frequency_relative_errors": frequency_relative_errors,
        "ab_reference_angular_frequency_rad_per_sample": float(
            ab_reference_frequency
        ),
        "a_b_frequency_relative_difference": float(
            ab_frequency_relative_difference
        ),
        "c_frequency_relative_difference_to_ab": float(
            c_frequency_relative_difference_to_ab
        ) if c_frequency_relative_difference_to_ab is not None else None,
        "c_model": diagnostics["c"].get("model", "independent-sine-line"),
        "period_validation_gates": period_gates,
        "coefficients": diagnostics,
        "validation_passed": True,
    }
    dtype = torch.float32
    return {
        "a_processed": torch.as_tensor(processed["a"], dtype=dtype),
        "b_processed": torch.as_tensor(processed["b"], dtype=dtype),
        "c_processed": torch.as_tensor(processed["c"], dtype=dtype),
        "c_raw": torch.as_tensor(c_raw, dtype=dtype),
        "c_raw_aligned": torch.as_tensor(c_aligned, dtype=dtype),
        "a_fit_input": torch.as_tensor(fit_inputs["a"], dtype=dtype),
        "b_fit_input": torch.as_tensor(fit_inputs["b"], dtype=dtype),
        "c_fit_input": torch.as_tensor(fit_inputs["c"], dtype=dtype),
        "fit_sample_mask": torch.as_tensor(fit_mask, dtype=torch.bool),
        "report": report,
    }
