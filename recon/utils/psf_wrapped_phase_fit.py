"""
psf_wrapped_phase_fit.py

Fit wrapped phase lines

    p(kx, r) = a(kx) r + c(kx)

independently for every readout-hybrid sample kx, for wrapped phase data p in
[-pi, pi) along one normalized encoding coordinate r, with optional magnitude
weighting and robust (Huber) iteratively reweighted least squares.

Typical use
-----------
from utils.psf_wrapped_phase_fit import fit_wrapped_phase_lines_2d

result = fit_wrapped_phase_lines_2d(psf_diff, hyb_nowave, y_norm, mask_mode="none")
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Union

import torch

TensorLike1D = Union[torch.Tensor, Sequence[float]]


def wrap_to_pi(x: torch.Tensor) -> torch.Tensor:
    """Wrap radians to [-pi, pi)."""
    return torch.atan2(torch.sin(x), torch.cos(x))


def _mean_over_dims(
    x: torch.Tensor,
    dims: Optional[Union[int, Sequence[int]]],
) -> torch.Tensor:
    """Mean over one or more dims, supporting negative dim indexing."""
    if dims is None:
        return x
    if isinstance(dims, int):
        dims = (dims,)
    else:
        dims = tuple(dims)
    ndim = x.ndim
    dims = tuple(d if d >= 0 else ndim + d for d in dims)
    for d in sorted(dims, reverse=True):
        x = x.mean(dim=d)
    return x


def _weighted_lstsq(
    A: torch.Tensor,
    b: torch.Tensor,
    w: torch.Tensor,
    ridge: float = 1e-8,
) -> torch.Tensor:
    """Solve min ||sqrt(w) * (A theta - b)||^2 with ridge stabilization."""
    w = torch.clamp(w, min=0.0)
    sqrt_w = torch.sqrt(w)
    Aw = A * sqrt_w[:, None]
    bw = b * sqrt_w

    ATA = Aw.T @ Aw
    ATb = Aw.T @ bw
    eye = torch.eye(ATA.shape[0], dtype=ATA.dtype, device=ATA.device)
    ridge_scaled = ridge * torch.trace(ATA).clamp_min(1.0) / ATA.shape[0]

    try:
        return torch.linalg.solve(ATA + ridge_scaled * eye, ATb)
    except RuntimeError:
        return torch.linalg.lstsq(Aw, bw).solution


def _run_wrapped_irls(
    A: torch.Tensor,
    p: torch.Tensor,
    w_base: torch.Tensor,
    theta0: torch.Tensor,
    *,
    n_irls: int,
    huber_delta: float,
    ridge: float,
    convergence_tol: float,
) -> torch.Tensor:
    """Robust model-based wrapped IRLS for one slice."""
    eps = torch.finfo(p.dtype).eps
    theta = theta0.clone()

    for _ in range(n_irls):
        pred = A @ theta
        r = wrap_to_pi(p - pred)
        p_model_unwrapped = pred + r

        abs_r = torch.abs(r)
        w_huber = torch.ones_like(abs_r)
        large = abs_r > huber_delta
        w_huber[large] = huber_delta / abs_r[large].clamp_min(eps)

        theta_new = _weighted_lstsq(A, p_model_unwrapped, w_base * w_huber, ridge=ridge)

        denom = torch.norm(theta).clamp_min(1.0)
        if torch.norm(theta_new - theta) < convergence_tol * denom:
            theta = theta_new
            break
        theta = theta_new

    return theta


# -----------------------------------------------------------------------------
# 2D phase-line fitting helper
# -----------------------------------------------------------------------------

def fit_wrapped_phase_lines_2d(
    psf_diff: torch.Tensor,
    hyb_nowave: torch.Tensor,
    y_norm: TensorLike1D,
    *,
    mag_threshold: float = 8e-10,
    max_masked_threshold: float = 0.94,
    n_irls: int = 10,
    huber_delta: float = 0.7,
    ridge: float = 1e-8,
    min_valid_pixels: int = 10,
    weight_clip: Optional[float] = 100.0,
    mag_average_dims: Optional[Union[int, Sequence[int]]] = -1,
    use_previous_fit_as_init: bool = True,
    initial_coefficients: Optional[torch.Tensor] = None,
    verbose: bool = True,
    return_wrapped_prediction: bool = True,
    mask_mode: str = "adaptive_quantile",
    mag_abs_floor: Optional[float] = None,
    mag_lower_quantile: float = 0.10,
    min_keep_fraction: float = 0.30,
    use_mag_weight: bool = True,
    convergence_tol: float = 1e-6,
) -> Dict[str, torch.Tensor]:
    """
    Fit wrapped phase lines independently for each readout-hybrid sample.

    This is the 2D counterpart of the 3D plane fitter used by Wave-MPRAGE.
    It models

        p(kx, y) = a(kx) * y + c(kx)

    and does not fit a z term.  This is appropriate for non-MB 2D multislice
    acquisitions where the last dimension is physical slice index, not a
    Fourier-encoded kz/partition dimension.

    Parameters
    ----------
    psf_diff : torch.Tensor
        Real wrapped phase difference, shape (Nx, Ny).
    hyb_nowave : torch.Tensor
        Hybrid-domain no-wave data whose first two dims are (Nx, Ny). Extra
        dims such as coil are averaged according to mag_average_dims.
    y_norm : sequence or torch.Tensor
        Normalized y coordinates, length Ny.

    Returns
    -------
    dict with keys:
        a_fit_all, c_fit_all, coefficients, psf_diff_pred, mask,
        wrapped_rms, valid_pixels, masked_ratio, skipped
    """
    if psf_diff.ndim != 2:
        raise ValueError(f"psf_diff must be [Nx, Ny], got {tuple(psf_diff.shape)}")
    if not torch.is_floating_point(psf_diff):
        raise TypeError("psf_diff must be a real floating-point tensor")
    if hyb_nowave.shape[:2] != psf_diff.shape:
        raise ValueError(
            f"hyb_nowave first 2 dims must match psf_diff. Got {tuple(hyb_nowave.shape)} vs {tuple(psf_diff.shape)}"
        )

    device = psf_diff.device
    dtype = psf_diff.dtype
    Nx, Ny = psf_diff.shape

    if mag_abs_floor is None:
        mag_abs_floor = mag_threshold if mask_mode == "absolute" else 0.0

    y = torch.as_tensor(y_norm, dtype=dtype, device=device)
    if y.numel() != Ny:
        raise ValueError(f"len(y_norm)={y.numel()} must match Ny={Ny}")

    A_full = torch.stack([y, torch.ones(Ny, dtype=dtype, device=device)], dim=1)

    mag2 = _mean_over_dims(hyb_nowave.abs().square(), mag_average_dims).to(dtype=dtype, device=device)
    if mag2.shape != psf_diff.shape:
        raise ValueError(
            f"After averaging, mag2 shape must equal psf_diff. Got {tuple(mag2.shape)}. "
            "Check mag_average_dims."
        )

    a_fit_all = torch.full((Nx,), torch.nan, dtype=dtype, device=device)
    c_fit_all = torch.full((Nx,), torch.nan, dtype=dtype, device=device)
    psf_diff_pred = torch.full_like(psf_diff, torch.nan)
    wrapped_rms = torch.full((Nx,), torch.nan, dtype=dtype, device=device)
    valid_pixels = torch.zeros((Nx,), dtype=torch.long, device=device)
    masked_ratio = torch.full((Nx,), torch.nan, dtype=dtype, device=device)
    skipped = torch.ones((Nx,), dtype=torch.bool, device=device)
    final_mask_all = torch.zeros((Nx, Ny), dtype=torch.bool, device=device)

    if initial_coefficients is None:
        prev_theta = None
    else:
        prev_theta = torch.as_tensor(initial_coefficients, dtype=dtype, device=device).reshape(2)

    eps = torch.finfo(dtype).eps

    with torch.no_grad():
        for k in range(Nx):
            phase_t = psf_diff[k]
            mag2_t = mag2[k]
            finite = torch.isfinite(phase_t) & torch.isfinite(mag2_t)

            if mask_mode == "none":
                mask = finite
            elif mask_mode == "absolute":
                mask = finite & (mag2_t > mag_threshold)
            elif mask_mode in ("adaptive_quantile", "combined"):
                valid_abs = finite & (mag2_t > float(mag_abs_floor))
                if int(valid_abs.sum().item()) == 0:
                    mask = valid_abs
                else:
                    q = torch.quantile(mag2_t[valid_abs], mag_lower_quantile)
                    mask = valid_abs & (mag2_t >= q)
                    min_count = int(round(min_keep_fraction * Ny))
                    if min_count > 0 and int(mask.sum().item()) < min_count:
                        valid_idx = torch.nonzero(valid_abs, as_tuple=False).flatten()
                        kk = min(min_count, valid_idx.numel())
                        _, order = torch.topk(mag2_t[valid_idx], k=kk, largest=True)
                        mask = torch.zeros_like(valid_abs, dtype=torch.bool)
                        mask[valid_idx[order]] = True
            else:
                raise ValueError("mask_mode must be one of: 'absolute', 'adaptive_quantile', 'combined', 'none'")

            total = mask.numel()
            n_valid = int(mask.sum().item())
            ratio = 1.0 - n_valid / total

            if n_valid < min_valid_pixels or ratio > max_masked_threshold:
                valid_pixels[k] = n_valid
                masked_ratio[k] = ratio
                if verbose:
                    print(f"Skipped {k}/{Nx}: masked_ratio={ratio:.3f}, valid={n_valid}")
                continue

            A = A_full[mask]
            p = phase_t[mask]

            if use_mag_weight:
                w_mag = mag2_t[mask]
                w_mag = w_mag / w_mag.median().clamp_min(eps)
                if weight_clip is not None:
                    w_mag = torch.clamp(w_mag, max=float(weight_clip))
            else:
                w_mag = torch.ones_like(p)

            if use_previous_fit_as_init and prev_theta is not None:
                theta0 = prev_theta.clone()
            elif initial_coefficients is not None:
                theta0 = torch.as_tensor(initial_coefficients, dtype=dtype, device=device).reshape(2).clone()
            else:
                theta0 = _weighted_lstsq(A, p, w_mag, ridge=ridge)

            theta = _run_wrapped_irls(
                A,
                p,
                w_mag,
                theta0,
                n_irls=n_irls,
                huber_delta=huber_delta,
                ridge=ridge,
                convergence_tol=convergence_tol,
            )

            final_mask_all[k] = mask
            a_fit_all[k], c_fit_all[k] = theta
            pred_full = A_full @ theta
            psf_diff_pred[k] = pred_full

            r_final = wrap_to_pi(p - A @ theta)
            wrapped_rms[k] = torch.sqrt(torch.mean(r_final.square()))
            valid_pixels[k] = n_valid
            masked_ratio[k] = ratio
            skipped[k] = False

            if use_previous_fit_as_init:
                prev_theta = theta.clone()

            if verbose:
                print(
                    f"Finished {k}/{Nx}: "
                    f"a={theta[0].item():.4g}, c={theta[1].item():.4g}, "
                    f"wrapped RMS={wrapped_rms[k].item():.4g} rad, valid={n_valid}/{total}"
                )

    coefficients = torch.stack([a_fit_all, c_fit_all], dim=1)
    result: Dict[str, torch.Tensor] = {
        "a_fit_all": a_fit_all,
        "c_fit_all": c_fit_all,
        "coefficients": coefficients,
        "psf_diff_pred": psf_diff_pred,
        "mask": final_mask_all,
        "wrapped_rms": wrapped_rms,
        "valid_pixels": valid_pixels,
        "masked_ratio": masked_ratio,
        "skipped": skipped,
    }
    if return_wrapped_prediction:
        result["psf_diff_pred_wrapped"] = wrap_to_pi(psf_diff_pred)
    return result
