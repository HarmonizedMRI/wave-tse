"""Two-dimensional MB2 Wave forward, adjoint, and PCG-SENSE operators.

For one SMS group, the forward model is

``M * sum_b D_b * F_y(P_b * F_x(S_b * x_b))``.

``P_b`` is a hybrid-space Wave PSF.  The module never transforms across the
physical-slice/band dimension; in particular, it performs no ``kz`` FFT.
The established no-Wave MB2 implementation remains in ``mb2_cg_sense.py``.
"""

import torch

from utils.wave_cg_sense_precondition import fftc_dim, ifftc_dim


def mb2_wave_forward(x, sens, caipi_phase, wave_psf, mask):
    """Apply the 2D MB2 Wave encoding operator.

    Parameters
    ----------
    x : torch.Tensor
        Two complex band images shaped ``(2, Nx_os, Ny)``.
    sens : torch.Tensor
        Receive sensitivities shaped ``(2, Ncoil, Nx_os, Ny)``.
    caipi_phase : torch.Tensor
        Unit-magnitude slice phase shaped ``(2, Ny)``.
    wave_psf : torch.Tensor
        Band-specific hybrid-space Wave PSF shaped ``(2, Nx_os, Ny)``.
    mask : torch.Tensor
        Sampling mask broadcastable to ``(Ncoil, Nx_os, Ny)``.

    Returns
    -------
    torch.Tensor
        Collapsed sampled coil k-space shaped ``(Ncoil, Nx_os, Ny)``.
    """
    _validate_shapes(x, sens, caipi_phase, wave_psf)
    image_coils = sens * x.unsqueeze(1)
    hybrid = fftc_dim(image_coils, dim=2)
    hybrid = hybrid * wave_psf[:, None, :, :]
    kspace_bands = fftc_dim(hybrid, dim=3)
    kspace_bands = kspace_bands * caipi_phase[:, None, None, :]
    return kspace_bands.sum(dim=0) * mask


def mb2_wave_adjoint(y, sens, caipi_phase, wave_psf, mask):
    """Apply the Hermitian adjoint of :func:`mb2_wave_forward`.

    Parameters
    ----------
    y : torch.Tensor
        Collapsed coil k-space shaped ``(Ncoil, Nx_os, Ny)``.
    sens, caipi_phase, wave_psf, mask : torch.Tensor
        Encoding terms with the shapes documented by
        :func:`mb2_wave_forward`.

    Returns
    -------
    torch.Tensor
        Two adjoint band images shaped ``(2, Nx_os, Ny)``.
    """
    if y.ndim != 3:
        raise ValueError(f"y must be (Ncoil, Nx_os, Ny), got {tuple(y.shape)}")
    _validate_encoding_shapes(sens, caipi_phase, wave_psf)
    if y.shape != sens.shape[1:]:
        raise ValueError(
            f"y shape {tuple(y.shape)} does not match sens {tuple(sens.shape[1:])}"
        )

    kspace_bands = (y * mask).unsqueeze(0)
    kspace_bands = kspace_bands * torch.conj(caipi_phase)[:, None, None, :]
    hybrid = ifftc_dim(kspace_bands, dim=3)
    hybrid = hybrid * torch.conj(wave_psf)[:, None, :, :]
    image_coils = ifftc_dim(hybrid, dim=2)
    return (torch.conj(sens) * image_coils).sum(dim=1)


def cg_sense_mb2_wave_2d(
    y,
    sens,
    caipi_phase,
    wave_psf,
    mask,
    n_iter=50,
    tol=1e-6,
    lambda_l2=0.0,
    use_preconditioner=True,
    verbose=True,
    return_info=False,
):
    """Solve the regularized MB2 Wave normal equations with PCG.

    Parameters
    ----------
    y, sens, caipi_phase, wave_psf, mask : torch.Tensor
        Measured data and encoding terms documented by
        :func:`mb2_wave_forward`.
    n_iter : int
        Maximum number of PCG iterations.
    tol : float
        Relative normal-equation residual tolerance.
    lambda_l2 : float
        Nonnegative L2 regularization weight.
    use_preconditioner : bool
        Apply the diagonal coil-power preconditioner when true.
    verbose : bool
        Print per-iteration convergence information when true.
    return_info : bool
        Return solver diagnostics together with the image when true.

    Returns
    -------
    torch.Tensor or tuple[torch.Tensor, dict]
        Reconstructed bands, optionally followed by convergence diagnostics.

    Raises
    ------
    ValueError
        If shapes or the regularization weight are invalid.
    RuntimeError
        If PCG encounters a nonpositive or nonfinite quadratic form.
    """
    if y.ndim != 3:
        raise ValueError(f"y must be (Ncoil, Nx_os, Ny), got {tuple(y.shape)}")
    _validate_encoding_shapes(sens, caipi_phase, wave_psf)
    if y.shape != sens.shape[1:]:
        raise ValueError(
            f"y shape {tuple(y.shape)} does not match sens {tuple(sens.shape[1:])}"
        )
    if lambda_l2 < 0:
        raise ValueError("lambda_l2 must be nonnegative.")

    def normal_operator(image):
        """Apply the regularized normal operator to two band images."""
        result = mb2_wave_adjoint(
            mb2_wave_forward(image, sens, caipi_phase, wave_psf, mask),
            sens,
            caipi_phase,
            wave_psf,
            mask,
        )
        if lambda_l2:
            result = result + lambda_l2 * image
        return result

    b = mb2_wave_adjoint(y, sens, caipi_phase, wave_psf, mask)
    x = torch.zeros_like(b)
    residual = b.clone()
    b_power = torch.vdot(b.reshape(-1), b.reshape(-1)).real + 1e-24

    sensitivity_power = torch.sum(torch.abs(sens) ** 2, dim=1).real
    stabilizer = 1e-8 * sensitivity_power.max()

    def precondition(value):
        """Apply the optional diagonal image-space preconditioner."""
        if not use_preconditioner:
            return value
        return value / (
            sensitivity_power + lambda_l2 + stabilizer + 1e-24
        )

    z = precondition(residual)
    direction = z.clone()
    residual_z = torch.vdot(residual.reshape(-1), z.reshape(-1)).real
    residual_history = [
        float(
            torch.sqrt(
                torch.vdot(residual.reshape(-1), residual.reshape(-1)).real
                / b_power
            )
        )
    ]

    converged = False
    iterations = 0
    for iteration in range(int(n_iter)):
        normal_direction = normal_operator(direction)
        curvature = torch.vdot(
            direction.reshape(-1), normal_direction.reshape(-1)
        ).real
        if not torch.isfinite(curvature) or curvature <= 0:
            raise RuntimeError(
                "PCG encountered invalid p^H A p="
                f"{curvature.item()} at iteration {iteration + 1}."
            )

        alpha = residual_z / curvature
        x = x + alpha * direction
        residual = residual - alpha * normal_direction
        residual_power = torch.vdot(
            residual.reshape(-1), residual.reshape(-1)
        ).real
        relative_residual = torch.sqrt(residual_power / b_power)
        residual_history.append(float(relative_residual))
        iterations = iteration + 1

        if verbose:
            print(
                f"MB2 Wave PCG {iterations:02d}/{n_iter}: relative residual "
                f"{relative_residual.item():.3e}"
            )
        if relative_residual < tol:
            converged = True
            break

        z = precondition(residual)
        residual_z_new = torch.vdot(
            residual.reshape(-1), z.reshape(-1)
        ).real
        beta = residual_z_new / (residual_z + 1e-24)
        direction = z + beta * direction
        residual_z = residual_z_new

    info = {
        "iterations": iterations,
        "converged": converged,
        "relative_residual": residual_history[-1],
        "residual_history": residual_history,
    }
    return (x, info) if return_info else x


def _validate_encoding_shapes(sens, caipi_phase, wave_psf):
    """Validate non-image encoding arrays for the MB2 Wave operator."""
    if sens.ndim != 4 or sens.shape[0] != 2:
        raise ValueError(
            f"sens must be (2, Ncoil, Nx_os, Ny), got {tuple(sens.shape)}"
        )
    expected_phase = (2, sens.shape[-1])
    if caipi_phase.shape != expected_phase:
        raise ValueError(
            f"caipi_phase must be {expected_phase}, got {tuple(caipi_phase.shape)}"
        )
    expected_psf = (2, sens.shape[-2], sens.shape[-1])
    if wave_psf.shape != expected_psf:
        raise ValueError(
            f"wave_psf must be {expected_psf}, got {tuple(wave_psf.shape)}"
        )


def _validate_shapes(x, sens, caipi_phase, wave_psf):
    """Validate image and encoding arrays for the MB2 Wave operator."""
    _validate_encoding_shapes(sens, caipi_phase, wave_psf)
    if x.ndim != 3 or x.shape[0] != 2:
        raise ValueError(f"x must be (2, Nx_os, Ny), got {tuple(x.shape)}")
    if x.shape[1:] != sens.shape[2:]:
        raise ValueError(
            f"x spatial shape {tuple(x.shape[1:])} does not match "
            f"sens {tuple(sens.shape[2:])}"
        )
