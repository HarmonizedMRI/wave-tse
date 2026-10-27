"""Center-pad slice-wise coil sensitivities onto the oversampled readout grid.

The function is copied unchanged from the research module
``experiment_recon/tse_2d_recon_update/mb2_reconstruction.py`` of
pulseq_basic (commit 6c01ffa).
"""

import torch


def build_oversampled_sens(csm, nx_os):
    """Center-pad physical-slice CSMs along readout.

    Parameters
    ----------
    csm : numpy.ndarray
        Complex CSM shaped ``(coil, x, y, physical_slice)``.
    nx_os : int
        Target oversampled readout size.

    Returns
    -------
    torch.Tensor
        Center-padded complex64 sensitivities.

    Raises
    ------
    ValueError
        If the target readout is smaller than the input readout.
    """
    ncc, nx, ny, nslices = csm.shape
    if nx_os < nx:
        raise ValueError(f"Nx_os={nx_os} is smaller than CSM Nx={nx}.")
    sens = torch.zeros((ncc, nx_os, ny, nslices), dtype=torch.complex64)
    x0 = (nx_os - nx) // 2
    sens[:, x0 : x0 + nx] = torch.from_numpy(csm).contiguous()
    return sens
