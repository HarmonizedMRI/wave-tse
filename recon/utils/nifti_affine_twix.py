"""Build a NIfTI affine from Siemens TWIX slice geometry.

The functions are copied unchanged from the research module
``experiment_recon/tse_2d_recon_update/nifti_conversion.py`` of pulseq_basic
(commit 6c01ffa).
"""

import numpy as np
import nibabel as nib
import mapvbvd


def _normalize(v, name="vector"):
    """Return a unit-length float vector.

    Parameters
    ----------
    v : array-like
        Vector to normalize.
    name : str
        Human-readable label used in validation errors.

    Returns
    -------
    numpy.ndarray
        Unit-length vector.

    Raises
    ------
    ValueError
        If the input vector has zero norm.
    """
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError(f"{name} has zero norm.")
    return v / n


def _twix_get(yaps, key, default=None):
    """
    Robust MeasYaps getter for mapVBVD-style nested tuple keys.
    """
    try:
        val = yaps.get(key, default)
    except Exception:
        try:
            val = yaps[key]
        except Exception:
            val = default

    if val is None:
        return default

    return val


def _get_twix_scan_and_hdr(twix_file=None, twix_obj=None, scan_index=-1):
    """
    Return scan, hdr, yaps from either a twix object or a twix filename.

    Requires Python mapvbvd only if twix_file is provided.
    """
    if twix_obj is None:
        if twix_file is None:
            raise ValueError("Either twix_file or twix_obj must be provided.")

        import mapvbvd
        twix_obj = mapvbvd.mapVBVD(str(twix_file))

    if isinstance(twix_obj, (list, tuple)):
        scan = twix_obj[scan_index]
    else:
        scan = twix_obj

    try:
        hdr = scan["hdr"]
    except Exception:
        hdr = scan.hdr

    yaps = hdr["MeasYaps"]

    return scan, hdr, yaps


def _sct_to_ras(v_sct, twix_coord_system="LPS"):
    """
    Convert Siemens Sag/Cor/Tra vector to NIfTI RAS.

    twix_coord_system="LPS" assumes:
        +Sag = Left
        +Cor = Posterior
        +Tra = Superior/Head

    Then:
        RAS = [-Sag, -Cor, +Tra]

    If your mapVBVD/YAPS convention has already been converted to RAS,
    use twix_coord_system="RAS".
    """
    v_sct = np.asarray(v_sct, dtype=float)

    if twix_coord_system.upper() == "LPS":
        return np.array([-v_sct[0], -v_sct[1], v_sct[2]], dtype=float)

    elif twix_coord_system.upper() == "RAS":
        return np.array([v_sct[0], v_sct[1], v_sct[2]], dtype=float)

    else:
        raise ValueError("twix_coord_system must be 'LPS' or 'RAS'.")


def _make_inplane_basis_from_normal_ras(
    normal_ras,
    inplane_rot=0.0,
    inplane_rot_sign=1.0,
):
    """
    Build readout, phase, slice directions from slice normal and in-plane rotation.

    Output is right-handed:
        readout x phase = slice

    This is a practical Siemens-style construction:
        - transverse/coronal-like planes use R as the default readout reference
        - sagittal-like planes use A as the default readout reference

    The sign of Siemens dInPlaneRot can differ by convention, so expose
    inplane_rot_sign. If DICOM comparison shows in-plane axes rotated the
    wrong way, change +1 to -1.
    """
    n = _normalize(normal_ras, "slice normal RAS")

    R = np.array([1.0, 0.0, 0.0])
    A = np.array([0.0, 1.0, 0.0])

    dominant_axis = int(np.argmax(np.abs(n)))

    # If normal is mostly R/L, use A as in-plane reference.
    # Otherwise use R as in-plane reference.
    if dominant_axis == 0:
        ref = A
    else:
        ref = R

    read0 = ref - n * np.dot(ref, n)

    if np.linalg.norm(read0) < 1e-8:
        ref = A if np.allclose(ref, R) else R
        read0 = ref - n * np.dot(ref, n)

    read0 = _normalize(read0, "unrotated readout direction")
    phase0 = _normalize(np.cross(n, read0), "unrotated phase direction")

    phi = float(inplane_rot_sign) * float(inplane_rot)

    read = np.cos(phi) * read0 + np.sin(phi) * phase0
    phase = -np.sin(phi) * read0 + np.cos(phi) * phase0

    read = _normalize(read, "readout direction")
    phase = _normalize(phase, "phase direction")

    return read, phase, n


def make_nifti_affine_from_twix(
    twix_file=None,
    twix_obj=None,
    scan_index=-1,
    slice_index=0,
    npy_shape=None,
    twix_array_axis_roles=("readout", "phase", "slice"),
    twix_array_axis_flips=(False, False, False),
    twix_coord_system="LPS",
    twix_inplane_rot_sign=1.0,
    twix_use_fov_for_voxel_size=True,
    voxel_size_mm=None,
    twix_fov_override=None,
    twix_center_sct_override=None,
):
    """
    Build a NIfTI RAS affine from Siemens Twix MeasYaps geometry.

    Parameters
    ----------
    twix_file / twix_obj
        Provide either a filename or an already-loaded mapVBVD object.

    npy_shape
        Shape of the image array being saved.

    twix_array_axis_roles
        Role of each npy array axis.

        Examples:
            If npy shape is [readout, phase, partition]:
                ("readout", "phase", "slice")

            If npy shape is [H/F, A/P, L/R], and for this acquisition
            H/F is partition, A/P is phase, L/R is readout:
                ("slice", "phase", "readout")

    twix_array_axis_flips
        Flip flags relative to the Twix physical directions.
        If L/R is visually flipped, usually toggle the readout-axis flip.

    twix_coord_system
        Usually "LPS" for Siemens Sag/Cor/Tra fields.

    twix_inplane_rot_sign
        Use +1 first. If comparison to DICOM shows in-plane rotation sign
        is wrong, use -1.

    twix_use_fov_for_voxel_size
        If True, voxel sizes are inferred from Twix FOV / npy_shape.
        If False, voxel_size_mm is used directly.

    Returns
    -------
    affine_ras : np.ndarray
    voxel_size_mm_out : tuple
    twix_info : dict
    """
    if npy_shape is None:
        raise ValueError("npy_shape must be provided.")

    if len(twix_array_axis_roles) != 3:
        raise ValueError("twix_array_axis_roles must have length 3.")

    if len(twix_array_axis_flips) != 3:
        raise ValueError("twix_array_axis_flips must have length 3.")

    valid_roles = {"readout", "phase", "slice"}
    for role in twix_array_axis_roles:
        if role not in valid_roles:
            raise ValueError(
                f"Invalid Twix axis role '{role}'. "
                "Use 'readout', 'phase', or 'slice'."
            )

    if len(set(twix_array_axis_roles)) != 3:
        raise ValueError(
            f"twix_array_axis_roles must contain each role once. "
            f"Got {twix_array_axis_roles}"
        )

    scan, hdr, yaps = _get_twix_scan_and_hdr(
        twix_file=twix_file,
        twix_obj=twix_obj,
        scan_index=scan_index,
    )

    Ns = str(slice_index)

    fov = {
        "readout": _twix_get(
            yaps,
            ("sSliceArray", "asSlice", Ns, "dReadoutFOV"),
            None,
        ),
        "phase": _twix_get(
            yaps,
            ("sSliceArray", "asSlice", Ns, "dPhaseFOV"),
            None,
        ),
        "slice": _twix_get(
            yaps,
            ("sSliceArray", "asSlice", Ns, "dThickness"),
            None,
        ),
    }

    if twix_fov_override is not None:
        fov.update(twix_fov_override)

    normal_sct = np.array([
        _twix_get(yaps, ("sSliceArray", "asSlice", Ns, "sNormal", "dSag"), 0.0),
        _twix_get(yaps, ("sSliceArray", "asSlice", Ns, "sNormal", "dCor"), 0.0),
        _twix_get(yaps, ("sSliceArray", "asSlice", Ns, "sNormal", "dTra"), 0.0),
    ], dtype=float)

    position_sct = np.array([
        _twix_get(yaps, ("sSliceArray", "asSlice", Ns, "sPosition", "dSag"), 0.0),
        _twix_get(yaps, ("sSliceArray", "asSlice", Ns, "sPosition", "dCor"), 0.0),
        _twix_get(yaps, ("sSliceArray", "asSlice", Ns, "sPosition", "dTra"), 0.0),
    ], dtype=float)
    header_position_sct = position_sct.copy()
    if twix_center_sct_override is not None:
        position_sct = np.asarray(twix_center_sct_override, dtype=float)
        if position_sct.shape != (3,):
            raise ValueError("twix_center_sct_override must contain Sag/Cor/Tra coordinates.")

    inplane_rot = (
        hdr.get("Meas", {}).get("dInPlaneRot")
        or hdr.get("Protocol", {}).get("dInPlaneRot")
        or _twix_get(yaps, ("sSliceArray", "asSlice", Ns, "dInPlaneRot"), 0.0)
        or 0.0
    )

    normal_ras = _normalize(
        _sct_to_ras(normal_sct, twix_coord_system=twix_coord_system),
        "Twix normal RAS",
    )

    center_ras = _sct_to_ras(
        position_sct,
        twix_coord_system=twix_coord_system,
    )

    read_dir, phase_dir, slice_dir = _make_inplane_basis_from_normal_ras(
        normal_ras,
        inplane_rot=float(inplane_rot),
        inplane_rot_sign=twix_inplane_rot_sign,
    )

    direction_by_role = {
        "readout": read_dir,
        "phase": phase_dir,
        "slice": slice_dir,
    }

    axis_for_role = {
        role: ax for ax, role in enumerate(twix_array_axis_roles)
    }

    if twix_use_fov_for_voxel_size:
        spacing_by_role = {}

        for role in ("readout", "phase", "slice"):
            if fov[role] is None:
                raise ValueError(
                    f"Twix FOV for role '{role}' is missing. "
                    "Either provide twix_fov_override, or set "
                    "twix_use_fov_for_voxel_size=False and provide voxel_size_mm."
                )

            ax = axis_for_role[role]
            spacing_by_role[role] = float(fov[role]) / float(npy_shape[ax])

    else:
        if voxel_size_mm is None:
            raise ValueError(
                "voxel_size_mm must be provided when "
                "twix_use_fov_for_voxel_size=False."
            )

        voxel_size_mm = tuple(float(x) for x in voxel_size_mm)

        spacing_by_role = {}
        for ax, role in enumerate(twix_array_axis_roles):
            spacing_by_role[role] = voxel_size_mm[ax]

    role_vectors = {
        role: direction_by_role[role] * spacing_by_role[role]
        for role in ("readout", "phase", "slice")
    }

    affine = np.eye(4, dtype=float)

    # Start from volume center, then go to voxel [0,0,0].
    origin_ras = center_ras.copy()

    for ax, role in enumerate(twix_array_axis_roles):
        v = role_vectors[role]
        origin_ras -= v * ((npy_shape[ax] - 1) / 2.0)

    affine[:3, 3] = origin_ras

    for ax, role in enumerate(twix_array_axis_roles):
        v = role_vectors[role].copy()

        if twix_array_axis_flips[ax]:
            affine[:3, 3] += v * (npy_shape[ax] - 1)
            v = -v

        affine[:3, ax] = v

    voxel_size_mm_out = tuple(float(x) for x in nib.affines.voxel_sizes(affine))

    twix_info = {
        "FOV": fov,
        "NormalSagCorTra": normal_sct.tolist(),
        "PositionSagCorTra": position_sct.tolist(),
        "HeaderPositionSagCorTra": header_position_sct.tolist(),
        "CenterSagCorTraOverrideApplied": twix_center_sct_override is not None,
        "NormalRAS": normal_ras.tolist(),
        "CenterRAS": center_ras.tolist(),
        "InPlaneRotationRad": float(inplane_rot),
        "InPlaneRotationSignUsed": float(twix_inplane_rot_sign),
        "TwixCoordinateSystemAssumption": twix_coord_system,
        "TwixArrayAxisRoles": list(twix_array_axis_roles),
        "TwixArrayAxisFlips": list(twix_array_axis_flips),
        "ReadoutDirectionRAS": read_dir.tolist(),
        "PhaseDirectionRAS": phase_dir.tolist(),
        "SliceDirectionRAS": slice_dir.tolist(),
    }

    return affine, voxel_size_mm_out, twix_info
