"""Standalone local Tycho-2 cone-search, proper-motion propagation, and admission gate.

Used to validate whether a locally-built Tycho-2 index (from the ESO Science Archive
mirror of CDS I/259) can replace astromon's live VizieR-based Tycho2 lookup in
:any:`astromon.cross_match`.

.. Important::

    This module is **not** wired into :any:`astromon.cross_match` or any other
    production code path. It exists purely to support the held-out regression test in
    :any:`astromon.tests.test_tycho2_local_index` and the companion fixture/validation
    scripts, while the live-VizieR-vs-local-index cutover decision is still pending (see
    the Part 2 validation campaign notes). Nothing here is imported by production code.

The built index itself (a parquet file of ~2.5 million parsed Tycho-2 + supplement
stars, built by a separate offline parsing script from the ESO flat files) is a few
hundred MB and is not committed to this repository; its location is given by the
``ASTROMON_TYCHO2_INDEX`` environment variable.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

ARCSEC_PER_DEG = 3600.0
MAS_PER_DEG = 3_600_000.0


def load_star_table(parquet_path: str) -> pd.DataFrame:
    """Load the parsed Tycho-2 star table (see the offline parsing script)."""
    return pd.read_parquet(parquet_path)


def build_kdtree(stars: pd.DataFrame) -> cKDTree:
    """Build a cKDTree over unit Cartesian vectors of each star's anchor position.

    The anchor position is ``(ra0, dec0)``: the mean J2000 position for ordinary
    main-catalog stars, or the observed/J1991.25 position for main-catalog
    ``pflag='X'`` stars and supplement-file stars (which have no mean J2000 position).
    """
    ra_rad = np.radians(stars["ra0"].to_numpy())
    dec_rad = np.radians(stars["dec0"].to_numpy())
    cos_dec = np.cos(dec_rad)
    xyz = np.column_stack(
        [
            cos_dec * np.cos(ra_rad),
            cos_dec * np.sin(ra_rad),
            np.sin(dec_rad),
        ]
    )
    return cKDTree(xyz)


def max_pm_mas_per_yr(stars: pd.DataFrame) -> float:
    """Maximum |proper motion| (mas/yr) over all stars that have one.

    Used to size the kd-tree query padding (see :any:`query_cone`).
    """
    proper_motion = np.hypot(stars["pmRA"].to_numpy(), stars["pmDE"].to_numpy())
    return float(np.nanmax(proper_motion))


def propagate_positions(
    stars: pd.DataFrame, frac_year: float
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate each star's anchor position to a fractional-year observation epoch.

    ``frac_year`` is a fractional Julian year (e.g. from ``CxoTime(...).frac_year``).
    Reproduces VizieR's
    ``_RAJ2000(J2000,Jepoch)`` / ``_DEJ2000(J2000,Jepoch)`` semantics as empirically
    calibrated against the stored production cross-matches (see the Part 2 validation
    report): stars with a published proper motion are propagated linearly from their
    anchor epoch using ``pmRA*`` (already ``cos(dec)``-scaled) and ``pmDE``, both in
    mas/yr; stars without a usable proper motion keep their anchor position unchanged at
    any query epoch.
    """
    ra0 = stars["ra0"].to_numpy()
    dec0 = stars["dec0"].to_numpy()
    has_pm = stars["has_pm"].to_numpy()
    dt_years = frac_year - stars["anchor_epoch"].to_numpy()

    dec_rad = np.radians(dec0)
    proper_motion_ra = stars["pmRA"].to_numpy()
    proper_motion_dec = stars["pmDE"].to_numpy()
    with np.errstate(invalid="ignore"):
        delta_ra = (proper_motion_ra * dt_years) / (MAS_PER_DEG * np.cos(dec_rad))
        delta_dec = (proper_motion_dec * dt_years) / MAS_PER_DEG

    ra = np.where(has_pm, ra0 + delta_ra, ra0)
    dec = np.where(has_pm, dec0 + delta_dec, dec0)
    return ra, dec


def angular_separation_arcsec(
    ra0: float, dec0: float, ra1: np.ndarray, dec1: np.ndarray
) -> np.ndarray:
    """Haversine angular separation in arcsec (accurate to sub-mas at these scales)."""
    ra0_rad, dec0_rad = np.radians(ra0), np.radians(dec0)
    ra1_rad, dec1_rad = np.radians(ra1), np.radians(dec1)
    delta_ra = ra1_rad - ra0_rad
    delta_dec = dec1_rad - dec0_rad
    haversine = (
        np.sin(delta_dec / 2) ** 2
        + np.cos(dec0_rad) * np.cos(dec1_rad) * np.sin(delta_ra / 2) ** 2
    )
    central_angle = 2 * np.arcsin(np.sqrt(np.clip(haversine, 0, 1)))
    return np.degrees(central_angle) * ARCSEC_PER_DEG


def query_cone(
    tree: cKDTree,
    stars: pd.DataFrame,
    ra: float,
    dec: float,
    frac_year: float,
    radius_arcsec: float = 3.0,
    global_max_pm_mas_per_yr: float | None = None,
) -> pd.DataFrame:
    """Cone search around ``(ra, dec)`` [deg] at observation epoch ``frac_year``.

    Returns all Tycho-2 stars (rows of ``stars``, plus ``ra_t``/``dec_t``/
    ``separation_arcsec`` columns) whose *propagated* position at ``frac_year`` is
    within ``radius_arcsec`` of ``(ra, dec)``.

    The kd-tree broad-phase search (on *anchor*-epoch positions) is padded to account
    for the maximum possible drift of any star between its anchor epoch and
    ``frac_year``, so high-proper-motion stars are not missed near the cone edge.
    """
    if global_max_pm_mas_per_yr is None:
        global_max_pm_mas_per_yr = max_pm_mas_per_yr(stars)
    global_max_pm_arcsec_per_yr = global_max_pm_mas_per_yr / 1000.0

    # anchor epochs are 2000.0, 1991.25, or an observed ~1990-1992 epoch; a conservative
    # max |dt| of 40 years safely covers the observation epoch range astromon cares
    # about (Chandra launched in 1999) regardless of the true anchor epoch distribution.
    max_dt_years = 40.0
    pad_arcsec = global_max_pm_arcsec_per_yr * max_dt_years
    broad_radius_arcsec = radius_arcsec + pad_arcsec

    ra_rad, dec_rad = np.radians(ra), np.radians(dec)
    cos_dec = np.cos(dec_rad)
    xyz = np.array(
        [cos_dec * np.cos(ra_rad), cos_dec * np.sin(ra_rad), np.sin(dec_rad)]
    )
    # chord length on the unit sphere corresponding to the padded angular radius
    chord = 2 * np.sin(np.radians(broad_radius_arcsec / ARCSEC_PER_DEG) / 2)
    idx = tree.query_ball_point(xyz, r=chord)
    if not idx:
        return stars.iloc[[]].assign(ra_t=[], dec_t=[], separation_arcsec=[])

    candidates = stars.iloc[idx].reset_index(drop=True)
    ra_t, dec_t = propagate_positions(candidates, frac_year)
    separation = angular_separation_arcsec(ra, dec, ra_t, dec_t)
    candidates = candidates.assign(ra_t=ra_t, dec_t=dec_t, separation_arcsec=separation)
    candidates = candidates[candidates["separation_arcsec"] <= radius_arcsec]
    return candidates.reset_index(drop=True)


def add_yagzag(candidates: pd.DataFrame, quat: Quat) -> pd.DataFrame:
    """Add ``y_angle``/``z_angle`` columns (tangent-plane projection) to candidates.

    Computed from the candidates' propagated ``ra_t``/``dec_t`` columns through
    ``quat``. This mirrors how production computes ``y_angle``/``z_angle`` for
    rough-match candidates (``radec_to_yagzag(ra, dec, quat)`` in
    ``astromon/scripts/get_cat_obs_data.py``).
    """
    candidates = candidates.copy()
    y_angle, z_angle = radec_to_yagzag(
        candidates["ra_t"].to_numpy(), candidates["dec_t"].to_numpy(), quat
    )
    candidates["y_angle"] = y_angle
    candidates["z_angle"] = z_angle
    return candidates


def admit_tycho2_match(
    candidates: pd.DataFrame,
    x_snr: float,
    x_near_neighbor_dist: float,
    x_r_angle: float,
    x_grating: str,
    x_y_angle: float,
    x_z_angle: float,
    *,
    snr: float = 3,
    r_angle: float = 120.0,
    r_angle_grating: float = 120.0,
    near_neighbor_dist: float = 6.0,
    dr: float = 3.0,
) -> tuple[pd.DataFrame, pd.Series | None]:
    """Apply the "tycho2" select_name admission gates to rough-match candidates.

    Mirrors ``simple_cross_match`` with ``catalogs=["Tycho2"]``, ``snr=3``,
    ``r_angle=120``, ``r_angle_grating=120``, ``near_neighbor_dist=6.0``, ``dr=3.0``
    (see :any:`astromon.cross_match`), applied to a set of rough-match candidates (as
    returned by :any:`query_cone`, with ``y_angle``/``z_angle`` already added via
    :any:`add_yagzag`) for a single x-ray source.

    Returns ``(candidates_with_dr, admitted_row_or_None)``.
    """
    candidates = candidates.copy()
    candidates["dy"] = x_y_angle - candidates["y_angle"]
    candidates["dz"] = x_z_angle - candidates["z_angle"]
    candidates["dr"] = np.sqrt(candidates["dy"] ** 2 + candidates["dz"] ** 2)

    r_angle_limit = r_angle_grating if x_grating != "NONE" else r_angle
    xray_source_ok = (
        (x_snr > snr)
        and (x_near_neighbor_dist > near_neighbor_dist)
        and (x_r_angle < r_angle_limit)
    )
    if not xray_source_ok:
        return candidates, None

    passing = candidates[candidates["dr"] < dr]
    if len(passing) == 0:
        return candidates, None

    best = passing.loc[passing["dr"].idxmin()]
    return candidates, best


def fit_pointing_quat(
    ra_nom: float,
    dec_nom: float,
    roll_nom: float,
    src_ra: np.ndarray,
    src_dec: np.ndarray,
    src_y_angle: np.ndarray,
    src_z_angle: np.ndarray,
) -> Quat:
    """Recover the pointing quaternion production used for one obsid's x-ray sources.

    ``astromon_obs`` only stores the *nominal* attitude (``ra``/``dec``/``roll``, i.e.
    ``ra_nom``/``dec_nom``/``roll_nom`` from the obspar), but production's
    ``radec_to_yagzag`` calls use the *pointing* quaternion (``ra_pnt``/``dec_pnt``/
    ``roll_pnt``), which can differ from nominal by a small dither/offset amount and is
    not stored in the production h5 tables. Rather than re-deriving it from scratch, fit
    a small ``(dra, ddec, droll)`` perturbation to the nominal attitude that reproduces
    the already-stored x-ray source ``y_angle``/``z_angle`` to high precision, using all
    x-ray sources in the obsid as constraints. With 2 or more sources this is well
    constrained (including roll); with exactly 1 source, roll is left at its nominal
    value (degenerate: a single 2-component constraint cannot separate roll from
    ra/dec).
    """
    n_sources = len(src_ra)
    fit_roll = n_sources >= 2

    def residuals(params: np.ndarray) -> np.ndarray:
        dra, ddec = params[0], params[1]
        droll = params[2] if fit_roll else 0.0
        quat = Quat(equatorial=(ra_nom + dra, dec_nom + ddec, roll_nom + droll))
        y_angle, z_angle = radec_to_yagzag(src_ra, src_dec, quat)
        return np.concatenate([y_angle - src_y_angle, z_angle - src_z_angle])

    x0 = np.zeros(3 if fit_roll else 2)
    method = "lm" if 2 * n_sources >= len(x0) else "trf"
    result = least_squares(residuals, x0, method=method, xtol=1e-14, ftol=1e-14)
    dra, ddec = result.x[0], result.x[1]
    droll = result.x[2] if fit_roll else 0.0
    return Quat(equatorial=(ra_nom + dra, dec_nom + ddec, roll_nom + droll))
