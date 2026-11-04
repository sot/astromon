"""Local Tycho-2 cone-search, proper-motion propagation, and admission gate.

Built from the ESO Science Archive mirror of CDS I/259, replacing astromon's former
live VizieR-based Tycho2 lookup (CDS Strasbourg/VizieR has been unavailable since
2026-10-28, no ETA): :any:`astromon.cross_match`'s ``_get_tycho2_local`` /
``_tycho2_cone_search`` call :func:`query_cone` and :func:`admit_tycho2_match` from
here directly, dispatched via ``LOCAL_CATALOG_GETTERS["Tycho2"]`` exactly like RFC/
ICRF3/ICRF2. :any:`astromon.tests.test_tycho2_local_index` (a held-out regression
against the production DB's prior VizieR-sourced cross-matches) and
:any:`astromon.tests.test_tycho2_hipparcos_crosscheck` (an offline cross-check against
HEASARC's independently-sourced Hipparcos catalog) are what validated this cutover; see
the Part 2 validation campaign notes.

:func:`get_tycho2_local_index` and :func:`get_hipparcos_snapshot` download and build
their respective index files on first use and cache them under
``$ASTROMON_DATA_DIR/tycho2_local_index/`` (default
``$SKA/data/astromon/tycho2_local_index/``), the same download-and-cache-on-first-use
convention :any:`astromon.cross_match` uses for Milliquas, ICRF3, etc. The
Tycho-2 index is a few hundred MB and takes a couple of minutes to build; see
``astromon/scripts/build_tycho2_local_index.py`` to pre-build it ahead of time rather
than paying that cost inside whichever process happens to need it first.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

from astromon import utils

logger = logging.getLogger("astromon")

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
    """Apply single-catalog Tycho2 admission gates (the retired "tycho2" selection).

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


# ---------------------------------------------------------------------------------
# Download-and-cache-on-first-use index building, mirroring astromon.cross_match's
# get_milliquas / get_icrf3 convention (pinned-release cache, rebuilt only when the
# recorded release no longer matches, not on a timer -- both source catalogs are
# static CDS releases).
# ---------------------------------------------------------------------------------

_TYCHO2_CACHE_DIR = utils.ASTROMON_DATA_DIR / "tycho2_local_index"
TYCHO2_INDEX_PATH = _TYCHO2_CACHE_DIR / "tyc2_stars.parquet"
HIPPARCOS_INDEX_PATH = _TYCHO2_CACHE_DIR / "hipparcos_stars.parquet"

_ESO_TYC2_BASE_URL = "https://archive.eso.org/ASTROM/TYC-2/data"
_TYCHO2_RAW_FILES = ("catalog.dat", "suppl_1.dat", "suppl_2.dat")
# CDS Tycho-2 catalog number; the catalog is a fixed historical dataset, so this only
# ever changes if a parsing bug forces a new build of this index.
_TYCHO2_VERSION = "I/259"

_HEASARC_TAP_URL = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
_HEASARC_HIP_COLUMNS = (
    "hip_number, ra_deg, dec_deg, pm_ra, pm_dec, vt_mag, bt_mag, ccdm_id, "
    "ccdm_n_comp, ccdm_history"
)
# HEASARC's TAP sync service silently truncates a single query at 100000 rows; the
# full catalog (118218 rows) is fetched in two halves split on this hip_number.
_HEASARC_HIP_SPLIT = 60000
# CDS Hipparcos catalog number; static like Tycho-2 above.
_HIPPARCOS_VERSION = "I/239"
_HIPPARCOS_EXPECTED_ROWS = 118218


def _release_marker(cache_path: Path) -> Path:
    """Path to the small sidecar file recording which release `cache_path` holds."""
    return cache_path.with_name(cache_path.name + ".release")


def _cache_needs_refresh(cache_path: Path, version: str) -> bool:
    """Whether a locally cached index file should be (re)built.

    Both source catalogs are pinned, static CDS releases, so a mismatch between the
    recorded release and the declared one is what invalidates a cache, not file age
    (see astromon.cross_match._cache_needs_refresh for the same convention applied to
    Milliquas/ICRF3/etc.). A cache with no recorded release is adopted as-is.
    """
    if not cache_path.exists():
        return True
    marker = _release_marker(cache_path)
    if not marker.exists():
        return False
    return marker.read_text().strip() != version


def _record_catalog_version(cache_path: Path, version: str) -> None:
    _release_marker(cache_path).write_text(version)


# --- Tycho-2: ESO Science Archive mirror of CDS I/259 -----------------------------

# (name, start_byte_1based, width, kind) kind in {"i", "f", "a"}; byte layout from the
# ESO ASTROM/TYC-2 format description (https://archive.eso.org/ASTROM/TYC-2/format.html),
# the same layout as CDS I/259.
_TYCHO2_RECORD_LEN = 206 + 2  # 206 data bytes + CRLF
_TYCHO2_CATALOG_FIELDS = [
    ("TYC1", 1, 4, "i"),
    ("TYC2", 6, 5, "i"),
    ("TYC3", 12, 1, "i"),
    ("pflag", 14, 1, "a"),
    ("mRAdeg", 16, 12, "f"),
    ("mDEdeg", 29, 12, "f"),
    ("pmRA", 42, 7, "f"),
    ("pmDE", 50, 7, "f"),
    ("BT", 111, 6, "f"),
    ("VT", 124, 6, "f"),
    ("HIP", 143, 6, "i"),
    ("RAdeg", 153, 12, "f"),
    ("DEdeg", 166, 12, "f"),
    ("epRA", 179, 4, "f"),
]
_TYCHO2_SUPPL_RECORD_LEN = 122 + 2
_TYCHO2_SUPPL_FIELDS = [
    ("TYC1", 1, 4, "i"),
    ("TYC2", 6, 5, "i"),
    ("TYC3", 12, 1, "i"),
    ("flag", 14, 1, "a"),
    ("RAdeg", 16, 12, "f"),
    ("DEdeg", 29, 12, "f"),
    ("pmRA", 42, 7, "f"),
    ("pmDE", 50, 7, "f"),
    ("BT", 84, 6, "f"),
    ("VT", 97, 6, "f"),
    ("HIP", 116, 6, "i"),
]


def _download_tycho2_raw_files(raw_dir: Path) -> None:
    """Download the ESO mirror of CDS I/259 into `raw_dir` (skips files already there).

    ESO independently mirrors the Tycho-2 flat files byte-for-byte identical to the
    CDS release, which matters while CDS Strasbourg/VizieR itself is down (as of
    2026-10-28, no ETA) -- see the Part 2 cutover validation report.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    for name in _TYCHO2_RAW_FILES:
        dest = raw_dir / name
        if dest.exists():
            continue
        url = f"{_ESO_TYC2_BASE_URL}/{name}"
        logger.info(f"Downloading {url} -> {dest}")
        resp = requests.get(url, timeout=600, stream=True)
        resp.raise_for_status()
        tmp = dest.with_name(dest.name + ".tmp")
        with tmp.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
        tmp.rename(dest)


def _read_tycho2_fixed_width(
    path: Path, record_len: int, fields: list[tuple[str, int, int, str]]
) -> pd.DataFrame:
    raw = np.fromfile(path, dtype=np.uint8)
    if raw.size % record_len != 0:
        raise ValueError(
            f"{path}: file size {raw.size} is not a multiple of record length "
            f"{record_len} -- corrupt or incomplete download"
        )
    n_records = raw.size // record_len
    raw = raw.reshape(n_records, record_len)

    columns: dict[str, np.ndarray] = {}
    for name, start, width, kind in fields:
        start0 = start - 1
        chunk = np.ascontiguousarray(raw[:, start0 : start0 + width])
        strs = np.char.decode(chunk.view(f"S{width}").reshape(n_records), "ascii")
        if kind == "a":
            columns[name] = np.char.strip(strs)
        else:
            stripped = np.char.strip(strs)
            values = np.full(n_records, np.nan, dtype=np.float64)
            non_blank = stripped != ""
            values[non_blank] = stripped[non_blank].astype(np.float64)
            columns[name] = values
    return pd.DataFrame(columns)


def _build_tycho2_star_table(raw_dir: Path) -> pd.DataFrame:
    """Parse the ESO/CDS I/259 flat files into one combined per-star table.

    Columns: ``tyc_name``, ``source`` ("main"/"suppl1"/"suppl2"), ``ra0``/``dec0``
    (anchor position, deg), ``anchor_epoch`` (fractional year), ``pmRA``/``pmDE``
    (mas/yr, NaN if unavailable), ``has_pm``, ``VT``/``BT``, ``HIP``.
    """
    catalog = _read_tycho2_fixed_width(
        raw_dir / "catalog.dat", _TYCHO2_RECORD_LEN, _TYCHO2_CATALOG_FIELDS
    )
    catalog["TYC1"] = catalog["TYC1"].astype(np.int32)
    catalog["TYC2"] = catalog["TYC2"].astype(np.int32)
    catalog["TYC3"] = catalog["TYC3"].astype(np.int16)

    # a blank pflag (the normal case: full astrometric + proper-motion solution) is
    # represented as a single space, matching the printed catalog's fixed-width column
    catalog["pflag"] = catalog["pflag"].replace("", " ")
    has_mean_pos = catalog["pflag"] != "X"
    ra0 = np.where(has_mean_pos, catalog["mRAdeg"], catalog["RAdeg"])
    dec0 = np.where(has_mean_pos, catalog["mDEdeg"], catalog["DEdeg"])
    anchor_epoch = np.where(
        has_mean_pos,
        2000.0,
        np.where(~np.isnan(catalog["epRA"]), 1990.0 + catalog["epRA"], 1990.0),
    )
    pmRA = np.where(has_mean_pos, catalog["pmRA"], np.nan)
    pmDE = np.where(has_mean_pos, catalog["pmDE"], np.nan)
    has_pm = has_mean_pos & ~np.isnan(catalog["pmRA"]) & ~np.isnan(catalog["pmDE"])

    tables = [
        pd.DataFrame(
            {
                "TYC1": catalog["TYC1"],
                "TYC2": catalog["TYC2"],
                "TYC3": catalog["TYC3"],
                "source": "main",
                "pflag": catalog["pflag"],
                "ra0": ra0,
                "dec0": dec0,
                "anchor_epoch": anchor_epoch,
                "pmRA": pmRA,
                "pmDE": pmDE,
                "has_pm": has_pm,
                "VT": catalog["VT"],
                "BT": catalog["BT"],
                "HIP": catalog["HIP"],
            }
        )
    ]
    for name, source_name in [("suppl_1.dat", "suppl1"), ("suppl_2.dat", "suppl2")]:
        suppl = _read_tycho2_fixed_width(
            raw_dir / name, _TYCHO2_SUPPL_RECORD_LEN, _TYCHO2_SUPPL_FIELDS
        )
        suppl["TYC1"] = suppl["TYC1"].astype(np.int32)
        suppl["TYC2"] = suppl["TYC2"].astype(np.int32)
        suppl["TYC3"] = suppl["TYC3"].astype(np.int16)
        has_pm_s = suppl["flag"] == "H"
        tables.append(
            pd.DataFrame(
                {
                    "TYC1": suppl["TYC1"],
                    "TYC2": suppl["TYC2"],
                    "TYC3": suppl["TYC3"],
                    "source": source_name,
                    "pflag": suppl["flag"],
                    "ra0": suppl["RAdeg"],
                    "dec0": suppl["DEdeg"],
                    # supplement positions are quoted at the fixed epoch J1991.25,
                    # unlike the main catalog's per-star mean epoch above.
                    "anchor_epoch": 1991.25,
                    "pmRA": np.where(has_pm_s, suppl["pmRA"], np.nan),
                    "pmDE": np.where(has_pm_s, suppl["pmDE"], np.nan),
                    "has_pm": (
                        has_pm_s & ~np.isnan(suppl["pmRA"]) & ~np.isnan(suppl["pmDE"])
                    ),
                    "VT": suppl["VT"],
                    "BT": suppl["BT"],
                    "HIP": suppl["HIP"],
                }
            )
        )

    combined = pd.concat(tables, ignore_index=True)
    no_position = np.isnan(combined["ra0"]) | np.isnan(combined["dec0"])
    if no_position.any():
        logger.warning(
            f"dropping {int(no_position.sum())} Tycho-2 rows with no usable position"
        )
        combined = combined[~no_position].reset_index(drop=True)
    combined["tyc_name"] = (
        combined["TYC1"].astype(str)
        + "-"
        + combined["TYC2"].astype(str)
        + "-"
        + combined["TYC3"].astype(str)
    )
    return combined


def _build_tycho2_index() -> None:
    import tempfile  # noqa: PLC0415 — only needed for this one-time build

    TYCHO2_INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="tycho2_raw_") as raw_dir_name:
        raw_dir = Path(raw_dir_name)
        _download_tycho2_raw_files(raw_dir)
        stars = _build_tycho2_star_table(raw_dir)
    stars.to_parquet(TYCHO2_INDEX_PATH, index=False)
    _record_catalog_version(TYCHO2_INDEX_PATH, _TYCHO2_VERSION)
    logger.info(f"Built {len(stars)}-row Tycho-2 local index at {TYCHO2_INDEX_PATH}")


def get_tycho2_local_index(force_rebuild: bool = False) -> Path:
    """Return the path to a built local Tycho-2 index, building it first if needed.

    Downloads the ESO mirror of CDS I/259 and parses it into
    ``$ASTROMON_DATA_DIR/tycho2_local_index/tyc2_stars.parquet`` on first use (a
    few hundred MB, taking a couple of minutes); later calls reuse that file until
    ``force_rebuild=True`` or the recorded release no longer matches
    :data:`_TYCHO2_VERSION`. Pass the returned path to :func:`load_star_table`.
    """
    if force_rebuild or _cache_needs_refresh(TYCHO2_INDEX_PATH, _TYCHO2_VERSION):
        _build_tycho2_index()
    return TYCHO2_INDEX_PATH


# --- Hipparcos: HEASARC's independently-sourced snapshot of CDS I/239 -------------


def _query_heasarc_hipparcos_half(hip_min: int | None, hip_max: int | None) -> str:
    where = []
    if hip_min is not None:
        where.append(f"hip_number > {hip_min}")
    if hip_max is not None:
        where.append(f"hip_number <= {hip_max}")
    query = f"SELECT {_HEASARC_HIP_COLUMNS} FROM hipparcos"
    if where:
        query += " WHERE " + " AND ".join(where)
    resp = requests.get(
        _HEASARC_TAP_URL,
        params={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "text", "QUERY": query},
        timeout=300,
    )
    resp.raise_for_status()
    if "OVERFLOW" in resp.text:
        raise RuntimeError(
            "HEASARC TAP query overflowed its row cap -- the hip_number split needs "
            "to be narrower"
        )
    return resp.text


def _parse_heasarc_hipparcos_text(text: str) -> pd.DataFrame:
    import io  # noqa: PLC0415 — only needed for this one-time parse

    lines = text.splitlines()
    # the data rows are followed by a blank line and then a "Number of rows: ..."
    # footer (an overflowed query instead gets a "---- Messages ----" block, but
    # that is already raised on by _query_heasarc_hipparcos_half before this is
    # ever called)
    end = next(i for i, line in enumerate(lines) if line.strip() == "")
    body = "\n".join(lines[:end])
    df = pd.read_csv(
        io.StringIO(body),
        sep="|",
        skiprows=1,
        engine="python",
        names=[
            "hip_number",
            "ra_deg",
            "dec_deg",
            "pm_ra",
            "pm_dec",
            "vt_mag",
            "bt_mag",
            "ccdm_id",
            "ccdm_n_comp",
            "ccdm_history",
        ],
    )
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip().replace("null", pd.NA)
    numeric_columns = [
        "hip_number",
        "ra_deg",
        "dec_deg",
        "pm_ra",
        "pm_dec",
        "vt_mag",
        "bt_mag",
        "ccdm_n_comp",
    ]
    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _build_hipparcos_snapshot() -> None:
    """Fetch and write a full snapshot of HEASARC's `hipparcos` table (CDS I/239).

    Independent of CDS/VizieR -- see the Part 2 cutover validation report -- and
    written to ``$ASTROMON_DATA_DIR/tycho2_local_index/hipparcos_stars.parquet``.
    """
    HIPPARCOS_INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Downloading Hipparcos Main Catalog snapshot from {_HEASARC_TAP_URL}")
    half1 = _parse_heasarc_hipparcos_text(
        _query_heasarc_hipparcos_half(None, _HEASARC_HIP_SPLIT)
    )
    half2 = _parse_heasarc_hipparcos_text(
        _query_heasarc_hipparcos_half(_HEASARC_HIP_SPLIT, None)
    )
    stars = pd.concat([half1, half2], ignore_index=True)
    if len(stars) != _HIPPARCOS_EXPECTED_ROWS:
        raise RuntimeError(
            f"expected {_HIPPARCOS_EXPECTED_ROWS} Hipparcos rows, got {len(stars)}"
        )

    stars["hip_number"] = stars["hip_number"].astype("int64")
    stars["ccdm_n_comp"] = stars["ccdm_n_comp"].fillna(1).astype("int64")
    stars = stars.rename(
        columns={
            "hip_number": "HIP",
            "ra_deg": "ra",
            "dec_deg": "dec",
            "pm_ra": "pmRA",
            "pm_dec": "pmDE",
            "vt_mag": "VT",
            "bt_mag": "BT",
        }
    )
    stars.to_parquet(HIPPARCOS_INDEX_PATH, index=False)
    _record_catalog_version(HIPPARCOS_INDEX_PATH, _HIPPARCOS_VERSION)
    logger.info(f"Built {len(stars)}-row Hipparcos snapshot at {HIPPARCOS_INDEX_PATH}")


def get_hipparcos_snapshot(force_rebuild: bool = False) -> Path:
    """Return the path to a built local Hipparcos snapshot, building it first if needed.

    See :func:`get_tycho2_local_index` for the caching convention; this one is a
    single bulk HEASARC query (a few tens of MB, well under a minute) rather than a
    multi-hundred-MB file parse.
    """
    if force_rebuild or _cache_needs_refresh(HIPPARCOS_INDEX_PATH, _HIPPARCOS_VERSION):
        _build_hipparcos_snapshot()
    return HIPPARCOS_INDEX_PATH
