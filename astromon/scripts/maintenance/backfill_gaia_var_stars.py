"""Backfill GaiaVarStar catalog matches for all obsids in astromon.h5.

Downloads the full Gaia DR3 vari_rotation_modulation catalog (RUWE-filtered) once,
then matches locally against existing xray sources with per-obsid proper motion
correction.  This avoids one TAP round-trip per obsid.

Usage:
    python -m astromon.scripts.maintenance.backfill_gaia_var_stars --db astromon.h5
    python -m astromon.scripts.maintenance.backfill_gaia_var_stars \\
        --db astromon.h5 --catalog-cache gaia_var_stars.ecsv
    python -m astromon.scripts.maintenance.backfill_gaia_var_stars --db astromon.h5 --dry-run
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from astropy import coordinates as coords
from astropy.table import Table, vstack
from astroquery.gaia import Gaia
from chandra_aca.transform import radec_to_yagzag
from cxotime import CxoTime
from Quaternion import Quat

from astromon import cross_match, db

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("backfill_gaia_var_stars")

# Generous radius for pre-filtering the catalog to each obsid's field.
# ACIS-I is ~16 arcmin across; HRC-I is ~30 arcmin. 0.35 deg covers both.
_FOV_RADIUS_DEG = 0.35

_TAP_QUERY = """
SELECT g.source_id, g.ra, g.dec, g.pmra, g.pmdec, g.phot_g_mean_mag
FROM gaiadr3.gaia_source g
INNER JOIN gaiadr3.vari_rotation_modulation v ON g.source_id = v.source_id
WHERE g.ruwe < {max_ruwe}
"""


def download_catalog(max_ruwe: float = 1.4, cache_path: Path | None = None) -> Table:
    """Download the full vari_rotation_modulation catalog from Gaia TAP.

    Parameters
    ----------
    max_ruwe
        RUWE threshold applied in the TAP query; sources above this are excluded.
    cache_path
        If given and the file exists, load from cache instead of querying.
        If given and the file does not exist, save the result there for next time.

    Returns
    -------
    Table with columns: source_id, ra, dec, pmra, pmdec, phot_g_mean_mag
    """
    if cache_path is not None and Path(cache_path).exists():
        logger.info(f"Loading cached catalog from {cache_path}")
        cat = Table.read(cache_path)
        logger.info(f"  {len(cat):,} sources")
        return cat

    logger.info(
        "Downloading full GaiaVarStar catalog via TAP (may take several minutes)…"
    )
    job = Gaia.launch_job_async(_TAP_QUERY.format(max_ruwe=max_ruwe), verbose=False)
    cat = Table(job.get_results())
    logger.info(f"  Downloaded {len(cat):,} sources")

    # Fill missing proper motions with 0 so downstream code needn't handle NaN.
    for col in ("pmra", "pmdec"):
        cat[col] = np.where(np.isfinite(cat[col]), cat[col], 0.0)

    if cache_path is not None:
        logger.info(f"Saving catalog cache to {cache_path}")
        cat.write(cache_path, overwrite=True)

    return cat


def build_varstar_cat_src(
    full_catalog: Table,
    new_xray: Table,
    new_obs: Table,
    existing_cat: Table,
) -> Table:
    """Build GaiaVarStar astromon_cat_src rows for all obsids.

    For each obsid:
    1. Coarse-filter the 130k-star catalog to sources within _FOV_RADIUS_DEG of
       the aimpoint (vectorized, no TAP query needed).
    2. Apply PM correction to the obs epoch.
    3. Match each candidate to its nearest xray source and record the pair.

    Parameters
    ----------
    full_catalog
        Full vari_rotation_modulation catalog from :func:`download_catalog`.
    new_xray
        Full astromon_xray_src table.
    new_obs
        Full astromon_obs table.
    existing_cat
        Existing astromon_cat_src rows (used to avoid id collisions).

    Returns
    -------
    Table compatible with ASTROMON_CAT_SRC_DTYPE.
    """
    # Use celldetect sources for x_id assignment; xcorr step re-matches per method.
    celldetect_xray = new_xray[new_xray["detect_method"] == "celldetect"]

    obs_by_id = {int(row["obsid"]): row for row in new_obs}

    # Max existing id per obsid so GaiaVarStar ids don't collide with Vizier/GaiaAGN ids.
    id_offset_by_obsid: dict[int, int] = {}
    for row in existing_cat:
        oid = int(row["obsid"])
        id_offset_by_obsid[oid] = (
            max(id_offset_by_obsid.get(oid, -1), int(row["id"])) + 1
        )

    cat_ra = np.asarray(full_catalog["ra"], dtype=float)
    cat_dec = np.asarray(full_catalog["dec"], dtype=float)
    cat_pmra = np.asarray(full_catalog["pmra"], dtype=float)
    cat_pmdec = np.asarray(full_catalog["pmdec"], dtype=float)

    obsids = np.unique(new_xray["obsid"])
    logger.info(f"Building GaiaVarStar cat_src for {len(obsids):,} obsids")

    all_rows: list[Table] = []

    for obsid in obsids:
        obsid_i = int(obsid)
        obs_row = obs_by_id.get(obsid_i)
        if obs_row is None:
            continue

        cd_sources = celldetect_xray[celldetect_xray["obsid"] == obsid_i]
        if len(cd_sources) == 0:
            logger.debug(f"OBSID={obsid_i} no celldetect sources, skipping")
            continue

        aimpoint_ra = float(obs_row["ra"])
        aimpoint_dec = float(obs_row["dec"])

        # Coarse cone filter in flat-sky approximation (~1 ms for 130k stars).
        cos_dec = np.cos(np.radians(aimpoint_dec))
        dra = (cat_ra - aimpoint_ra) * cos_dec
        ddec = cat_dec - aimpoint_dec
        in_fov = np.sqrt(dra**2 + ddec**2) < _FOV_RADIUS_DEG

        if not np.any(in_fov):
            continue

        nearby = full_catalog[in_fov]
        obs_epoch_yr = CxoTime(str(obs_row["date_obs"])).jyear
        ra_corr, dec_corr = cross_match._apply_proper_motion(
            cat_ra[in_fov],
            cat_dec[in_fov],
            cat_pmra[in_fov],
            cat_pmdec[in_fov],
            obs_epoch_yr,
        )

        cat_sc = coords.SkyCoord(ra_corr, dec_corr, unit="deg")
        xray_sc = coords.SkyCoord(
            np.asarray(cd_sources["ra"], dtype=float),
            np.asarray(cd_sources["dec"], dtype=float),
            unit="deg",
        )
        xray_idx, sep2d, _ = cat_sc.match_to_catalog_sky(xray_sc)

        q = Quat(equatorial=(aimpoint_ra, aimpoint_dec, float(obs_row["roll"])))
        y_angle, z_angle = radec_to_yagzag(ra_corr, dec_corr, q)

        id_offset = id_offset_by_obsid.get(obsid_i, 0)
        id_offset_by_obsid[obsid_i] = id_offset + len(nearby)

        t = Table(
            {
                "obsid": np.full(len(nearby), obsid_i, dtype=np.int32),
                "id": np.arange(len(nearby), dtype=np.int32) + id_offset,
                "x_id": cd_sources["id"][xray_idx].astype(np.int32),
                "catalog": np.full(len(nearby), "GaiaVarStar"),
                "name": [f"GaiaVarStar-{sid}" for sid in nearby["source_id"]],
                "ra": ra_corr,
                "dec": dec_corr,
                "separation": sep2d.arcsec.astype(np.float32),
                "mag": np.asarray(nearby["phot_g_mean_mag"], dtype=np.float32),
                "y_angle": y_angle.astype(np.float32),
                "z_angle": z_angle.astype(np.float32),
            }
        )
        all_rows.append(t)

    if not all_rows:
        return Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)
    return vstack(all_rows)


def run_xcorr_for_obsids(
    obsids: np.ndarray,
    varstar_cat: Table,
    new_xray: Table,
    new_obs: Table,
) -> Table:
    """Run compute_cross_matches('gaia_var_star') for every (obsid, detect_method) pair.

    For each detect_method, re-matches x_id so the join in simple_cross_match is correct.

    Returns
    -------
    Table of xcorr rows with detect_method column.
    """
    xcorr_cols = list(db.ASTROMON_XCORR_DTYPE.names)
    obs_by_id = {int(row["obsid"]): Table([row]) for row in new_obs}
    all_xcorr: list[Table] = []

    for obsid in obsids:
        obsid_i = int(obsid)
        cat_for_obsid = varstar_cat[varstar_cat["obsid"] == obsid_i]
        if len(cat_for_obsid) == 0:
            continue

        obspar = obs_by_id.get(obsid_i)
        if obspar is None:
            continue

        obsid_xray = new_xray[new_xray["obsid"] == obsid_i]
        detect_methods = np.unique(obsid_xray["detect_method"])
        cat_sc = coords.SkyCoord(cat_for_obsid["ra"], cat_for_obsid["dec"], unit="deg")

        for detect_method in detect_methods:
            sources = obsid_xray[obsid_xray["detect_method"] == detect_method]
            xray_sc = coords.SkyCoord(
                np.asarray(sources["ra"], dtype=float),
                np.asarray(sources["dec"], dtype=float),
                unit="deg",
            )
            xray_idx, _, _ = cat_sc.match_to_catalog_sky(xray_sc)

            candidates = cat_for_obsid.copy()
            candidates["x_id"] = sources["id"][xray_idx].astype(np.int32)

            try:
                matches = cross_match.compute_cross_matches(
                    "gaia_var_star",
                    astromon_obs=obspar,
                    astromon_xray_src=sources,
                    astromon_cat_src=candidates,
                )
            except Exception as exc:
                logger.warning(f"OBSID={obsid_i} {detect_method} xcorr failed: {exc}")
                continue

            if len(matches) == 0:
                continue

            matches["detect_method"] = detect_method
            all_xcorr.append(matches[xcorr_cols])
            logger.debug(
                f"OBSID={obsid_i} {detect_method}: {len(matches)} gaia_var_star match(es)"
            )

    if not all_xcorr:
        return Table(names=xcorr_cols)
    return vstack(all_xcorr)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path, help="astromon.h5 to update")
    parser.add_argument(
        "--catalog-cache",
        type=Path,
        default=None,
        help="Path to cache the downloaded Gaia catalog (ECSV or FITS). "
        "Loaded on subsequent runs instead of re-querying TAP.",
    )
    parser.add_argument(
        "--max-ruwe",
        type=float,
        default=1.4,
        help="Maximum RUWE for catalog download (default 1.4).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Build cat_src and xcorr but do not write to the database.",
    )
    args = parser.parse_args()

    cat = download_catalog(max_ruwe=args.max_ruwe, cache_path=args.catalog_cache)

    logger.info(f"Loading tables from {args.db}")
    existing_cat = db.get_table("astromon_cat_src", dbfile=args.db)
    new_xray = db.get_table("astromon_xray_src", dbfile=args.db)
    new_obs = db.get_table("astromon_obs", dbfile=args.db)

    existing_var = existing_cat[existing_cat["catalog"] == "GaiaVarStar"]
    if len(existing_var):
        logger.info(
            f"Removing {len(existing_var)} existing GaiaVarStar cat_src rows "
            f"({len(np.unique(existing_var['obsid']))} obsids) before rebuild"
        )
        existing_cat = existing_cat[existing_cat["catalog"] != "GaiaVarStar"]

    # Pass 1: build cat_src for all obsids.
    logger.info("Pass 1: building GaiaVarStar cat_src")
    varstar_cat = build_varstar_cat_src(cat, new_xray, new_obs, existing_cat)
    n_obsids = len(np.unique(varstar_cat["obsid"])) if len(varstar_cat) else 0
    logger.info(f"  {len(varstar_cat):,} candidate rows across {n_obsids:,} obsids")

    if len(varstar_cat) == 0:
        logger.warning(
            "No GaiaVarStar candidates found — check catalog and xray sources"
        )
        return

    if not args.dry_run:
        db.save(
            "astromon_cat_src", varstar_cat, dbfile=args.db, replace_keys=("catalog",)
        )
        logger.info("  Saved cat_src (GaiaVarStar rows replaced, others untouched)")

    # Pass 2: compute xcorr for all obsids × detect_methods.
    var_obsids = np.unique(varstar_cat["obsid"])
    logger.info(
        f"Pass 2: computing gaia_var_star xcorr for {len(var_obsids):,} obsids "
        f"× all detect_methods"
    )
    xcorr_result = run_xcorr_for_obsids(var_obsids, varstar_cat, new_xray, new_obs)
    logger.info(f"  {len(xcorr_result):,} gaia_var_star xcorr rows total")

    if len(xcorr_result) and not args.dry_run:
        from collections import Counter

        by_method = Counter(str(m) for m in xcorr_result["detect_method"])
        for method, count in sorted(by_method.items()):
            logger.info(f"    {method}: {count}")
        db.save("astromon_xcorr", xcorr_result, dbfile=args.db, select_name_key=True)
        logger.info("  Saved xcorr")
    elif len(xcorr_result) and args.dry_run:
        logger.info("  (dry-run: not writing to database)")
        from collections import Counter

        by_method = Counter(str(m) for m in xcorr_result["detect_method"])
        for method, count in sorted(by_method.items()):
            logger.info(f"    {method}: {count}")
    else:
        logger.warning(
            "No gaia_var_star xcorr rows produced — "
            "check GaiaVarStar cat_src and xray sources"
        )


if __name__ == "__main__":
    main()
