"""Backfill ICRF3 catalog matches for astromon_21 in astromon.h5.

ICRF3 S/X (Charlot et al. 2020, J/A+A/644/A159/table10) supersedes ICRF2
(I/323) with 4,588 sources vs 3,414 and improved astrometry.

Steps:
  1. Download the full ICRF3 S/X catalog from VizieR (~4,588 rows — trivial).
  2. Replace all catalog=="ICRS" rows in astromon_cat_src with ICRF3 matches,
     cone-filtered to each obsid's FOV and matched to xray sources.
  3. Re-run astromon_21 xcorr for affected obsids using the full updated
     cat_src (ICRS + existing Tycho2 entries).  Other select_names (gaia_agn,
     gaia_var_star) are preserved.

Usage:
    python -m astromon.scripts.maintenance.backfill_icrf3 --db astromon_merged.h5 [--dry-run]
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from astropy import coordinates as coords
from astropy.table import Table, vstack
from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat

from astromon import cross_match, db

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("backfill_icrf3")

# ICRF3 S/X band — Charlot et al. 2020, A&A 644, A159

# Cone pre-filter radius: ACIS-I ~16 arcmin, HRC-I ~30 arcmin; 0.35 deg covers both.
_FOV_RADIUS_DEG = 0.35


def build_icrs_cat_src(
    icrf3: Table,
    new_xray: Table,
    new_obs: Table,
    existing_cat: Table,
) -> Table:
    """Build ICRS astromon_cat_src rows for all obsids using ICRF3 positions.

    Parameters
    ----------
    icrf3
        Full ICRF3 catalog from :func:`astromon.cross_match.get_icrf3`.
    new_xray
        Full astromon_xray_src table.
    new_obs
        Full astromon_obs table.
    existing_cat
        Existing astromon_cat_src rows (non-ICRS, used to avoid id collisions).

    Returns
    -------
    Table compatible with ASTROMON_CAT_SRC_DTYPE.
    """
    celldetect_xray = new_xray[new_xray["detect_method"] == "celldetect"]
    obs_by_id = {int(row["obsid"]): row for row in new_obs}

    # Max existing non-ICRS id per obsid to avoid collisions.
    id_offset_by_obsid: dict[int, int] = {}
    for row in existing_cat:
        oid = int(row["obsid"])
        id_offset_by_obsid[oid] = (
            max(id_offset_by_obsid.get(oid, -1), int(row["id"])) + 1
        )

    cat_ra = np.asarray(icrf3["ra"], dtype=float)
    cat_dec = np.asarray(icrf3["dec"], dtype=float)

    obsids = np.unique(new_xray["obsid"])
    logger.info(f"Building ICRS cat_src for {len(obsids):,} obsids")

    all_rows: list[Table] = []

    for obsid in obsids:
        obsid_i = int(obsid)
        obs_row = obs_by_id.get(obsid_i)
        if obs_row is None:
            continue

        cd_sources = celldetect_xray[celldetect_xray["obsid"] == obsid_i]
        if len(cd_sources) == 0:
            continue

        aimpoint_ra = float(obs_row["ra"])
        aimpoint_dec = float(obs_row["dec"])

        cos_dec = np.cos(np.radians(aimpoint_dec))
        dra = (cat_ra - aimpoint_ra) * cos_dec
        ddec = cat_dec - aimpoint_dec
        in_fov = np.sqrt(dra**2 + ddec**2) < _FOV_RADIUS_DEG

        if not np.any(in_fov):
            continue

        nearby = icrf3[in_fov]
        cat_sc = coords.SkyCoord(
            np.asarray(nearby["ra"], dtype=float),
            np.asarray(nearby["dec"], dtype=float),
            unit="deg",
        )
        xray_sc = coords.SkyCoord(
            np.asarray(cd_sources["ra"], dtype=float),
            np.asarray(cd_sources["dec"], dtype=float),
            unit="deg",
        )
        xray_idx, sep2d, _ = cat_sc.match_to_catalog_sky(xray_sc)

        q = Quat(equatorial=(aimpoint_ra, aimpoint_dec, float(obs_row["roll"])))
        y_angle, z_angle = radec_to_yagzag(
            np.asarray(nearby["ra"], dtype=float),
            np.asarray(nearby["dec"], dtype=float),
            q,
        )

        id_offset = id_offset_by_obsid.get(obsid_i, 0)
        id_offset_by_obsid[obsid_i] = id_offset + len(nearby)

        t = Table(
            {
                "obsid": np.full(len(nearby), obsid_i, dtype=np.int32),
                "id": np.arange(len(nearby), dtype=np.int32) + id_offset,
                "x_id": cd_sources["id"][xray_idx].astype(np.int32),
                "catalog": np.full(len(nearby), "ICRS"),
                "name": list(nearby["name"]),
                "ra": np.asarray(nearby["ra"], dtype=np.float64),
                "dec": np.asarray(nearby["dec"], dtype=np.float64),
                "separation": sep2d.arcsec.astype(np.float32),
                "mag": np.full(len(nearby), np.nan, dtype=np.float32),
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
    full_cat: Table,
    new_xray: Table,
    new_obs: Table,
) -> Table:
    """Re-run astromon_21 xcorr for the given obsids using the full cat_src.

    Returns
    -------
    Table of astromon_21 xcorr rows with detect_method column.
    """
    xcorr_cols = list(db.ASTROMON_XCORR_DTYPE.names)
    obs_by_id = {int(row["obsid"]): Table([row]) for row in new_obs}
    all_xcorr: list[Table] = []

    for obsid in obsids:
        obsid_i = int(obsid)
        cat_for_obsid = full_cat[full_cat["obsid"] == obsid_i]
        if len(cat_for_obsid) == 0:
            continue

        obspar = obs_by_id.get(obsid_i)
        if obspar is None:
            continue

        obsid_xray = new_xray[new_xray["obsid"] == obsid_i]
        for detect_method in np.unique(obsid_xray["detect_method"]):
            sources = obsid_xray[obsid_xray["detect_method"] == detect_method]

            cat_sc = coords.SkyCoord(
                cat_for_obsid["ra"], cat_for_obsid["dec"], unit="deg"
            )
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
                    "astromon_21",
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
                f"OBSID={obsid_i} {detect_method}: {len(matches)} astromon_21 match(es)"
            )

    if not all_xcorr:
        return Table(names=xcorr_cols)
    return vstack(all_xcorr)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    icrf3 = cross_match.get_icrf3()

    logger.info(f"Loading tables from {args.db}")
    existing_cat = db.get_table("astromon_cat_src", dbfile=args.db)
    new_xray = db.get_table("astromon_xray_src", dbfile=args.db)
    new_obs = db.get_table("astromon_obs", dbfile=args.db)

    n_old_icrs = np.sum(existing_cat["catalog"] == "ICRS")
    logger.info(f"  Existing ICRS cat_src rows: {n_old_icrs:,}")

    non_icrs_cat = existing_cat[existing_cat["catalog"] != "ICRS"]

    # Pass 1: build new ICRS cat_src from ICRF3.
    logger.info("Pass 1: building ICRS cat_src from ICRF3")
    new_icrs_cat = build_icrs_cat_src(icrf3, new_xray, new_obs, non_icrs_cat)
    n_obsids = len(np.unique(new_icrs_cat["obsid"])) if len(new_icrs_cat) else 0
    logger.info(f"  {len(new_icrs_cat):,} ICRS cat_src rows across {n_obsids:,} obsids")

    merged_cat = vstack([non_icrs_cat, new_icrs_cat], metadata_conflicts="silent")

    # Pass 2: re-run astromon_21 xcorr for obsids with ICRS candidates.
    icrs_obsids = (
        np.unique(new_icrs_cat["obsid"]) if len(new_icrs_cat) else np.array([])
    )
    logger.info(f"Pass 2: re-running astromon_21 xcorr for {len(icrs_obsids):,} obsids")
    new_xcorr = run_xcorr_for_obsids(icrs_obsids, merged_cat, new_xray, new_obs)
    logger.info(f"  {len(new_xcorr):,} astromon_21 xcorr rows")

    if args.dry_run:
        logger.info("  (dry-run: not writing to database)")
        return

    if not len(new_icrs_cat):
        logger.warning("No ICRS cat_src rows produced — nothing to save")
        return

    db.save("astromon_cat_src", new_icrs_cat, dbfile=args.db, replace_keys=("catalog",))
    logger.info("  Saved cat_src (ICRS rows replaced, other catalogs untouched)")

    if len(new_xcorr):
        # select_name_key=True ensures only astromon_21 rows for the affected
        # obsids are replaced; gaia_agn and gaia_var_star rows are untouched.
        db.save("astromon_xcorr", new_xcorr, dbfile=args.db, select_name_key=True)
        logger.info("  Saved xcorr (astromon_21 updated, other select_names preserved)")
    else:
        logger.info("  No new xcorr rows; xcorr table unchanged")


if __name__ == "__main__":
    main()
