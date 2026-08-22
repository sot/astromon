"""Backfill missing ICRF3 catalog sources and icrf3 xcorr rows.

Between 2026-08-04 and 2026-08-10, get_cat_obs_data.py used only RFC in its
rough_match catalog list (ICRF3 was temporarily removed).  Observations processed
during that window have RFC entries in astromon_cat_src but no ICRF3 entries, so
the ``icrf3`` cross-match finds zero rows for them.

This script is purely additive:
  1. Identifies obsids that have other catalog rows in astromon_cat_src but are
     missing ICRF3.
  2. Looks up ICRF3 sources within the FOV of each such obsid using the local
     cached catalog (no VizieR query needed).
  3. Appends the new ICRF3 astromon_cat_src rows to the existing table.
  4. Runs the ``icrf3`` cross-match for all detect_methods present in
     astromon_xray_src for the affected obsids.
  5. Saves the new ``icrf3`` xcorr rows (select_name_key=True so other
     select_names are untouched).

Usage::

    python -m astromon.scripts.maintenance.backfill_icrf3_addonly --db /Volumes/Black/data/astromon/astromon.h5
    python -m astromon.scripts.maintenance.backfill_icrf3_addonly --db /path/to/astromon.h5 --dry-run
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
from astromon.cross_match import get_icrf3

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("backfill_icrf3_addonly")

# Pre-filter radius: ACIS-I ~16', HRC-I ~30'; 0.35° generously covers both.
_FOV_RADIUS_DEG = 0.35


def build_icrf3_cat_src(
    icrf3: Table,
    celldetect_xray: Table,
    obspar: Table,
    target_obsids: np.ndarray,
    existing_cat: Table,
) -> Table:
    """Build ICRF3 astromon_cat_src rows for *target_obsids*.

    Parameters
    ----------
    icrf3
        Full ICRF3 catalog from :func:`get_icrf3`.
    celldetect_xray
        celldetect rows from astromon_xray_src.
    obspar
        Full astromon_obs table.
    target_obsids
        Obsids to process (already filtered to those missing ICRF3 in cat_src).
    existing_cat
        Existing astromon_cat_src rows; used to compute id offsets so new rows
        do not collide with existing ones.

    Returns
    -------
    Table in astromon_cat_src format containing only the new ICRF3 rows.
    """
    obs_by_id = {int(row["obsid"]): row for row in obspar}

    # Max existing id per obsid, so new ICRF3 ids start above them.
    id_offset: dict[int, int] = {}
    for row in existing_cat:
        oid = int(row["obsid"])
        id_offset[oid] = max(id_offset.get(oid, -1), int(row["id"])) + 1

    cat_ra = np.asarray(icrf3["ra"], dtype=float)
    cat_dec = np.asarray(icrf3["dec"], dtype=float)

    all_rows: list[Table] = []
    n_no_xray = 0
    n_no_fov = 0

    for obsid in target_obsids:
        obsid_i = int(obsid)
        obs_row = obs_by_id.get(obsid_i)
        if obs_row is None:
            continue

        cd_sources = celldetect_xray[celldetect_xray["obsid"] == obsid_i]
        if len(cd_sources) == 0:
            n_no_xray += 1
            continue

        aimpoint_ra = float(obs_row["ra"])
        aimpoint_dec = float(obs_row["dec"])

        cos_dec = np.cos(np.radians(aimpoint_dec))
        dra = (cat_ra - aimpoint_ra) * cos_dec
        ddec = cat_dec - aimpoint_dec
        in_fov = np.sqrt(dra**2 + ddec**2) < _FOV_RADIUS_DEG

        if not np.any(in_fov):
            n_no_fov += 1
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

        base_id = id_offset.get(obsid_i, 0)
        id_offset[obsid_i] = base_id + len(nearby)

        all_rows.append(
            Table(
                {
                    "obsid": np.full(len(nearby), obsid_i, dtype=np.int32),
                    "id": np.arange(len(nearby), dtype=np.int32) + base_id,
                    "x_id": cd_sources["id"][xray_idx].astype(np.int32),
                    "catalog": np.full(len(nearby), "ICRF3"),
                    "name": list(nearby["name"]),
                    "ra": np.asarray(nearby["ra"], dtype=np.float64),
                    "dec": np.asarray(nearby["dec"], dtype=np.float64),
                    "separation": sep2d.arcsec.astype(np.float32),
                    "mag": np.full(len(nearby), np.nan, dtype=np.float32),
                    "y_angle": y_angle.astype(np.float32),
                    "z_angle": z_angle.astype(np.float32),
                }
            )
        )

    logger.debug(
        f"  {n_no_xray} obsids skipped (no celldetect sources), "
        f"{n_no_fov} skipped (no ICRF3 in FOV)"
    )

    if not all_rows:
        return Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)
    return vstack(all_rows)


def run_icrf3_xcorr(
    obsids: np.ndarray,
    full_cat: Table,
    xray: Table,
    obspar: Table,
) -> Table:
    """Run the ``icrf3`` cross-match for *obsids* across all detect_methods.

    Parameters
    ----------
    obsids
        Obsids to process (those that gained new ICRF3 cat_src entries).
    full_cat
        Full astromon_cat_src, including newly appended ICRF3 rows.
    xray
        Full astromon_xray_src.
    obspar
        Full astromon_obs.

    Returns
    -------
    Table of ``icrf3`` xcorr rows with detect_method column.
    """
    xcorr_cols = list(db.ASTROMON_XCORR_DTYPE.names)
    obs_by_id = {int(row["obsid"]): Table([row]) for row in obspar}
    all_xcorr: list[Table] = []

    for obsid in obsids:
        obsid_i = int(obsid)
        cat_for_obsid = full_cat[full_cat["obsid"] == obsid_i]
        if len(cat_for_obsid) == 0:
            continue

        obs_row = obs_by_id.get(obsid_i)
        if obs_row is None:
            continue

        obsid_xray = xray[xray["obsid"] == obsid_i]
        for detect_method in np.unique(obsid_xray["detect_method"]):
            sources = obsid_xray[obsid_xray["detect_method"] == detect_method]

            # Remap x_id to the nearest source in this detect_method version.
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
                    "icrf3",
                    astromon_obs=obs_row,
                    astromon_xray_src=sources,
                    astromon_cat_src=candidates,
                )
            except Exception as exc:
                logger.warning(
                    f"OBSID={obsid_i} {detect_method} icrf3 xcorr failed: {exc}"
                )
                continue

            if not len(matches):
                continue

            matches["detect_method"] = detect_method
            all_xcorr.append(matches[xcorr_cols])
            logger.debug(
                f"OBSID={obsid_i} {detect_method}: {len(matches)} icrf3 match(es)"
            )

    if not all_xcorr:
        return Table(names=xcorr_cols)
    return vstack(all_xcorr)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path, help="Path to astromon.h5")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without writing to the database.",
    )
    args = parser.parse_args()

    logger.info("Loading ICRF3 catalog from local cache")
    icrf3 = get_icrf3()
    logger.info(f"  {len(icrf3):,} ICRF3 sources")

    logger.info(f"Loading tables from {args.db}")
    existing_cat = db.get_table("astromon_cat_src", dbfile=args.db)
    xray = db.get_table("astromon_xray_src", dbfile=args.db)
    obspar = db.get_table("astromon_obs", dbfile=args.db)

    # Obsids with any cat_src rows but no ICRF3 — the target for this backfill.
    obsids_with_cat = np.unique(existing_cat["obsid"])
    obsids_with_icrf3 = np.unique(
        existing_cat[existing_cat["catalog"] == "ICRF3"]["obsid"]
    )
    target_obsids = obsids_with_cat[~np.isin(obsids_with_cat, obsids_with_icrf3)]
    logger.info(
        f"  {len(obsids_with_cat):,} obsids in cat_src, "
        f"{len(obsids_with_icrf3):,} already have ICRF3, "
        f"{len(target_obsids):,} need ICRF3 backfill"
    )

    celldetect_xray = xray[xray["detect_method"] == "celldetect"]

    # Pass 1: build ICRF3 cat_src rows for target obsids.
    logger.info(f"Pass 1: building ICRF3 cat_src for {len(target_obsids):,} obsids")
    new_icrf3_cat = build_icrf3_cat_src(
        icrf3, celldetect_xray, obspar, target_obsids, existing_cat
    )
    if len(new_icrf3_cat):
        n_icrf3_obsids = len(np.unique(new_icrf3_cat["obsid"]))
    else:
        n_icrf3_obsids = 0
    logger.info(
        f"  {len(new_icrf3_cat):,} new ICRF3 cat_src rows "
        f"across {n_icrf3_obsids:,} obsids with ICRF3 in FOV"
    )

    if len(new_icrf3_cat) == 0:
        logger.info("No ICRF3 sources found in FOV of any target obsid — nothing to do")
        return

    merged_cat = vstack([existing_cat, new_icrf3_cat], metadata_conflicts="silent")
    icrf3_obsids = np.unique(new_icrf3_cat["obsid"])

    # Pass 2: run icrf3 xcorr for the obsids that got new ICRF3 cat_src entries.
    logger.info(f"Pass 2: running icrf3 xcorr for {len(icrf3_obsids):,} obsids")
    new_xcorr = run_icrf3_xcorr(icrf3_obsids, merged_cat, xray, obspar)
    n_xcorr_obsids = len(np.unique(new_xcorr["obsid"])) if len(new_xcorr) else 0
    logger.info(
        f"  {len(new_xcorr):,} icrf3 xcorr rows across {n_xcorr_obsids:,} obsids"
    )

    if args.dry_run:
        logger.info("(dry-run: not writing to database)")
        return

    logger.info("Saving new ICRF3 cat_src rows (replace_keys=('obsid', 'catalog'))")
    db.save(
        "astromon_cat_src",
        new_icrf3_cat,
        dbfile=args.db,
        replace_keys=("obsid", "catalog"),
    )
    logger.info(
        f"  cat_src saved: {len(new_icrf3_cat):,} new ICRF3 rows added to "
        f"{len(existing_cat):,} existing"
    )

    if len(new_xcorr):
        logger.info(
            "Saving icrf3 xcorr "
            "(select_name_key=True: replaces only icrf3 rows for affected obsids)"
        )
        db.save("astromon_xcorr", new_xcorr, dbfile=args.db, select_name_key=True)
        logger.info(f"  {len(new_xcorr):,} icrf3 xcorr rows saved")
    else:
        logger.info("No icrf3 xcorr rows to save")

    logger.info("Done")


if __name__ == "__main__":
    main()
