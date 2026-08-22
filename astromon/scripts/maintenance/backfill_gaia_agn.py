"""Merge GaiaAGN catalog data from an old astromon.h5 into a new one, then
backfill gaia_agn xcorr rows for all detect_methods in the new h5.

Usage:
    python -m astromon.scripts.maintenance.backfill_gaia_agn --old-db ~/git/astromon/astromon.h5 --new-db astromon.h5

The old h5 has GaiaAGN cat_src and xcorr rows but no detect_method column.
The new h5 has celldetect + gaussian_detect + peak_gaussian_detect but no GaiaAGN.

Three passes:
  1. Rebuild GaiaAGN astromon_cat_src using positions from the old h5, re-matched
     to the new h5's xray sources.
  2. Run compute_cross_matches("gaia_agn") for every (obsid, detect_method)
     combination that has GaiaAGN candidates, collecting xcorr rows.
  3. Save cat_src (merged with existing Vizier entries) and xcorr to the new h5.
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
logger = logging.getLogger("backfill_gaia_agn")


def build_gaia_cat_src(
    old_gaia: Table,
    new_xray: Table,
    new_obs: Table,
    existing_cat: Table,
) -> Table:
    """Rebuild GaiaAGN cat_src rows with x_id matched to new h5's celldetect sources.

    Parameters
    ----------
    old_gaia
        GaiaAGN rows from old astromon_cat_src (catalog == "GaiaAGN").
    new_xray
        Full astromon_xray_src from new h5.
    new_obs
        Full astromon_obs from new h5.
    existing_cat
        Non-GaiaAGN astromon_cat_src from new h5 (Vizier entries).

    Returns
    -------
    Table with ASTROMON_CAT_SRC_DTYPE-compatible columns for all GaiaAGN sources.
    """
    celldetect_xray = new_xray[new_xray["detect_method"] == "celldetect"]

    # Pre-build a lookup: obsid -> (obs_row, celldetect_sources)
    obs_by_id = {int(row["obsid"]): row for row in new_obs}

    # Max existing id per obsid so GaiaAGN ids don't collide with Vizier ids
    id_offset_by_obsid: dict[int, int] = {}
    for row in existing_cat:
        oid = int(row["obsid"])
        id_offset_by_obsid[oid] = (
            max(id_offset_by_obsid.get(oid, -1), int(row["id"])) + 1
        )

    obsids_with_gaia = np.unique(old_gaia["obsid"])
    logger.info(f"Building GaiaAGN cat_src for {len(obsids_with_gaia)} obsids")

    all_rows: list[Table] = []

    for obsid in obsids_with_gaia:
        obsid_i = int(obsid)
        gaia_rows = old_gaia[old_gaia["obsid"] == obsid_i]
        cd_sources = celldetect_xray[celldetect_xray["obsid"] == obsid_i]

        if len(cd_sources) == 0:
            logger.debug(f"OBSID={obsid_i} no celldetect sources in new h5, skipping")
            continue

        obs_row = obs_by_id.get(obsid_i)
        if obs_row is None:
            logger.debug(f"OBSID={obsid_i} not in new obs table, skipping")
            continue

        q = Quat(
            equatorial=(
                float(obs_row["ra"]),
                float(obs_row["dec"]),
                float(obs_row["roll"]),
            )
        )

        agn_sc = coords.SkyCoord(gaia_rows["ra"], gaia_rows["dec"], unit="deg")
        xray_sc = coords.SkyCoord(cd_sources["ra"], cd_sources["dec"], unit="deg")
        xray_idx, sep2d, _ = agn_sc.match_to_catalog_sky(xray_sc)

        id_offset = id_offset_by_obsid.get(obsid_i, 0)
        y_angle, z_angle = radec_to_yagzag(
            np.array(gaia_rows["ra"], dtype=float),
            np.array(gaia_rows["dec"], dtype=float),
            q,
        )

        t = Table(
            {
                "obsid": np.full(len(gaia_rows), obsid_i, dtype=np.int32),
                "id": np.arange(len(gaia_rows), dtype=np.int32) + id_offset,
                "x_id": cd_sources["id"][xray_idx].astype(np.int32),
                "catalog": gaia_rows["catalog"],
                "name": gaia_rows["name"],
                "ra": np.array(gaia_rows["ra"], dtype=np.float64),
                "dec": np.array(gaia_rows["dec"], dtype=np.float64),
                "separation": sep2d.arcsec.astype(np.float32),
                "mag": np.array(gaia_rows["mag"], dtype=np.float32),
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
    gaia_cat: Table,
    new_xray: Table,
    new_obs: Table,
) -> Table:
    """Run compute_cross_matches('gaia_agn') for every (obsid, detect_method) pair.

    For each detect_method, GaiaAGN x_id is re-matched to that method's xray source
    IDs so the join in simple_cross_match works correctly.

    Returns
    -------
    Table of xcorr rows with detect_method column.
    """
    xcorr_cols = list(db.ASTROMON_XCORR_DTYPE.names)
    obs_by_id = {int(row["obsid"]): Table([row]) for row in new_obs}
    all_xcorr: list[Table] = []

    for obsid in obsids:
        obsid_i = int(obsid)
        gaia_for_obsid = gaia_cat[gaia_cat["obsid"] == obsid_i]
        if len(gaia_for_obsid) == 0:
            continue

        obspar = obs_by_id.get(obsid_i)
        if obspar is None:
            continue

        obsid_xray = new_xray[new_xray["obsid"] == obsid_i]
        detect_methods = np.unique(obsid_xray["detect_method"])

        agn_sc = coords.SkyCoord(
            gaia_for_obsid["ra"], gaia_for_obsid["dec"], unit="deg"
        )

        for detect_method in detect_methods:
            sources = obsid_xray[obsid_xray["detect_method"] == detect_method]

            xray_sc = coords.SkyCoord(sources["ra"], sources["dec"], unit="deg")
            xray_idx, _, _ = agn_sc.match_to_catalog_sky(xray_sc)

            candidates = gaia_for_obsid.copy()
            candidates["x_id"] = sources["id"][xray_idx].astype(np.int32)

            try:
                gaia_matches = cross_match.compute_cross_matches(
                    "gaia_agn",
                    astromon_obs=obspar,
                    astromon_xray_src=sources,
                    astromon_cat_src=candidates,
                )
            except Exception as exc:
                logger.warning(f"OBSID={obsid_i} {detect_method} xcorr failed: {exc}")
                continue

            if len(gaia_matches) == 0:
                continue

            gaia_matches["detect_method"] = detect_method
            all_xcorr.append(gaia_matches[xcorr_cols])
            logger.debug(
                f"OBSID={obsid_i} {detect_method}: {len(gaia_matches)} gaia_agn match(es)"
            )

    if not all_xcorr:
        return Table(names=xcorr_cols)
    return vstack(all_xcorr)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--old-db", required=True, type=Path, help="Old astromon.h5 with GaiaAGN"
    )
    parser.add_argument(
        "--new-db", required=True, type=Path, help="New astromon.h5 to update"
    )
    args = parser.parse_args()

    logger.info(f"Loading tables from old h5: {args.old_db}")
    old_cat = db.get_table("astromon_cat_src", dbfile=args.old_db)

    logger.info(f"Loading tables from new h5: {args.new_db}")
    new_cat = db.get_table("astromon_cat_src", dbfile=args.new_db)
    new_xray = db.get_table("astromon_xray_src", dbfile=args.new_db)
    new_obs = db.get_table("astromon_obs", dbfile=args.new_db)

    old_gaia = old_cat[old_cat["catalog"] == "GaiaAGN"]
    logger.info(
        f"Old h5 GaiaAGN: {len(old_gaia)} cat_src rows, "
        f"{len(np.unique(old_gaia['obsid']))} obsids"
    )

    # ── Pass 1: rebuild GaiaAGN cat_src with new x_ids ───────────────────────
    logger.info("Pass 1: rebuilding GaiaAGN cat_src")
    new_gaia_cat = build_gaia_cat_src(old_gaia, new_xray, new_obs, new_cat)
    logger.info(f"  Built {len(new_gaia_cat)} GaiaAGN cat_src rows")

    merged_cat = vstack([new_cat, new_gaia_cat], metadata_conflicts="silent")
    logger.info(f"  Merged cat_src: {len(merged_cat)} total rows")

    db.save("astromon_cat_src", merged_cat, dbfile=args.new_db, ignore_obsid=True)
    logger.info("  Saved merged cat_src")

    # ── Pass 2 + 3: compute xcorr for all obsids with GaiaAGN ────────────────
    gaia_obsids = np.unique(new_gaia_cat["obsid"])
    logger.info(
        f"Pass 2/3: computing gaia_agn xcorr for {len(gaia_obsids)} obsids "
        f"× all detect_methods"
    )
    xcorr_result = run_xcorr_for_obsids(gaia_obsids, new_gaia_cat, new_xray, new_obs)
    logger.info(f"  {len(xcorr_result)} gaia_agn xcorr rows total")

    if len(xcorr_result):
        from collections import Counter

        by_method = Counter(str(m) for m in xcorr_result["detect_method"])
        for method, count in sorted(by_method.items()):
            logger.info(f"    {method}: {count}")
        # This backfill computes only gaia_agn rows, so key on
        # (obsid, detect_method, select_name) to avoid deleting other
        # select_names already present for the same obsid+detect_method.
        db.save(
            "astromon_xcorr",
            xcorr_result,
            dbfile=args.new_db,
            select_name_key=True,
        )
        logger.info("  Saved xcorr")
    else:
        logger.warning(
            "  No gaia_agn xcorr rows produced — check GaiaAGN cat_src and xray sources"
        )


if __name__ == "__main__":
    main()
