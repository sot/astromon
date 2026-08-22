"""Add GaiaAGN matches for the 107 reprocessed obsids via Gaia TAP.

The migration script (backfill_gaia_agn.py) already copied GaiaAGN entries from the
old h5 and recomputed xcorr for all detect_methods. This script is a second pass that
queries Gaia TAP directly using the new detection positions — which may be slightly
different from the old celldetect positions — and adds any Gaia AGN sources that were
not already in the catalog for each obsid.

Usage:
    SKA=/Users/jean/ska python -m astromon.scripts.maintenance.tap_backfill_107 \\
        [--db astromon_merged.h5]
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from astropy import coordinates as coords
from astropy import units as u
from astropy.table import Table, vstack
from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat

from astromon import cross_match, db

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("tap_backfill_107")


def get_new_gaia_for_obsid(
    obsid: int,
    xray_sources: Table,
    existing_cat: Table,
    obs_row: Table,
    radius: u.Quantity = 3 * u.arcsec,
) -> Table:
    """Query Gaia TAP for AGN near xray_sources and return only sources not already in cat_src.

    Parameters
    ----------
    obsid
        The obsid to process.
    xray_sources
        All xray sources for this obsid (any detect_method).
    existing_cat
        Current astromon_cat_src rows for this obsid.
    obs_row
        Single-row Table from astromon_obs for this obsid.
    radius
        TAP search radius.

    Returns
    -------
    Table with ASTROMON_CAT_SRC_DTYPE-compatible columns for new GaiaAGN sources only,
    or empty table if nothing new is found.
    """
    pos = coords.SkyCoord(
        np.array(xray_sources["ra"], dtype=float),
        np.array(xray_sources["dec"], dtype=float),
        unit="deg",
    )
    gaia_results = cross_match.get_gaia_agn(
        pos, radius=radius, logging_tag=f"OBSID={obsid}"
    )

    if len(gaia_results) == 0:
        return Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)

    existing_names = set(
        existing_cat[existing_cat["catalog"].astype(str) == "GaiaAGN"]["name"].astype(
            str
        )
    )
    new_mask = np.array([str(n) not in existing_names for n in gaia_results["name"]])
    new_gaia = gaia_results[new_mask]

    if len(new_gaia) == 0:
        logger.debug(
            f"OBSID={obsid}: {len(gaia_results)} Gaia AGN found, all already in catalog"
        )
        return Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)

    logger.info(
        f"OBSID={obsid}: {len(new_gaia)} NEW Gaia AGN (of {len(gaia_results)} found)"
    )

    # Match each new Gaia source to nearest celldetect xray source for x_id
    celldetect_src = xray_sources[
        xray_sources["detect_method"].astype(str) == "celldetect"
    ]
    if len(celldetect_src) == 0:
        celldetect_src = xray_sources  # fall back to whatever is available

    agn_sc = coords.SkyCoord(
        np.array(new_gaia["ra"], dtype=float),
        np.array(new_gaia["dec"], dtype=float),
        unit="deg",
    )
    xray_sc = coords.SkyCoord(
        np.array(celldetect_src["ra"], dtype=float),
        np.array(celldetect_src["dec"], dtype=float),
        unit="deg",
    )
    xray_idx, sep2d, _ = agn_sc.match_to_catalog_sky(xray_sc)

    q = Quat(
        equatorial=(float(obs_row["ra"]), float(obs_row["dec"]), float(obs_row["roll"]))
    )
    y_angle, z_angle = radec_to_yagzag(
        np.array(new_gaia["ra"], dtype=float),
        np.array(new_gaia["dec"], dtype=float),
        q,
    )

    # Assign ids above the existing max for this obsid
    existing_max_id = int(np.max(existing_cat["id"])) if len(existing_cat) else -1
    new_ids = np.arange(len(new_gaia), dtype=np.int32) + existing_max_id + 1

    return Table(
        {
            "obsid": np.full(len(new_gaia), obsid, dtype=np.int32),
            "id": new_ids,
            "x_id": celldetect_src["id"][xray_idx].astype(np.int32),
            "catalog": new_gaia["catalog"],
            "name": new_gaia["name"],
            "ra": np.array(new_gaia["ra"], dtype=np.float64),
            "dec": np.array(new_gaia["dec"], dtype=np.float64),
            "separation": sep2d.arcsec.astype(np.float32),
            "mag": np.array(new_gaia["mag"], dtype=np.float32),
            "y_angle": y_angle.astype(np.float32),
            "z_angle": z_angle.astype(np.float32),
        }
    )


def main() -> None:  # noqa: PLR0915
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="astromon_merged.h5", type=Path)
    args = parser.parse_args()

    logger.info(f"Loading tables from {args.db}")
    cat = db.get_table("astromon_cat_src", dbfile=args.db)
    xray = db.get_table("astromon_xray_src", dbfile=args.db)
    obs = db.get_table("astromon_obs", dbfile=args.db)

    new_detect_xray = xray[xray["detect_method"].astype(str) != "celldetect"]
    target_obsids = np.unique(new_detect_xray["obsid"])
    logger.info(
        f"Processing {len(target_obsids)} obsids with gaussian/peak_gaussian_detect sources"
    )

    obs_by_id = {int(row["obsid"]): Table([row]) for row in obs}
    cat_by_obsid = {
        int(oid): cat[cat["obsid"].astype(int) == int(oid)] for oid in target_obsids
    }
    xray_by_obsid = {
        int(oid): xray[xray["obsid"].astype(int) == int(oid)] for oid in target_obsids
    }

    new_cat_rows: list[Table] = []

    for obsid in target_obsids:
        obsid_i = int(obsid)
        obs_row = obs_by_id.get(obsid_i)
        if obs_row is None:
            logger.warning(f"OBSID={obsid_i} not in obs table, skipping")
            continue

        new_rows = get_new_gaia_for_obsid(
            obsid_i,
            xray_sources=xray_by_obsid[obsid_i],
            existing_cat=cat_by_obsid[obsid_i],
            obs_row=obs_row[0],
        )
        if len(new_rows):
            new_cat_rows.append(new_rows)

    if not new_cat_rows:
        logger.info(
            "No new GaiaAGN sources found — catalog already complete for these obsids"
        )
        return

    new_gaia_combined = vstack(new_cat_rows)
    logger.info(
        f"Found {len(new_gaia_combined)} new GaiaAGN cat_src rows across "
        f"{len(np.unique(new_gaia_combined['obsid']))} obsids"
    )

    merged_cat = vstack([cat, new_gaia_combined], metadata_conflicts="silent")

    db.save("astromon_cat_src", merged_cat, dbfile=args.db, ignore_obsid=True)
    logger.info("Saved updated cat_src")

    # Reload and compute xcorr for new entries
    updated_cat = db.get_table("astromon_cat_src", dbfile=args.db)

    xcorr_cols = list(db.ASTROMON_XCORR_DTYPE.names)
    all_xcorr: list[Table] = []

    new_obsids_with_gaia = np.unique(new_gaia_combined["obsid"])
    for obsid in new_obsids_with_gaia:
        obsid_i = int(obsid)
        obspar = obs_by_id.get(obsid_i)
        if obspar is None:
            continue

        gaia_for_obsid = updated_cat[
            (updated_cat["obsid"].astype(int) == obsid_i)
            & (updated_cat["catalog"].astype(str) == "GaiaAGN")
        ]
        obsid_xray = xray_by_obsid[obsid_i]

        agn_sc = coords.SkyCoord(
            np.array(gaia_for_obsid["ra"], dtype=float),
            np.array(gaia_for_obsid["dec"], dtype=float),
            unit="deg",
        )

        for detect_method in np.unique(obsid_xray["detect_method"].astype(str)):
            sources = obsid_xray[
                obsid_xray["detect_method"].astype(str) == detect_method
            ]
            xray_sc = coords.SkyCoord(
                np.array(sources["ra"], dtype=float),
                np.array(sources["dec"], dtype=float),
                unit="deg",
            )
            xray_idx, _, _ = agn_sc.match_to_catalog_sky(xray_sc)

            candidates = gaia_for_obsid.copy()
            candidates["x_id"] = sources["id"][xray_idx].astype(np.int32)

            try:
                matches = cross_match.compute_cross_matches(
                    "gaia_agn",
                    astromon_obs=obspar,
                    astromon_xray_src=sources,
                    astromon_cat_src=candidates,
                )
            except Exception as exc:
                logger.warning(f"OBSID={obsid_i} {detect_method} xcorr failed: {exc}")
                continue

            if len(matches):
                matches["detect_method"] = detect_method
                all_xcorr.append(matches[xcorr_cols])

    if all_xcorr:
        xcorr_result = vstack(all_xcorr)
        logger.info(f"Saving {len(xcorr_result)} new gaia_agn xcorr rows")
        # This incremental pass adds only gaia_agn rows, so preserve
        # unrelated select_names for the same obsid+detect_method.
        db.save(
            "astromon_xcorr",
            xcorr_result,
            dbfile=args.db,
            select_name_key=True,
        )
    else:
        logger.info("No new xcorr rows to add")


if __name__ == "__main__":
    main()
