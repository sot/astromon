"""Backfill Radio Fundamental Catalog (RFC) matches in astromon.h5.

RFC (Petrov & Kovalev 2025, ApJS 276:38) provides VLBI-measured positions for
22,839 compact radio sources with sub-milliarcsecond accuracy.  It is a strict
superset of ICRF3 and adds ~18k additional calibrators.

RFC replaces ICRF3 as the source of catalog=="RFC" entries in astromon_cat_src
(previously stored as catalog=="ICRS" before the rename).  Matches are re-run
under select_name="astromon_21", so the result folds cleanly into the existing
astromon_21 pipeline without deduplication complexity.

Usage:
    python -m astromon.scripts.maintenance.backfill_rfc \\
        --db astromon_merged.h5 --catalog-cache rfc_2026b_cat.txt
    python -m astromon.scripts.maintenance.backfill_rfc \\
        --db astromon_merged.h5 --catalog-cache rfc_2026b_cat.txt --dry-run
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
from astromon.cross_match import get_rfc

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("backfill_rfc")

# Cone pre-filter radius: ACIS-I ~16 arcmin, HRC-I ~30 arcmin; 0.35 deg covers both.
_FOV_RADIUS_DEG = 0.35


def build_rfc_cat_src(
    rfc: Table,
    new_xray: Table,
    new_obs: Table,
    non_rfc_cat: Table,
) -> Table:
    """Build RFC astromon_cat_src rows for all obsids.

    Parameters
    ----------
    rfc
        Full RFC catalog from :func:`download_rfc`.
    new_xray
        Full astromon_xray_src table.
    new_obs
        Full astromon_obs table.
    non_rfc_cat
        Existing astromon_cat_src rows with catalog not in ("ICRS", "RFC"),
        used for id-offset calculation to avoid collisions.
    """
    celldetect_xray = new_xray[new_xray["detect_method"] == "celldetect"]
    obs_by_id = {int(row["obsid"]): row for row in new_obs}

    id_offset_by_obsid: dict[int, int] = {}
    for row in non_rfc_cat:
        oid = int(row["obsid"])
        id_offset_by_obsid[oid] = (
            max(id_offset_by_obsid.get(oid, -1), int(row["id"])) + 1
        )

    cat_ra = np.asarray(rfc["ra"], dtype=float)
    cat_dec = np.asarray(rfc["dec"], dtype=float)

    obsids = np.unique(new_xray["obsid"])
    logger.info(f"Building RFC cat_src for {len(obsids):,} obsids")

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

        nearby = rfc[in_fov]
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

        all_rows.append(
            Table(
                {
                    "obsid": np.full(len(nearby), obsid_i, dtype=np.int32),
                    "id": np.arange(len(nearby), dtype=np.int32) + id_offset,
                    "x_id": cd_sources["id"][xray_idx].astype(np.int32),
                    "catalog": np.full(len(nearby), "RFC"),
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

    if not all_rows:
        return Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)
    return vstack(all_rows)


def run_xcorr_for_obsids(
    select_name: str,
    obsids: np.ndarray,
    full_cat: Table,
    new_xray: Table,
    new_obs: Table,
) -> Table:
    """Re-run cross-matches for *select_name* for the given obsids using the full cat_src.

    Parameters
    ----------
    select_name:
        Cross-match set to run, e.g. ``'astromon_21'`` or ``'astromon_23'``.
        Controls which catalog priority order and match parameters are applied.
    obsids:
        Obsids to process.
    full_cat:
        Full astromon_cat_src including newly added RFC rows.
    new_xray:
        Full astromon_xray_src.
    new_obs:
        Full astromon_obs.

    Returns
    -------
    Table of xcorr rows with detect_method column.
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
                    select_name,
                    astromon_obs=obspar,
                    astromon_xray_src=sources,
                    astromon_cat_src=candidates,
                )
            except Exception as exc:
                logger.warning(
                    f"OBSID={obsid_i} {detect_method} {select_name} xcorr failed: {exc}"
                )
                continue

            if len(matches) == 0:
                continue

            matches["detect_method"] = detect_method
            all_xcorr.append(matches[xcorr_cols])

    if not all_xcorr:
        return Table(names=xcorr_cols)
    return vstack(all_xcorr)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument(
        "--catalog-cache",
        required=True,
        type=Path,
        metavar="PATH",
        help="Path to cache the RFC ASCII file. Downloaded on first run, reused thereafter.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    rfc = get_rfc(cache_path=args.catalog_cache, max_age_days=None)
    logger.info(f"  {len(rfc):,} RFC sources")

    logger.info(f"Loading tables from {args.db}")
    existing_cat = db.get_table("astromon_cat_src", dbfile=args.db)
    new_xray = db.get_table("astromon_xray_src", dbfile=args.db)
    new_obs = db.get_table("astromon_obs", dbfile=args.db)

    # Strip both "ICRS" (old name) and "RFC" (new name) for idempotency on re-runs.
    n_old = np.sum(np.isin(existing_cat["catalog"], ["ICRS", "RFC"]))
    logger.info(
        f"  Removing {n_old:,} existing ICRS/RFC cat_src rows before RFC rebuild"
    )
    non_rfc_cat = existing_cat[~np.isin(existing_cat["catalog"], ["ICRS", "RFC"])]

    logger.info("Pass 1: building RFC cat_src")
    new_rfc_cat = build_rfc_cat_src(rfc, new_xray, new_obs, non_rfc_cat)
    n_obsids = len(np.unique(new_rfc_cat["obsid"])) if len(new_rfc_cat) else 0
    logger.info(f"  {len(new_rfc_cat):,} RFC cat_src rows across {n_obsids:,} obsids")

    merged_cat = vstack([non_rfc_cat, new_rfc_cat], metadata_conflicts="silent")

    rfc_obsids = np.unique(new_rfc_cat["obsid"]) if len(new_rfc_cat) else np.array([])

    # Run xcorr for every select_name that includes RFC in its catalog list.
    # astromon_21 is the traditional RFC+Tycho2 set; astromon_23 is the full
    # new catalog set where RFC is the highest-priority source.  Both need to be
    # updated so that get_cross_matches() for either select_name sees RFC.
    # Filter by the "name" key to avoid running aliases (e.g. "default" → "astromon_21")
    # twice.  Use the args-dict "name" as the canonical select_name.
    seen_names: set[str] = set()
    select_names_with_rfc = []
    for sn, args_dict in cross_match.CROSS_MATCHES_ARGS.items():
        canonical = args_dict.get("name", sn)
        if "RFC" in args_dict.get("catalogs", []) and canonical not in seen_names:
            seen_names.add(canonical)
            select_names_with_rfc.append(canonical)
    logger.info(
        f"Pass 2: re-running xcorr for {len(rfc_obsids):,} obsids × "
        f"{select_names_with_rfc} (all select_names that include RFC)"
    )

    all_new_xcorr: list[Table] = []
    for sn in select_names_with_rfc:
        logger.info(f"  Running {sn}...")
        xcorr_sn = run_xcorr_for_obsids(sn, rfc_obsids, merged_cat, new_xray, new_obs)
        logger.info(f"    {len(xcorr_sn):,} rows")
        if len(xcorr_sn):
            all_new_xcorr.append(xcorr_sn)

    new_xcorr = (
        vstack(all_new_xcorr)
        if all_new_xcorr
        else Table(
            names=[
                "select_name",
                "obsid",
                "c_id",
                "x_id",
                "dy",
                "dz",
                "dr",
                "detect_method",
            ]
        )
    )
    logger.info(f"  {len(new_xcorr):,} total xcorr rows across all select_names")

    if args.dry_run:
        logger.info("  (dry-run: not writing to database)")
        return

    if not len(new_rfc_cat):
        logger.warning("No RFC cat_src rows produced — nothing to save")
        return

    db.save("astromon_cat_src", merged_cat, dbfile=args.db, ignore_obsid=True)
    logger.info(f"  Saved merged cat_src ({len(new_rfc_cat):,} RFC rows)")

    if len(new_xcorr):
        # select_name_key=True replaces only the rows for the affected select_names
        # and obsids; gaia_agn, gaia_var_star, and milliquas_gaia rows are untouched.
        db.save("astromon_xcorr", new_xcorr, dbfile=args.db, select_name_key=True)
        sn_counts = {
            sn: int(np.sum(np.array(new_xcorr["select_name"]) == sn))
            for sn in select_names_with_rfc
        }
        logger.info(f"  Saved xcorr: {sn_counts}")
    else:
        logger.info("  No new xcorr rows; xcorr table unchanged")


if __name__ == "__main__":
    main()
