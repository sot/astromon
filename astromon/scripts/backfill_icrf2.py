"""Backfill ICRF2 as a reference-baseline catalog, across every obsid already in the DB.

ICRF2 (Ma et al. 2009, IERS TN 35; VizieR I/323) was fully retired from the pipeline
when ICRF3/RFC replaced it -- see :any:`astromon.cross_match.get_icrf2`, which brings
it back purely as a reference baseline for measuring how many additional matches the
new catalogs (RFC, ICRF3, GaiaAGN, ...) get us relative to the pre-this-branch catalog
set. This is *not* part of any combined hierarchy (astromon_23-27), so those do not
need reprocessing -- only the standalone "icrf2" single-catalog selection is backfilled
here.

This works entirely from tables already in the DB (astromon_obs, astromon_xray_src) --
no live observation reprocessing, no CIAO. For each obsid:

  1. Rough-match ICRF2 against that obsid's celldetect sources (mirrors the RFC/ICRF3
     rough_match call in astromon/scripts/get_cat_obs_data.py).
  2. For each detect_method present, remap the rough-matched candidates' x_id to that
     method's own source ids (:any:`astromon.cross_match.remap_x_id_to_sources`) before
     cross-matching -- required because celldetect and gaussian_detect each number
     their sources independently per obsid, so a candidate's x_id from step 1 (matched
     to celldetect ids) does not carry over to gaussian_detect's ids. See the
     remap_x_id_to_sources docstring: skipping this silently produces wrong or missing
     matches for every detect_method beyond the one used to build the rough match.
  3. Cross-match with :any:`astromon.cross_match.compute_cross_matches` (name="icrf2").

Usage
-----
    python -m astromon.scripts.backfill_icrf2 --dbfile PATH [--dry-run] [--limit N]
"""

import argparse
import logging

import numpy as np
from astropy.table import Table, vstack
from chandra_aca.transform import radec_to_yagzag
from cxotime import CxoTime
from Quaternion import Quat

from astromon import cross_match, db

logger = logging.getLogger("backfill_icrf2")


def _build_icrf2_candidates(obsid: int, obs_row, celldetect_sources: Table) -> Table:
    """Rough-match ICRF2 against one obsid's celldetect sources.

    Returns an astromon_cat_src-shaped Table (catalog="ICRF2"), x_id matched to
    `celldetect_sources`; empty (but correctly shaped) if there are no candidates.
    """
    obs_time = CxoTime(str(obs_row["date_obs"]))
    candidates = cross_match.rough_match(
        celldetect_sources, obs_time, catalogs=("ICRF2",)
    )
    if len(candidates) == 0:
        return candidates

    q = Quat(
        equatorial=(float(obs_row["ra"]), float(obs_row["dec"]), float(obs_row["roll"]))
    )
    candidates["obsid"] = obsid
    candidates["id"] = np.arange(len(candidates))
    candidates["y_angle"], candidates["z_angle"] = radec_to_yagzag(
        candidates["ra"], candidates["dec"], q
    )
    return candidates


def backfill(dbfile, limit: int | None = None) -> tuple[Table, Table]:
    """Build ICRF2 cat_src and icrf2 xcorr rows for every obsid already in `dbfile`.

    Returns
    -------
    (cat_src, xcorr): the new rows to save, not yet written to the DB.
    """
    astromon_obs = db.get_table("astromon_obs", dbfile)
    astromon_xray_src = db.get_table("astromon_xray_src", dbfile)
    obs_by_id = {int(row["obsid"]): row for row in astromon_obs}

    obsids = np.unique(astromon_xray_src["obsid"])
    if limit is not None:
        obsids = obsids[:limit]

    xcorr_cols = list(db.ASTROMON_XCORR_DTYPE.names)
    all_cat_src: list[Table] = []
    all_xcorr: list[Table] = []
    n_obsids_with_candidates = 0

    for i, raw_obsid in enumerate(obsids):
        obsid = int(raw_obsid)
        obs_row = obs_by_id.get(obsid)
        if obs_row is None:
            logger.warning(f"OBSID={obsid} missing from astromon_obs, skipping")
            continue

        obsid_xray = astromon_xray_src[astromon_xray_src["obsid"] == obsid]
        celldetect_sources = obsid_xray[obsid_xray["detect_method"] == "celldetect"]
        if len(celldetect_sources) == 0:
            continue

        candidates = _build_icrf2_candidates(obsid, obs_row, celldetect_sources)
        if len(candidates) == 0:
            continue
        all_cat_src.append(candidates)
        n_obsids_with_candidates += 1

        obspar = Table([obs_row])
        for detect_method in np.unique(obsid_xray["detect_method"]):
            sources = obsid_xray[obsid_xray["detect_method"] == detect_method]
            candidates_for_version = cross_match.remap_x_id_to_sources(
                candidates, sources
            )
            try:
                matches = cross_match.compute_cross_matches(
                    "icrf2",
                    astromon_obs=obspar,
                    astromon_xray_src=sources,
                    astromon_cat_src=candidates_for_version,
                    logging_tag=f"OBSID={obsid}",
                )
            except Exception as exc:
                logger.warning(
                    f"OBSID={obsid} {detect_method} icrf2 xcorr failed: {exc}"
                )
                continue
            if len(matches) == 0:
                continue
            matches["detect_method"] = detect_method
            all_xcorr.append(matches[xcorr_cols])

        if (i + 1) % 200 == 0:
            logger.info(f"...{i + 1}/{len(obsids)} obsids processed")

    cat_src = (
        vstack(all_cat_src) if all_cat_src else Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)
    )
    xcorr = vstack(all_xcorr) if all_xcorr else Table(names=xcorr_cols)
    logger.info(
        f"{n_obsids_with_candidates:,}/{len(obsids):,} obsids had an ICRF2 candidate"
    )
    return cat_src, xcorr


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dbfile", required=True, help="Path to the astromon DB file")
    parser.add_argument(
        "--dry-run", action="store_true", help="Build the rows but do not save them"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N obsids (for testing)",
    )
    args = parser.parse_args()

    cat_src, xcorr = backfill(args.dbfile, limit=args.limit)
    n_obsids_cat = len(np.unique(cat_src["obsid"])) if len(cat_src) else 0
    n_obsids_xcorr = len(np.unique(xcorr["obsid"])) if len(xcorr) else 0
    logger.info(
        f"Built {len(cat_src):,} ICRF2 cat_src rows across {n_obsids_cat:,} obsids"
    )
    logger.info(
        f"Built {len(xcorr):,} icrf2 xcorr rows across {n_obsids_xcorr:,} obsids"
    )

    if args.dry_run:
        logger.info("(dry-run: not writing to database)")
        return

    if len(cat_src):
        db.save(
            "astromon_cat_src", cat_src, dbfile=args.dbfile, replace_keys=("catalog",)
        )
        logger.info("Saved cat_src (ICRF2 rows replaced, other catalogs untouched)")
    else:
        logger.warning("No ICRF2 cat_src rows produced -- nothing to save")

    if len(xcorr):
        db.save("astromon_xcorr", xcorr, dbfile=args.dbfile, select_name_key=True)
        logger.info("Saved xcorr (icrf2 select_name updated, others untouched)")
    else:
        logger.warning("No icrf2 xcorr rows produced -- xcorr table unchanged")


if __name__ == "__main__":
    main()
