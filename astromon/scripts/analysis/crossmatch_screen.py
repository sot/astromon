"""Bulk-screen xray sources against new Gaia-based catalogs to prioritise obsids for reprocessing.

Screens astromon_xray_src rows from the original database against GaiaAGN,
GaiaVarStar, and Milliquas v8 to identify which obsids have at least one
candidate match with the new catalogs added since the original pipeline ran.

All sources for all target obsids are queried in one batched pass per
catalog rather than obsid-by-obsid:

  - GaiaAGN / GaiaVarStar: batched Gaia TAP queries (~4 k round-trips at
    the default batch_size=50; ~450 at batch_size=500).
  - Milliquas v8: downloaded once from CDS, cached locally as a FITS file,
    then queried via astropy search_around_sky — no network after the first
    run.

Matches are mapped back to obsids via search_around_sky spatial joins, so
the result is a simple obsid-level flag: "has at least one new-catalog match
within <radius>".

Output files (written to --output-dir):
  screen_results.csv   — one row per obsid: obsid, n_xray_src, n_gaia_agn,
                         n_gaia_var_star, n_milliquas
  priority_obsids.txt  — obsids with ≥1 match, sorted by total match count
  no_match_obsids.txt  — obsids with 0 matches (lower priority for reprocessing)

Usage::

    python -m astromon.scripts.analysis.crossmatch_screen \\
        --obsid-list obsids_to_screen.txt \\
        --original-db /path/to/original/astromon.h5 \\
        --output-dir /path/to/output \\
        [--radius 3.0] \\
        [--gaia-batch-size 500]
"""

import argparse
import logging
import time
from pathlib import Path

import numpy as np
from astropy import coordinates as coords
from astropy import table
from astropy import units as u
from astropy.coordinates import search_around_sky
from ska_helpers.logging import basic_logger

from astromon import db
from astromon.cross_match import get_gaia_agn, get_gaia_var_stars, get_milliquas


def load_xray_sources(
    db_path: Path,
    obsids: list[int],
) -> table.Table:
    """Load astromon_xray_src rows for *obsids* from the original database.

    Parameters
    ----------
    db_path
        Path to the original astromon HDF5 database.
    obsids
        Target obsids; rows for any obsid not in the database are silently
        omitted.

    Returns
    -------
    table.Table
        Subset of astromon_xray_src with at least columns obsid, ra, dec.
    """
    with db.connect(db_path, mode="r") as con:
        xray = db.get_table("astromon_xray_src", con)
    mask = np.isin(xray["obsid"], obsids)
    return xray[mask]


def map_hits_to_obsids(
    all_pos: coords.SkyCoord,
    all_obsids: np.ndarray,
    hit_table: table.Table,
    radius: u.Quantity,
) -> dict[int, int]:
    """Return a count of catalog hits per obsid.

    For each row in *hit_table* (which has 'ra' and 'dec' columns), find
    all xray sources in *all_pos* within *radius* and credit their obsid.
    One catalog source can credit multiple obsids when the same sky position
    appears in several observations, which is intentional.

    Parameters
    ----------
    all_pos
        SkyCoord for all xray sources in the screening set.
    all_obsids
        Parallel array of obsids, one per element of *all_pos*.
    hit_table
        Catalog sources returned by a bulk query function.
    radius
        Search radius.

    Returns
    -------
    dict[int, int]
        Mapping of obsid → number of catalog sources within *radius* of any
        of that obsid's xray sources.  Obsids with no hits are absent.
    """
    if len(hit_table) == 0:
        return {}

    hit_pos = coords.SkyCoord(hit_table["ra"], hit_table["dec"], unit="deg")
    # search_around_sky(coords1, coords2, seplimit) → idx1, idx2, d2d, _
    # idx1[i] indexes into coords1 (all_pos); idx2[i] into coords2 (hit_pos)
    xray_idx, _hit_idx, _d2d, _ = search_around_sky(all_pos, hit_pos, radius)

    counts: dict[int, int] = {}
    for xi in xray_idx:
        obsid = int(all_obsids[xi])
        counts[obsid] = counts.get(obsid, 0) + 1
    return counts


def screen(
    xray_src: table.Table,
    radius: u.Quantity = 3 * u.arcsec,
    gaia_batch_size: int = 50,
    milliquas_cache: Path | None = None,
) -> table.Table:
    """Run the three-catalog bulk screen and return a per-obsid results table.

    Parameters
    ----------
    xray_src
        Combined astromon_xray_src rows for all target obsids.
    radius
        Match radius for all catalogs.  Default 3 arcsec.
    gaia_batch_size
        Positions per TAP CONTAINS-OR batch for Gaia queries.  Default 50
        (conservative).  500 gives ~10× fewer round-trips; increase if the
        TAP server accepts larger queries without timing out.
    milliquas_cache
        Override the default Milliquas cache path.

    Returns
    -------
    table.Table
        One row per unique obsid in *xray_src* with columns:
        obsid, n_xray_src, n_gaia_agn, n_gaia_var_star, n_milliquas.
    """
    logger = logging.getLogger("astromon")

    unique_obsids, counts_per_obs = np.unique(xray_src["obsid"], return_counts=True)
    n_src = len(xray_src)
    n_obs = len(unique_obsids)
    logger.info(
        f"Screening {n_src} xray sources across {n_obs} obsids "
        f"(radius={radius}, gaia_batch_size={gaia_batch_size})"
    )

    all_pos = coords.SkyCoord(xray_src["ra"], xray_src["dec"], unit="deg")
    all_obsids = np.asarray(xray_src["obsid"])

    # --- GaiaAGN (local cached catalog, no TAP) ---------------------------
    t0 = time.time()
    logger.info(f"GaiaAGN: searching {n_src} positions in the cached catalog...")
    agn_hits = get_gaia_agn(all_pos, radius=radius, logging_tag="screen")
    agn_counts = map_hits_to_obsids(all_pos, all_obsids, agn_hits, radius)
    logger.info(
        f"GaiaAGN done in {time.time() - t0:.0f}s: "
        f"{len(agn_hits)} hits across {len(agn_counts)} obsids"
    )

    # --- GaiaVarStar -------------------------------------------------------
    # Pass pos without obstime so the function skips PM correction (fine for
    # screening; full reprocessing applies per-obsid correction).
    t0 = time.time()
    logger.info(
        f"GaiaVarStar: querying {n_src} positions "
        f"(~{n_src // gaia_batch_size + 1} TAP batches)..."
    )
    var_hits = get_gaia_var_stars(
        all_pos, radius=radius, batch_size=gaia_batch_size, logging_tag="screen"
    )
    var_counts = map_hits_to_obsids(all_pos, all_obsids, var_hits, radius)
    logger.info(
        f"GaiaVarStar done in {time.time() - t0:.0f}s: "
        f"{len(var_hits)} hits across {len(var_counts)} obsids"
    )

    # --- Milliquas (local cache + search_around_sky) -----------------------
    t0 = time.time()
    logger.info("Milliquas: loading local catalog...")
    mq = get_milliquas(cache_path=milliquas_cache)
    mq_pos = coords.SkyCoord(mq["ra"], mq["dec"], unit="deg")
    logger.info(
        f"Milliquas: spatial join against {len(mq)} sources "
        f"and {n_src} xray positions..."
    )
    # search_around_sky(coords1, coords2, seplimit): find pairs within seplimit.
    # We want xray sources near Milliquas sources, so coords1=all_pos, coords2=mq_pos.
    xray_idx_mq, _mq_idx, _d2d_mq, _ = search_around_sky(all_pos, mq_pos, radius)
    mq_counts: dict[int, int] = {}
    for xi in xray_idx_mq:
        obsid = int(all_obsids[xi])
        mq_counts[obsid] = mq_counts.get(obsid, 0) + 1
    logger.info(
        f"Milliquas done in {time.time() - t0:.0f}s: "
        f"{len(xray_idx_mq)} hits across {len(mq_counts)} obsids"
    )

    # --- Assemble per-obsid result table -----------------------------------
    obs_n_src = dict(zip(unique_obsids.tolist(), counts_per_obs.tolist(), strict=True))
    rows = [
        {
            "obsid": int(obsid),
            "n_xray_src": obs_n_src.get(int(obsid), 0),
            "n_gaia_agn": agn_counts.get(int(obsid), 0),
            "n_gaia_var_star": var_counts.get(int(obsid), 0),
            "n_milliquas": mq_counts.get(int(obsid), 0),
        }
        for obsid in unique_obsids
    ]
    col_names = ["obsid", "n_xray_src", "n_gaia_agn", "n_gaia_var_star", "n_milliquas"]
    result = table.Table(rows=rows, names=col_names)
    result["n_total"] = (
        result["n_gaia_agn"] + result["n_gaia_var_star"] + result["n_milliquas"]
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--obsid-list",
        required=True,
        type=Path,
        help="Text file with one obsid per line (obsids to screen).",
    )
    parser.add_argument(
        "--original-db",
        required=True,
        type=Path,
        dest="original_db",
        help="Path to the original astromon HDF5 database (source of xray positions).",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        dest="output_dir",
        help="Directory for output files (created if absent).",
    )
    parser.add_argument(
        "--radius",
        default=3.0,
        type=float,
        help="Match radius in arcsec (default: 3.0).",
    )
    parser.add_argument(
        "--gaia-batch-size",
        default=50,
        type=int,
        dest="gaia_batch_size",
        help=(
            "Positions per Gaia TAP CONTAINS-OR query (default: 50). "
            "Increase to 200-500 to reduce round-trips; lower if TAP times out."
        ),
    )
    parser.add_argument(
        "--milliquas-cache",
        default=None,
        type=Path,
        dest="milliquas_cache",
        help="Override the default Milliquas FITS cache path.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        dest="log_level",
        help="Logging level (default: INFO).",
    )
    args = parser.parse_args()

    basic_logger(
        "astromon",
        level=args.log_level,
        format="%(asctime)s %(levelname)-8s %(message)s",
    )
    logger = logging.getLogger("astromon")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    obsids = [
        int(line.strip())
        for line in args.obsid_list.read_text().splitlines()
        if line.strip()
    ]
    logger.info(f"Loaded {len(obsids)} obsids from {args.obsid_list}")

    logger.info(f"Loading xray sources from {args.original_db}...")
    xray_src = load_xray_sources(args.original_db, obsids)
    loaded_obsids = {int(o) for o in np.unique(xray_src["obsid"])}
    missing = set(obsids) - loaded_obsids
    if missing:
        logger.warning(
            f"{len(missing)} obsids from list not found in original DB: "
            f"{sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}"
        )
    logger.info(f"Loaded {len(xray_src)} xray sources for {len(loaded_obsids)} obsids")

    result = screen(
        xray_src,
        radius=args.radius * u.arcsec,
        gaia_batch_size=args.gaia_batch_size,
        milliquas_cache=args.milliquas_cache,
    )

    # Write CSV with all results.
    csv_path = args.output_dir / "screen_results.csv"
    result.write(csv_path, format="ascii.csv", overwrite=True)
    logger.info(f"Wrote {len(result)} rows to {csv_path}")

    # Priority obsids: any new-catalog match, sorted by total matches desc.
    has_match = result["n_total"] > 0
    priority = result[has_match]
    priority.sort("n_total", reverse=True)
    priority_path = args.output_dir / "priority_obsids.txt"
    priority_path.write_text("\n".join(str(row["obsid"]) for row in priority) + "\n")
    logger.info(f"{len(priority)} obsids with ≥1 new-catalog match → {priority_path}")

    # Remaining obsids: no matches.
    no_match = result[~has_match]
    no_match_path = args.output_dir / "no_match_obsids.txt"
    no_match_path.write_text("\n".join(str(row["obsid"]) for row in no_match) + "\n")
    logger.info(f"{len(no_match)} obsids with no new-catalog match → {no_match_path}")

    n_agn = int(np.sum(result["n_gaia_agn"] > 0))
    n_var = int(np.sum(result["n_gaia_var_star"] > 0))
    n_mq = int(np.sum(result["n_milliquas"] > 0))
    logger.info(
        f"Summary: GaiaAGN hits in {n_agn} obsids, "
        f"GaiaVarStar in {n_var}, Milliquas in {n_mq}, "
        f"any-match in {len(priority)}"
    )


if __name__ == "__main__":
    main()
