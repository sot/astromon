"""Build per-catalog cross-match residuals independent of the hierarchy.

The standard astromon cross-match assigns each X-ray source to the
highest-priority catalog that has a nearby counterpart.  This makes it
impossible to compare catalog quality fairly: RFC wins contested sources,
Tycho2 only sees the leftovers, etc.

This script instead matches every catalog independently: for each X-ray
source (celldetect only, since astromon_cat_src.x_id references celldetect
sources) it finds the nearest counterpart in *each* catalog and records the
residuals.  The result is a CSV with one row per (obsid, x_id, catalog) so
you can compare RFC vs Tycho2 vs GaiaAGN on the same footing.

Usage
-----
    python -m astromon.scripts.analysis.solo_catalog_matches \\
        --db ~/git/astromon_hrc/astromon_part.h5 \\
        --out solo_catalog_matches.csv \\
        [--sep 3.0] [--snr 3.0]
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from astropy.table import Table, join
from cxotime import CxoTime

from astromon import db
from astromon.cross_match import get_bad_source_mask

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("solo_catalog_matches")


def build_solo_matches(
    cat_src: Table,
    xray_src: Table,
    obs: Table,
    max_sep_arcsec: float = 3.0,
    min_snr: float = 3.0,
) -> Table:
    """Compute per-catalog residuals without the priority hierarchy.

    Parameters
    ----------
    cat_src
        Full ``astromon_cat_src`` table.
    xray_src
        Full ``astromon_xray_src`` table (celldetect rows only are used, since
        ``cat_src.x_id`` references celldetect source ids).
    obs
        Full ``astromon_obs`` table.
    max_sep_arcsec
        Maximum angular separation (arcsec) to include a match.
    min_snr
        Minimum X-ray source SNR.

    Returns
    -------
    Table with columns:
        obsid, x_id, catalog, cat_name, sep_arcsec, dy, dz, dr, snr,
        psfratio, concentration_ratio, near_neighbor_dist, pileup, r_angle,
        target, tstart, year, is_outlier, is_bad_source
    """
    # Restrict cat_src to matches within the search radius.
    cat_sep_ok = cat_src[np.array(cat_src["separation"]) <= max_sep_arcsec]
    logger.info(
        f'cat_src rows within {max_sep_arcsec}" sep: {len(cat_sep_ok):,} '
        f"(of {len(cat_src):,})"
    )

    # Keep only celldetect xray sources.
    cd_xray = xray_src[np.asarray(xray_src["detect_method"]) == "celldetect"]
    logger.info(f"celldetect xray sources: {len(cd_xray):,}")

    # Build fast lookup: (obsid, x_id) -> xray row index.
    cd_obsids = np.array(cd_xray["obsid"], dtype=np.int32)
    cd_ids = np.array(cd_xray["id"], dtype=np.int32)
    xray_index: dict[tuple[int, int], int] = {
        (int(o), int(i)): idx
        for idx, (o, i) in enumerate(zip(cd_obsids, cd_ids, strict=True))
    }

    # For each cat_src row, look up the matched xray source and compute dy/dz.
    cat_obsids = np.array(cat_sep_ok["obsid"], dtype=np.int32)
    cat_x_ids = np.array(cat_sep_ok["x_id"], dtype=np.int32)
    cat_y = np.array(cat_sep_ok["y_angle"], dtype=np.float32)
    cat_z = np.array(cat_sep_ok["z_angle"], dtype=np.float32)
    cat_catalogs = np.asarray(cat_sep_ok["catalog"])
    cat_names = np.asarray(cat_sep_ok["name"])

    matched_rows: list[dict] = []
    n_nomatch = 0

    for i in range(len(cat_sep_ok)):
        key = (int(cat_obsids[i]), int(cat_x_ids[i]))
        xray_idx = xray_index.get(key)
        if xray_idx is None:
            n_nomatch += 1
            continue

        xray_row = cd_xray[xray_idx]
        snr = float(xray_row["snr"])
        if snr < min_snr:
            continue

        dy = float(xray_row["y_angle"]) - float(cat_y[i])
        dz = float(xray_row["z_angle"]) - float(cat_z[i])
        dr = float(np.hypot(dy, dz))

        matched_rows.append(
            {
                "obsid": int(cat_obsids[i]),
                "x_id": int(cat_x_ids[i]),
                "catalog": cat_catalogs[i],
                "cat_name": cat_names[i],
                "sep_arcsec": float(cat_sep_ok["separation"][i]),
                "dy": dy,
                "dz": dz,
                "dr": dr,
                "snr": snr,
                "psfratio": float(xray_row["psfratio"]),
                # concentration_ratio: fraction of counts within 2" vs 10".
                # Near 1 = point source; near 0 = extended.  Available for all
                # celldetect rows.
                "concentration_ratio": float(xray_row["concentration_ratio"]),
                # near_neighbor_dist: arcsec to nearest other detected source.
                # Low values flag crowded/confused detections.
                "near_neighbor_dist": float(xray_row["near_neighbor_dist"]),
                # pileup: per-source pile-up fraction estimate.
                "pileup": float(xray_row["pileup"]),
                # r_angle: off-axis angle in arcmin; PSF degrades with radius.
                "r_angle": float(xray_row["r_angle"]),
            }
        )

    logger.info(
        f"Matched {len(matched_rows):,} pairs; "
        f"{n_nomatch:,} cat_src rows had no celldetect xray source"
    )

    if not matched_rows:
        raise RuntimeError("No matches found — check db path and filter thresholds")

    matches = Table(rows=matched_rows)

    # For each (obsid, x_id, catalog) keep only the nearest match.
    # (Multiple catalog sources in the same catalog can share an x_id.)
    matches.sort("dr")
    seen: set[tuple[int, int, str]] = set()
    keep = np.zeros(len(matches), dtype=bool)
    for i, row in enumerate(matches):
        key = (int(row["obsid"]), int(row["x_id"]), str(row["catalog"]))
        if key not in seen:
            seen.add(key)
            keep[i] = True
    matches = matches[keep]
    logger.info(
        f"After deduplication (nearest per obsid/x_id/catalog): {len(matches):,} rows"
    )

    # Join obs info: target, tstart, year.
    obs_cols = Table(
        {
            "obsid": np.array(obs["obsid"], dtype=np.int32),
            "target": np.asarray(obs["target"]),
            "tstart": np.array(obs["tstart"], dtype=np.float64),
        }
    )
    matches = join(matches, obs_cols, keys="obsid", join_type="left")
    matches["year"] = CxoTime(matches["tstart"]).decimalyear

    # Bad-source flag (uses target-prefix, cat_source, obsid_src exclusions).
    # get_bad_source_mask expects a "name" column holding the individual catalog
    # source name (e.g. "MCG 7-41-003"), not the catalog type string.
    matches_for_mask = matches.copy()
    matches_for_mask["name"] = matches_for_mask["cat_name"]
    matches["is_bad_source"] = get_bad_source_mask(matches_for_mask)

    # Outlier flag consistent with the 1-arcsec threshold used in the main plots.
    matches["is_outlier"] = matches["dr"] > 1.0

    return matches


def main() -> None:
    default_db = Path.home() / "ska" / "data" / "astromon" / "astromon.h5"

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--db", type=Path, default=default_db)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("solo_catalog_matches.csv"),
        help="Output CSV path (default: solo_catalog_matches.csv)",
    )
    parser.add_argument(
        "--sep",
        type=float,
        default=3.0,
        metavar="ARCSEC",
        help="Maximum separation to include (default: 3.0 arcsec)",
    )
    parser.add_argument(
        "--snr",
        type=float,
        default=3.0,
        help="Minimum X-ray source SNR (default: 3.0)",
    )
    args = parser.parse_args()

    if not args.db.exists():
        parser.error(f"DB not found: {args.db}")

    logger.info(f"Loading tables from {args.db}")
    cat_src = db.get_table("astromon_cat_src", dbfile=args.db)
    xray_src = db.get_table("astromon_xray_src", dbfile=args.db)
    obs = db.get_table("astromon_obs", dbfile=args.db)
    logger.info(
        f"  cat_src: {len(cat_src):,}  xray_src: {len(xray_src):,}  obs: {len(obs):,}"
    )

    matches = build_solo_matches(
        cat_src,
        xray_src,
        obs,
        max_sep_arcsec=args.sep,
        min_snr=args.snr,
    )

    # Print per-catalog summary.
    catalogs = sorted(set(matches["catalog"]))
    logger.info("Per-catalog summary (excluding bad sources):")
    logger.info(
        f"  {'Catalog':<20} {'n':>6}  {'MAD_dy':>8}  {'MAD_dz':>8}  {'outlier%':>9}"
    )
    for cat in catalogs:
        ok = (matches["catalog"] == cat) & ~matches["is_bad_source"]
        sub = matches[ok]
        if len(sub) == 0:
            continue
        mad_dy = float(np.median(np.abs(sub["dy"] - np.median(sub["dy"]))))
        mad_dz = float(np.median(np.abs(sub["dz"] - np.median(sub["dz"]))))
        outlier_pct = 100.0 * float(np.mean(sub["is_outlier"]))
        logger.info(
            f"  {cat:<20} {len(sub):>6}  {mad_dy:>8.3f}  {mad_dz:>8.3f}  {outlier_pct:>8.1f}%"
        )

    matches.write(args.out, format="csv", overwrite=True)
    logger.info(f"Wrote {len(matches):,} rows to {args.out}")


if __name__ == "__main__":
    main()
