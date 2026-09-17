"""Extract small, shareable samples of cross-matches for the new catalogs.

The validation notebook (``catalog_matches_validation.ipynb``) is meant to run without
a live connection to the production astromon DB. This script pulls two things from a
live DB and writes them to compact ECSV files that can be hosted for the notebook to
fetch:

1. A capped random sample of the precomputed single-catalog cross-matches for each new
   catalog added in this branch (RFC, ICRF3, GaiaAGN, GaiaQSO, Quaia, MilliquasGaia,
   DESIV161, GaiaVarStar), plus the pre-existing Tycho2 selection as a familiar
   baseline (``extract_sample`` / ``--out``).
2. A "paired" sample: for each catalog, the subset of catalog counterparts (c_id) that
   were matched independently by *both* a celldetect and a gaussian_detect X-ray
   source in the same obsid. Comparing dr between the two rows in a pair -- same
   physical source, same counterpart, different detection/centroiding algorithm --
   isolates which detection method centroids closer to the true position
   (``extract_paired_sample`` / ``--paired-out``).

The main sample also carries ``dy_rebased``/``dz_rebased``: dy/dz rebased to the single
most-recently-dated CALALIGN alignment matrix per detector, so a time series doesn't mix
real signal with step changes caused purely by CALALIGN updates (see
:any:`astromon.utils.get_rebased_offsets` and ``--calalign-dir``).

The astromon_xcorr rows for each of these ``select_name`` values are already the
*filtered* matches (dr, snr, r_angle, near_neighbor_dist cuts baked in at the time
:any:`astromon.cross_match.compute_cross_matches` populated the DB) -- see
CROSS_MATCHES_ARGS in astromon/cross_match.py. This script does not re-filter; it only
samples down to a shareable size.

Usage
-----
    python notebooks/extract_catalog_validation_sample.py [--n-per-catalog 400]
        [--n-pairs-per-catalog 400] [--out PATH] [--paired-out PATH]

Requires ASTROMON_FILE (or SKA) to point at a DB with these select_names already
populated.
"""

import argparse
from collections import defaultdict

import numpy as np
from astropy import table

import astromon
from astromon import utils

# New catalogs introduced in this branch, plus Tycho2 as a familiar baseline. Ordered to
# match cross_match._MATCH_HIERARCHY exactly: ICRF3 before RFC, since ICRF3 is a subset of
# RFC and has to be tried first in a combined selection or it would never win a match.
# ICRF2 is not a new catalog -- it's the pre-this-branch reference catalog, kept only as
# the "old baseline" comparison point (see compute_match_gains_summary).
SELECT_NAMES = [
    "icrf3",
    "rfc",
    "gaia_agn",
    "quaia",
    "desi_v161",
    "milliquas_gaia",
    "gaia_qso",
    "gaia_var_star",
    "tycho2",
    "icrf2",
]

# Columns to keep in the sample -- enough to judge match quality (dr, snr, r_angle,
# near_neighbor_dist, mag) and to group/label plots (catalog, obsid, target, category,
# detect_method), without carrying the full joined row (~45 columns).
SAMPLE_COLUMNS = [
    "select_name",
    "catalog",
    "name",
    "obsid",
    "target",
    "category",
    "detect_method",
    "detector",
    "grating",
    "date_obs",
    "tstart",
    "x_id",
    "c_id",
    "dr",
    "dy",
    "dz",
    "snr",
    "r_angle",
    "near_neighbor_dist",
    "net_counts",
    "mag",
    "x_ra",
    "x_dec",
    "c_ra",
    "c_dec",
    "pileup",
    "psfratio",
    "concentration_ratio",
    "peak_offset",
    "dy_rebased",
    "dz_rebased",
]

RANDOM_SEED = 20260825


# Per-source columns kept for each side (celldetect/gaussian_detect) of a pair, before
# the _cell/_gauss suffix is appended.
PAIR_VALUE_COLUMNS = [
    "dr",
    "dy",
    "dz",
    "snr",
    "r_angle",
    "near_neighbor_dist",
    "net_counts",
    "x_ra",
    "x_dec",
]


def _build_pairs_for_catalog(matches: table.Table) -> table.Table:
    """Pivot one catalog's matches into celldetect/gaussian_detect pairs on (obsid, c_id).

    Only (obsid, c_id) counterparts matched by both detect methods are kept -- these
    are the same physical X-ray source and the same catalog counterpart, detected
    independently by each method, so the two dr values are directly comparable.
    """
    by_key = defaultdict(dict)
    for row in matches:
        by_key[(int(row["obsid"]), int(row["c_id"]))][row["detect_method"]] = row

    pair_rows = []
    for (obsid, c_id), by_method in by_key.items():
        if "celldetect" not in by_method or "gaussian_detect" not in by_method:
            continue
        cell, gauss = by_method["celldetect"], by_method["gaussian_detect"]
        pair_rows.append(
            (
                cell["select_name"],
                cell["catalog"],
                obsid,
                c_id,
                cell["target"],
                cell["tstart"],
                cell["mag"],
                cell["c_ra"],
                cell["c_dec"],
                *[cell[col] for col in PAIR_VALUE_COLUMNS],
                *[gauss[col] for col in PAIR_VALUE_COLUMNS],
            )
        )

    names = [
        "select_name",
        "catalog",
        "obsid",
        "c_id",
        "target",
        "tstart",
        "mag",
        "c_ra",
        "c_dec",
        *[f"{col}_cell" for col in PAIR_VALUE_COLUMNS],
        *[f"{col}_gauss" for col in PAIR_VALUE_COLUMNS],
    ]
    return (
        table.Table(rows=pair_rows, names=names)
        if pair_rows
        else table.Table(names=names)
    )


def extract_paired_sample(n_pairs_per_catalog: int, dbfile=None) -> table.Table:
    """Return a capped random sample of celldetect/gaussian_detect pairs per catalog."""
    rng = np.random.default_rng(RANDOM_SEED)
    samples = []
    for select_name in SELECT_NAMES:
        matches = astromon.get_cross_matches(name=select_name, dbfile=dbfile)
        pairs = _build_pairs_for_catalog(matches)
        if len(pairs) > n_pairs_per_catalog:
            keep = rng.choice(len(pairs), size=n_pairs_per_catalog, replace=False)
            keep.sort()
            pairs = pairs[keep]
        samples.append(pairs)
        print(f"{select_name:16s} kept {len(pairs)} celldetect/gaussian_detect pairs")
    return table.vstack(samples)


def extract_sample(n_per_catalog: int, calalign_dir, dbfile=None) -> table.Table:
    """Return a capped random sample of matches for each catalog in SELECT_NAMES."""
    rng = np.random.default_rng(RANDOM_SEED)
    samples = []
    for select_name in SELECT_NAMES:
        matches = astromon.get_cross_matches(name=select_name, dbfile=dbfile)
        matches = utils.get_rebased_offsets(matches, calalign_dir=calalign_dir)
        matches = matches[SAMPLE_COLUMNS]
        if len(matches) > n_per_catalog:
            keep = rng.choice(len(matches), size=n_per_catalog, replace=False)
            keep.sort()
            matches = matches[keep]
        samples.append(matches)
        print(f"{select_name:16s} kept {len(matches)} of the precomputed matches")
    return table.vstack(samples)


# The 8 catalogs added in this branch, each as its own independent single-catalog
# selection -- deliberately not astromon_23 (the combined, priority-ordered
# hierarchy): a plain union across independent selections says "does this source
# have a counterpart in any new catalog" without needing to reason about which
# catalog wins when several would match the same source.
NEW_CATALOG_SELECT_NAMES = [
    "icrf3",
    "rfc",
    "gaia_agn",
    "quaia",
    "desi_v161",
    "milliquas_gaia",
    "gaia_qso",
    "gaia_var_star",
]


def compute_match_gains_summary(dbfile=None) -> table.Table:
    """How many X-ray sources gain a counterpart that had none under the old baseline.

    "Old baseline" is the pre-this-branch reference catalog set: ICRF2 union Tycho2
    (matches the historical astromon_21 catalogs=["ICRS", "Tycho2"] -- see
    docs/quick.rst's pre-branch example -- with "ICRS" being what ICRF2 was called
    there). "New" is the union of the 8 independent single-catalog selections for
    the catalogs added in this branch (NEW_CATALOG_SELECT_NAMES) -- not astromon_23,
    so this doesn't depend on the combined hierarchy's catalog priority order at all.

    A source counts as "gained" only if it has *no* counterpart at all under the old
    baseline -- not merely a different one -- so this measures new astrometric
    coverage, not just a catalog-label change for a source that was already usable.
    Uses the full DB (not a capped sample): these are small summary counts, not
    per-match rows, so there is no size/shareability reason to cap them.

    Returns
    -------
    table.Table
        One row per new catalog plus a "TOTAL" row, with columns catalog, n_gained,
        n_new_total, pct_of_new. Per-catalog rows can double-count: a source with no
        old-baseline counterpart that matches two new catalogs at once (e.g. both
        GaiaAGN and Quaia independently) is counted in both those rows, so they need
        not sum to TOTAL -- TOTAL is the distinct count, not a per-catalog sum.
    """
    icrf2 = astromon.get_cross_matches(name="icrf2", dbfile=dbfile)
    tycho2 = astromon.get_cross_matches(name="tycho2", dbfile=dbfile)

    def _match_keys(matches):
        return set(
            zip(
                matches["obsid"], matches["x_id"], matches["detect_method"], strict=True
            )
        )

    old_ids = _match_keys(icrf2) | _match_keys(tycho2)

    new_matches_by_catalog = {
        name: astromon.get_cross_matches(name=name, dbfile=dbfile)
        for name in NEW_CATALOG_SELECT_NAMES
    }
    new_ids_by_catalog = {
        name: _match_keys(m) for name, m in new_matches_by_catalog.items()
    }
    new_ids = set().union(*new_ids_by_catalog.values())
    gained_ids = new_ids - old_ids

    rows = [
        (m["catalog"][0], len(ids - old_ids))
        for name, ids in new_ids_by_catalog.items()
        if len(ids) > 0
        for m in [new_matches_by_catalog[name]]
    ]
    rows.sort(key=lambda cn: -cn[1])
    rows.append(("TOTAL (distinct sources)", len(gained_ids)))

    summary = table.Table(
        rows=[(c, n, len(new_ids), 100 * n / len(new_ids)) for c, n in rows],
        names=["catalog", "n_gained", "n_new_total", "pct_of_new"],
    )
    summary["pct_of_new"].format = "%.1f"
    print(
        f"Old baseline (ICRF2 union Tycho2): {len(old_ids):,} matched sources\n"
        f"New (union of the 8 new-catalog selections): {len(new_ids):,} matched sources\n"
        f"Genuinely new (no old-baseline counterpart at all): {len(gained_ids):,}"
        f" ({100 * len(gained_ids) / len(new_ids):.1f}% of the new-catalog union)"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n-per-catalog",
        type=int,
        default=400,
        help="Maximum rows to keep per catalog (default: 400)",
    )
    parser.add_argument(
        "--out",
        default="notebooks/data/catalog_match_validation_sample.ecsv",
        help="Output ECSV path",
    )
    parser.add_argument(
        "--n-pairs-per-catalog",
        type=int,
        default=400,
        help="Maximum celldetect/gaussian_detect pairs to keep per catalog (default: 400)",
    )
    parser.add_argument(
        "--paired-out",
        default="notebooks/data/catalog_match_paired_sample.ecsv",
        help="Output ECSV path for the celldetect/gaussian_detect paired sample",
    )
    parser.add_argument(
        "--dbfile",
        default=None,
        help="Path to the astromon DB file (default: $ASTROMON_FILE or"
        " $SKA/data/astromon/astromon.h5)",
    )
    parser.add_argument(
        "--calalign-dir",
        default=None,
        help="Directory of CALALIGN FITS files, for rebasing dy/dz to the latest"
        " alignment epoch (default: /data/caldb/data/chandra/pcad/align, per"
        " astromon.utils.calalign_from_files)",
    )
    parser.add_argument(
        "--gains-out",
        default="notebooks/data/catalog_match_gains_summary.ecsv",
        help="Output ECSV path for the new-match-gains summary",
    )
    args = parser.parse_args()

    sample = extract_sample(args.n_per_catalog, args.calalign_dir, dbfile=args.dbfile)
    sample.write(args.out, overwrite=True)
    print(f"\nWrote {len(sample)} rows to {args.out}")

    pairs = extract_paired_sample(args.n_pairs_per_catalog, dbfile=args.dbfile)
    pairs.write(args.paired_out, overwrite=True)
    print(f"\nWrote {len(pairs)} rows to {args.paired_out}")

    gains = compute_match_gains_summary(dbfile=args.dbfile)
    gains.write(args.gains_out, overwrite=True)
    print(f"\nWrote {len(gains)} rows to {args.gains_out}")


if __name__ == "__main__":
    main()
