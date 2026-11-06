"""Recompute y_angle/z_angle/r_angle that disagree with the obsid's evt2 pointing.

Why this is needed
------------------
Source angles have to come from the pointing in the evt2 file's own header
(RA_PNT/DEC_PNT/ROLL_PNT). They used to come from mica's quick-look obspar,
which differs from it by up to ~40" (see ``reprocess_pointing.py``). That was
fixed, and the database was recomputed, but ``Observation.get_sources`` is a
``@stored_result`` whose cache key is only its arguments. A ``sources__*.fits.gz``
cache written with the old pointing is therefore still a valid hit, and every
later run that reused it saved the old angles back to the database. RA/Dec are
not affected: they come from the detection's own WCS.

The result is rows whose RA/Dec are right and whose y/z are off by a constant
for the whole (obsid, detect_method), so ``simple_cross_match``, which pairs on
y/z, can match a catalog source to the wrong X-ray source or miss the right one.

What it does
------------
For every obsid in ``astromon_xray_src`` and ``astromon_cat_src`` (or the
``--obsid-file`` subset), reads the pointing from the obsid's local evt2 header
and finds the groups -- (obsid, detect_method) for X-ray sources, obsid for
catalog sources, i.e. the keys ``db.save`` replaces by -- with any row whose
stored y/z differ by more than ``--tolerance`` from the ones its RA/Dec give.
Only those groups are rewritten, with y/z/r_angle recomputed from RA/Dec. Every
other column, ids included, is kept, so astromon_cat_src.celldetect_x_id and
astromon_xcorr's x_id/c_id still point at the same rows.

The pointing is read straight from the header (``<root>/obsNN/<obsid>/primary``,
the archive ``*_evt2.fits*`` first, then filter_events' ``*_evt2_filtered``
copy, which keeps the header). It does not go through
``requery_cat_src.evt2_pointing``: that builds an ``Observation``, which writes
an evt2_info cache into the work dir, may try to download, and does not read
the filtered copy -- which is why ``fix_cat_src_pointing.py`` skipped work dirs
that only have one. Obsids with no local evt2 are reported and left alone.

What it does not do
-------------------
- Rebuild astromon_xcorr. Its dy/dz/dr, and which pairs were matched at all,
  came from the stale angles. Run ``rebuild_xcorr.py --obsids <touched>`` on
  the obsids this writes to ``--touched-out``.
- Re-apply ``get_sources``' r_angle < 180" cut. A row that now lands at or past
  180" stays, and is reported; a source the stale pointing put past 180" was
  never stored and is not recovered here.

Usage
-----
::

    python -m astromon.scripts.maintenance.fix_pointing_angles \\
        --db astromon_copy.h5 --root /Volumes/Black/data/astromon/work \\
        --touched-out touched_obsids.txt [--dry-run]
"""

import argparse
import logging
import zlib
from collections.abc import Iterator, Sequence
from pathlib import Path

import numpy as np
from astropy.table import Table
from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat

from astromon import db
from astromon.scripts.maintenance.backfill_obs_gti import obsid_work_dir
from astromon.utils import events_header

logger = logging.getLogger("fix_pointing_angles")

#: (RA_PNT, DEC_PNT, ROLL_PNT) in degrees.
Pointing = tuple[float, float, float]

#: Files in ``primary/`` to read the pointing from, in order of preference. The
#: first is what Observation.get_evt2_info reads; the filtered copy keeps its header.
EVT2_PATTERNS = ("*_evt2.fits*", "*_evt2_filtered.fits*")

#: Tables repaired, each with the key db.save replaces its rows by.
TABLE_KEYS = {
    "astromon_xray_src": ("obsid", "detect_method"),
    "astromon_cat_src": ("obsid",),
}

#: Stored angles are float32, so recomputing a correct row reproduces it to ~1e-4".
DEFAULT_TOLERANCE = 0.1  # arcsec


def local_evt2_pointing(obsid: int, roots: Sequence[Path]) -> Pointing | None:
    """`obsid`'s pointing from the first local evt2 found under `roots`, or None.

    Example::

        work = Path("/Volumes/Black/data/astromon/work")
        ra_pnt, dec_pnt, roll_pnt = local_evt2_pointing(27059, [work])
    """
    for pattern in EVT2_PATTERNS:
        for root in roots:
            primary = obsid_work_dir(root, obsid) / "primary"
            files = sorted(primary.glob(pattern)) if primary.is_dir() else []
            if not files:
                continue
            try:
                header = events_header(files[0])
                return (
                    float(header["RA_PNT"]),
                    float(header["DEC_PNT"]),
                    float(header["ROLL_PNT"]),
                )
            # EOFError and zlib.error: a gzip cut off or corrupted mid-download
            except (OSError, EOFError, zlib.error, ValueError, KeyError) as exc:
                logger.warning(
                    f"OBSID={obsid} cannot read pointing from {files[0]}: {exc}"
                )
    return None


def _rows_by_obsid(obsids: np.ndarray) -> Iterator[tuple[int, np.ndarray]]:
    """(obsid, row indices) for each obsid, without a full scan per obsid."""
    order = np.argsort(obsids, kind="stable")
    unique, starts = np.unique(obsids[order], return_index=True)
    stops = np.append(starts[1:], len(order))
    for obsid, start, stop in zip(unique, starts, stops, strict=True):
        yield int(obsid), order[start:stop]


def _angles_from_radec(
    table: Table, pointings: dict[int, Pointing]
) -> tuple[np.ndarray, np.ndarray]:
    """y/z from each row's ra/dec and its obsid's pointing; NaN where there is none."""
    y_angle = np.full(len(table), np.nan)
    z_angle = np.full(len(table), np.nan)
    ra = np.asarray(table["ra"], dtype=float)
    dec = np.asarray(table["dec"], dtype=float)
    for obsid, rows in _rows_by_obsid(np.asarray(table["obsid"])):
        pointing = pointings.get(obsid)
        if pointing is None:
            continue
        y_angle[rows], z_angle[rows] = radec_to_yagzag(
            ra[rows], dec[rows], Quat(equatorial=pointing)
        )
    return y_angle, z_angle


def angle_errors(table: Table, pointings: dict[int, Pointing]) -> np.ndarray:
    """Per row, how far (arcsec) the stored y/z are from what ra/dec give; NaN if no pointing."""
    y_angle, z_angle = _angles_from_radec(table, pointings)
    return np.hypot(
        np.asarray(table["y_angle"], dtype=float) - y_angle,
        np.asarray(table["z_angle"], dtype=float) - z_angle,
    )


def recompute_angles(table: Table, pointings: dict[int, Pointing]) -> Table:
    """A copy of `table` with y_angle/z_angle (and r_angle, if present) from ra/dec.

    Rows whose obsid has no pointing keep their stored values.
    """
    fixed = table.copy()
    y_angle, z_angle = _angles_from_radec(table, pointings)
    has_pointing = np.isfinite(y_angle)
    fixed["y_angle"][has_pointing] = y_angle[has_pointing]
    fixed["z_angle"][has_pointing] = z_angle[has_pointing]
    if "r_angle" in fixed.colnames:
        fixed["r_angle"][has_pointing] = np.hypot(
            y_angle[has_pointing], z_angle[has_pointing]
        )
    return fixed


def _key_values(table: Table, keys: Sequence[str]) -> list[tuple]:
    """Each row's key as a tuple of plain Python values."""
    columns = []
    for key in keys:
        values = np.asarray(table[key])
        columns.append(values.astype(str) if values.dtype.kind in "SU" else values)
    return [tuple(v.item() for v in row) for row in zip(*columns, strict=True)]


def stale_groups(
    table: Table,
    pointings: dict[int, Pointing],
    keys: Sequence[str],
    tolerance: float = DEFAULT_TOLERANCE,
) -> list[tuple]:
    """Sorted `keys` values of the groups with a row off by more than `tolerance` arcsec."""
    stale = angle_errors(table, pointings) > tolerance  # NaN compares False
    return sorted(set(_key_values(table[stale], keys)))


def repair(
    dbfile: Path,
    roots: Sequence[Path],
    obsids: Sequence[int] | None = None,
    tolerance: float = DEFAULT_TOLERANCE,
    dry_run: bool = False,
) -> dict:
    """Rewrite the stale groups of astromon_xray_src and astromon_cat_src in `dbfile`.

    Returns a summary: for each table, the groups rewritten; ``no_pointing``, the
    obsids with no local evt2; ``past_180``, (obsid, detect_method, id) of X-ray
    rows the repair moved from inside 180" to at or past it; ``touched_obsids``,
    every obsid rewritten in either table, i.e. the ones whose xcorr to rebuild.
    """
    tables = {name: db.get_table(name, dbfile) for name in TABLE_KEYS}
    if obsids is not None:
        tables = {
            name: table[np.isin(np.asarray(table["obsid"]), list(obsids))]
            for name, table in tables.items()
        }
    all_obsids = sorted(
        {int(o) for table in tables.values() for o in np.unique(table["obsid"])}
    )
    logger.info(f"reading the evt2 pointing of {len(all_obsids):,} obsid(s)")
    pointings = {}
    no_pointing = []
    for obsid in all_obsids:
        pointing = local_evt2_pointing(obsid, roots)
        if pointing is None:
            no_pointing.append(obsid)
        else:
            pointings[obsid] = pointing

    summary = {"no_pointing": no_pointing, "past_180": []}
    touched = set()
    for name, keys in TABLE_KEYS.items():
        table = tables[name]
        groups = stale_groups(table, pointings, keys, tolerance)
        summary[name] = groups
        if not groups:
            logger.info(f"{name}: nothing to repair")
            continue
        stale_keys = set(groups)
        in_group = np.array([key in stale_keys for key in _key_values(table, keys)])
        fixed = recompute_angles(table[in_group], pointings)
        if "r_angle" in fixed.colnames:
            crossed = (np.asarray(table[in_group]["r_angle"]) < 180) & (
                np.asarray(fixed["r_angle"]) >= 180
            )
            summary["past_180"] += _key_values(fixed[crossed], [*keys, "id"])
        touched.update(group[0] for group in groups)
        logger.info(
            f"{name}: {len(fixed):,} row(s) in {len(groups):,} group(s)"
            f" over {len({g[0] for g in groups}):,} obsid(s)"
        )
        if not dry_run:
            db.save(name, fixed, dbfile, expect_existing=True)
    summary["touched_obsids"] = sorted(touched)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--db", required=True, type=Path, dest="db_file")
    parser.add_argument(
        "--root",
        required=True,
        type=Path,
        action="append",
        dest="roots",
        help="directory laid out as obsNN/<obsid>/primary; repeat to search several",
    )
    parser.add_argument(
        "--obsid-file", type=Path, help="only these obsids (one per line)"
    )
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    parser.add_argument(
        "--touched-out",
        type=Path,
        help="write the rewritten obsids here, for rebuild_xcorr.py --obsids",
    )
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    if not args.db_file.exists():
        parser.error(f"db does not exist: {args.db_file}")

    obsids = (
        [int(x) for x in args.obsid_file.read_text().split()]
        if args.obsid_file
        else None
    )
    summary = repair(
        args.db_file, args.roots, obsids, tolerance=args.tolerance, dry_run=args.dry_run
    )
    logger.info(f"obsids with no local evt2: {summary['no_pointing']}")
    logger.info(f"X-ray rows moved to r_angle >= 180: {summary['past_180']}")
    logger.info(f"obsids rewritten: {len(summary['touched_obsids']):,}")
    if args.touched_out:
        args.touched_out.write_text(
            "".join(f"{o}\n" for o in summary["touched_obsids"])
        )
        logger.info(f"wrote {args.touched_out}; rebuild their xcorr next")
    if args.dry_run:
        logger.info("dry run -- nothing written")


if __name__ == "__main__":
    main()
