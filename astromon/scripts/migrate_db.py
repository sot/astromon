"""
Migrate an astromon HDF5 database to the current schema.

Adds any columns present in the target dtypes that are missing from the stored
tables and fills them with appropriate defaults:

astromon_xray_src:
    psfratio           -> nan        (not computed for legacy rows)
    concentration_ratio -> nan       (not computed for legacy rows)
    detect_method      -> "celldetect"  (the only method used before this change)

astromon_xcorr:
    detect_method      -> "celldetect"  (the only method used before this change)

Usage
-----
    python -m astromon.scripts.migrate_db astromon.h5 [--dry-run]
"""

import argparse
import shutil
from pathlib import Path

import numpy as np
from astropy.table import Table

from astromon import db


XRAY_SRC_DEFAULTS = {
    "psfratio": float("nan"),
    "concentration_ratio": float("nan"),
    "detect_method": "celldetect",
}

XCORR_DEFAULTS = {
    "detect_method": "celldetect",
}


def _migrate_table(
    h5,
    node_name: str,
    target_dtype: np.dtype,
    defaults: dict,
    dry_run: bool,
) -> bool:
    """Migrate one table node. Returns True if changes were needed."""
    raw = getattr(h5.root, node_name)[:]
    tbl = Table(raw)
    tbl.convert_bytestring_to_unicode()

    missing = [n for n in target_dtype.names if n not in tbl.colnames]
    if not missing:
        print(f"{node_name}: already up to date.")
        return False

    print(f"{node_name}: columns to add: {missing}")
    for col in missing:
        if col not in defaults:
            raise ValueError(
                f"No default defined for new column '{col}' in {node_name}. "
                "Add it to the defaults dict in this script."
            )
        default = defaults[col]
        fill = np.full(len(tbl), default, dtype=target_dtype[col])
        tbl[col] = fill
        print(f"  {col}: filled {len(tbl)} rows with {default!r}")

    if not dry_run:
        db.save(node_name, tbl, h5, ignore_obsid=True)

    return True


def migrate(dbfile: Path, dry_run: bool = False) -> None:
    with db.connect(dbfile) as h5:
        raw_xray = h5.root.astromon_xray_src[:]
        raw_xcorr = h5.root.astromon_xcorr[:]

    xray_src = Table(raw_xray)
    xray_src.convert_bytestring_to_unicode()
    xcorr = Table(raw_xcorr)
    xcorr.convert_bytestring_to_unicode()

    xray_missing = [n for n in db.ASTROMON_XRAY_SRC_DTYPE.names if n not in xray_src.colnames]
    xcorr_missing = [n for n in db.ASTROMON_XCORR_DTYPE.names if n not in xcorr.colnames]

    if not xray_missing and not xcorr_missing:
        print("Both tables are already up to date — nothing to do.")
        return

    if xray_missing:
        print(f"astromon_xray_src: columns to add: {xray_missing}")
        for col in xray_missing:
            if col not in XRAY_SRC_DEFAULTS:
                raise ValueError(
                    f"No default defined for new column '{col}' in astromon_xray_src. "
                    "Add it to XRAY_SRC_DEFAULTS in this script."
                )
            default = XRAY_SRC_DEFAULTS[col]
            fill = np.full(len(xray_src), default, dtype=db.ASTROMON_XRAY_SRC_DTYPE[col])
            xray_src[col] = fill
            print(f"  {col}: filled {len(xray_src)} rows with {default!r}")

    if xcorr_missing:
        print(f"astromon_xcorr: columns to add: {xcorr_missing}")
        for col in xcorr_missing:
            if col not in XCORR_DEFAULTS:
                raise ValueError(
                    f"No default defined for new column '{col}' in astromon_xcorr. "
                    "Add it to XCORR_DEFAULTS in this script."
                )
            default = XCORR_DEFAULTS[col]
            fill = np.full(len(xcorr), default, dtype=db.ASTROMON_XCORR_DTYPE[col])
            xcorr[col] = fill
            print(f"  {col}: filled {len(xcorr)} rows with {default!r}")

    if dry_run:
        print("Dry run — no changes written.")
        return

    backup = dbfile.with_suffix(".h5.pre_migrate")
    print(f"Backing up to {backup}")
    shutil.copy2(dbfile, backup)

    if xray_missing:
        db.save("astromon_xray_src", xray_src, dbfile, ignore_obsid=True)
        print("astromon_xray_src migration complete.")
    if xcorr_missing:
        db.save("astromon_xcorr", xcorr, dbfile, ignore_obsid=True)
        print("astromon_xcorr migration complete.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dbfile", type=Path, help="Path to astromon.h5")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report what would change without writing anything",
    )
    args = parser.parse_args()
    migrate(args.dbfile, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
