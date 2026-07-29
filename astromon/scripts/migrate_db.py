"""
Migrate an astromon HDF5 database to the current schema.

Adds any columns present in ASTROMON_XRAY_SRC_DTYPE that are missing from the
stored table and fills them with appropriate defaults:

    psfratio           -> nan   (not computed for legacy rows)
    concentration_ratio -> nan  (not computed for legacy rows)
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


COLUMN_DEFAULTS = {
    "psfratio": float("nan"),
    "concentration_ratio": float("nan"),
    "detect_method": "celldetect",
}


def migrate(dbfile: Path, dry_run: bool = False) -> None:
    # Read the raw stored array without enforcing the current dtype, so we can
    # handle files written before the new columns existed.
    with db.connect(dbfile) as h5:
        raw = h5.root.astromon_xray_src[:]
    xray_src = Table(raw)
    xray_src.convert_bytestring_to_unicode()

    target_names = db.ASTROMON_XRAY_SRC_DTYPE.names
    missing = [n for n in target_names if n not in xray_src.colnames]

    if not missing:
        print("astromon_xray_src is already up to date — nothing to do.")
        return

    print(f"Columns to add: {missing}")

    for col in missing:
        if col not in COLUMN_DEFAULTS:
            raise ValueError(
                f"No default defined for new column '{col}'. "
                "Add it to COLUMN_DEFAULTS in this script."
            )
        default = COLUMN_DEFAULTS[col]
        dtype = db.ASTROMON_XRAY_SRC_DTYPE[col]
        fill = np.full(len(xray_src), default, dtype=dtype)
        xray_src[col] = fill
        print(f"  {col}: filled {len(xray_src)} rows with {default!r}")

    if dry_run:
        print("Dry run — no changes written.")
        return

    backup = dbfile.with_suffix(".h5.pre_migrate")
    print(f"Backing up to {backup}")
    shutil.copy2(dbfile, backup)

    db.save("astromon_xray_src", xray_src, dbfile, ignore_obsid=True)
    print("Migration complete.")


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
