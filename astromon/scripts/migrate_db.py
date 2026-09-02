"""
Migrate an astromon HDF5 database to the current schema.

Adds any columns present in the target dtypes that are missing from the stored
tables, filling them with appropriate defaults:

astromon_xray_src:
    psfratio           -> nan        (not computed for legacy rows)
    concentration_ratio -> nan       (not computed for legacy rows)
    peak_offset        -> nan       (not computed for legacy rows)
    detect_method      -> "celldetect"  (the only method used before this change)

astromon_xcorr:
    detect_method      -> "celldetect"  (the only method used before this change)

Drops columns that have been retired from the schema (none currently; the
production and dev databases have neither had one).

A column found in the database but neither in the target dtype nor in the
retired list is an error, so a schema mistake cannot silently delete data.

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

DEFAULTS = {
    "astromon_xray_src": {
        "psfratio": float("nan"),
        "concentration_ratio": float("nan"),
        "peak_offset": float("nan"),
        "detect_method": "celldetect",
    },
    "astromon_xcorr": {
        "detect_method": "celldetect",
    },
}

# Columns that used to be in the schema and may still exist in older database
# files. These are dropped on migration. Currently empty; nothing has been
# retired that isn't already gone from the production and dev databases.
RETIRED_COLUMNS = {
    "astromon_xray_src": set(),
    "astromon_xcorr": set(),
}


def _prepare_table(table_name: str, raw: np.ndarray) -> tuple[Table, bool]:
    """Return (migrated table, changes_needed) for one stored table.

    Adds missing target columns with defaults and reports retired columns to
    drop. The actual drop happens in db.save(), which casts to the target
    dtype and ignores columns not in it.
    """
    target_dtype = db.DTYPES[table_name]
    defaults = DEFAULTS[table_name]
    retired = RETIRED_COLUMNS[table_name]

    tbl = Table(raw)
    tbl.convert_bytestring_to_unicode()

    missing = [n for n in target_dtype.names if n not in tbl.colnames]
    extra = [n for n in tbl.colnames if n not in target_dtype.names]

    unexpected = [n for n in extra if n not in retired]
    if unexpected:
        raise ValueError(
            f"{table_name}: columns {unexpected} are in the database but neither in "
            "the target dtype nor in RETIRED_COLUMNS. Refusing to drop them; add "
            "them to RETIRED_COLUMNS in this script if the removal is intended."
        )

    if missing:
        print(f"{table_name}: columns to add: {missing}")
        for col in missing:
            if col not in defaults:
                raise ValueError(
                    f"No default defined for new column '{col}' in {table_name}. "
                    "Add it to DEFAULTS in this script."
                )
            default = defaults[col]
            tbl[col] = np.full(len(tbl), default, dtype=target_dtype[col])
            print(f"  {col}: filled {len(tbl)} rows with {default!r}")

    if extra:
        print(f"{table_name}: retired columns to drop: {extra}")

    return tbl, bool(missing or extra)


def migrate(dbfile: Path, dry_run: bool = False) -> None:
    table_names = ["astromon_xray_src", "astromon_xcorr"]
    with db.connect(dbfile) as h5:
        raw_tables = {name: getattr(h5.root, name)[:] for name in table_names}

    migrated = {}
    for name in table_names:
        tbl, changed = _prepare_table(name, raw_tables[name])
        if changed:
            migrated[name] = tbl
        else:
            print(f"{name}: already up to date.")

    if not migrated:
        print("Both tables are already up to date — nothing to do.")
        return

    if dry_run:
        print("Dry run — no changes written.")
        return

    backup = dbfile.with_suffix(".h5.pre_migrate")
    print(f"Backing up to {backup}")
    shutil.copy2(dbfile, backup)

    for name, tbl in migrated.items():
        db.save(name, tbl, dbfile, ignore_obsid=True)
        print(f"{name} migration complete.")


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
