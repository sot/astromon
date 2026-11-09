#!/usr/bin/env python
"""Pre-build the local DESIV161 (DESI DR1 QSO/GALAXY) catalog.

Production cross-matching gets its DESIV161 data from a local cache built from
DESI's own DR1 data release (astromon.cross_match's "DESIV161" catalog), not
live VizieR -- see astromon.cross_match._build_desi_v161_catalog for why (CDS/
VizieR has been unavailable since 2026-10-28, no ETA; the same situation that
moved Tycho2 off VizieR, see astromon.tycho2_local_index). get_desi_v161_catalog()
already downloads and builds this file automatically on first use (the same
download-and-cache-on-first-use convention astromon.cross_match uses for RFC,
ICRF3, etc.), so running this script is never required. It exists to pay the
one-time cost (~21 GB transient download, several minutes to filter and
deduplicate) up front, in one place, rather than inside whichever of several
parallel per-obsid worker processes happens to need it first.

The ~21 GB raw download needs a temporary directory with that much free space;
set TMPDIR first if the default location does not have it.

Usage
-----
::

    python -m astromon.scripts.build_desi_v161_local_index
    python -m astromon.scripts.build_desi_v161_local_index --force
"""

import argparse

from ska_helpers.logging import basic_logger

from astromon.cross_match import _DESI_V161_CACHE_PATH, get_desi_v161_catalog


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Rebuild even if a cached copy with the current release already exists.",
    )
    parser.add_argument(
        "--log-level",
        choices=["debug", "info", "warning", "error", "fatal"],
        default="info",
    )
    return parser


def main() -> None:
    args = get_parser().parse_args()
    logger = basic_logger("astromon", level=args.log_level.upper())

    logger.info(f"DESIV161 catalog -> {_DESI_V161_CACHE_PATH}")
    catalog = get_desi_v161_catalog(force_rebuild=args.force)
    logger.info(f"{len(catalog)} DESIV161 sources")


if __name__ == "__main__":
    main()
