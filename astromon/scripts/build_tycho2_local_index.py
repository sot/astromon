#!/usr/bin/env python
"""Pre-build the local Tycho-2 + Hipparcos validation index.

astromon.tycho2_local_index.get_tycho2_local_index() / get_hipparcos_snapshot()
already download and build these files automatically on first use (the same
download-and-cache-on-first-use convention astromon.cross_match uses for Milliquas,
Quaia, ICRF3, etc.), so running this script is never required. It exists to pay the
one-time cost (the Tycho-2 index is a few hundred MB and takes a couple of minutes to
parse) up front, in one place, rather than inside whichever of several parallel
per-obsid worker processes happens to need it first.

This module is purely a validation tool for a potential VizieR cutover (see
astromon.tycho2_local_index); nothing here is called by production cross-matching.

Usage
-----
::

    python -m astromon.scripts.build_tycho2_local_index
    python -m astromon.scripts.build_tycho2_local_index --force
"""

import argparse

from ska_helpers.logging import basic_logger

from astromon.tycho2_local_index import (
    HIPPARCOS_INDEX_PATH,
    TYCHO2_INDEX_PATH,
    get_hipparcos_snapshot,
    get_tycho2_local_index,
)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Rebuild both files even if a cached copy with the current release "
        "already exists.",
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

    logger.info(f"Tycho-2 index -> {TYCHO2_INDEX_PATH}")
    get_tycho2_local_index(force_rebuild=args.force)

    logger.info(f"Hipparcos snapshot -> {HIPPARCOS_INDEX_PATH}")
    get_hipparcos_snapshot(force_rebuild=args.force)


if __name__ == "__main__":
    main()
