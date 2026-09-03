#!/usr/bin/env python
"""
Confirm that one obsid can be downloaded and processed via both data-source paths.

The two paths are also checked for agreement on where they place celldetect
sources. This is meant to run on a CXC Linux system, since arc5gl needs the
CXC network and both paths need CIAO. For each of "archive" (CIAO's
download_chandra_obsid against the public Chandra archive) and "arc5gl", this
downloads the obsid into its own temporary working directory, runs
Observation.process(), and reads back the celldetect source list. If both
paths succeed, celldetect source positions are compared pairwise (nearest
neighbor) and flagged if they differ by more than --tolerance arcsec.

Nothing downloaded is kept -- each path gets a fresh temporary directory that
is removed when the script exits.

Usage
-----
    python -m astromon.scripts.check_obsid_sources OBSID
    python -m astromon.scripts.check_obsid_sources OBSID --source archive
"""

import argparse
import sys
import tempfile
import traceback
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.table import Table
from ska_helpers.logging import basic_logger

from astromon.observation import Observation

logger = basic_logger(__name__, level="INFO")

ALL_SOURCES = ("archive", "arc5gl")

# Sources within this separation (arcsec) are considered the same source when
# comparing celldetect positions between the two paths.
DEFAULT_TOLERANCE_ARCSEC = 1.0


@dataclass
class PathResult:
    source: str
    ok: bool
    msg: str = ""
    sources: Table | None = None


def _run_path(
    obsid: str, source: str, workdir: Path, ciao_prefix: Path | None
) -> PathResult:
    """Download and process `obsid` via `source`; return a PathResult."""
    try:
        observation = Observation(
            obsid, workdir=workdir, source=source, ciao_prefix=ciao_prefix
        )
        rv = observation.process()
        if not rv:
            return PathResult(source, ok=False, msg=f"process() failed: {rv.msg}")
        sources = observation.get_sources(version="celldetect")
        if len(sources) == 0:
            return PathResult(source, ok=False, msg="no celldetect sources found")
        return PathResult(source, ok=True, sources=sources)
    except Exception as e:  # noqa: BLE001 -- report any failure, keep checking the other path
        logger.error(f"{source}: {e}")
        traceback.print_exc()
        return PathResult(source, ok=False, msg=str(e))


def _compare_positions(
    a: PathResult, b: PathResult, tolerance_arcsec: float
) -> list[str]:
    """Compare celldetect source positions between two successful PathResults.

    Each source in `a` is matched to its nearest neighbor in `b`. Returns a
    list of human-readable problems; empty if every source matches within
    `tolerance_arcsec` and both paths found the same number of sources.
    """
    problems = []
    if len(a.sources) != len(b.sources):
        problems.append(
            f"source count differs: {a.source}={len(a.sources)} vs "
            f"{b.source}={len(b.sources)}"
        )

    ra_a, dec_a = a.sources["ra"], a.sources["dec"]
    ra_b, dec_b = b.sources["ra"], b.sources["dec"]
    cos_dec = np.cos(np.radians(dec_a[:, None]))
    # small-angle approximation is fine at celldetect-scale separations
    sep = 3600 * np.sqrt(
        ((ra_a[:, None] - ra_b[None, :]) * cos_dec) ** 2
        + (dec_a[:, None] - dec_b[None, :]) ** 2
    )
    nearest = np.argmin(sep, axis=1)
    for i, j in enumerate(nearest):
        if sep[i, j] > tolerance_arcsec:
            problems.append(
                f"source {i}: {a.source} ({ra_a[i]:.5f}, {dec_a[i]:.5f}) has no "
                f"match in {b.source} within {tolerance_arcsec} arcsec "
                f"(nearest is {sep[i, j]:.2f} arcsec away)"
            )
    return problems


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("obsid", help="OBSID to check")
    parser.add_argument(
        "--source",
        action="append",
        choices=ALL_SOURCES,
        dest="sources",
        help="Data source to check. Repeat to check more than one. Default: both.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE_ARCSEC,
        help=(
            "Maximum separation, in arcsec, for two celldetect sources from "
            f"different paths to count as the same source. Default: "
            f"{DEFAULT_TOLERANCE_ARCSEC}."
        ),
    )
    parser.add_argument("--ciao-prefix", help="CIAO prefix", default=None, type=Path)
    parser.add_argument(
        "--log-level",
        choices=["debug", "info", "warning"],
        default="info",
    )
    return parser


def main():
    args = get_parser().parse_args()
    logger.setLevel(args.log_level.upper())
    sources = args.sources or list(ALL_SOURCES)

    results = {}
    for source in sources:
        with tempfile.TemporaryDirectory() as workdir:
            logger.info(f"obsid {args.obsid}: checking source={source!r}")
            results[source] = _run_path(
                args.obsid, source, Path(workdir), args.ciao_prefix
            )

    ok = True
    for source, result in results.items():
        if result.ok:
            logger.info(
                f"obsid {args.obsid}: source={source!r} OK "
                f"({len(result.sources)} celldetect sources)"
            )
        else:
            ok = False
            logger.error(f"obsid {args.obsid}: source={source!r} FAILED: {result.msg}")

    if ok and len(results) > 1:
        source_a, source_b = sources[:2]
        problems = _compare_positions(
            results[source_a], results[source_b], args.tolerance
        )
        if problems:
            ok = False
            logger.error(
                f"obsid {args.obsid}: celldetect positions disagree between "
                f"{source_a!r} and {source_b!r}:"
            )
            for problem in problems:
                logger.error(f"  {problem}")
        else:
            logger.info(
                f"obsid {args.obsid}: celldetect positions agree between "
                f"{source_a!r} and {source_b!r} within {args.tolerance} arcsec"
            )

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
