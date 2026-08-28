#!/usr/bin/env python
"""
Script to find x-ray sources in observations, and a tentative set of optical/radio counterparts.
"""

import argparse
import logging
import os
import re
import shutil
import sys
import tempfile
import traceback
from multiprocessing import Pool, set_start_method
from pathlib import Path

import numpy as np
import stk
from astropy import coordinates as coords
from astropy.table import Table, vstack
from cxotime import CxoTime
from cxotime import units as u
from Quaternion import Quat
from Ska.arc5gl import Arc5gl
from ska_helpers.logging import basic_logger

from astromon import db, utils
from astromon.cross_match import (
    assign_cat_src_ids,
    compute_cross_matches,
    get_desi_v161_candidates,
    get_gaia_agn,
    get_gaia_qso_candidates,
    get_gaia_var_stars,
    get_milliquas_gaia,
    get_quaia_candidates,
    rough_match,
)
from astromon.observation import Observation, ReturnCode
from astromon.task import run_tasks


class Skipped(Exception):
    """
    Exception class used to abort and silently skip processing an observation.
    """


class SkippedWithWarning(Exception):
    """
    Exception class used to abort, issue a warning, and skip processing an observation.
    """


def get_obsids(tstart, tstop):
    tstart = CxoTime(tstart)
    tstop = CxoTime(tstop)
    with tempfile.TemporaryDirectory() as td, utils.chdir(td):
        path = Path(td)
        arc5gl = Arc5gl()
        arc5gl.sendline(f"tstart={tstart.date}")
        arc5gl.sendline(f"tstop={tstop.date}")
        arc5gl.sendline("get obspar")
        names = [str(p) for p in path.glob("*")]
        return [re.search("axaff([0-9]+)_", name).group(1) for name in names]


def save(data, db_file):  # noqa: PLR0912
    logger = logging.getLogger(name="astromon")

    errors = {}
    skipped = {}
    n_skip = 0
    for d in [d for d in data if d["msg"]]:
        if d["ok"]:
            obs = skipped.get(d["msg"], [])
            skipped[d["msg"]] = obs + [str(d["obsid"])]
        else:
            obs = errors.get(d["msg"], [])
            errors[d["msg"]] = obs + [str(d["obsid"])]
        n_skip += 1
    n = len(data)
    if n_skip:
        logger.warning(f"WARNING: {n} observations , {n_skip} skipped.")
    else:
        logger.info(f"{n} observations , {n_skip} skipped.")
    for key, val in skipped.items():
        logger.warning(f"WARNING: {len(val)} {key} ({', '.join(val)})")
    for key, val in errors.items():
        logger.warning(f"ERROR: {len(val)} {key} ({', '.join(val)})")

    failed_obsids = [str(d["obsid"]) for d in data if not d["ok"]]
    failures_file = Path(str(db_file)).with_suffix(".failures")
    if failed_obsids:
        failures_file.write_text("\n".join(failed_obsids) + "\n")
        logger.warning(
            f"{len(failed_obsids)} failed obsids written to {failures_file} "
            f"-- retry with: astromon-cross-match --obsid {failures_file}"
        )
    elif failures_file.exists():
        failures_file.unlink()

    data = [d for d in data if "astromon_obs" in d and len(d["astromon_obs"]) > 0]
    if len(data) == 0:
        logger.info("Nothing to save")
        return

    names = ["astromon_obs", "astromon_xray_src", "astromon_cat_src", "astromon_xcorr"]
    # these following baroque lines are here because there are some columns we cannot vstack
    # so I decided to only vstack the columns that are in the dtype,
    # but it turns out that not all columns in the dtype are actually in the data...
    # and the data itself could be an empty array...
    data = {name: [d[name] for d in data if len(d[name])] for name in names}
    for name in names:
        if len(data[name]):
            cols = [
                col for col in db.DTYPES[name].names if col in data[name][0].colnames
            ]
            data[name] = [d[cols] for d in data[name]]
            data[name] = vstack(data[name], metadata_conflicts="silent")

    logger.debug(
        f"About to write {len(data['astromon_obs'])} observations to {db_file}"
    )
    with db.connect(db_file, mode="r+") as con:
        for name in names:
            if len(data[name]):
                db.save(name, data[name], con)


# Detection versions that are run for a side effect rather than saved as their own
# detect_method. peak_gaussian_detect writes the .src file that get_sources reads to
# compute peak_offset, but its fitted positions are not persisted: they never won on
# match count, and what they carry is the offset itself. The .src stays on disk, so
# the positions remain recoverable.
DIAGNOSTIC_VERSIONS = ("peak_gaussian_detect",)


def _split_versions(versions):
    """Split requested versions into ``(saved, diagnostic)``, preserving order.

    Diagnostic versions must run before any ``get_sources`` call, because that is
    what reads their output to compute ``peak_offset`` -- and its result is cached
    by ``@stored_result``, so a version run afterwards would leave the column NaN
    with no second chance.
    """
    saved = [v for v in versions if v not in DIAGNOSTIC_VERSIONS]
    diagnostic = [v for v in versions if v in DIAGNOSTIC_VERSIONS]
    return saved, diagnostic


def process_obsid(  # noqa: PLR0915, PLR0912, PLR0917
    obsid,
    workdir,
    archive_dir=None,
    versions=("celldetect",),
    source: str = "arc5gl",
    ciao_prefix: str | None = None,
    cleanup: bool = False,
):
    """
    Core processing logic for a single obsid. Raises on failure.

    Returns a dict with keys: obsid, astromon_obs, astromon_xray_src,
    astromon_cat_src, astromon_xcorr.  Raises Skipped or SkippedWithWarning
    for expected non-fatal cases, and Exception for real failures.

    This is the single source of truth for per-obsid processing. Both the
    astromon-cross-match CLI (via process()) and process_one_obsid.py call
    this function; neither duplicates the logic.

    Parameters
    ----------
    obsid : int
        Chandra observation ID.
    workdir : str or Path
        Top-level working directory.
    archive_dir : str or Path, optional
        Pre-populated archive directory; files are symlinked from here
        rather than downloaded.
    versions : tuple of str
        Detection versions to run (e.g. ``("celldetect", "gaussian_detect")``).
    source : str
        Data download method: ``"arc5gl"`` (default, CXC-internal) or
        ``"archive"`` (public Chandra archive via ``download_chandra_obsid``).
    ciao_prefix : str or None
        Path to the CIAO installation (e.g. ``"/soft/ciao"``).  ``None``
        falls back to the ``Ciao`` class default (``/soft/ciao``).
    cleanup : bool
        If True, delete downloaded files that are not needed after source
        detection (e.g. level-1 events, VV report, preview JPEGs).  See
        ``Observation.cleanup_downloads`` for the full list.  Default False.
    """
    logger = logging.getLogger("astromon")

    exceptions = {
        ReturnCode.SKIP: Skipped,
        ReturnCode.WARNING: SkippedWithWarning,
        ReturnCode.ERROR: Exception,
    }

    logger.info(f"OBSID={obsid} *** Processing OBSID {obsid} ***")
    observation = Observation(
        obsid,
        workdir,
        archive_dir=archive_dir,
        source=source,
        ciao_prefix=ciao_prefix,
    )

    if archive_dir:
        observation.create_archive_symlinks()

    if not (rv := observation.process()):
        raise exceptions.get(rv.return_code, Exception)(rv.msg)

    # Archive celldetect intermediate files before running additional versions.
    if archive_dir:
        observation.archive()

    saved_versions, diagnostic_versions = _split_versions(versions)

    # Before any get_sources call: it reads the diagnostic versions' .src output to
    # compute peak_offset, and caches its result.
    for version in diagnostic_versions:
        try:
            rv = run_tasks(
                observation, requested_files=[f"sources/{obsid}_{version}.src"]
            )
        except Exception as exc:
            logger.warning(f"OBSID={obsid} {version} raised exception, skipping: {exc}")
            continue
        failed = [
            v for v in rv.values() if v.return_code.value >= ReturnCode.ERROR.value
        ]
        if failed:
            logger.warning(f"OBSID={obsid} {version} failed: {failed[0].msg}")

    obspar = Table([observation.get_info()])

    # Rough-match catalog sources once from celldetect positions, which are always
    # available after observation.process(). The positions are nearly identical across
    # detection methods so this candidate set is valid for all of them.
    celldetect_sources = observation.get_sources(version="celldetect")
    if len(celldetect_sources) == 0:
        raise SkippedWithWarning("No x-ray sources found")

    q = Quat(
        equatorial=(
            obspar["ra_pnt"][0],
            obspar["dec_pnt"][0],
            obspar["roll_pnt"][0],
        )
    )
    obs_time = CxoTime(obspar["date_obs"][0])
    # Include RFC and ICRF3 (both local-cache, no network) so astromon_23 can match
    # against radio sources and icrf3 can label ICRF3-specific matches separately.
    # Tycho2/USNO-B1.0/2MASS/SDSS are VizieR queries.
    match_candidates = rough_match(
        celldetect_sources,
        obs_time,
        catalogs=("RFC", "ICRF3", "Tycho2", "USNO-B1.0", "2MASS", "SDSS"),
    )
    # rough_match already sets one row per (source, candidate) pair, so
    # assign_cat_src_ids only stamps the remaining bookkeeping fields.
    assign_cat_src_ids(match_candidates, obsid, celldetect_sources, q)

    # Build a time-stamped SkyCoord for all Gaia catalog queries.
    # obstime is required by get_gaia_var_stars for proper motion correction.
    xray_pos = coords.SkyCoord(
        celldetect_sources["ra"],
        celldetect_sources["dec"],
        unit="deg",
        obstime=obs_time,
    )

    def _add_obsid_and_anchor(candidates, id_offset):
        """Stamp obsid, sequential ids, y/z angles and the celldetect anchor.

        The anchor is the nearest celldetect source: provenance for why this
        catalogue row was fetched, and what `separation` measures from. It is not a
        match key -- the pairing is decided per method by position and the dr cut in
        simple_cross_match.
        """
        assign_cat_src_ids(candidates, obsid, celldetect_sources, q, id_offset)

    gaia_agn_candidates = get_gaia_agn(
        xray_pos, radius=3 * u.arcsec, logging_tag=f"OBSID={obsid}"
    )
    if len(gaia_agn_candidates):
        _add_obsid_and_anchor(gaia_agn_candidates, id_offset=len(match_candidates))
        logger.info(f"OBSID={obsid} GaiaAGN: {len(gaia_agn_candidates)} candidate(s)")
    else:
        logger.debug(f"OBSID={obsid} GaiaAGN: no candidates")

    # Gaia DR3 qso_candidates -- extends agn_cross_id to ~6.6M sources (G up to ~21.5).
    gaia_qso_candidates = get_gaia_qso_candidates(
        xray_pos, radius=3 * u.arcsec, logging_tag=f"OBSID={obsid}"
    )
    if len(gaia_qso_candidates):
        id_offset = len(match_candidates) + len(gaia_agn_candidates)
        _add_obsid_and_anchor(gaia_qso_candidates, id_offset=id_offset)
        logger.info(f"OBSID={obsid} GaiaQSO: {len(gaia_qso_candidates)} candidate(s)")
    else:
        logger.debug(f"OBSID={obsid} GaiaQSO: no candidates")

    # Milliquas v8 + Gaia DR3 candidates -- VizieR cone search + Gaia position upgrade.
    milliquas_candidates = get_milliquas_gaia(
        xray_pos, radius=3 * u.arcsec, logging_tag=f"OBSID={obsid}"
    )
    if len(milliquas_candidates):
        id_offset = (
            len(match_candidates) + len(gaia_agn_candidates) + len(gaia_qso_candidates)
        )
        _add_obsid_and_anchor(milliquas_candidates, id_offset=id_offset)
        logger.info(
            f"OBSID={obsid} MilliquasGaia: {len(milliquas_candidates)} candidate(s)"
        )
    else:
        logger.debug(f"OBSID={obsid} MilliquasGaia: no candidates")

    # Quaia (Gaia DR3 + unWISE) candidates -- local cache, bounding-cone filter.
    quaia_candidates = get_quaia_candidates(
        xray_pos, radius=3 * u.arcsec, logging_tag=f"OBSID={obsid}"
    )
    if len(quaia_candidates):
        id_offset = (
            len(match_candidates)
            + len(gaia_agn_candidates)
            + len(gaia_qso_candidates)
            + len(milliquas_candidates)
        )
        _add_obsid_and_anchor(quaia_candidates, id_offset=id_offset)
        logger.info(f"OBSID={obsid} Quaia: {len(quaia_candidates)} candidate(s)")
    else:
        logger.debug(f"OBSID={obsid} Quaia: no candidates")

    # DESI EDR (V/161) candidates -- VizieR cone search, QSO+GALAXY, ZWARN==0.
    desi_candidates = get_desi_v161_candidates(
        xray_pos, radius=3 * u.arcsec, logging_tag=f"OBSID={obsid}"
    )
    if len(desi_candidates):
        id_offset = (
            len(match_candidates)
            + len(gaia_agn_candidates)
            + len(gaia_qso_candidates)
            + len(milliquas_candidates)
            + len(quaia_candidates)
        )
        _add_obsid_and_anchor(desi_candidates, id_offset=id_offset)
        logger.info(f"OBSID={obsid} DESIV161: {len(desi_candidates)} candidate(s)")
    else:
        logger.debug(f"OBSID={obsid} DESIV161: no candidates")

    # GaiaVarStar catalog candidates -- proper motion corrected to obs epoch.
    var_star_candidates = get_gaia_var_stars(
        xray_pos, radius=3 * u.arcsec, logging_tag=f"OBSID={obsid}"
    )
    if len(var_star_candidates):
        id_offset = (
            len(match_candidates)
            + len(gaia_agn_candidates)
            + len(gaia_qso_candidates)
            + len(milliquas_candidates)
            + len(quaia_candidates)
            + len(desi_candidates)
        )
        _add_obsid_and_anchor(var_star_candidates, id_offset=id_offset)
        logger.info(
            f"OBSID={obsid} GaiaVarStar: {len(var_star_candidates)} candidate(s)"
        )
    else:
        logger.debug(f"OBSID={obsid} GaiaVarStar: no candidates")

    xcorr_cols = list(db.ASTROMON_XCORR_DTYPE.names)
    all_sources = []
    all_xcorr = []

    for version in saved_versions:
        if version != "celldetect":
            try:
                rv = run_tasks(
                    observation,
                    requested_files=[f"sources/{obsid}_{version}.src"],
                )
            except Exception as exc:
                logger.warning(
                    f"OBSID={obsid} {version} raised exception, skipping version: {exc}"
                )
                continue
            failed = [
                v for v in rv.values() if v.return_code.value >= ReturnCode.ERROR.value
            ]
            if failed:
                logger.warning(f"OBSID={obsid} {version} failed: {failed[0].msg}")
                continue

        sources = observation.get_sources(version=version)
        if len(sources) == 0:
            logger.warning(f"OBSID={obsid} no sources for {version}")
            continue
        all_sources.append(sources)

        logger.debug(f"OBSID={obsid} About to cross-match ({version})")
        matches = compute_cross_matches(
            "astromon_21",
            astromon_obs=obspar,
            astromon_xray_src=sources,
            astromon_cat_src=match_candidates,
            logging_tag=f"OBSID={obsid}",
        )
        if len(matches):
            all_xcorr.append(matches[xcorr_cols])

        icrf3_matches = compute_cross_matches(
            "icrf3",
            astromon_obs=obspar,
            astromon_xray_src=sources,
            astromon_cat_src=match_candidates,
            logging_tag=f"OBSID={obsid}",
        )
        if len(icrf3_matches):
            all_xcorr.append(icrf3_matches[xcorr_cols])

        rfc_matches = compute_cross_matches(
            "rfc",
            astromon_obs=obspar,
            astromon_xray_src=sources,
            astromon_cat_src=match_candidates,
            logging_tag=f"OBSID={obsid}",
        )
        if len(rfc_matches):
            all_xcorr.append(rfc_matches[xcorr_cols])

        tycho2_matches = compute_cross_matches(
            "tycho2",
            astromon_obs=obspar,
            astromon_xray_src=sources,
            astromon_cat_src=match_candidates,
            logging_tag=f"OBSID={obsid}",
        )
        if len(tycho2_matches):
            all_xcorr.append(tycho2_matches[xcorr_cols])

        if len(gaia_agn_candidates):
            gaia_matches = compute_cross_matches(
                "gaia_agn",
                astromon_obs=obspar,
                astromon_xray_src=sources,
                astromon_cat_src=gaia_agn_candidates,
                logging_tag=f"OBSID={obsid}",
            )
            if len(gaia_matches):
                all_xcorr.append(gaia_matches[xcorr_cols])

        if len(gaia_qso_candidates):
            gaia_qso_matches = compute_cross_matches(
                "gaia_qso",
                astromon_obs=obspar,
                astromon_xray_src=sources,
                astromon_cat_src=gaia_qso_candidates,
                logging_tag=f"OBSID={obsid}",
            )
            if len(gaia_qso_matches):
                all_xcorr.append(gaia_qso_matches[xcorr_cols])

        if len(milliquas_candidates):
            mq_matches = compute_cross_matches(
                "milliquas_gaia",
                astromon_obs=obspar,
                astromon_xray_src=sources,
                astromon_cat_src=milliquas_candidates,
                logging_tag=f"OBSID={obsid}",
            )
            if len(mq_matches):
                all_xcorr.append(mq_matches[xcorr_cols])

        if len(quaia_candidates):
            quaia_matches = compute_cross_matches(
                "quaia",
                astromon_obs=obspar,
                astromon_xray_src=sources,
                astromon_cat_src=quaia_candidates,
                logging_tag=f"OBSID={obsid}",
            )
            if len(quaia_matches):
                all_xcorr.append(quaia_matches[xcorr_cols])

        if len(desi_candidates):
            desi_matches = compute_cross_matches(
                "desi_v161",
                astromon_obs=obspar,
                astromon_xray_src=sources,
                astromon_cat_src=desi_candidates,
                logging_tag=f"OBSID={obsid}",
            )
            if len(desi_matches):
                all_xcorr.append(desi_matches[xcorr_cols])

        if len(var_star_candidates):
            vs_matches = compute_cross_matches(
                "gaia_var_star",
                astromon_obs=obspar,
                astromon_xray_src=sources,
                astromon_cat_src=var_star_candidates,
                logging_tag=f"OBSID={obsid}",
            )
            if len(vs_matches):
                all_xcorr.append(vs_matches[xcorr_cols])

        # Combined hierarchical matches using all available catalogs.
        # astromon_23: RFC > GaiaAGN > GaiaQSO > Quaia > MilliquasGaia > DESIV161
        #              > GaiaVarStar > Tycho2.
        # astromon_24: same + ICRF3 prepended, celldetect only.
        # astromon_25: same + ICRF3 prepended, gaussian_detect only.
        all_candidates = [
            t
            for t in [
                match_candidates,
                gaia_agn_candidates,
                gaia_qso_candidates,
                quaia_candidates,
                milliquas_candidates,
                desi_candidates,
                var_star_candidates,
            ]
            if len(t)
        ]
        if all_candidates:
            combined_for_version = vstack(all_candidates, metadata_conflicts="silent")
            for select_name in (
                "astromon_23",
                "astromon_24",
                "astromon_25",
                "astromon_26",
                "astromon_27",
            ):
                combined_matches = compute_cross_matches(
                    select_name,
                    astromon_obs=obspar,
                    astromon_xray_src=sources,
                    astromon_cat_src=combined_for_version,
                    logging_tag=f"OBSID={obsid}",
                )
                if len(combined_matches):
                    all_xcorr.append(combined_matches[xcorr_cols])

    if not all_sources:
        raise SkippedWithWarning("No x-ray sources found for any detection version")

    sources = vstack(all_sources)
    matches = vstack(all_xcorr) if all_xcorr else Table(names=xcorr_cols)

    # Archive any new detection-version output files produced above.
    if archive_dir and len(versions) > 1:
        observation.archive()

    _cat_src_parts = [
        t
        for t in [
            match_candidates,
            gaia_agn_candidates,
            gaia_qso_candidates,
            quaia_candidates,
            milliquas_candidates,
            desi_candidates,
            var_star_candidates,
        ]
        if len(t)
    ]
    # match_candidates is always defined with the right schema (even when empty),
    # so it serves as the zero-row fallback when no catalogs have candidates.
    _all_cat_src = (
        vstack(_cat_src_parts, metadata_conflicts="silent")
        if _cat_src_parts
        else match_candidates
    )

    result = {
        "obsid": obsid,
        "astromon_obs": obspar,
        "astromon_xray_src": sources,
        "astromon_cat_src": _all_cat_src,
        "astromon_xcorr": matches,
    }

    (observation.workdir / "results").mkdir(exist_ok=True)
    # unicode characters in observer names and other fields require this
    result["astromon_obs"].convert_unicode_to_bytestring()
    for name in [
        "astromon_obs",
        "astromon_xray_src",
        "astromon_cat_src",
        "astromon_xcorr",
    ]:
        fn = observation.workdir / "results" / f"{name}.fits"
        fn.unlink(missing_ok=True)
        if len(result[name]) > 0:
            result[name].write(fn)

    if archive_dir:
        observation.archive(
            "astromon_obs.fits",
            "astromon_xray_src.fits",
            "astromon_cat_src.fits",
            "astromon_xcorr.fits",
        )

    if cleanup:
        observation.cleanup_downloads()

    return result


def process(  # noqa: PLR0917
    obsid,
    workdir,
    log_level,
    archive_dir,
    versions=("celldetect",),
    source: str = "arc5gl",
    ciao_prefix: str | None = None,
):
    """
    Catch exceptions from process_obsid() and return a status dict for the Pool.

    Suitable for use with the multiprocessing Pool in the astromon-cross-match CLI.
    """
    logger = basic_logger(
        "astromon", level=log_level, format="%(asctime)s %(funcName)-25s: %(message)s"
    )
    try:
        result = process_obsid(
            obsid,
            workdir,
            archive_dir=archive_dir,
            versions=versions,
            source=source,
            ciao_prefix=ciao_prefix,
        )
        return {"ok": True, "msg": "", **result}
    except Skipped as e:
        logger.info(f"OBSID={obsid} skipped")
        return {
            "ok": True,
            "msg": f"skipped: {e}",
            "obsid": obsid,
            "astromon_obs": [],
            "astromon_xray_src": [],
            "astromon_cat_src": [],
            "astromon_xcorr": [],
        }
    except SkippedWithWarning as e:
        logger.warning(f"OBSID={obsid} WARNING - skipped: {e}")
        return {
            "ok": True,
            "msg": f"skipped: {e}",
            "obsid": obsid,
            "astromon_obs": [],
            "astromon_xray_src": [],
            "astromon_cat_src": [],
            "astromon_xcorr": [],
        }
    except Exception as e:
        logger.error(f"OBSID={obsid} failed: {e}")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        trace = traceback.extract_tb(exc_traceback)
        for step in trace:
            logger.error(
                f"OBSID={obsid}            in {step.filename}:{step.lineno}/{step.name}:"
            )
            logger.error(f"OBSID={obsid}            {step.line}")
        return {
            "ok": False,
            "msg": f"error: {e}",
            "obsid": obsid,
            "astromon_obs": [],
            "astromon_xray_src": [],
            "astromon_cat_src": [],
            "astromon_xcorr": [],
        }


def get_parser():
    """
    Get the argument parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db-file",
        help="SQLite file where data is saved",
        default="astromon.h5",
        type=Path,
    )
    parser.add_argument(
        "--workdir",
        help="Working directory. A temp directory is created for each observation.",
        default=None,
    )
    parser.add_argument("--archive-dir", help="Archive directory.", default=None)
    parser.add_argument(
        "--start",
        help="Start of the time range to process. Default: stop - 30 days.",
    )
    parser.add_argument(
        "--stop",
        help="End of the time range to process. Default: now.",
    )
    parser.add_argument(
        "--obsid",
        help="An OBSID to process or a file with a list of OBSIDs",
    )
    parser.add_argument(
        "--no-exception",
        help="Do not skip observations already in file",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--threads",
        help="The number of processes to use.",
        # '9' because the archive limits number of concurrent connections
        # from same IP address (well, it did when it was 'ftp').
        default=9,
        type=int,
    )
    parser.add_argument(
        "--log-level",
        help="Logging severity level",
        choices=["debug", "info", "warning", "error", "fatal"],
        default="debug",
    )
    parser.add_argument(
        "--versions",
        nargs="+",
        default=["celldetect"],
        help=(
            "Detection versions to run and store. Default: celldetect. "
            "Adding gaussian_detect or peak_gaussian_detect runs those methods "
            "and saves additional rows keyed by detect_method."
        ),
    )
    parser.add_argument(
        "--source",
        choices=["arc5gl", "archive"],
        default="arc5gl",
        help=(
            "Data download method. 'arc5gl' (default) uses the CXC-internal "
            "arc5gl tool. 'archive' uses CIAO's download_chandra_obsid to pull "
            "from the public Chandra Data Archive — suitable for running outside "
            "the CXC network."
        ),
    )
    parser.add_argument(
        "--ciao-prefix",
        default=None,
        dest="ciao_prefix",
        help=(
            "Path to the CIAO installation (e.g. /soft/ciao). "
            "Required when CIAO is not at the default /soft/ciao location."
        ),
    )
    return parser


def main():  # noqa: PLR0912, PLR0915
    """
    Main routine that deals with processing arguments, output and such.
    """
    set_start_method("spawn")

    parser = get_parser()
    args = parser.parse_args()

    logger = basic_logger(
        "astromon",
        level=args.log_level.upper(),
        format="%(asctime)s %(funcName)-25s: %(message)s",
    )

    if args.workdir:
        workdir = Path(args.workdir)
        workdir.mkdir(exist_ok=True, parents=True)
    else:
        # if args.workdir is not given, create one for the result file
        # do not pass this downstream, because it can fill the tmp directory
        tmpdir = (
            tempfile.TemporaryDirectory()
        )  # this variable should live until the end
        workdir = Path(tmpdir.name)

    # all changes go in this file which is copied back to its final destination at the end.
    db_file = workdir / args.db_file.name
    db_file.parent.mkdir(exist_ok=True, parents=True)
    if args.db_file.exists():
        shutil.copy(args.db_file, db_file)
    else:
        logger.info(f"File does not exist: {args.db_file}. Will create")
        for name in db.DTYPES:
            tab = db.create_table(name)
            db.save(name, tab, dbfile=db_file)

    if args.obsid is not None:
        if os.path.exists(args.obsid) and os.path.isfile(args.obsid):
            # '@-' builds stack but does not include path name
            obsids = stk.build("@-" + args.obsid)
        else:
            obsids = stk.build(args.obsid)
    else:
        args.stop = CxoTime(args.stop) if args.stop is not None else CxoTime()
        args.start = (
            CxoTime(args.start) if args.start is not None else args.stop - 30 * u.day
        )
        obsids = get_obsids(args.start, args.stop)

    if len(obsids) == 0:
        logger.info("No OBSIDs found as specified")
        return

    if db_file.exists() and not args.no_exception:
        # "no_exception" means all OBSIDs are processes,
        # not no_exception means we need to look for existing OBSIDs and skip them
        obsids = np.array(obsids, dtype=int)
        try:
            astromon_obs = db.get_table("astromon_obs", dbfile=db_file)
            exceptions = np.isin(obsids, astromon_obs["obsid"])
            n_exceptions = np.sum(exceptions)
            if n_exceptions:
                exceptions_str = str(obsids[exceptions])[1:-1]
                logger.info(
                    f"skipping {n_exceptions} OBSID{'s' if n_exceptions > 1 else ''} "
                    f"that {'are' if n_exceptions > 1 else 'is'} on file already: {exceptions_str}"
                )
            obsids = obsids[~exceptions].astype(str)
        except utils.MissingTableException:
            pass

    if len(obsids) == 0:
        logger.info("All OBSIDs already processed")
        return

    logger.info(f"will process the following obsids: {', '.join(obsids)}")
    versions = tuple(args.versions)
    task_args = [
        (
            int(obsid),
            args.workdir,
            args.log_level.upper(),
            args.archive_dir,
            versions,
            args.source,
            args.ciao_prefix,
        )
        for obsid in obsids
    ]

    with Pool(processes=args.threads) as pool:
        results = pool.starmap(process, task_args)

    save(results, db_file)

    args.db_file.parent.mkdir(exist_ok=True, parents=True)
    shutil.copy(db_file, args.db_file)


if __name__ == "__main__":
    main()
