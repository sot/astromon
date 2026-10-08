"""Find completed observations more recent than anything astromon has.

Why this exists
----------------
The normal new-observation discovery path (``get_cat_obs_data.get_obsids``)
queries the archive through ``arc5gl``, which only works from a machine on
the CXC-internal network. This Mac processes obsids through the public
Chandra Data Archive instead (``--source archive``), so discovery needs a
public-archive-compatible source too: ``mica.archive.cda.get_ocat_local``
reads a local, cron-synced copy of the Ocat target table -- no network
round-trip, and no arc5gl dependency.

This only finds candidates; it does not decide whether an obsid is suitable
for astromon processing. ``Observation.is_selected`` already does that
(excludes the engineering/cal-replay obsid block, non-pointing observations,
ACIS non-TE/duty-cycled readouts, etc.), and process_one_obsid.py already
handles an obsid that turns out to be embargoed or not yet public
(``ObsidNotPubliclyAvailable`` -> ``status="skipped_not_public"``). Reusing
that gate instead of re-deriving a "science observation" filter from ocat's
own type/category fields keeps eligibility decided in exactly one place.

Usage
-----
::

    python -m astromon.scripts.maintenance.find_new_observations \\
        --db-file /Users/jean/ska/data/astromon/astromon_gaussian_redo.h5 \\
        --obsid-list-out /tmp/new_obsids.txt
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from cxotime import CxoTime

from astromon import db

logger = logging.getLogger("find_new_observations")

COMPLETED_STATUSES = ("archived", "observed")


def find_new_candidate_obsids(
    ocat, known_obsids: set[int], since: CxoTime
) -> list[int]:
    """Obsids in `ocat` observed strictly after `since`, not already known.

    `ocat` is a table with at least `obsid`, `status` and `start_date`
    columns, as returned by `mica.archive.cda.get_ocat_local`. Only rows
    whose status marks the observation as having actually happened
    (`COMPLETED_STATUSES`) are considered -- scheduled, canceled,
    unobserved, untriggered and discarded rows are not candidates, whatever
    their start_date.
    """
    status = np.asarray(ocat["status"]).astype(str)
    raw_start_date = np.asarray(ocat["start_date"]).astype(str)

    # A row with no real start_date yet carries the literal string "0"
    # instead of a date ocat could parse -- status alone does not predict
    # this: every "archived" row has a real date, but most "observed" rows
    # (ocat's engineering/replay block) do not. Without a real date there is
    # no way to tell whether a row is new, so it is dropped before CxoTime
    # ever sees it rather than guessed at.
    has_date = raw_start_date != "0"
    usable = np.isin(status, COMPLETED_STATUSES) & has_date

    obsid = np.asarray(ocat["obsid"]).astype(int)[usable]
    start_date = CxoTime(raw_start_date[usable])

    is_new = start_date.secs > since.secs
    candidates = [int(o) for o in obsid[is_new] if int(o) not in known_obsids]
    return sorted(candidates)


def main() -> None:
    from mica.archive import cda

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-file", required=True, type=Path)
    parser.add_argument(
        "--since",
        default=None,
        help="Only consider obsids observed after this date. "
        "Default: the latest date_obs already in astromon_obs.",
    )
    parser.add_argument(
        "--obsid-list-out",
        default=None,
        type=Path,
        help="Write the candidate obsids here, one per line, suitable for "
        "run_all.py's --obsid-list. Default: print them instead.",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    obs = db.get_table("astromon_obs", args.db_file)
    known_obsids = {int(o) for o in np.asarray(obs["obsid"])}

    if args.since is not None:
        since = CxoTime(args.since)
    else:
        date_obs = np.asarray(obs["date_obs"]).astype(str)
        since = CxoTime(date_obs[np.argmax(CxoTime(date_obs).secs)])
    logger.info(f"latest date_obs already known: {since.date}")

    ocat = cda.get_ocat_local()
    candidates = find_new_candidate_obsids(ocat, known_obsids, since)
    logger.info(f"candidate obsids observed since {since.date}: {len(candidates)}")

    if args.obsid_list_out is not None:
        args.obsid_list_out.write_text("\n".join(str(o) for o in candidates) + "\n")
        logger.info(f"written to {args.obsid_list_out}")
    else:
        for o in candidates:
            print(o)


if __name__ == "__main__":
    main()
