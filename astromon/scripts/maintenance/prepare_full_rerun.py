"""Generate an obsid list for the full pipeline rerun from the DB.

Writes all obsids in astromon_obs to a file, then prints the run_all.py
command to run against preserve-workdir.

Usage::

    python -m astromon.scripts.maintenance.prepare_full_rerun [--db PATH] [--out obsids_rerun.txt]
"""

import argparse
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("prepare_full_rerun")

_DEFAULT_DB = Path("/Volumes/Black/data/astromon/astromon.h5")
_DEFAULT_OUT = Path("/Volumes/Black/data/astromon/obsids_full_rerun.txt")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        type=Path,
        default=_DEFAULT_DB,
        help=f"Path to astromon.h5 (default: {_DEFAULT_DB})",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=_DEFAULT_OUT,
        help=f"Output obsid list file (default: {_DEFAULT_OUT})",
    )
    args = parser.parse_args()

    import sys  # noqa: PLC0415

    sys.path.insert(0, str(Path(__file__).parent))
    from astromon import db  # noqa: PLC0415

    logger.info(f"Loading obsids from {args.db}")
    obs = db.get_table("astromon_obs", dbfile=args.db)
    obsids = sorted(int(x) for x in obs["obsid"])
    logger.info(f"Found {len(obsids)} obsids in astromon_obs")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(str(o) for o in obsids) + "\n")
    logger.info(f"Wrote {len(obsids)} obsids to {args.out}")

    workdir = "/Volumes/Black/data/astromon/work"
    db_file = str(args.db)
    tracking = "/Volumes/Black/data/astromon/tracking_full_rerun.csv"
    log_dir = "/Volumes/Black/data/astromon/logs_full_rerun"

    print("\n# --- Full pipeline rerun command ---")
    print(
        f"python -m astromon.scripts.maintenance.run_all \\\n"
        f"    --db-file {db_file} \\\n"
        f"    --workdir {workdir} \\\n"
        f"    --obsid-list {args.out} \\\n"
        f"    --tracking-csv {tracking} \\\n"
        f"    --log-dir {log_dir} \\\n"
        f"    --parallel 4 \\\n"
        f"    --source archive \\\n"
        f"    --versions celldetect gaussian_detect peak_gaussian_detect"
    )
    print(
        "\n# Notes:"
        "\n#  - 1,785 obsids with invalidated make_images cache will re-run fluximage + acis_streak_map(ssigma=3)"
        "\n#  - All other preserve-workdir obsids skip CIAO tasks (cached) and only re-run catalog+xcorr stage"
        "\n#  - Pool F retry obsids not yet in preserve-workdir will re-download from public archive"
        "\n#  - resume-safe: re-run the same command if interrupted (tracking-csv tracks done obsids)"
    )


if __name__ == "__main__":
    main()
