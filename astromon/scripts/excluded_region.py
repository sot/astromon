#!/usr/bin/env python
"""
Script to add excluded regions to astromon.db
"""

import argparse
import getpass
import logging

from astropy.coordinates import SkyCoord
from ska_helpers.logging import basic_logger

from astromon import db

logger = basic_logger("regions", level="INFO")


class ArgumentError(Exception):
    pass


def get_parser():
    config_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    config_levels += [x.lower() for x in config_levels]

    parser = argparse.ArgumentParser(
        description="Add excluded regions to astromon.db",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "action",
        choices=list(ACTIONS),
        help="Action to perform",
    )
    parser.add_argument(
        "--ra",
        help="Right Ascension of the center of the excluded region (in degrees)",
    )
    parser.add_argument(
        "--dec",
        help="Declination of the center of the excluded region (in degrees)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=5.0,
        help="Radius of the excluded region (in arcseconds)",
    )
    parser.add_argument(
        "--username",
        help="Username (default: current user)",
        default=getpass.getuser(),
    )
    parser.add_argument(
        "--comment",
        help="Comment for the excluded region (truncated to 200 characters)",
        default="",
    )
    parser.add_argument(
        "--obsid",
        type=int,
        default=0,
        help="OBSID associated with the excluded region (default is to apply to all OBSIDs)",
    )
    parser.add_argument(
        "--region-id",
        type=int,
        help="Region ID to list, modify or remove (ignored when adding a region)",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        choices=config_levels,
        default="INFO",
    )
    return parser


def main(args=None):
    parser = get_parser()
    args = parser.parse_args(args)
    args.log_level = args.log_level.upper()

    logger.setLevel(args.log_level)
    # this downgrades astromon's INFO to WARNING to avoid too much output
    logging.getLogger("astromon").setLevel(
        args.log_level if args.log_level != "INFO" else "WARNING"
    )

    try:
        ACTIONS[args.action](args)
    except ArgumentError as e:
        parser.error(str(e))
    except Exception as e:
        logger.error(f"Error performing action {args.action}: {e}")
        raise


def add_region(args):
    if args.ra is None or args.dec is None:
        raise ArgumentError("RA and DEC must be specified when adding a region")
    try:
        # if both can be converted to float, we assume they are correct ra/dec in degrees
        args.ra = float(args.ra)
        args.dec = float(args.dec)
    except ValueError:
        try:
            c = SkyCoord(args.ra, args.dec)
            args.ra = c.ra.deg
            args.dec = c.dec.deg
        except ValueError:
            get_parser().error(
                "RA and DEC must be specified in degrees"
                " or in a format understood by astropy.coordinates.SkyCoord"
            )

    result = db.add_regions(
        regions=[
            {
                "ra": args.ra,
                "dec": args.dec,
                "radius": args.radius,
                "user": args.username,
                "comments": args.comment[:200],
                "obsid": args.obsid,
            }
        ]
    )

    region = dict(result[0])
    logger.info(OUT.format(**region))


def remove_region(args):
    if args.region_id is None:
        raise ValueError("region_id is required to remove a region")
    db.remove_regions([args.region_id])


def list_regions(args):
    regions = db.get_regions(args.obsid if args.obsid > 0 else None)
    if args.region_id is not None:
        regions = regions[regions["region_id"] == args.region_id]
    regions.pprint(max_width=-1, max_lines=-1)


OUT = """Added region:
    ID:       {region_id}
    RA:       {ra}
    DEC:      {dec}
    Radius:   {radius}
    OBSID:    {obsid}
    User:     {user}
    Comment:  {comments}
"""

ACTIONS = {
    "add": add_region,
    "remove": remove_region,
    "rm": remove_region,
    "list": list_regions,
}

if __name__ == "__main__":
    main()
