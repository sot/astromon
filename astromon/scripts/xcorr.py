#!/usr/bin/env python
"""
"""

import logging
import argparse
import pyyaks
from pathlib import Path

from astromon.cross_match import cross_match
from astromon import db


logger = logging.getLogger('astromon')


def get_parser():
    """
    Get the argument parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--db-file',
        help='SQLite file where data is saved',
        required=True,
        type=Path
    )
    parser.add_argument(
        '--selection',
        help='Name of a standard selection for cross-matching',
        default=None,
    )
    parser.add_argument(
        '--log-level',
        help='Logging severity level',
        choices=['debug', 'info', 'warning', 'error', 'fatal']
    )
    return parser


def main():
    """
    Main routine that deals with processing arguments, output and such.
    """
    args = get_parser().parse_args()

    pyyaks.logger.get_logger(
        name='astromon',
        level=args.log_level.upper(),
        format="%(asctime)s %(message)s"
    )

    assert args.db_file.exists(), f'File does not exist: {args.db_file}'

    x_match = cross_match(args.selection)

    db.save('astromon_xcorr', x_match)


if __name__ == '__main__':
    main()
