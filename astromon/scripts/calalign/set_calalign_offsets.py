#!/usr/bin/env python

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from cxotime import CxoTime

A2R = np.pi / 180 / 3600  # conversion from arcsec to rad


def update_tstop(in_file, out_file, tstop, date, clobber=False):
    hdus = fits.open(in_file)
    for hdu in hdus:
        hdu.header["DATE"] = date.isot.split(".")[0]
        hdu.header["TSTOP"] = np.floor(tstop.secs - 1)
    if clobber and out_file.exists():
        os.remove(out_file)
    hdus.writeto(out_file)


def apply_calalign_shift(tstart, tstop, in_file, out_file, dy, dz, date, clobber=False):
    hdus = fits.open(in_file)
    for hdu in hdus:
        hdu.header["DATE"] = date.isot.split(".")[0]
        hdu.header["TSTART"] = np.floor(tstart.secs)
        hdu.header["TSTOP"] = np.floor(tstop.secs - 1)
        hdu.header["CVSD0001"] = tstart.isot.split(".")[0]
        hdu.header["CVST0001"] = tstart.isot.split("T")[1].split(".")[0]

    hdus[0].header["HISTORY"] = "CALALIGN updated"

    assert hdus[1].name == "CALALIGN"
    for row in hdus[1].data:
        sy = np.sin(dy * A2R)
        cy = np.cos(dy * A2R)
        sz = np.sin(dz * A2R)
        cz = np.cos(dz * A2R)

        Rz = np.array([[cz, 0.0, -sz], [0.0, 1.0, 0.0], [sz, 0.0, cz]])
        Ry = np.array([[cy, -sy, 0.0], [sy, cy, 0.0], [0.0, 0.0, 1.0]])
        R = Rz @ Ry

        in_aca = row["ACA_MISALIGN"]
        in_fts = row["FTS_MISALIGN"]

        # in the following, we use the fact that the transpose of these matrices is their inverse:

        # For FTS align drift
        # row['ACA_MISALIGN'] = in_aca
        # row['FTS_MISALIGN'] = R.T @ in_fts

        # For ACA misalign drift
        row["ACA_MISALIGN"] = R @ in_aca
        row["FTS_MISALIGN"] = in_aca.T @ R.T @ in_aca @ in_fts

        msg = (
            f"{row['INSTR_ID'].strip()} alignment shift (dy, dz) = ({dy:.2f}, {dz:.2f})"
        )
        # print(msg)
        hdus[0].header["HISTORY"] = msg
    if clobber and out_file.exists():
        os.remove(out_file)
    out_file.parent.mkdir(exist_ok=True, parents=True)
    hdus.writeto(out_file)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--offsets",
        default=Path("offsets.json"),
        type=Path,
        help="JSON file with offsets. The output of calc_median.py. Default: offsets.json",
    )
    parser.add_argument(
        "--clobber",
        default=False,
        action="store_true",
        help="Overwrite existing files.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("files"),
        help="Output directory. Default: files",
    )
    parser.add_argument(
        "--version", default="N0010", help="The new file version. Default: N0010."
    )
    parser.add_argument(
        "--base",
        type=Path,
        default=Path("pcadD2013-01-19alignN0009.fits"),
        help=(
            "Align file to be used as the basis for all output files. "
            "Default: pcadD2013-01-19alignN0009.fits"
        ),
    )
    parser.add_argument(
        "--calalign-dir",
        type=Path,
        default=Path("/data/caldb/data/chandra/pcad/align"),
        help=(
            "Directory where CALALIGN files are located. "
            "Default: /data/caldb/data/chandra/pcad/align"
        ),
    )
    return parser


def main():
    args = get_parser().parse_args()

    args.base = args.calalign_dir / args.base

    date = CxoTime(CxoTime().isot, scale="tt")

    with open(args.offsets) as fh:
        offsets = json.load(fh)

    for sh in offsets:
        start = CxoTime(sh["tstart"], scale="tt").iso.split()[0]
        out_file = args.out_dir / f"pcadD{start}align{args.version}.fits"
        apply_calalign_shift(
            tstart=CxoTime(sh["tstart"], scale="tt"),
            tstop=CxoTime(sh["tstop"], scale="tt"),
            in_file=args.base,
            out_file=out_file,
            # note the minus sign, because the CALALIGN offset corrects for the measured offset
            dy=-sh["dy"],
            dz=-sh["dz"],
            date=date,
            clobber=args.clobber,
        )

    m = re.match(
        r"pcad(?P<date>D[0-9\-]+)align(?P<version>N[0-9]+).fits", args.base.name
    )
    update_tstop(
        in_file=args.base,
        tstop=CxoTime(offsets[0]["tstart"], scale="tt"),
        out_file=args.out_dir
        / args.base.name.replace(m.groupdict()["version"], args.version),
        date=date,
        clobber=args.clobber,
    )


if __name__ == "__main__":
    main()
