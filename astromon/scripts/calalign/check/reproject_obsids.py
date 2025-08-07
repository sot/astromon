#!/usr/bin/env python

import subprocess
import sys
import traceback
from pathlib import Path

import obsids_to_check
import pyyaks.logger
import Ska.arc5gl
from astropy.io import ascii
from Ska.File import chdir

from astromon.utils import Ciao


def celldetect(evt, src, si, asol):
    CIAO = Ciao(
        workdir=evt.parent, logger=pyyaks.logger.get_logger("ciao", level="WARNING")
    )

    pixel = 1 if si == "hrc" else 0.5
    band = "wide" if si == "hrc" else "broad"
    radius = 180
    workdir = evt.parent
    root = evt.name.replace("_evt2.fits", "")
    evt2 = workdir / f"{root}_filtered_evt2.fits"
    CIAO("dmkeypar", evt, "RA_PNT")
    ra = CIAO.pget("dmkeypar")
    CIAO("dmkeypar", evt, "DEC_PNT")
    dec = CIAO.pget("dmkeypar")

    CIAO("dmcoords", evt, op="cel", celfmt="deg", ra=ra, dec=dec)
    x = CIAO.pget("dmcoords", "x")
    y = CIAO.pget("dmcoords", "y")

    CIAO(
        "dmcopy",
        f"{evt}[(x,y)=circle({x},{y},{radius / pixel})]",
        evt2,
        clobber="yes",
    )

    CIAO(
        "fluximage",
        infile=evt2,
        outroot=str(workdir / root),
        bands=band,
        binsize="1",
        psfecf="0.9",
        background="none",
        asol=asol,
    )
    if not (workdir / f"{root}_{band}_thresh.img").exists():
        raise Exception(f"{workdir}/{root}_{band}_thresh.img does not exist")
    if not (workdir / f"{root}_{band}_thresh.psfmap").exists():
        raise Exception(f"{workdir}/{root}_{band}_thresh.psfmap does not exist")

    CIAO(
        "celldetect",
        workdir / f"{root}_{band}_thresh.img",
        src,
        psffile=workdir / f"{root}_{band}_thresh.psfmap",
        thresh=3.0,
        maxlogicalwindow=2048,
        clobber="yes",
    )


def run_cmd(*cmd):
    print(" ".join(cmd))
    subprocess.run(cmd, check=False)


def main():  # noqa: PLR0915
    RUNASP = "/proj/sot/ska/jgonzalez/aca_cal_align/update_2022-feb/runasp/runasp.py"

    for _, obsid, i in obsids_to_check.OBS:
        try:
            obsdir = (
                Path("reprodata").absolute()
                / f"obs{obsid // 1000:02d}"
                / f"{obsid:05d}"
            )
            obsdir.mkdir(parents=True, exist_ok=True)
            CIAO = Ciao(
                workdir=obsdir, logger=pyyaks.logger.get_logger("ciao", level="INFO")
            )

            with chdir(obsdir):
                # make an aspect solution using the new alignment
                new_asols = list(obsdir.glob("archive/ASP*/out1/*asol*"))
                if len(new_asols) < 1:
                    print(f"running aspect pipeline on {obsid}")
                    align = obsids_to_check.ALIGN_DIR / obsids_to_check.ALIGN[i]
                    run_cmd(
                        "python",
                        RUNASP,
                        "--obsid",
                        str(obsid),
                        "--param",
                        f"align={align}",
                        "--dir",
                        "archive",
                    )
                    new_asols = list(obsdir.glob("archive/ASP*/out1/*asol*"))
                new_asol = new_asols[0]

                # make an aspect solution using the current pipeline and alignment
                old_asols = list(obsdir.glob("archive/uncorr/out1/*asol*"))
                if len(old_asols) < 1:
                    print("running aspect pipeline on {} uncorrected".format(obsid))
                    run_cmd(
                        "python",
                        RUNASP,
                        "--obsid",
                        str(obsid),
                        "--label",
                        "uncorr",
                        "--dir",
                        "archive",
                    )
                    old_asols = list(obsdir.glob("archive/uncorr/out1/*asol*"))
                old_asol = old_asols[0]

            obsid_info = {
                r[0]: r[3]
                for r in ascii.read(list(obsdir.glob("archive/obspar/*obs0a.par"))[0])
            }
            si = obsid_info["instrume"].lower()

            (obsdir / "reproject").mkdir(exist_ok=True)

            evt = list(obsdir.glob(f"reproject/{si}f*evt2*"))
            if not evt:
                with chdir(obsdir / "reproject"):
                    print(f"Fetching archive OBSID {obsid} files")
                    arc = Ska.arc5gl.Arc5gl()
                    arc.sendline("obsid={}".format(obsid))
                    arc.sendline(f"get {si}2{{evt2}}")
                    arc.sendline("get asp1{aspsol}")
                    arc.sendline("get obspar")
                    arc.sendline(f"get {si}1{{msk}}")
                    arc.sendline(f"get {si}1{{fov}}")
                    arc.sendline(f"get {si}1[*bpix*]")
                    del arc
                evt = list(obsdir.glob(f"reproject/{si}f*evt2*"))
            evt = evt[0]

            corr_evt2 = obsdir / "reproject" / f"{si}_corr_evt2.fits"
            if not corr_evt2.exists():
                print(f"reprojecting {obsid} using the new solution")
                CIAO("reproject_events", evt, corr_evt2, match="none", aspect=new_asol)

            uncorr_evt2 = obsdir / "reproject" / f"{si}_uncorr_evt2.fits"
            if not uncorr_evt2.exists():
                print(f"reprojecting {obsid} using the old solution")
                CIAO(
                    "reproject_events", evt, uncorr_evt2, match="none", aspect=old_asol
                )

            new_src2 = obsdir / "reproject" / "new_src2.fits"
            if not new_src2.exists():
                print(f"making new src list for {obsid}")
                celldetect(evt=corr_evt2, src=new_src2, si=si, asol=new_asol)
                if not new_src2.exists():
                    raise Exception("Failed creating new source list")

            old_src2 = obsdir / "reproject" / "old_src2.fits"
            if not old_src2.exists():
                print(f"making old src list for {obsid}")
                celldetect(evt=uncorr_evt2, src=old_src2, si=si, asol=old_asol)
                if not old_src2.exists():
                    raise Exception("Failed creating old source list")
        except Exception as e:
            print(f"failed {obsid}: {e}")
            exc_type, exc_value, exc_traceback = sys.exc_info()
            trace = traceback.extract_tb(exc_traceback)
            for step in trace:
                print(f"  in {step.filename}:{step.lineno}/{step.name}:")
                print(f"    {step.line}")


if __name__ == "__main__":
    main()
