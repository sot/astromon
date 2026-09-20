"""Remediate the historical obspar-vs-evt2 pointing mismatch in the astromon DB.

STATUS: executed against /Volumes/Black/data/astromon/astromon.h5 on 2026-09-04
(backup: astromon.h5.backup-before-pointing-fix). Updated all 13378 astromon_obs
rows and recomputed y_angle/z_angle/r_angle for all 197662 astromon_xray_src and
94449 astromon_cat_src rows. Verified: r_angle >= 180" dropped from 9916/193030
(5.1%) to 2/197662 (0.00%) dataset-wide. Kept here for the record and in case a
similar remediation is needed again; DRY_RUN below defaults to True, as it did
for the original DRY_RUN pass before this was run for real.

Background
----------
astromon_obs.ra/dec/roll, and astromon_xray_src/astromon_cat_src's y_angle/z_angle/
r_angle, were computed (historically) from Observation.get_obspar()'s ra_pnt/dec_pnt/
roll_pnt. mica's local obspar archive holds only quick-look ("ql" mode) obspars for
this entire dataset (confirmed: 114131 of 114131 local records) -- never the final,
standard-processing product. The evt2 file's own header (RA_PNT/DEC_PNT/ROLL_PNT,
read by Observation.get_evt2_info(), added this session) is the correct, measured
pointing -- confirmed against the actual aspect solution for one example (obsid
16608: evt2 RA_PNT matched the mean of 135,186 real aspect-solution samples to 6
decimal places; the quick-look obspar did not, by 38.6 arcsec).

Across the 13378 obsids in astromon_obs, 12456 (93%) show a pointing mismatch
> 10 arcsec between the stored value and the evt2 header, typically ~40 arcsec.
See /tmp/pointing_mismatch_gt10arcsec.csv for the full list (obsid, separation,
stored ra/dec, evt2 ra_pnt/dec_pnt), built from local evt2 files (12396 obsids,
paths in /tmp/evt2_paths.txt) plus freshly-downloaded fov files for the remaining
982 (under .../scratchpad/fov_missing/<obsid>/primary/*fov*, fov's RA_PNT/DEC_PNT
confirmed identical to evt2's for a spot-checked example).

What this script does
----------------------
1. astromon_obs.ra/dec/roll -> evt2/fov RA_PNT/DEC_PNT/ROLL_PNT, for every obsid.
2. astromon_xray_src.y_angle/z_angle/r_angle -> recomputed from each source's
   existing RA/DEC (unaffected by the bug -- celldetect/gaussian_detect fit sky
   pixels via the image's own WCS, not obspar) rotated through the corrected
   pointing quaternion for that source's obsid.
3. astromon_cat_src.y_angle/z_angle -> same recomputation, for consistency (its
   own r_angle is not stored in this table, so there's no admission-gate
   consequence on this side -- see the NOT-done item below).

What this script deliberately does NOT do
-------------------------------------------
- Does not touch dy/dz/dr in any already-persisted match table: dy = x_y_angle -
  c_y_angle is a difference of two angles that were both rotated by the SAME
  (wrong) quaternion for a given obsid, so a common pointing error cancels to
  first order in that difference. Verified: both sides used the same quaternion
  source historically (get_obspar(), single cached call per obsid) and now use
  the same evt2_info source -- there was never a split between the two sides.
- Does not re-run detection (celldetect/gaussian_detect/peak_gaussian_detect):
  RA/DEC per source is independent of this bug.
- Does not re-run the r_angle admission gate in cross-matching
  (cross_match.py:~2433-2436, `matches["r_angle"] < r_angle` at the default
  120"/24"-grating threshold). That gate uses the x-ray source's ABSOLUTE
  r_angle, not a difference, so it is NOT protected by the cancellation
  argument above: a source whose true r_angle sits within ~40" of that
  boundary could have been wrongly admitted or excluded. Fixing this properly
  means re-running rough_match for boundary-region sources against the
  relevant catalogs (astromon_cat_src may not even have an entry for a source
  that was previously excluded before ever being queried) -- a bigger, separate
  task, not just an angle recomputation. Left as a TODO; do not assume it is
  covered by this script.

Usage
-----
Review DRY_RUN below and the printed summary before ever setting it to False.
"""

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / "git" / "astromon_hrc"))
import os

os.environ.setdefault("SKA", str(Path.home() / "ska"))

import numpy as np
from astropy.io import fits
from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat

from astromon import db

DBFILE = "/Volumes/Black/data/astromon/astromon.h5"

# Toggle to False only after reviewing the printed summary from a DRY_RUN pass.
DRY_RUN = True

FOV_ROOT = Path(
    "/private/tmp/claude-502/-Users-jean-git-astromon-project--claude-worktrees"
    "-pytest-hanging-astromon-0a071a/51a68aa6-6287-4012-9f75-7e644f34b7b2"
    "/scratchpad/fov_missing"
)


def backup_dbfile():
    """Copy DBFILE to a sibling ``.backup-before-pointing-fix``.

    Matches the naming convention already used elsewhere in the archive (e.g.
    astromon.h5.backup-before-status-backfill). Refuses to overwrite an existing
    backup -- if one is already there, either it's from this same remediation
    (nothing to redo) or from something else, and clobbering it silently would be
    exactly the kind of mistake this backup is meant to protect against.
    """
    dbfile = Path(DBFILE)
    backup = dbfile.with_name(dbfile.name + ".backup-before-pointing-fix")
    if backup.exists():
        raise FileExistsError(
            f"{backup} already exists -- not overwriting. Remove it first if you "
            "really want a fresh backup, or restore from it if this remediation "
            "already ran."
        )
    print(f"backing up {dbfile} -> {backup}")
    shutil.copy2(dbfile, backup)
    return backup


def get_corrected_pointing():
    """Return {obsid: (ra_pnt, dec_pnt, roll_pnt)} from evt2 (preferred) or fov headers.

    Uses the evt2 paths already gathered in /tmp/evt2_paths.txt (12396 obsids) and the
    fov files fetched into FOV_ROOT for the remaining 982. Re-run the gathering step in
    this session's transcript if those files are no longer present (e.g. a fresh
    worktree/session) -- this function does not re-download anything itself.
    """
    evt2_paths = {}
    with open("/tmp/evt2_paths.txt") as f:
        for line in f:
            oid, p = line.strip().split("\t")
            evt2_paths[int(oid)] = p

    pointing = {}
    n_missing = 0
    for oid, path in evt2_paths.items():
        h = fits.getheader(path, 1)
        pointing[oid] = (float(h["RA_PNT"]), float(h["DEC_PNT"]), float(h["ROLL_PNT"]))

    for d in FOV_ROOT.glob("*/primary"):
        oid = int(d.parent.name)
        if oid in pointing:
            continue
        cand = list(d.glob("*fov*"))
        if not cand:
            n_missing += 1
            continue
        h = fits.getheader(str(cand[0]), 1)
        pointing[oid] = (float(h["RA_PNT"]), float(h["DEC_PNT"]), float(h["ROLL_PNT"]))

    print(f"corrected pointing for {len(pointing)} obsids, {n_missing} still missing")
    return pointing


def update_astromon_obs(pointing):
    obs = db.get_table("astromon_obs", DBFILE)
    obs = obs.copy()
    n_updated = 0
    for i, oid in enumerate(np.asarray(obs["obsid"])):
        p = pointing.get(int(oid))
        if p is None:
            continue
        obs["ra"][i], obs["dec"][i], obs["roll"][i] = p
        n_updated += 1
    print(f"astromon_obs: {n_updated} of {len(obs)} rows updated")

    if not DRY_RUN:
        db.save("astromon_obs", obs, DBFILE, ignore_obsid=True)
    return obs


def recompute_angles(table_name, pointing):
    """Recompute y_angle/z_angle (and r_angle, if present) in-place from RA/DEC."""
    t = db.get_table(table_name, DBFILE)
    t = t.copy()
    obsids = np.asarray(t["obsid"])
    ra = np.asarray(t["RA" if "RA" in t.colnames else "ra"], dtype=float)
    dec = np.asarray(t["DEC" if "DEC" in t.colnames else "dec"], dtype=float)

    # Start from the stored values, not zero: an obsid missing from `pointing`
    # (get_corrected_pointing logs these as n_missing) must keep its existing
    # angles rather than have them overwritten with 0.0, which reads as "source
    # located exactly at the aimpoint" instead of "not recomputed".
    y_angle = np.asarray(t["y_angle"], dtype=np.float32).copy()
    z_angle = np.asarray(t["z_angle"], dtype=np.float32).copy()
    updated = np.zeros(len(t), dtype=bool)
    for oid in np.unique(obsids):
        p = pointing.get(int(oid))
        if p is None:
            continue
        sel = obsids == oid
        q = Quat(equatorial=p)
        y_angle[sel], z_angle[sel] = radec_to_yagzag(ra[sel], dec[sel], q)
        updated |= sel

    t["y_angle"] = y_angle
    t["z_angle"] = z_angle
    if "r_angle" in t.colnames:
        r_angle = np.asarray(t["r_angle"], dtype=np.float32).copy()
        r_angle[updated] = np.sqrt(y_angle[updated] ** 2 + z_angle[updated] ** 2)
        t["r_angle"] = r_angle
    print(f"{table_name}: {int(updated.sum())} of {len(t)} rows recomputed")

    if not DRY_RUN:
        db.save(table_name, t, DBFILE, ignore_obsid=True)
    return t


def main():
    if not DRY_RUN:
        backup_dbfile()
    pointing = get_corrected_pointing()
    update_astromon_obs(pointing)
    recompute_angles("astromon_xray_src", pointing)
    recompute_angles("astromon_cat_src", pointing)
    if DRY_RUN:
        print(
            "\nDRY_RUN is True -- nothing was written. Review the summary above, "
            "then set DRY_RUN = False to actually save."
        )


if __name__ == "__main__":
    main()
