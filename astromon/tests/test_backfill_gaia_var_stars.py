"""Tests for astromon.scripts.maintenance.backfill_gaia_var_stars."""

import numpy as np
from astropy.table import Table

from astromon import db
from astromon.scripts.maintenance import backfill_gaia_var_stars as bf

OBSID = 9001


def _obs_row(obsid=OBSID, ra=150.0, dec=2.0):
    row = np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE)
    row["obsid"] = obsid
    row["ra"] = ra
    row["dec"] = dec
    row["date_obs"] = "2020-01-01T00:00:00"
    return Table(row)


def _xray_row(obsid, source_id, ra, dec):
    row = np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = source_id
    row["ra"] = ra
    row["dec"] = dec
    row["detect_method"] = "celldetect"
    return Table(row)


def _varstar_catalog(ra, dec, source_id=1):
    return Table(
        {
            "source_id": [source_id],
            "ra": [ra],
            "dec": [dec],
            "pmra": [0.0],
            "pmdec": [0.0],
            "phot_g_mean_mag": [15.0],
        }
    )


def test_build_varstar_cat_src_matches_across_ra_meridian():
    """A real GaiaVarStar match just across RA=0/360 must not be dropped.

    The cone pre-filter computed dra = (cat_ra - aimpoint_ra) * cos_dec with no
    wraparound handling: an obsid at ra=0.05 and a catalog source at ra=359.98
    are ~0.07 deg apart on the sky, but the unwrapped difference is ~-359.93,
    so the cone check silently rejected the match.
    """
    aimpoint_ra, dec = 0.05, 10.0
    source_ra = 359.98

    new_obs = _obs_row(ra=aimpoint_ra, dec=dec)
    new_xray = _xray_row(OBSID, 1, aimpoint_ra, dec)
    full_catalog = _varstar_catalog(source_ra, dec)
    existing_cat = Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)

    new_cat = bf.build_varstar_cat_src(full_catalog, new_xray, new_obs, existing_cat)

    assert len(new_cat) == 1, (
        "a source ~0.07 deg away across the RA=0/360 meridian must be matched"
    )
