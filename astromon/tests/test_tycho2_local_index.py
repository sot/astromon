"""
Regression test for the local Tycho-2 index (astromon/tycho2_local_index.py).

This checks that a locally-built Tycho-2 spatial index (parsed from the ESO Science
Archive mirror of CDS I/259, see the Part 2 Tycho-2 cutover validation report) continues
to reproduce the "tycho2" cross-match results currently stored in the production
database (astromon_gaussian_redo.h5), without needing access to that database or to any
live network service at test time.

The held-out expected values come from a small fixture file
(astromon/tests/data/tycho2_local_index_fixture.ecsv) generated ahead of time by
scratchpad/scripts/regen_tycho2_fixture.py (not part of this repository) from the
production DB. To regenerate the fixture or re-run the full ~2000-row validation, see
that script's docstring.

The built index itself is a few hundred MB (parsed from ~2.5 million Tycho-2 + Tycho-2
supplement stars) and is therefore NOT committed to this repository. This test is
skipped cleanly if it can't find one:

    - set the ASTROMON_TYCHO2_INDEX environment variable to the built
      "tyc2_stars.parquet" file, or
    - place it at $SKA/data/astromon/tycho2_local_index/tyc2_stars.parquet (the
      documented default location).

A second, separate test (test_known_cone_edge_cases) pins four specific x-ray sources
found during the validation campaign where the local index legitimately disagrees with
the production DB for understood reasons (documented inline) -- not local-index bugs.
"""

import os
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table
from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat

from astromon.tycho2_local_index import (
    add_yagzag,
    admit_tycho2_match,
    build_kdtree,
    load_star_table,
    max_pm_mas_per_yr,
    query_cone,
)

DATA_DIR = Path(__file__).parent / "data"
FIXTURE_PATH = DATA_DIR / "tycho2_local_index_fixture.ecsv"

ROUGH_MATCH_RADIUS_ARCSEC = 3.0
SIMPLE_XMATCH_ARGS = {
    "snr": 3,
    "r_angle": 120.0,
    "r_angle_grating": 120.0,
    "near_neighbor_dist": 6.0,
    "dr": 3.0,
}

# generous relative to the sub-0.1 arcsec / exact-magnitude agreement actually observed
# in the validation campaign (see the Part 2 report); this only needs to catch real
# regressions, not pin floating-point noise.
POSITION_TOLERANCE_ARCSEC = 0.5
MAG_TOLERANCE = 0.02

assert FIXTURE_PATH.exists(), f"missing fixture file {FIXTURE_PATH}"
FIXTURE_TABLE = Table.read(FIXTURE_PATH)


def _default_index_path() -> Path | None:
    env_path = os.environ.get("ASTROMON_TYCHO2_INDEX")
    if env_path:
        return Path(env_path)
    if "SKA" in os.environ:
        default = (
            Path(os.environ["SKA"])
            / "data"
            / "astromon"
            / "tycho2_local_index"
            / "tyc2_stars.parquet"
        )
        if default.exists():
            return default
    return None


@pytest.fixture(scope="module")
def tycho2_index():
    path = _default_index_path()
    if path is None or not path.exists():
        pytest.skip(
            "local Tycho-2 index not found; set the ASTROMON_TYCHO2_INDEX environment "
            "variable to a built tyc2_stars.parquet file (see "
            "astromon/tycho2_local_index.py) to run this test"
        )
    stars = load_star_table(str(path))
    tree = build_kdtree(stars)
    global_max_pm = max_pm_mas_per_yr(stars)
    return stars, tree, global_max_pm


def _recompute_match(
    stars,
    tree,
    global_max_pm,
    *,
    x_ra,
    x_dec,
    x_snr,
    x_near_neighbor_dist,
    x_r_angle,
    x_grating,
    quat,
    frac_year,
):
    candidates = query_cone(
        tree,
        stars,
        x_ra,
        x_dec,
        frac_year,
        radius_arcsec=ROUGH_MATCH_RADIUS_ARCSEC,
        global_max_pm_mas_per_yr=global_max_pm,
    )
    candidates = add_yagzag(candidates, quat)
    x_y_angle, x_z_angle = radec_to_yagzag(x_ra, x_dec, quat)
    _, admitted = admit_tycho2_match(
        candidates,
        x_snr,
        x_near_neighbor_dist,
        x_r_angle,
        x_grating,
        x_y_angle,
        x_z_angle,
        **SIMPLE_XMATCH_ARGS,
    )
    return admitted


@pytest.mark.parametrize(
    "row",
    list(FIXTURE_TABLE),
    ids=[f"obsid{row['obsid']}_x{row['x_id']}" for row in FIXTURE_TABLE],
)
def test_matches_stored_cross_match(tycho2_index, row):
    stars, tree, global_max_pm = tycho2_index

    quat = Quat(
        equatorial=(
            row["quat_ra_pnt_fit"],
            row["quat_dec_pnt_fit"],
            row["quat_roll_pnt_fit"],
        )
    )
    admitted = _recompute_match(
        stars,
        tree,
        global_max_pm,
        x_ra=row["x_ra"],
        x_dec=row["x_dec"],
        x_snr=row["x_snr"],
        x_near_neighbor_dist=row["x_near_neighbor_dist"],
        x_r_angle=row["x_r_angle"],
        x_grating=row["x_grating"],
        quat=quat,
        frac_year=row["frac_year"],
    )

    assert (admitted is not None) == bool(row["expected_has_match"])
    if not row["expected_has_match"]:
        return

    assert admitted["tyc_name"] == row["expected_tyc_name"]
    sep_arcsec = 3600 * np.hypot(
        (admitted["ra_t"] - row["expected_ra"])
        * np.cos(np.radians(row["expected_dec"])),
        admitted["dec_t"] - row["expected_dec"],
    )
    assert sep_arcsec < POSITION_TOLERANCE_ARCSEC
    assert abs(admitted["VT"] - row["expected_mag"]) < MAG_TOLERANCE


# Four x-ray sources (found while sampling obsids for the fixture above) where the
# local-index recompute legitimately disagrees with the stored production DB. Kept out
# of the generic fixture (and the generic tolerance checks above) because the
# "difference" here is understood and expected, not a bug:
#
#   - obsid 2569 x_id 4, obsid 7326 x_id 4: the matched star is present only in the
#     Tycho-2 *supplement* files (suppl_1.dat, flag='T', no proper motion), not in the
#     main catalog. Production's live VizieR query uses catalog identifier
#     "I/259/tyc2", which is the main Tycho-2 table only -- it never queries the
#     supplement tables, so these candidates were never seen by production at all.
#   - obsid 5392 x_id 2: same, but the supplement star has flag='H' (has a proper
#     motion, propagated from J1991.25).
#   - obsid 12344 x_id 1: a genuine main-catalog star (TYC 2141-972-1, HIP 98505) with a
#     large declination proper motion (pmDE=-250.7 mas/yr). At its *mean* J2000 epoch
#     (anchor epoch), it is 3.21 arcsec from the x-ray source -- just outside the 3
#     arcsec rough_match cone -- but after propagating to the observation epoch
#     (frac_year=2011.54, ~11.5 years after J2000) it is only 0.72 arcsec away. This
#     implies VizieR's live cone-search cut is applied at the catalog's native (mean)
#     epoch, before propagating to the requested output epoch, so it misses
#     high-proper-motion stars that only enter the search radius after propagation.
#     That makes a (correctly epoch-padded) local index *more* complete than the live
#     VizieR path for this class of star, not less.
KNOWN_CONE_EDGE_CASES = [
    pytest.param(
        250.35776015035617,
        -48.75398404842161,
        29.897384643554688,
        18.833053588867188,
        80.95392608642578,
        "HETG",
        (250.32557587954113, -48.76143388168115, 312.6440510420793),
        2002.4855409056318,
        "8334-1053-1",
        id="obsid2569_x4_supplement_star",
    ),
    pytest.param(
        36.49821797372745,
        -15.339781564222859,
        122.96605682373047,
        12.168434143066406,
        26.616840362548828,
        "NONE",
        (36.49903535567111, -15.332430139150654, 135.93652259257559),
        2005.4246866438357,
        "5860-2591-1",
        id="obsid5392_x2_supplement_star",
    ),
    pytest.param(
        179.9500464268888,
        -78.19586752067015,
        16.760520935058594,
        53.080116271972656,
        36.916290283203125,
        "HETG",
        (179.99659498197553, -78.19206585817356, 273.6919318276069),
        2006.4472078576864,
        "9415-2387-1",
        id="obsid7326_x4_supplement_star",
    ),
    pytest.param(
        300.18192270865745,
        22.709984979184103,
        804.7455444335938,
        11.337236404418945,
        9.23698616027832,
        "NONE",
        (300.17916358953653, 22.71030966027409, 173.7158581514334),
        2011.538559899797,
        "2141-972-1",
        id="obsid12344_x1_high_pm_cone_edge",
    ),
]


@pytest.mark.parametrize(
    "x_ra,x_dec,x_snr,x_near_neighbor_dist,x_r_angle,x_grating,quat_equatorial,"
    "frac_year,expected_tyc_name",
    KNOWN_CONE_EDGE_CASES,
)
def test_known_cone_edge_cases(
    tycho2_index,
    x_ra,
    x_dec,
    x_snr,
    x_near_neighbor_dist,
    x_r_angle,
    x_grating,
    quat_equatorial,
    frac_year,
    expected_tyc_name,
):
    stars, tree, global_max_pm = tycho2_index
    quat = Quat(equatorial=quat_equatorial)
    admitted = _recompute_match(
        stars,
        tree,
        global_max_pm,
        x_ra=x_ra,
        x_dec=x_dec,
        x_snr=x_snr,
        x_near_neighbor_dist=x_near_neighbor_dist,
        x_r_angle=x_r_angle,
        x_grating=x_grating,
        quat=quat,
        frac_year=frac_year,
    )
    assert admitted is not None
    assert admitted["tyc_name"] == expected_tyc_name
