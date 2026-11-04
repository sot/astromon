"""The default cross-match selection, computed from the stored tables on demand.

astromon stores one named selection, astromon_21 (flight's). Everything else is
the default selection: the absolute-astrometry calibration's own cuts, computed
from astromon_obs, astromon_xray_src and astromon_cat_src whenever it is asked
for, so it can never be stale or missing for an obsid whose candidates exist.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table, vstack

from astromon import cross_match, db
from astromon.scripts import get_cat_obs_data
from astromon.tests.utils import (
    minimal_cat_src,
    minimal_obs,
    minimal_xray_src,
    two_method_xray_src,
)

CALIBRATION_CATALOGS = [
    "ICRF3",
    "RFC",
    "GaiaAGN",
    "DESIV161",
    "MilliquasGaia",
    "GaiaQSO",
    "Tycho2",
]


def _candidates(*catalogs: str, obsid: int = 99900) -> Table:
    """One candidate per catalog, all at the X-ray source's position."""
    rows = vstack(
        [
            minimal_cat_src(obsid, catalog=catalog, name=f"{catalog}-src")
            for catalog in catalogs
        ]
    )
    rows["id"] = np.arange(len(rows))
    return rows


def _default(obs: Table, xray: Table, cat: Table) -> Table:
    return cross_match.compute_cross_matches(
        astromon_obs=obs, astromon_xray_src=xray, astromon_cat_src=cat
    )


def test_only_astromon_21_and_the_default_are_defined():
    assert set(cross_match.CROSS_MATCHES_ARGS) == {"astromon_21", "default"}


def test_the_default_is_the_calibration_selection():
    args = cross_match.CROSS_MATCHES_ARGS["default"]
    assert args["catalogs"] == CALIBRATION_CATALOGS
    assert args["detect_method_filter"] == "celldetect"
    assert args["snr"] == 3
    assert args["r_angle"] == 160.0
    assert args["r_angle_grating"] == 160.0
    assert args["dr"] == 3.0
    assert args["near_neighbor_dist"] == 6.0


def test_the_default_prefers_the_earlier_catalog():
    matches = _default(
        minimal_obs(), minimal_xray_src(), _candidates("Tycho2", "GaiaQSO", "ICRF3")
    )
    assert list(matches["catalog"]) == ["ICRF3"]


def test_the_default_ignores_catalogs_outside_it():
    matches = _default(
        minimal_obs(), minimal_xray_src(), _candidates("Quaia", "GaiaVarStar", "SDSS")
    )
    assert len(matches) == 0


def test_the_default_uses_celldetect_only():
    matches = _default(minimal_obs(), two_method_xray_src(), _candidates("Tycho2"))
    assert list(matches["detect_method"]) == ["celldetect"]


@pytest.mark.parametrize(
    ("r_angle", "grating", "kept"),
    [
        (150.0, "NONE", True),
        (165.0, "NONE", False),
        (150.0, "HETG", True),
        (165.0, "HETG", False),
    ],
)
def test_the_default_off_axis_limit_is_160_arcsec_with_or_without_gratings(
    r_angle, grating, kept
):
    obs = minimal_obs()
    obs["grating"] = grating
    xray = minimal_xray_src()
    xray["r_angle"] = r_angle
    assert (len(_default(obs, xray, _candidates("Tycho2"))) == 1) is kept


def test_compute_cross_matches_without_arguments_is_the_default():
    obs, xray, cat = minimal_obs(), minimal_xray_src(), _candidates("GaiaQSO")
    implicit = cross_match.compute_cross_matches(
        astromon_obs=obs, astromon_xray_src=xray, astromon_cat_src=cat
    )
    explicit = cross_match.compute_cross_matches(
        "default", astromon_obs=obs, astromon_xray_src=xray, astromon_cat_src=cat
    )
    assert list(implicit["select_name"]) == ["default"]
    assert list(implicit["c_id"]) == list(explicit["c_id"])


def _seeded_db(dbfile: Path, *, region_at_source: bool = False) -> None:
    """A database with one matchable source and no stored cross-matches at all."""
    obs = Table(np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE))
    template = minimal_obs()
    for name in template.colnames:
        if name in obs.colnames:
            obs[name] = template[name]
    obs["detector"] = "ACIS-S"
    obs["sim_z"] = cross_match.SIM_Z["ACIS-S"]
    db.save("astromon_obs", obs, dbfile)

    xray = Table(np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    for name, value in zip(
        minimal_xray_src().colnames, minimal_xray_src()[0], strict=True
    ):
        xray[name] = value
    db.save("astromon_xray_src", xray, dbfile)

    cat = Table(np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE))
    for name, value in zip(
        _candidates("Tycho2").colnames, _candidates("Tycho2")[0], strict=True
    ):
        if name in cat.colnames:
            cat[name] = value
    db.save("astromon_cat_src", cat, dbfile)

    db.save("astromon_xcorr", db.create_table("astromon_xcorr"), dbfile)
    regions = Table(np.zeros(0, dtype=db.ASTROMON_REGION_DTYPE))
    db.save("astromon_regions", regions, dbfile)
    if region_at_source:
        db.add_regions(
            Table(
                {
                    "ra": [float(xray["ra"][0])],
                    "dec": [float(xray["dec"][0])],
                    "radius": [3.0],
                    "obsid": [0],
                    "comments": ["test region"],
                }
            ),
            dbfile,
        )


def test_get_cross_matches_without_a_name_computes_the_default():
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "default.h5"
        _seeded_db(dbfile)

        computed = db.get_cross_matches(dbfile=dbfile)

    # astromon_xcorr is empty, so this can only have come from the stored candidates
    assert len(computed) == 1
    assert list(computed["select_name"]) == ["default"]


def test_get_cross_matches_applies_filters_to_the_default():
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "default.h5"
        _seeded_db(dbfile)

        assert len(db.get_cross_matches(dbfile=dbfile, snr=5)) == 1
        assert len(db.get_cross_matches(dbfile=dbfile, snr=50)) == 0


def test_exclude_regions_reads_the_regions_of_the_given_database():
    """Not the default database's: the filter used to ignore dbfile."""
    with tempfile.TemporaryDirectory() as tmpdir:
        with_region = Path(tmpdir) / "with_region.h5"
        without_region = Path(tmpdir) / "without_region.h5"
        _seeded_db(with_region, region_at_source=True)
        _seeded_db(without_region)

        assert len(db.get_cross_matches(dbfile=with_region, exclude_regions=True)) == 0
        assert (
            len(db.get_cross_matches(dbfile=without_region, exclude_regions=True)) == 1
        )


def test_the_pipeline_stores_only_astromon_21():
    obs, xray = minimal_obs(), two_method_xray_src()
    candidates = _candidates("RFC", "Tycho2", "GaiaQSO")
    stored = get_cat_obs_data.cross_matches_to_store(obs, xray, candidates)
    assert set(stored.colnames) == set(db.ASTROMON_XCORR_DTYPE.names)
    assert set(np.asarray(stored["select_name"]).astype(str)) == {"astromon_21"}
    assert len(stored) == 2  # one per detect_method, RFC beats Tycho2
