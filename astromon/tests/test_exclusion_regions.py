"""The absolute-astrometry calibration's exclusions, as astromon_regions rows.

The list is kept in astromon/data/excluded_sources.ecsv so it is versioned with the
code; load_exclusion_regions places each entry on the sky and adds it to a database's
astromon_regions, which db.get_cross_matches(exclude_regions=True) applies.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from astromon import cross_match, db
from astromon.scripts.maintenance import load_exclusion_regions as loader
from astromon.tests.utils import minimal_cat_src, minimal_obs, minimal_xray_src


def _exclusions(rows: list[tuple[str, str, int, str]]) -> Table:
    return Table(rows=rows, names=["catalog", "name", "obsid", "reason"])


def _cat_src(rows: list[tuple[int, str, str, float, float]]) -> Table:
    return Table(rows=rows, names=["obsid", "catalog", "name", "ra", "dec"])


def test_the_list_holds_the_calibrations_exclusions():
    exclusions = loader.read_exclusions()
    pairs = {
        (str(catalog), str(name))
        for catalog, name in zip(exclusions["catalog"], exclusions["name"], strict=True)
        if catalog
    }
    obsids = {int(obsid) for obsid in exclusions["obsid"] if obsid > 0}

    assert len(pairs) == 29  # 27 catalog sources and the stars of zeta Tau, HD 97950 C
    assert ("MilliquasGaia", "ESO 253-3") in pairs
    assert ("Tycho2", "1310-2697-1") in pairs  # zeta Tau
    assert ("Tycho2", "8959-1959-1") in pairs  # HD 97950 C
    assert ("Tycho2", "2389-955-1") in pairs  # RW Aur
    assert ("Tycho2", "4774-931-1") in pairs  # Orion's Trapezium star, not Orion whole
    assert obsids == {9687, 11963, 11967}
    assert all(str(reason).strip() for reason in exclusions["reason"])


def test_a_catalog_source_is_a_circle_at_its_position_in_every_observation():
    regions, missing = loader.exclusion_regions(
        _exclusions([("Tycho2", "star-a", 0, "blended")]),
        _cat_src(
            [
                (1, "Tycho2", "star-a", 10.0, 20.0),
                (2, "Tycho2", "star-a", 10.0002, 20.0002),
                (2, "RFC", "star-a", 50.0, 50.0),  # same name, other catalog
            ]
        ),
        obs=Table({"obsid": [1, 2], "ra": [10.0, 10.0], "dec": [20.0, 20.0]}),
    )

    assert missing == []
    assert len(regions) == 1
    assert regions["ra"][0] == pytest.approx(10.0001, abs=1e-6)
    assert regions["dec"][0] == pytest.approx(20.0001, abs=1e-6)
    assert regions["radius"][0] == loader.CATALOG_SOURCE_RADIUS_ARCSEC == 3.0
    assert regions["obsid"][0] == 0
    assert regions["comments"][0].startswith("exclusion Tycho2 star-a: blended")


def test_a_catalog_source_across_ra_zero_is_centred_there():
    regions, _ = loader.exclusion_regions(
        _exclusions([("RFC", "wrap", 0, "test")]),
        _cat_src([(1, "RFC", "wrap", 359.9999, 0.0), (2, "RFC", "wrap", 0.0001, 0.0)]),
        obs=Table({"obsid": [1, 2], "ra": [0.0, 0.0], "dec": [0.0, 0.0]}),
    )
    assert min(regions["ra"][0], 360 - regions["ra"][0]) < 1e-6


def test_an_excluded_observation_is_a_circle_over_its_whole_field():
    regions, missing = loader.exclusion_regions(
        _exclusions([("", "", 9687, "wrong alignment")]),
        _cat_src([]),
        obs=Table({"obsid": [9687], "ra": [333.0], "dec": [45.7]}),
    )

    assert missing == []
    assert list(regions["obsid"]) == [9687]
    assert regions["ra"][0] == 333.0
    assert regions["dec"][0] == 45.7
    assert regions["radius"][0] == loader.WHOLE_OBSERVATION_RADIUS_ARCSEC
    assert regions["comments"][0].startswith("exclusion obsid 9687: wrong alignment")


def test_an_entry_the_database_cannot_place_is_reported_not_guessed():
    regions, missing = loader.exclusion_regions(
        _exclusions([("Tycho2", "absent", 0, "x"), ("", "", 4242, "y")]),
        _cat_src([]),
        obs=Table({"obsid": [1], "ra": [0.0], "dec": [0.0]}),
    )
    assert len(regions) == 0
    assert missing == ["Tycho2 absent", "obsid 4242"]


def _seeded_db(dbfile: Path) -> None:
    """One observation with one Tycho2 match, and no regions yet."""
    obs = Table(np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE))
    for name, value in zip(minimal_obs().colnames, minimal_obs()[0], strict=True):
        if name in obs.colnames:
            obs[name] = value
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
        minimal_cat_src(name="star-a").colnames,
        minimal_cat_src(name="star-a")[0],
        strict=True,
    ):
        if name in cat.colnames:
            cat[name] = value
    db.save("astromon_cat_src", cat, dbfile)
    db.save("astromon_xcorr", db.create_table("astromon_xcorr"), dbfile)
    db.save(
        "astromon_regions",
        Table(np.zeros(0, dtype=db.ASTROMON_REGION_DTYPE)),
        dbfile,
    )


def test_loading_excludes_the_match_and_is_idempotent(tmp_path):
    exclusions_file = tmp_path / "excluded.ecsv"
    _exclusions([("Tycho2", "star-a", 0, "blended")]).write(
        exclusions_file, format="ascii.ecsv"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seeded_db(dbfile)
        assert len(db.get_cross_matches(dbfile=dbfile, exclude_regions=True)) == 1

        first = loader.load(dbfile, exclusions_file=exclusions_file)
        second = loader.load(dbfile, exclusions_file=exclusions_file)

        regions = db.get_table("astromon_regions", dbfile)
        assert (first["added"], second["added"]) == (1, 0)
        assert len(regions) == 1
        assert len(db.get_cross_matches(dbfile=dbfile, exclude_regions=True)) == 0
        assert len(db.get_cross_matches(dbfile=dbfile)) == 1


def test_a_dry_run_writes_nothing(tmp_path):
    exclusions_file = tmp_path / "excluded.ecsv"
    _exclusions([("Tycho2", "star-a", 0, "blended")]).write(
        exclusions_file, format="ascii.ecsv"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seeded_db(dbfile)
        summary = loader.load(dbfile, exclusions_file=exclusions_file, dry_run=True)
        assert summary["added"] == 1
        assert len(db.get_table("astromon_regions", dbfile)) == 0


def test_orion_is_no_longer_excluded_wholesale():
    """Only its Trapezium star is, as a region; the field's other sources stay."""
    matches = Table(
        {
            "target": ["Orion Nebula Cluster", "RW Aur"],
            "name": ["4774-933-1", "2389-955-1"],
            "obsid": [1, 2],
            "x_id": [1, 1],
        }
    )
    assert list(cross_match.get_bad_source_mask(matches)) == [False, True]
