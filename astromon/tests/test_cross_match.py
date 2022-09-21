from pathlib import Path

import numpy as np
from astropy.table import Table

from astromon import cross_match, db

DATA_DIR = Path(__file__).parent / "data"


def _get_table(name, *args, **kawargs):
    filename = DATA_DIR / f"{name}.ecsv"
    assert filename.exists()
    return Table.read(filename)


def test_db_agreement(monkeypatch):
    monkeypatch.setattr(db, "get_table", _get_table)

    matches = cross_match.compute_cross_matches("astromon_21")
    matches2 = db.get_cross_matches()
    # the cross_match method can potentially return more columns, that's ok.
    assert not [m for m in matches2.colnames if m not in matches.colnames]
    assert len(matches) > 0
    assert len(matches) == len(matches2)
    assert np.all(matches2["category"] == matches["category"])


def test_default_cross_match(monkeypatch):
    monkeypatch.setattr(db, "get_table", _get_table)

    # cross_match without arguments uses the default standard cross_match
    matches2 = cross_match.compute_cross_matches()
    matches = cross_match.compute_cross_matches("astromon_21")

    assert matches2.colnames == matches.colnames
    assert len(matches) > 0
    assert len(matches) == len(matches2)
    assert np.all(matches2["category"] == matches["category"])


def test_custom_args(monkeypatch):
    monkeypatch.setattr(db, "get_table", _get_table)

    # cross_match with custom arguments uses the "simple" algorithm
    matches = cross_match.compute_cross_matches(
        catalogs=["Tycho2"],
        snr=5,
        r_angle=120.0,
        r_angle_grating=0.0,
        near_neighbor_dist=6.0,
        dr=3.0,
    )
    assert len(matches) == 5
