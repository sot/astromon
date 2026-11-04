"""drop_retired_rows removes what astromon no longer computes or queries."""

import numpy as np
from astropy.table import Table, vstack

from astromon import db
from astromon.scripts.maintenance import drop_retired_rows


def _cat_src(rows: list[tuple[int, int, str]]) -> Table:
    table = Table(np.zeros(len(rows), dtype=db.ASTROMON_CAT_SRC_DTYPE))
    for index, (obsid, cat_id, catalog) in enumerate(rows):
        table["obsid"][index] = obsid
        table["id"][index] = cat_id
        table["catalog"][index] = catalog
        table["name"][index] = f"{catalog}-{cat_id}"
    return table


def _xcorr(rows: list[tuple[int, int, str]]) -> Table:
    table = Table(np.zeros(len(rows), dtype=db.ASTROMON_XCORR_DTYPE))
    for index, (obsid, c_id, select_name) in enumerate(rows):
        table["obsid"][index] = obsid
        table["c_id"][index] = c_id
        table["x_id"][index] = 1
        table["select_name"][index] = select_name
        table["detect_method"][index] = "celldetect"
    return table


def _seed(dbfile) -> None:
    db.create_empty_tables(dbfile)
    db.save(
        "astromon_cat_src",
        vstack(
            [
                _cat_src([(1, 0, "RFC"), (1, 1, "Quaia"), (1, 2, "Tycho2")]),
                _cat_src([(2, 0, "GaiaVarStar"), (2, 1, "GaiaQSO")]),
            ]
        ),
        dbfile,
        ignore_obsid=True,
        expect_existing=True,
    )
    db.save(
        "astromon_xcorr",
        _xcorr(
            [
                (1, 0, "astromon_21"),
                (1, 1, "astromon_26"),
                (1, 1, "quaia"),
                (2, 0, "gaia_var_star"),
                (2, 1, "astromon_24"),
            ]
        ),
        dbfile,
        ignore_obsid=True,
        expect_existing=True,
    )


def test_drops_retired_catalogs_and_selections(tmp_path):
    dbfile = tmp_path / "retire.h5"
    _seed(dbfile)

    summary = drop_retired_rows.drop_retired(dbfile)

    cat_src = db.get_table("astromon_cat_src", dbfile)
    xcorr = db.get_table("astromon_xcorr", dbfile)
    assert sorted(np.asarray(cat_src["catalog"]).astype(str)) == [
        "GaiaQSO",
        "RFC",
        "Tycho2",
    ]
    # ids are kept, so stored astromon_21 rows still point at the same sources
    assert sorted(
        zip(cat_src["obsid"].tolist(), cat_src["id"].tolist(), strict=True)
    ) == [
        (1, 0),
        (1, 2),
        (2, 1),
    ]
    assert list(np.asarray(xcorr["select_name"]).astype(str)) == ["astromon_21"]
    assert summary == {
        "cat_src_dropped": {"GaiaVarStar": 1, "Quaia": 1},
        "xcorr_dropped": {
            "astromon_24": 1,
            "astromon_26": 1,
            "gaia_var_star": 1,
            "quaia": 1,
        },
    }


def test_a_dry_run_writes_nothing(tmp_path):
    dbfile = tmp_path / "retire.h5"
    _seed(dbfile)

    summary = drop_retired_rows.drop_retired(dbfile, dry_run=True)

    assert summary["cat_src_dropped"] == {"GaiaVarStar": 1, "Quaia": 1}
    assert len(db.get_table("astromon_cat_src", dbfile)) == 5
    assert len(db.get_table("astromon_xcorr", dbfile)) == 5


def test_refuses_to_leave_a_stored_match_pointing_at_a_dropped_source(tmp_path):
    dbfile = tmp_path / "retire.h5"
    _seed(dbfile)
    db.save(
        "astromon_xcorr",
        _xcorr([(1, 1, "astromon_21")]),  # an astromon_21 row on the Quaia source
        dbfile,
        ignore_obsid=True,
        expect_existing=True,
    )

    try:
        drop_retired_rows.drop_retired(dbfile)
    except RuntimeError as error:
        assert "astromon_21" in str(error)
    else:
        raise AssertionError("expected a refusal")
    assert len(db.get_table("astromon_cat_src", dbfile)) == 5
