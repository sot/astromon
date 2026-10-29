"""Tests for astromon_obs_gti: the time span each observation's event data cover.

The span is the first GTI START and last GTI STOP over every GTI extension of the
archive evt2 (ACIS has one per CCD). The fixtures are tiny real evt2 files: an
EVENTS table and a few GTI extensions.
"""

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from astromon import db, observation, utils

# obsid 25741's real span, which float32 would round by tens of seconds.
FIRST_START = 758669665.0724581
LAST_STOP = 758678754.2725282

# Two CCDs: the earliest START is on the second, the latest STOP on the first.
ACIS_GTIS = [
    ([758669700.5, 758675000.25], [758674000.75, LAST_STOP]),
    ([FIRST_START], [758678700.125]),
]


def write_evt2(
    path: Path,
    gtis: list[tuple[list[float], list[float]]],
    tstart: float = 758668597.12245,
    tstop: float = 758682267.0358,
) -> Path:
    """A minimal evt2: EVENTS, then one GTI extension per (starts, stops) pair.

    TSTART/TSTOP default to obsid 25741's header values, which run past its GTIs.
    """
    events = fits.BinTableHDU.from_columns(
        [fits.Column(name="TIME", format="D", array=np.array([tstart + 100.0]))],
        name="EVENTS",
    )
    events.header["TSTART"] = tstart
    events.header["TSTOP"] = tstop
    hdus = [fits.PrimaryHDU(), events]
    for ccd_id, (starts, stops) in enumerate(gtis):
        gti = fits.BinTableHDU.from_columns(
            [
                fits.Column(name="START", format="D", array=np.array(starts)),
                fits.Column(name="STOP", format="D", array=np.array(stops)),
            ],
            name="GTI",
        )
        gti.header["CCD_ID"] = ccd_id
        hdus.append(gti)
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(path)
    return path


def gti_row(obsid: int, gti_start: float, gti_stop: float) -> Table:
    row = np.zeros(1, dtype=db.ASTROMON_OBS_GTI_DTYPE)
    row["obsid"] = obsid
    row["gti_start"] = gti_start
    row["gti_stop"] = gti_stop
    return Table(row)


def create_pre_gti_db(dbfile: Path) -> None:
    """A database as it was before astromon_obs_gti existed."""
    db.create_empty_tables(
        dbfile, table_names=[name for name in db.DTYPES if name != "astromon_obs_gti"]
    )


def make_observation(tmp_path: Path, obsid: int = 25741) -> observation.Observation:
    obs = observation.Observation(
        obsid,
        workdir=tmp_path / "work",
        archive_dir=tmp_path / "archive",
        use_ciao=False,
    )
    # The files are on disk already; never reach for arc5gl or the archive.
    obs.download = lambda *args, **kwargs: None
    return obs


# ---- reading the span from an evt2 file ----


def test_read_gti_span_is_min_start_and_max_stop_over_all_gti_extensions(tmp_path):
    evt2 = write_evt2(tmp_path / "acisf25741N002_evt2.fits", ACIS_GTIS)

    gti_start, gti_stop = utils.read_gti_span(evt2)

    # Exact: the span must be the stored float64 values, not the header's
    # TSTART/TSTOP and not a rounded copy.
    assert gti_start == FIRST_START
    assert gti_stop == LAST_STOP


def test_read_gti_span_reads_a_gzipped_evt2(tmp_path):
    evt2 = write_evt2(tmp_path / "acisf25741N002_evt2.fits.gz", ACIS_GTIS)

    assert utils.read_gti_span(evt2) == (FIRST_START, LAST_STOP)


def test_read_gti_span_with_a_single_gti_extension(tmp_path):
    """HRC evt2 files have one GTI extension, with many intervals."""
    starts = [829007582.6328558, 829010000.0, 829015000.0]
    stops = [829009000.0, 829014000.0, 829021485.9900084]
    evt2 = write_evt2(tmp_path / "hrcf28864N001_evt2.fits", [(starts, stops)])

    assert utils.read_gti_span(evt2) == (829007582.6328558, 829021485.9900084)


def test_read_gti_span_raises_without_a_gti_extension(tmp_path):
    evt2 = write_evt2(tmp_path / "acisf1N001_evt2.fits", [])

    with pytest.raises(ValueError, match="no GTI"):
        utils.read_gti_span(evt2)


def test_read_gti_span_raises_when_every_gti_extension_is_empty(tmp_path):
    evt2 = write_evt2(tmp_path / "acisf1N001_evt2.fits", [([], []), ([], [])])

    with pytest.raises(ValueError, match="no GTI"):
        utils.read_gti_span(evt2)


# ---- choosing the archive evt2 ----


def test_find_archive_evt2_skips_the_filtered_copy(tmp_path):
    """primary/ also holds <obsid>_evt2_filtered.fits.gz, which sorts first."""
    paths = [
        tmp_path / "25741_evt2_filtered.fits.gz",
        tmp_path / "acisf25741N002_evt2.fits.gz",
        tmp_path / "acisf25741_000N002_fov1.fits.gz",
    ]

    assert utils.find_archive_evt2(paths) == tmp_path / "acisf25741N002_evt2.fits.gz"


def test_find_archive_evt2_is_none_when_only_the_filtered_copy_is_left(tmp_path):
    """cleanup_downloads keeps the filtered copy and deletes the archive evt2."""
    paths = [tmp_path / "25741_evt2_filtered.fits.gz"]

    assert utils.find_archive_evt2(paths) is None


def test_find_archive_evt2_refuses_to_pick_between_two_revisions(tmp_path):
    paths = [
        tmp_path / "acisf25741N001_evt2.fits.gz",
        tmp_path / "acisf25741N002_evt2.fits.gz",
    ]

    with pytest.raises(ValueError, match="N001.*N002"):
        utils.find_archive_evt2(paths)


# ---- Observation.get_gti_span ----


def test_get_gti_span_reads_the_archive_evt2_not_the_filtered_copy(tmp_path):
    obs = make_observation(tmp_path)
    primary = obs.workdir / "primary"
    write_evt2(primary / "acisf25741N002_evt2.fits.gz", ACIS_GTIS)
    # A filtered copy with a different span: reading it would give a wrong answer.
    write_evt2(
        primary / "25741_evt2_filtered.fits.gz", [([758670000.0], [758678000.0])]
    )

    span = obs.get_gti_span()

    assert span == {"gti_start": FIRST_START, "gti_stop": LAST_STOP}


def test_get_gti_span_raises_without_an_archive_evt2(tmp_path):
    obs = make_observation(tmp_path)
    write_evt2(
        obs.workdir / "primary" / "25741_evt2_filtered.fits.gz",
        [([758670000.0], [758678000.0])],
    )

    with pytest.raises(FileNotFoundError, match="archive evt2"):
        obs.get_gti_span()


# ---- the astromon_obs_gti table ----


def test_obs_gti_schema():
    dtype = db.DTYPES["astromon_obs_gti"]

    assert dtype.names == ("obsid", "gti_start", "gti_stop")
    assert dtype["obsid"] == np.int32
    assert dtype["gti_start"] == np.float64
    assert dtype["gti_stop"] == np.float64


def test_get_table_reads_an_absent_obs_gti_as_empty(tmp_path):
    """A database from before astromon_obs_gti existed still reads it."""
    dbfile = tmp_path / "pre_gti.h5"
    create_pre_gti_db(dbfile)

    obs_gti = db.get_table("astromon_obs_gti", dbfile)

    assert len(obs_gti) == 0
    assert obs_gti.dtype == db.ASTROMON_OBS_GTI_DTYPE


def test_get_table_still_raises_for_other_absent_tables(tmp_path):
    """Only astromon_obs_gti reads as empty when absent; the rest stay strict."""
    dbfile = tmp_path / "obs_only.h5"
    db.create_empty_tables(dbfile, table_names=["astromon_obs"])

    with pytest.raises(utils.MissingTableException, match="astromon_xcorr"):
        db.get_table("astromon_xcorr", dbfile)


def test_create_empty_tables_adds_obs_gti_to_a_pre_gti_db(tmp_path):
    dbfile = tmp_path / "pre_gti.h5"
    create_pre_gti_db(dbfile)
    obs = Table(np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE))
    obs["obsid"] = 25741
    db.save("astromon_obs", obs, dbfile)

    created = db.create_empty_tables(dbfile)

    assert created == ["astromon_obs_gti"]
    assert list(db.get_table("astromon_obs", dbfile)["obsid"]) == [25741]


@pytest.mark.parametrize("expect_existing", [False, True])
def test_save_creates_obs_gti_in_a_pre_gti_db(tmp_path, expect_existing):
    """The pipeline saves with expect_existing=True, and every existing DB predates
    this table: its first write must create it rather than call it corruption."""
    dbfile = tmp_path / "pre_gti.h5"
    create_pre_gti_db(dbfile)

    db.save(
        "astromon_obs_gti",
        gti_row(25741, FIRST_START, LAST_STOP),
        dbfile,
        expect_existing=expect_existing,
    )

    obs_gti = db.get_table("astromon_obs_gti", dbfile)
    assert list(obs_gti["obsid"]) == [25741]


def test_save_keeps_the_span_at_float64_precision(tmp_path):
    dbfile = tmp_path / "gti.h5"

    db.save("astromon_obs_gti", gti_row(25741, FIRST_START, LAST_STOP), dbfile)

    obs_gti = db.get_table("astromon_obs_gti", dbfile)
    assert obs_gti["gti_start"][0] == FIRST_START
    assert obs_gti["gti_stop"][0] == LAST_STOP


def test_save_replaces_obs_gti_rows_per_obsid(tmp_path):
    dbfile = tmp_path / "gti.h5"
    db.save("astromon_obs_gti", gti_row(1, 10.0, 20.0), dbfile)
    db.save("astromon_obs_gti", gti_row(2, 30.0, 40.0), dbfile)

    db.save("astromon_obs_gti", gti_row(2, 31.0, 41.0), dbfile)

    obs_gti = db.get_table("astromon_obs_gti", dbfile)
    spans = {int(row["obsid"]): (row["gti_start"], row["gti_stop"]) for row in obs_gti}
    assert spans == {1: (10.0, 20.0), 2: (31.0, 41.0)}


def test_save_expect_existing_still_catches_an_interrupted_obs_gti_swap(tmp_path):
    """A save killed between removing the table and renaming its replacement leaves
    the data under the temporary name. That is damage, not a first write."""
    dbfile = tmp_path / "gti.h5"
    db.save("astromon_obs_gti", gti_row(1, 10.0, 20.0), dbfile)
    with db.connect(dbfile, mode="r+") as h5:
        h5.rename_node("/astromon_obs_gti", "__astromon_obs_gti_tmp")

    with pytest.raises(utils.MissingTableException, match="astromon_obs_gti"):
        db.save(
            "astromon_obs_gti", gti_row(2, 30.0, 40.0), dbfile, expect_existing=True
        )


# ---- filling it during processing ----


def test_obs_gti_table_is_one_row_from_the_archive_evt2(tmp_path):
    from astromon.scripts.get_cat_obs_data import obs_gti_table

    obs = make_observation(tmp_path)
    write_evt2(obs.workdir / "primary" / "acisf25741N002_evt2.fits.gz", ACIS_GTIS)

    obs_gti = obs_gti_table(obs)

    assert obs_gti.dtype == db.ASTROMON_OBS_GTI_DTYPE
    assert list(obs_gti["obsid"]) == [25741]
    assert obs_gti["gti_start"][0] == FIRST_START
    assert obs_gti["gti_stop"][0] == LAST_STOP


def test_obs_gti_table_is_empty_when_the_span_cannot_be_read(tmp_path):
    """No archive evt2 (e.g. a cleaned work dir) leaves the row out; it does not
    fail an obsid whose other results are good."""
    from astromon.scripts.get_cat_obs_data import obs_gti_table

    obs = make_observation(tmp_path)
    write_evt2(
        obs.workdir / "primary" / "25741_evt2_filtered.fits.gz",
        [([758670000.0], [758678000.0])],
    )

    obs_gti = obs_gti_table(obs)

    assert len(obs_gti) == 0
    assert obs_gti.dtype == db.ASTROMON_OBS_GTI_DTYPE


def test_get_cat_obs_data_save_writes_obs_gti_into_a_pre_gti_db(tmp_path):
    """Results without the key (skips, failures, older callers) are fine too."""
    from astromon.scripts.get_cat_obs_data import save

    dbfile = tmp_path / "pipeline.h5"
    create_pre_gti_db(dbfile)
    obs = Table(np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE))
    obs["obsid"] = 25741
    obs["ascdsver"] = "10.14.4.1"
    results = [
        {
            "ok": True,
            "msg": "",
            "obsid": 25741,
            "astromon_obs": obs,
            "astromon_xray_src": db.create_table("astromon_xray_src"),
            "astromon_cat_src": db.create_table("astromon_cat_src"),
            "astromon_xcorr": db.create_table("astromon_xcorr"),
            "astromon_obs_gti": gti_row(25741, FIRST_START, LAST_STOP),
        },
        {
            "ok": False,
            "msg": "error: CIAO crashed",
            "obsid": 25742,
            "astromon_obs": [],
            "astromon_xray_src": [],
            "astromon_cat_src": [],
            "astromon_xcorr": [],
        },
    ]

    save(results, dbfile)

    obs_gti = db.get_table("astromon_obs_gti", dbfile)
    assert list(obs_gti["obsid"]) == [25741]
    assert obs_gti["gti_start"][0] == FIRST_START


def test_save_with_lock_writes_obs_gti_into_a_pre_gti_db(tmp_path):
    from astromon.scripts.maintenance.process_one_obsid import save_with_lock

    dbfile = tmp_path / "pipeline.h5"
    create_pre_gti_db(dbfile)

    save_with_lock(
        dbfile,
        {"astromon_obs_gti": gti_row(25741, FIRST_START, LAST_STOP)},
        obsid=25741,
    )

    assert list(db.get_table("astromon_obs_gti", dbfile)["obsid"]) == [25741]


def test_merge_scratch_into_main_carries_obs_gti_into_a_pre_gti_main_db(tmp_path):
    from astromon.scripts.maintenance import run_all

    main_db = tmp_path / "main.h5"
    scratch_db = tmp_path / "scratch.h5"
    create_pre_gti_db(main_db)
    db.create_empty_tables(scratch_db)
    db.save("astromon_obs_gti", gti_row(25741, FIRST_START, LAST_STOP), scratch_db)

    run_all.merge_scratch_into_main(scratch_db, main_db)

    obs_gti = db.get_table("astromon_obs_gti", main_db)
    assert list(obs_gti["obsid"]) == [25741]
    assert obs_gti["gti_stop"][0] == LAST_STOP
