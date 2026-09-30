"""Tests for astromon.utils helpers that are not tied to a single pipeline stage."""

from unittest.mock import patch

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from cxotime import CxoTime
from Quaternion import Quat

from astromon import utils


def _fake_ciao_env(prefix):
    """A CIAO environment delta as Ska.Shell.getenv would return it."""
    return {"ASCDS_INSTALL": str(prefix), "PATH": "/usr/bin:/bin"}


@pytest.fixture
def clean_ciao_env_cache():
    """Isolate the module-level CIAO_ENV cache from other tests."""
    saved = dict(utils.CIAO_ENV)
    utils.CIAO_ENV.clear()
    yield utils.CIAO_ENV
    utils.CIAO_ENV.clear()
    utils.CIAO_ENV.update(saved)


def test_ciao_env_cache_not_polluted_by_workdir(tmp_path, clean_ciao_env_cache):
    """The cached CIAO environment holds no instance-specific parameter paths.

    Ciao caches the expensive `source ciao.sh` result per prefix. Storing the live
    instance dict let the per-observation ASCDS_WORK_PATH and PFILES leak into the
    cache, so a later Ciao(prefix) with no workdir inherited a param directory that
    Observation.get_ciao had already removed.
    """
    prefix = tmp_path / "ciao"
    (prefix / "param").mkdir(parents=True)
    workdir_a = tmp_path / "obs_a" / "param"

    with patch.object(utils.Ska.Shell, "getenv", return_value=_fake_ciao_env(prefix)):
        ciao_a = utils.Ciao(prefix=prefix, workdir=workdir_a, logger="astromon")

        assert ciao_a.env["ASCDS_WORK_PATH"] == str(workdir_a)
        cached = clean_ciao_env_cache[prefix]
        assert "ASCDS_WORK_PATH" not in cached
        assert "PFILES" not in cached

        # A later instance with no workdir must not inherit obs_a's param path.
        ciao_b = utils.Ciao(prefix=prefix, logger="astromon")

    assert "ASCDS_WORK_PATH" not in ciao_b.env
    assert "PFILES" not in ciao_b.env


def test_ciao_env_cache_avoids_repeat_getenv_call(tmp_path, clean_ciao_env_cache):
    """A cached prefix must not re-invoke the expensive `source ciao.sh` call.

    CIAO_ENV.get(prefix, Ska.Shell.getenv(...)) evaluates the default argument
    eagerly, so the subprocess call ran on every Ciao() construction regardless
    of whether prefix was already cached.
    """
    prefix = tmp_path / "ciao"
    (prefix / "param").mkdir(parents=True)

    with patch.object(
        utils.Ska.Shell, "getenv", return_value=_fake_ciao_env(prefix)
    ) as mock_getenv:
        utils.Ciao(prefix=prefix, logger="astromon")
        assert mock_getenv.call_count == 1

        utils.Ciao(prefix=prefix, logger="astromon")
        assert mock_getenv.call_count == 1, (
            "a cached prefix must not re-invoke Ska.Shell.getenv"
        )


def _fake_calalign_table():
    """A minimal two-version calalign table, as calalign_from_files would return."""
    aca_misalign = np.tile(np.eye(3), (2, 1, 1))
    fts_misalign = np.tile(np.eye(3), (2, 1, 1))
    dy, dz = utils.get_offsets(aca_misalign)
    return Table(
        {
            "start": CxoTime(["1999:001:00:00:00", "2010:001:00:00:00"]),
            "stop": CxoTime(["2010:001:00:00:00", "2050:001:00:00:00"]),
            "detector": ["ACIS-S", "ACIS-S"],
            "caldb_version": ["4.0.0", "4.10.0"],
            "since": CxoTime(["1999:001:00:00:00", "2010:001:00:00:00"]),
            "aca_misalign": aca_misalign,
            "fts_misalign": fts_misalign,
            "dy": dy,
            "dz": dz,
        }
    )


def _write_calalign_file(path, start, offsets_by_detector):
    """Write a minimal CALALIGN file, shaped like the CALDB's pcad/align files.

    ``offsets_by_detector`` maps each INSTR_ID to the (dy, dz) arcsec offsets
    its ACA_MISALIGN matrix encodes, as get_offsets reads them back.
    """
    detectors = list(offsets_by_detector)
    aca_misalign = np.array(
        [
            Quat(equatorial=[dy / 3600, dz / 3600, 0]).transform
            for dy, dz in offsets_by_detector.values()
        ]
    )
    table = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="INSTR_ID", format="9A", array=np.array(detectors)),
            fits.Column(
                name="ACA_MISALIGN", format="9D", dim="(3,3)", array=aca_misalign
            ),
            fits.Column(
                name="FTS_MISALIGN",
                format="9D",
                dim="(3,3)",
                array=np.tile(np.eye(3), (len(detectors), 1, 1)),
            ),
        ],
        name="CALALIGN",
    )
    table.header["CVSD0001"] = start
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(path)


def test_get_calalign_offsets_n0010_arrived_in_caldb_4_9_8(tmp_path):
    """An observation processed with CALDB 4.9.8 got the N0010 matrix for its date.

    CALDB 4.9.8's pcad index (caldbN0405.indx) is the first to list the N0010
    files, and it retires pcadD2013-01-19alignN0009.fits. With N0010 dated to
    4.10.0 instead, a 4.9.8 observation after 2013-01-19 is reconstructed with
    the N0009 matrix it was never processed with.
    """
    _write_calalign_file(
        tmp_path / "pcadD2012-09-13alignN0009.fits",
        "2012-09-13T00:00:00",
        {"ACIS-S": (10.0, -5.0)},
    )
    _write_calalign_file(
        tmp_path / "pcadD2013-01-19alignN0010.fits",
        "2013-01-19T00:00:00",
        {"ACIS-S": (20.0, 3.0)},
    )
    all_matches = Table(
        {
            "obsid": [1, 2],
            "x_id": [1, 1],
            "detector": ["ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00"] * 2),
            "caldb_version": ["4.9.8", "4.9.7"],
        }
    )

    result = utils.get_calalign_offsets(all_matches, calalign_dir=tmp_path)

    # 4.9.8 applied N0010's matrix; 4.9.7 predates N0010 and applied N0009's
    np.testing.assert_allclose(result["calalign_dy"], [20.0, 10.0])
    assert list(result["calalign_version"]) == ["4.9.8", "4.6.2"]


def test_get_calalign_offsets_row_order():
    """A shuffled input row order must come back out in that same order.

    join() sorts by the join keys internally and gives no guarantee that the
    input row order survives. get_calalign_offsets used to raise whenever the
    join's output order didn't already match all_matches' input order -- an
    unnecessary restriction on a valid, arbitrarily-ordered input table, since
    match_id_keys (obsid, x_id[, detect_method]) already uniquely identify
    each row. It must instead restore all_matches' own order rather than
    reject it.
    """
    all_matches = Table(
        {
            "obsid": [2, 1, 3],
            "x_id": [1, 1, 1],
            "detector": ["ACIS-S", "ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00"] * 3),
            # one row uses a differently-shaped version string so this column stays
            # ragged (dtype=object) rather than collapsing into a uniform 2D array,
            # matching how get_calalign_offsets expects mixed-length caldb versions.
            "caldb_version": ["4.10.0", "4.10.0.0", "4.10.0"],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table()
    ):
        result = utils.get_calalign_offsets(all_matches)

    assert list(result["obsid"]) == [2, 1, 3]


def test_get_calalign_offsets_raises_on_malformed_version_string():
    """A non-numeric version segment must raise rather than being silently dropped.

    get_calalign_offsets used to build the caldb_version/calalign_version
    comparison tuples with ``if f.isdigit()``, which drops any non-numeric
    dot-separated segment instead of raising. A truncated tuple then compares
    incorrectly (lexicographically shorter-vs-longer) against a full-length
    one, silently picking the wrong CalDB row as "actual" or "reference"
    instead of failing loudly on the unexpected input.
    """
    all_matches = Table(
        {
            "obsid": [1],
            "x_id": [1],
            "detector": ["ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00"]),
            "caldb_version": ["4.N0.0"],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table()
    ):
        with pytest.raises(ValueError):
            utils.get_calalign_offsets(all_matches)


def test_get_calalign_offsets_disambiguates_x_id_by_detect_method():
    """x_id numbering restarts per detect_method, so obsid+x_id alone can collide.

    celldetect and gaussian_detect each number their sources from 1 for a given
    obsid, so a matches table spanning both methods can have two physically
    different sources sharing (obsid, x_id). Without 'detect_method' in the
    grouping key, get_calalign_offsets collapses/mismatches those rows and used
    to raise RuntimeError("len(all_matches) != len(actual)").
    """
    all_matches = Table(
        {
            "obsid": [1, 1],
            "x_id": [1, 1],
            "detect_method": ["celldetect", "gaussian_detect"],
            "detector": ["ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00", "2015:001:00:00:00"]),
            "caldb_version": ["4.10.0", "4.10.0"],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table()
    ):
        result = utils.get_calalign_offsets(all_matches)

    assert len(result) == 2
    assert list(result["detect_method"]) == ["celldetect", "gaussian_detect"]


def test_get_calalign_offsets_without_detect_method_column():
    """Tables without a 'detect_method' column (e.g. from older DBs) still work."""
    all_matches = Table(
        {
            "obsid": [1, 2],
            "x_id": [1, 1],
            "detector": ["ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00", "2015:001:00:00:00"]),
            "caldb_version": ["4.10.0", "4.10.0"],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table()
    ):
        result = utils.get_calalign_offsets(all_matches)

    assert list(result["obsid"]) == [1, 2]
    assert "detect_method" not in result.colnames


def _fake_calalign_table_multi_matrix():
    """Three CALALIGN entries for one detector, the last two sharing a version label.

    Mirrors CALDB's real pattern: many distinct, periodically-updated alignment
    matrices (tracking periscope drift) can share one version label -- e.g. every
    file from 2013-01-19 through 2021-07-02 in a real CALDB checkout is tagged
    "N0010". get_latest_calalign_matrix must pick the row by date (start), not by
    that shared version label.
    """
    n = 3
    aca_misalign = np.tile(np.eye(3), (n, 1, 1))
    fts_misalign = np.tile(np.eye(3), (n, 1, 1))
    return Table(
        {
            "start": CxoTime(
                ["1999:001:00:00:00", "2013:001:00:00:00", "2020:001:00:00:00"]
            ),
            "stop": CxoTime(
                ["2013:001:00:00:00", "2050:001:00:00:00", "2050:001:00:00:00"]
            ),
            "detector": ["ACIS-S"] * n,
            "caldb_version": ["4.4.4", "4.10.0", "4.10.0"],
            "since": CxoTime(["1999:001:00:00:00", "2013:001:00:00:00"] * n)[:n],
            "aca_misalign": aca_misalign,
            "fts_misalign": fts_misalign,
            "dy": np.array([0.0, 1.0, 2.0]),
            "dz": np.array([0.0, -1.0, -2.0]),
        }
    )


def test_get_latest_calalign_matrix_picks_by_date_not_version():
    """The most-recently-dated matrix wins even when an older matrix shares its
    version label with a still-older one -- see _fake_calalign_table_multi_matrix.
    """
    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table_multi_matrix()
    ):
        latest = utils.get_latest_calalign_matrix()

    assert latest["ACIS-S"] == (2.0, -2.0)


def test_get_rebased_offsets():
    """dy_rebased undoes the as-processed CALALIGN and reapplies the latest matrix.

    One match processed under the middle (dy=1.0) matrix should end up shifted to
    read as if the latest (dy=2.0) matrix had been used instead: dy_rebased =
    dy - (calalign_dy - ref_dy) = 5.0 - (1.0 - 2.0) = 6.0. A second match with
    caldb_version "0.0" (no real CalDB version recorded) is preserved in the output
    but gets NaN instead of a computed value.
    """
    all_matches = Table(
        {
            "obsid": [1, 2],
            "x_id": [1, 1],
            "detector": ["ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00", "2015:001:00:00:00"]),
            "caldb_version": ["4.10.0", "0.0"],
            "dy": [5.0, 5.0],
            "dz": [5.0, 5.0],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table_multi_matrix()
    ):
        result = utils.get_rebased_offsets(all_matches)

    assert list(result["obsid"]) == [1, 2]
    assert result["dy_rebased"][0] == pytest.approx(6.0)
    assert result["dz_rebased"][0] == pytest.approx(4.0)
    assert np.isnan(result["dy_rebased"][1])
    assert np.isnan(result["dz_rebased"][1])


def test_get_rebased_offsets_preserves_row_order_with_interleaved_no_caldb():
    """The returned table's row order must match all_matches, even when
    caldb_version == "0.0" rows are interleaved with real ones.

    get_rebased_offsets used to split all_matches into a has-caldb table and a
    no-caldb table, compute dy_rebased/dz_rebased on the first, and vstack the
    two back together -- always has-caldb rows first, then no-caldb rows. Given
    obsids [2, 1, 3] with the middle one ("1") lacking a real caldb_version, the
    output row order came back as [2, 3, 1] instead of [2, 1, 3], even though
    the docstring promises "a copy of all_matches" (same row order, extra
    columns). Values stayed attached to the right obsid either way -- the bug
    was purely in row order, which is what this checks.
    """
    all_matches = Table(
        {
            "obsid": [2, 1, 3],
            "x_id": [1, 1, 1],
            "detector": ["ACIS-S", "ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00"] * 3),
            "caldb_version": ["4.10.0", "0.0", "4.10.0"],
            "dy": [5.0, 5.0, 5.0],
            "dz": [5.0, 5.0, 5.0],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table_multi_matrix()
    ):
        result = utils.get_rebased_offsets(all_matches)

    assert list(result["obsid"]) == [2, 1, 3]
    assert np.isnan(result["dy_rebased"][1])
    assert not np.isnan(result["dy_rebased"][0])
    assert not np.isnan(result["dy_rebased"][2])


def _write_fits_table(path, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.table_to_hdu(Table(columns)).writeto(path)


def _caldb_tree(root):
    """A CALDB-shaped tree: release table plus historical pcad indexes.

    4.9.3 and 4.9.6 share an index that still ships
    pcadD2013-01-19alignN0009.fits; 4.9.8's index retires it for N0010, as the
    real caldbN0404/caldbN0405 do. The index also lists a bad-quality ALIGN
    file and a non-ALIGN file, neither of which a directory needs.
    """
    caldb = root / "CALDB"
    _write_fits_table(
        caldb / "docs/chandra/caldb_version/caldb_version.fits",
        {
            "CALDB_VER": ["4.9.3", "4.9.6", "4.9.8"],
            "PCAD_INDEX": ["caldbN0404.indx", "caldbN0404.indx", "caldbN0405.indx"],
        },
    )
    _write_fits_table(
        caldb / "data/chandra/pcad/index/caldbN0404.indx",
        {
            "CAL_FILE": [
                "pcadD2012-09-13alignN0009.fits",
                "pcadD2013-01-19alignN0009.fits",
                "pcadD2003-06-09alignN0008.fits",
                "acapD2012-01-01darkN0002.fits",
            ],
            "CAL_CNAM": ["ALIGN", "ALIGN", "ALIGN", "DARK_CURR"],
            "CAL_QUAL": [0, 0, 5, 0],
        },
    )
    _write_fits_table(
        caldb / "data/chandra/pcad/index/caldbN0405.indx",
        {
            "CAL_FILE": [
                "pcadD2012-09-13alignN0009.fits",
                "pcadD2013-01-19alignN0010.fits",
            ],
            "CAL_CNAM": ["ALIGN", "ALIGN"],
            "CAL_QUAL": [0, 0],
        },
    )
    return caldb


def _calalign_dir(root, names):
    calalign = root / "calalign"
    calalign.mkdir()
    for name in names:
        (calalign / name).touch()
    return calalign


def test_missing_calalign_files_names_a_retired_file_the_directory_lacks(tmp_path):
    """A copy of today's CALDB lacks the file 4.6.2-4.9.7 processing applied."""
    caldb = _caldb_tree(tmp_path)
    calalign = _calalign_dir(
        tmp_path, ["pcadD2012-09-13alignN0009.fits", "pcadD2013-01-19alignN0010.fits"]
    )

    missing = utils.missing_calalign_files(["4.9.3", "4.9.6", "4.9.8"], calalign, caldb)

    assert missing == {
        "4.9.3": ["pcadD2013-01-19alignN0009.fits"],
        "4.9.6": ["pcadD2013-01-19alignN0009.fits"],
    }


def test_missing_calalign_files_is_empty_for_a_complete_directory(tmp_path):
    caldb = _caldb_tree(tmp_path)
    calalign = _calalign_dir(
        tmp_path,
        [
            "pcadD2012-09-13alignN0009.fits",
            "pcadD2013-01-19alignN0009.fits",
            "pcadD2013-01-19alignN0010.fits",
        ],
    )

    assert utils.missing_calalign_files(["4.9.3", "4.9.8"], calalign, caldb) == {}


def test_missing_calalign_files_ignores_a_trailing_dot(tmp_path):
    caldb = _caldb_tree(tmp_path)
    calalign = _calalign_dir(tmp_path, ["pcadD2012-09-13alignN0009.fits"])

    missing = utils.missing_calalign_files(["4.9.6."], calalign, caldb)

    assert missing == {"4.9.6": ["pcadD2013-01-19alignN0009.fits"]}


def test_missing_calalign_files_rejects_a_version_the_release_table_lacks(tmp_path):
    """A CIAO CALDB's table stops at its public release; SDP runs newer ones."""
    caldb = _caldb_tree(tmp_path)
    calalign = _calalign_dir(tmp_path, ["pcadD2012-09-13alignN0009.fits"])

    with pytest.raises(ValueError, match=r"4\.12\.6"):
        utils.missing_calalign_files(["4.9.3", "4.12.6"], calalign, caldb)


def test_get_calalign_offsets_with_caldb_dir_refuses_an_incomplete_directory(
    tmp_path,
):
    """Given caldb_dir, a directory missing a shipped file fails loudly.

    Without the retired pcadD2013-01-19alignN0009.fits, a 4.9.6 observation
    from 2015 is silently reconstructed with the 2012-09-13 N0009 matrix.
    """
    caldb = _caldb_tree(tmp_path)
    calalign = tmp_path / "calalign"
    calalign.mkdir()
    _write_calalign_file(
        calalign / "pcadD2012-09-13alignN0009.fits",
        "2012-09-13T00:00:00",
        {"ACIS-S": (10.0, -5.0)},
    )
    _write_calalign_file(
        calalign / "pcadD2013-01-19alignN0010.fits",
        "2013-01-19T00:00:00",
        {"ACIS-S": (20.0, 3.0)},
    )
    all_matches = Table(
        {
            "obsid": [1],
            "x_id": [1],
            "detector": ["ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00"]),
            "caldb_version": ["4.9.6"],
        }
    )

    with pytest.raises(ValueError, match="pcadD2013-01-19alignN0009.fits"):
        utils.get_calalign_offsets(all_matches, calalign_dir=calalign, caldb_dir=caldb)
    # without caldb_dir nothing is checked, as before
    result = utils.get_calalign_offsets(all_matches, calalign_dir=calalign)
    np.testing.assert_allclose(result["calalign_dy"], [10.0])


def test_missing_calalign_files_accepts_the_database_bytes_versions(tmp_path):
    """astromon's database stores caldb_version as bytes (S10); it must still match."""
    caldb = _caldb_tree(tmp_path)
    calalign = _calalign_dir(tmp_path, ["pcadD2012-09-13alignN0009.fits"])
    versions = np.unique(np.array([b"4.9.6", b"4.9.6"]))

    missing = utils.missing_calalign_files(versions, calalign, caldb)

    assert missing == {"4.9.6": ["pcadD2013-01-19alignN0009.fits"]}


def _write_acal_file(path, caldb_version, offsets, obsid=4686):
    """A minimal acal1 file: one row of applied matrices, CALDBVER and OBS_ID.

    ``offsets`` is the (dy, dz) arcsec its ACA_MISALIGN encodes.
    """
    dy, dz = offsets
    columns = {
        "aca_align": np.eye(3),
        "aca_misalign": Quat(equatorial=[dy / 3600, dz / 3600, 0]).transform,
        "fts_misalign": np.eye(3),
    }
    hdu = fits.BinTableHDU.from_columns(
        [
            fits.Column(name=name, format="9D", dim="(3,3)", array=matrix[np.newaxis])
            for name, matrix in columns.items()
        ]
    )
    hdu.header["CALDBVER"] = caldb_version
    hdu.header["OBS_ID"] = str(obsid)
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def test_read_acal_files_returns_the_applied_matrix_and_version(tmp_path):
    path = tmp_path / "pcadf192702697N004_acal1.fits.gz"
    _write_acal_file(path, "4.9.2.", (84.15, 51.25))

    acal = utils.read_acal_files([path])

    assert acal["caldb_version"] == "4.9.2"
    dy, dz = utils.get_offsets(np.array([acal["aca_misalign"]]))
    np.testing.assert_allclose([dy[0], dz[0]], [84.15, 51.25])


def test_read_acal_files_keeps_a_matrix_every_file_agrees_on(tmp_path):
    paths = [tmp_path / f"pcadf{n}N001_acal1.fits.gz" for n in (1, 2)]
    for path in paths:
        _write_acal_file(path, "4.9.4", (84.15, 51.25))

    acal = utils.read_acal_files(paths)

    assert acal["caldb_version"] == "4.9.4"
    assert acal["aca_misalign"] is not None


def test_read_acal_files_drops_a_matrix_the_files_disagree_on(tmp_path):
    """OBIs observed across an alignment change (obsid 108: 2000 and 2001).

    No one matrix was applied to the whole observation, so none is returned;
    the CALDB version they share still is.
    """
    paths = [tmp_path / f"pcadf{n}N005_acal1.fits.gz" for n in (1, 2)]
    _write_acal_file(paths[0], "4.9.4", (84.15, 51.25))
    _write_acal_file(paths[1], "4.9.4", (84.40, 51.10))

    acal = utils.read_acal_files(paths)

    assert acal["aca_misalign"] is None
    assert acal["caldb_version"] == "4.9.4"


def test_read_acal_files_drops_a_version_the_files_disagree_on(tmp_path):
    paths = [tmp_path / f"pcadf{n}N006_acal1.fits.gz" for n in (1, 2)]
    _write_acal_file(paths[0], "4.9.4", (84.15, 51.25))
    _write_acal_file(paths[1], "4.9.5", (84.15, 51.25))

    acal = utils.read_acal_files(paths)

    assert acal["caldb_version"] is None
    assert acal["aca_misalign"] is not None


def _two_epoch_calalign_dir(root):
    """CALALIGN files for ACIS-I and ACIS-S: a 2012 N0009 and the 2021 reference."""
    calalign = root / "calalign"
    calalign.mkdir()
    _write_calalign_file(
        calalign / "pcadD2012-09-13alignN0009.fits",
        "2012-09-13T00:00:00",
        {"ACIS-I": (10.0, -5.0), "ACIS-S": (30.0, 7.0)},
    )
    _write_calalign_file(
        calalign / "pcadD2021-07-02alignN0010.fits",
        "2021-07-02T12:00:00",
        {"ACIS-I": (12.0, -4.0), "ACIS-S": (32.0, 8.0)},
    )
    return calalign


def _acis_i_matches(acal_offsets):
    """ACIS-I matches from 2015 (CALDB 4.9.4), one per (acal_dy, acal_dz) pair."""
    n = len(acal_offsets)
    return Table(
        {
            "obsid": np.arange(1, n + 1),
            "x_id": np.ones(n, dtype=int),
            "detector": ["ACIS-I"] * n,
            "time": CxoTime(["2015:001:00:00:00"] * n),
            "caldb_version": ["4.9.4"] * n,
            "dy": np.full(n, 100.0),
            "dz": np.full(n, 50.0),
            "acal_dy": [dy for dy, _ in acal_offsets],
            "acal_dz": [dz for _, dz in acal_offsets],
        }
    )


def test_get_rebased_offsets_prefers_the_applied_matrix_acal1_recorded(tmp_path):
    """The applied matrix, where recorded, replaces the reconstruction.

    Obsid 9687 in miniature: labeled ACIS-I, but its aspect run applied the
    ACIS-S entry (30, 7), which the reconstruction cannot know; the second
    match has no acal1 record and is reconstructed from its CALDB version.
    Both are rebased onto ACIS-I's latest-dated matrix (12, -4).
    """
    calalign = _two_epoch_calalign_dir(tmp_path)
    matches = _acis_i_matches([(30.0, 7.0), (np.nan, np.nan)])

    result = utils.get_rebased_offsets(matches, calalign_dir=calalign)

    np.testing.assert_allclose(result["dy_rebased"], [100 - (30 - 12), 100 - (10 - 12)])
    np.testing.assert_allclose(result["dz_rebased"], [50 - (7 + 4), 50 - (-5 + 4)])
    assert list(result["calalign_source"]) == ["acal1", "reconstructed"]


def test_get_rebased_offsets_checks_only_the_rows_it_reconstructs(tmp_path):
    """An incomplete CALALIGN directory only matters for reconstructed rows."""
    caldb = _caldb_tree(tmp_path)
    calalign = _two_epoch_calalign_dir(tmp_path)  # lacks the retired 2013-01-19 N0009
    recorded = _acis_i_matches([(30.0, 7.0)])
    recorded["caldb_version"] = ["4.9.6"]
    reconstructed = _acis_i_matches([(np.nan, np.nan)])
    reconstructed["caldb_version"] = ["4.9.6"]

    result = utils.get_rebased_offsets(recorded, calalign_dir=calalign, caldb_dir=caldb)
    assert list(result["calalign_source"]) == ["acal1"]
    with pytest.raises(ValueError, match="pcadD2013-01-19alignN0009.fits"):
        utils.get_rebased_offsets(reconstructed, calalign_dir=calalign, caldb_dir=caldb)


def test_get_rebased_offsets_without_acal_columns_reconstructs_every_row(tmp_path):
    """Tables from databases without the acal columns behave as before."""
    calalign = _two_epoch_calalign_dir(tmp_path)
    matches = _acis_i_matches([(np.nan, np.nan)])
    matches.remove_columns(["acal_dy", "acal_dz"])

    result = utils.get_rebased_offsets(matches, calalign_dir=calalign)

    np.testing.assert_allclose(result["dy_rebased"], [100 - (10 - 12)])
    assert list(result["calalign_source"]) == ["reconstructed"]


def test_get_rebased_offsets_leaves_a_bytes_no_version_row_unrebased(tmp_path):
    """The database's caldb_version is bytes: b"0.0" still means no version."""
    calalign = _two_epoch_calalign_dir(tmp_path)
    matches = _acis_i_matches([(np.nan, np.nan), (np.nan, np.nan)])
    matches["caldb_version"] = np.array([b"4.9.4", b"0.0"])

    result = utils.get_rebased_offsets(matches, calalign_dir=calalign)

    assert np.isfinite(result["dy_rebased"][0])
    assert np.isnan(result["dy_rebased"][1])
    assert list(result["calalign_source"]) == ["reconstructed", ""]
