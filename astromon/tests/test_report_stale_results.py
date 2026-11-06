"""Tests for the report of task and cached results a version bump makes stale."""

import json
from pathlib import Path

from astropy.table import Table

from astromon.scripts.maintenance import report_stale_results as report
from astromon.stored_result import FORMATS

#: images -> detect -> fit, like make_images -> celldetect -> gaussian_detect
DEPENDENTS = {"images": {"detect", "fit"}, "detect": {"fit"}, "fit": set()}


def obsid_dir(root: Path, obsid: int) -> Path:
    return root / f"obs{obsid // 1000:02d}" / str(obsid)


def write_task_result(root: Path, obsid: int, name: str, version: int | None) -> None:
    result = {
        "return_code": "OK",
        "start": "2026-01-01T00:00:00",
        "stop": "2026-01-01T00:00:01",
    }
    if version is not None:
        result["version"] = version
    path = obsid_dir(root, obsid) / "cache" / f"task_{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result))


def write_cached(root: Path, obsid: int, name: str, fmt: str, version: int) -> None:
    value = Table({"x": [1.0]}) if fmt == "table" else {"x": 1}
    filename = FORMATS[fmt].sanitize_filename(f"{name}__0123abcd")
    FORMATS[fmt].save(
        value, obsid_dir(root, obsid) / "cache" / filename, force=True, version=version
    )


def write_filtered_evt2(root: Path, obsid: int) -> None:
    path = obsid_dir(root, obsid) / "primary" / f"{obsid}_evt2_filtered.fits.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")


def _seed(tmp_path: Path) -> tuple[Path, Path]:
    """7001: all at version 0, with inputs. 7002: detect at version 1, no filtered evt2.

    7003 has only an archive copy, made at version 0.
    """
    work, archive = tmp_path / "work", tmp_path / "archive"
    for name in ("images", "detect", "fit"):
        write_task_result(work, 7001, name, None)
        write_task_result(work, 7002, name, 1 if name == "detect" else 0)
        write_task_result(archive, 7003, name, 0)
    write_cached(work, 7001, "sources", "table", 0)
    write_cached(work, 7002, "sources", "table", 2)
    write_cached(work, 7001, "evt2_info", "json", 0)
    write_filtered_evt2(work, 7001)
    return work, archive


def test_stored_versions_are_read_from_each_format(tmp_path):
    work, _ = _seed(tmp_path)
    cache = obsid_dir(work, 7002) / "cache"

    assert report.stored_task_version(cache / "task_detect.json") == 1
    assert (
        report.stored_task_version(obsid_dir(work, 7001) / "cache" / "task_detect.json")
        == 0
    )
    assert report.stored_cache_version(cache / "sources__0123abcd.fits.gz") == 2
    assert (
        report.stored_cache_version(
            obsid_dir(work, 7001) / "cache" / "evt2_info__0123abcd.json"
        )
        == 0
    )


def test_results_from_another_version_are_listed(tmp_path):
    work, archive = _seed(tmp_path)

    stale = report.stale_results(
        [work, archive],
        task_versions={"images": 0, "detect": 1, "fit": 0},
        cache_versions={"sources": 0, "evt2_info": 0},
        dependents=DEPENDENTS,
    )

    # 7001 and 7003 ran detect at version 0, the code is at 1; 7002's sources table was
    # cached at version 2, the code is at 0
    assert sorted((r["obsid"], r["kind"], r["name"]) for r in stale) == [
        (7001, "task", "detect"),
        (7001, "task", "fit"),
        (7002, "cache", "sources"),
        (7003, "task", "detect"),
        (7003, "task", "fit"),
    ]


def test_a_bump_marks_the_task_and_everything_downstream(tmp_path):
    work, archive = _seed(tmp_path)

    stale = report.stale_results(
        [work, archive],
        task_versions={"images": 1, "detect": 0, "fit": 0},
        cache_versions={"sources": 2, "evt2_info": 0},
        dependents=DEPENDENTS,
    )

    by_obsid = {}
    for row in stale:
        by_obsid.setdefault(row["obsid"], set()).add((row["name"], row["reason"]))
    assert by_obsid[7001] == {
        ("images", "version 0, code 1"),
        ("detect", "downstream of images"),
        ("fit", "downstream of images"),
        ("sources", "version 0, code 2"),
    }
    # 7002's detect is stale in its own right (run at 1, code at 0): that reason wins
    assert by_obsid[7002] == {
        ("images", "version 0, code 1"),
        ("detect", "version 1, code 0"),
        ("fit", "downstream of detect, images"),
    }


def test_rows_say_whether_a_rerun_has_its_inputs(tmp_path):
    work, archive = _seed(tmp_path)

    stale = report.stale_results(
        [work, archive],
        task_versions={"images": 1, "detect": 0, "fit": 0},
        cache_versions={},
        dependents=DEPENDENTS,
    )

    has_inputs = {row["obsid"]: row["filtered_evt2"] for row in stale}
    assert has_inputs == {7001: True, 7002: False, 7003: False}


def test_bumped_versions_add_one_to_the_named_results():
    versions = report.bumped(
        {"images": 0, "detect": 3}, {"sources": 0}, bump=["detect", "sources"]
    )

    assert versions == ({"images": 0, "detect": 4}, {"sources": 1})


def test_bumping_an_unknown_name_is_an_error():
    import pytest  # noqa: PLC0415

    with pytest.raises(ValueError, match="nothing named"):
        report.bumped({"images": 0}, {"sources": 0}, bump=["typo"])


def test_an_unreadable_result_is_reported_as_stale_not_fatal(tmp_path):
    """The pipeline reruns a task whose result file it cannot read; so does the report."""
    work, archive = _seed(tmp_path)
    (obsid_dir(work, 7001) / "cache" / "task_detect.json").write_bytes(bytes(100))
    (obsid_dir(work, 7001) / "cache" / "sources__0123abcd.fits.gz").write_bytes(
        b"<html>"
    )

    stale = report.stale_results(
        [work, archive],
        task_versions={"images": 0, "detect": 0, "fit": 0},
        cache_versions={"sources": 0, "evt2_info": 0},
        dependents=DEPENDENTS,
    )

    rows = {(r["obsid"], r["name"]): r["reason"] for r in stale}
    assert rows[(7001, "detect")].startswith("unreadable")
    assert rows[(7001, "fit")] == "downstream of detect"
    assert rows[(7001, "sources")].startswith("unreadable")
