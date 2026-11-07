"""Tests for the arc5gl download path's "already downloaded" check.

arc5gl itself and the ocat lookup are the external boundaries and are replaced;
the files in the work directory are real.
"""

from pathlib import Path

from astromon import observation


def _obs(tmp_path: Path, monkeypatch, obsid: int = 8007) -> observation.Observation:
    """An arc5gl Observation with a private workdir, as an ACIS observation."""
    obs = observation.Observation(
        obsid,
        workdir=tmp_path / "work",
        archive_dir=tmp_path / "archive",
        source="arc5gl",
        use_ciao=False,
    )
    monkeypatch.setattr(obs, "_get_sequence_summary", lambda: {"instr": "ACIS-S"})
    return obs


def _stub_arc5gl(monkeypatch) -> list[str]:
    """Replace Ska.arc5gl.Arc5gl; return the list of lines sent to it."""
    sent = []

    class FakeArc5gl:
        def sendline(self, line: str) -> None:
            sent.append(line)

    monkeypatch.setattr(observation.Ska.arc5gl, "Arc5gl", FakeArc5gl)
    return sent


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")


def _write_fits(path: Path) -> None:
    """A complete FITS file: what an earlier, good download left."""
    from astropy.io import fits  # noqa: PLC0415

    path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU().writeto(path, overwrite=True)


def _write_truncated_fits_gz(path: Path) -> None:
    """A .fits.gz cut off halfway, as a download that stopped mid-file leaves it."""
    import numpy as np  # noqa: PLC0415
    from astropy.io import fits  # noqa: PLC0415
    from astropy.table import Table  # noqa: PLC0415

    path.parent.mkdir(parents=True, exist_ok=True)
    values = np.random.default_rng(0).random(50_000)
    fits.BinTableHDU(Table({"time": values})).writeto(path, overwrite=True)
    complete = path.read_bytes()
    path.write_bytes(complete[: len(complete) // 2])


def test_evt2_is_downloaded_when_only_the_filtered_copy_is_left(tmp_path, monkeypatch):
    """cleanup_downloads keeps <obsid>_evt2_filtered.fits.gz and deletes the archive
    evt2. The filtered copy matches primary/*evt2*, so it used to count as the evt2
    being downloaded, and anything needing the archive file then failed."""
    obs = _obs(tmp_path, monkeypatch)
    _touch(obs.workdir / "primary" / "8007_evt2_filtered.fits.gz")
    sent = _stub_arc5gl(monkeypatch)

    obs._download_arc5gl(["evt2"])

    assert sent == ["obsid=8007", "get acis2{evt2}"]


def test_evt2_is_not_downloaded_again_when_the_archive_file_is_there(
    tmp_path, monkeypatch
):
    obs = _obs(tmp_path, monkeypatch)
    _touch(obs.workdir / "primary" / "8007_evt2_filtered.fits.gz")
    _write_fits(obs.workdir / "primary" / "acisf08007N002_evt2.fits.gz")
    sent = _stub_arc5gl(monkeypatch)

    obs._download_arc5gl(["evt2"])

    assert sent == []


def test_other_file_types_still_skip_when_present(tmp_path, monkeypatch):
    obs = _obs(tmp_path, monkeypatch)
    _write_fits(obs.workdir / "primary" / "acisf08007_000N002_fov1.fits.gz")
    sent = _stub_arc5gl(monkeypatch)

    obs._download_arc5gl(["fov"])

    assert sent == []


def test_a_file_cut_short_in_the_work_dir_is_fetched_again_by_arc5gl(
    tmp_path, monkeypatch
):
    """A present file only means the type is downloaded if it is whole; arc5gl skips a
    type while any of its files is there, so all of that type's files are moved aside."""
    obs = _obs(tmp_path, monkeypatch)
    secondary = obs.workdir / "secondary"
    _write_truncated_fits_gz(secondary / "pcadf08007_000N001_asol1.fits.gz")
    _write_fits(secondary / "pcadf08007_001N001_asol1.fits.gz")
    gets = []

    class FakeArc5gl:
        def sendline(self, line: str) -> None:
            if line.startswith("get"):
                gets.append(line)
                _write_fits(secondary / "pcadf08007_000N001_asol1.fits.gz")
                _write_fits(secondary / "pcadf08007_001N001_asol1.fits.gz")

    monkeypatch.setattr(observation.Ska.arc5gl, "Arc5gl", FakeArc5gl)

    obs._download_arc5gl(["asol"])

    assert len(gets) == 1
    rejected = obs.workdir / "rejected_downloads" / "secondary"
    assert sorted(path.name for path in rejected.iterdir()) == [
        "pcadf08007_000N001_asol1.fits.gz",
        "pcadf08007_001N001_asol1.fits.gz",
    ]
    assert observation.utils.is_fits_file(
        secondary / "pcadf08007_000N001_asol1.fits.gz"
    )


def test_an_arc5gl_download_that_is_not_fits_is_moved_aside(tmp_path, monkeypatch):
    """arc5gl saves whatever the server sends; an error page must not pass as FITS."""
    import pytest  # noqa: PLC0415

    obs = _obs(tmp_path, monkeypatch)
    primary = obs.workdir / "primary"

    class FakeArc5gl:
        def sendline(self, line: str) -> None:
            if line.startswith("get"):
                primary.mkdir(parents=True, exist_ok=True)
                (primary / "acisf08007_000N002_fov1.fits.gz").write_bytes(
                    b"<html>503</html>"
                )

    monkeypatch.setattr(observation.Ska.arc5gl, "Arc5gl", FakeArc5gl)

    with pytest.raises(RuntimeError, match="not FITS"):
        obs._download_arc5gl(["fov"])

    assert (
        obs.workdir
        / "rejected_downloads"
        / "primary"
        / "acisf08007_000N002_fov1.fits.gz"
    ).exists()
    assert not (primary / "acisf08007_000N002_fov1.fits.gz").exists()


def test_a_retry_after_a_partly_bad_arc5gl_get_fetches_again(tmp_path, monkeypatch):
    """One bad asol segment among good ones: moving only it aside left the good ones,
    which made the retry skip the download and lose the bad segment for good."""
    import pytest  # noqa: PLC0415
    from astropy.io import fits  # noqa: PLC0415

    obs = _obs(tmp_path, monkeypatch)
    secondary = obs.workdir / "secondary"
    gets = []

    class FakeArc5gl:
        def sendline(self, line: str) -> None:
            if line.startswith("get"):
                gets.append(line)
                secondary.mkdir(parents=True, exist_ok=True)
                fits.PrimaryHDU().writeto(
                    secondary / "pcadf08007_000N001_asol1.fits", overwrite=True
                )
                second = secondary / "pcadf08007_001N001_asol1.fits"
                if len(gets) == 1:
                    second.write_bytes(b"<html>503</html>")
                else:
                    fits.PrimaryHDU().writeto(second, overwrite=True)

    monkeypatch.setattr(observation.Ska.arc5gl, "Arc5gl", FakeArc5gl)

    with pytest.raises(RuntimeError, match="not FITS"):
        obs._download_arc5gl(["asol"])
    obs._download_arc5gl(["asol"])

    assert len(gets) == 2
    assert (secondary / "pcadf08007_001N001_asol1.fits").exists()


def test_a_broken_symlink_in_the_destination_does_not_stop_a_download(
    tmp_path, monkeypatch
):
    import os  # noqa: PLC0415

    obs = _obs(tmp_path, monkeypatch)
    secondary = obs.workdir / "secondary"
    secondary.mkdir(parents=True)
    os.symlink(tmp_path / "gone" / "old_asol1.fits", secondary / "stale_link.fits")
    sent = _stub_arc5gl(monkeypatch)

    obs._download_arc5gl(["msk"])

    assert sent[-1].startswith("get")
