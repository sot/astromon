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
    _touch(obs.workdir / "primary" / "acisf08007N002_evt2.fits.gz")
    sent = _stub_arc5gl(monkeypatch)

    obs._download_arc5gl(["evt2"])

    assert sent == []


def test_other_file_types_still_skip_when_present(tmp_path, monkeypatch):
    obs = _obs(tmp_path, monkeypatch)
    _touch(obs.workdir / "primary" / "acisf08007_000N002_fov1.fits.gz")
    sent = _stub_arc5gl(monkeypatch)

    obs._download_arc5gl(["fov"])

    assert sent == []
