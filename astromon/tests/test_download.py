"""Tests for the public-archive download path (``source="archive"``).

This path exists so the pipeline can run without arc5gl or CXC-internal
products, which means it is exactly the path that cannot be covered by a
HEAD-network integration test.  Everything here runs offline: the release-date
check is pure logic over an ocat row, the layout fix-ups are filesystem work,
and download_chandra_obsid itself is replaced by a fake subprocess.
"""

import gzip
import subprocess
from datetime import datetime, timedelta
from types import SimpleNamespace

import numpy as np
import pytest

from astromon import observation
from astromon.observation import ObsidNotPubliclyAvailable


def _obs(tmp_path, obsid=8007, source="archive"):
    """An Observation with a private workdir and no CIAO requirement."""
    return observation.Observation(
        obsid,
        workdir=tmp_path / "work",
        archive_dir=tmp_path / "archive",
        source=source,
        use_ciao=False,
    )


def _stub_ciao(monkeypatch, env=None):
    """Give Observation.ciao a stand-in so _download_archive can read .env.

    Observation.ciao is ``property(get_ciao)``, which captured the original
    function at class-creation time -- patching get_ciao does not affect it, so
    the property itself has to be replaced.
    """
    stand_in = SimpleNamespace(env=env if env is not None else {"PATH": "/usr/bin"})
    monkeypatch.setattr(
        observation.Observation, "ciao", property(lambda self: stand_in)
    )


def _stub_ocat(monkeypatch, public_avail):
    """Make the release-date check see a given public_avail value."""
    monkeypatch.setattr(
        observation.cda, "get_ocat_local", lambda obsid: {"public_avail": public_avail}
    )


class FakeProcess:
    """Stand-in for the download_chandra_obsid subprocess."""

    def __init__(self, stdout=b"", returncode=0, timeout=False, on_communicate=None):
        self._stdout = stdout
        self.returncode = returncode
        self._timeout = timeout
        self._on_communicate = on_communicate
        self.killed = False
        self.communicate_calls = 0

    def communicate(self, timeout=None):
        self.communicate_calls += 1
        if self._timeout and self.communicate_calls == 1:
            raise subprocess.TimeoutExpired(
                cmd="download_chandra_obsid", timeout=timeout
            )
        if self._on_communicate is not None:
            self._on_communicate()
        return self._stdout, b""

    def kill(self):
        self.killed = True


def _stub_popen(monkeypatch, proc):
    """Replace subprocess.Popen and record the argv it was called with."""
    seen = {}

    def fake_popen(cmd, **kwargs):
        seen["cmd"] = cmd
        seen["kwargs"] = kwargs
        return proc

    monkeypatch.setattr(observation.subprocess, "Popen", fake_popen)
    return seen


# --- release-date gate ------------------------------------------------------
#
# ObsidNotPubliclyAvailable is what run_all turns into a permanent "skip"
# rather than a retryable failure, so each of these branches decides whether an
# obsid is dropped for good or retried forever.


def test_no_release_date_raises_not_publicly_available(tmp_path, monkeypatch):
    """public_avail='0' means embargoed or permanently proprietary."""
    _stub_ocat(monkeypatch, "0")
    obs = _obs(tmp_path)
    with pytest.raises(ObsidNotPubliclyAvailable, match="no public release date"):
        obs._download_archive(["evt2"])


def test_masked_release_date_raises_not_publicly_available(tmp_path, monkeypatch):
    """get_ocat_web returns np.ma.masked for the same case get_ocat_local calls '0'."""
    _stub_ocat(monkeypatch, np.ma.masked)
    obs = _obs(tmp_path)
    with pytest.raises(ObsidNotPubliclyAvailable, match="no public release date"):
        obs._download_archive(["evt2"])


def test_future_release_date_raises_not_publicly_available(tmp_path, monkeypatch):
    """An obsid whose embargo has not expired is a skip, not a failure."""
    future = (datetime.now() + timedelta(days=30)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, future)
    obs = _obs(tmp_path)
    with pytest.raises(ObsidNotPubliclyAvailable, match="not yet publicly available"):
        obs._download_archive(["evt2"])


def test_past_release_date_proceeds_to_download(tmp_path, monkeypatch):
    """A released obsid gets as far as invoking download_chandra_obsid."""
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    obs = _obs(tmp_path)
    seen = _stub_popen(monkeypatch, FakeProcess(stdout=b"done\n", returncode=0))
    obs._download_archive(["evt2"])
    assert seen["cmd"][0] == "download_chandra_obsid"
    assert seen["cmd"][1] == "8007"
    assert seen["cmd"][2] == observation.Observation._ARCHIVE_FILETYPES


# --- mirror ------------------------------------------------------------------
#
# download_chandra_obsid does not fall back to the primary CDA site if the
# mirror does not have the obsid -- it is skipped, same as a genuine "not
# found". So mirror=None (the default) has to mean "use the primary CDA site",
# not "use some hardcoded mirror" -- these confirm --mirror is only added to
# the command when actually requested.


def test_mirror_is_passed_to_download_chandra_obsid_when_set(tmp_path, monkeypatch):
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    obs = observation.Observation(
        8007,
        workdir=tmp_path / "work",
        archive_dir=tmp_path / "archive",
        source="archive",
        use_ciao=False,
        mirror="https://heasarc.gsfc.nasa.gov/FTP/chandra/data/",
    )
    seen = _stub_popen(monkeypatch, FakeProcess(stdout=b"done\n", returncode=0))
    obs._download_archive(["evt2"])
    assert seen["cmd"][0] == "download_chandra_obsid"
    assert seen["cmd"][1] == "--mirror"
    assert seen["cmd"][2] == "https://heasarc.gsfc.nasa.gov/FTP/chandra/data/"
    assert seen["cmd"][3] == "8007"


def test_no_mirror_by_default(tmp_path, monkeypatch):
    """The default (mirror=None) uses the primary CDA site, not some fixed mirror."""
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    obs = _obs(tmp_path)
    seen = _stub_popen(monkeypatch, FakeProcess(stdout=b"done\n", returncode=0))
    obs._download_archive(["evt2"])
    assert "--mirror" not in seen["cmd"]
    assert seen["cmd"][0] == "download_chandra_obsid"
    assert seen["cmd"][1] == "8007"


def test_ocat_lookup_failure_is_a_warning_not_a_skip(tmp_path, monkeypatch):
    """If the ocat cannot be consulted at all, the download still goes ahead.

    An unreachable ocat says nothing about whether the obsid is public, so
    treating it as ObsidNotPubliclyAvailable would permanently drop obsids over
    a transient problem.
    """

    def boom(*args, **kwargs):
        raise OSError("ocat unavailable")

    monkeypatch.setattr(observation.cda, "get_ocat_local", boom)
    monkeypatch.setattr(observation.cda, "get_ocat_web", boom)
    _stub_ciao(monkeypatch)
    obs = _obs(tmp_path)
    seen = _stub_popen(monkeypatch, FakeProcess(stdout=b"ok\n", returncode=0))
    obs._download_archive(["evt2"])
    assert seen["cmd"][0] == "download_chandra_obsid"


# --- download_chandra_obsid outcomes ---------------------------------------


def test_not_found_on_archive_site_raises_retryable_runtime_error(
    tmp_path, monkeypatch
):
    """download_chandra_obsid exits 0 when the obsid is missing, so parse stdout.

    This must NOT be ObsidNotPubliclyAvailable: the same message appears on
    transient network failures, and only the ocat check can declare an obsid
    permanently unavailable.
    """
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    obs = _obs(tmp_path)
    _stub_popen(
        monkeypatch,
        FakeProcess(
            stdout=b"Obsid 8007 was not found on the archive site\n", returncode=0
        ),
    )
    with pytest.raises(RuntimeError, match="not found on the CDA archive site") as exc:
        obs._download_archive(["evt2"])
    assert not isinstance(exc.value, ObsidNotPubliclyAvailable)


def test_nonzero_exit_code_raises(tmp_path, monkeypatch):
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    obs = _obs(tmp_path)
    _stub_popen(monkeypatch, FakeProcess(stdout=b"partial\n", returncode=3))
    with pytest.raises(Exception, match="exit code 3"):
        obs._download_archive(["evt2"])


def test_timeout_kills_the_process_and_raises(tmp_path, monkeypatch):
    """A hung connection must not stall the run forever, and must be reaped."""
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    obs = _obs(tmp_path)
    proc = FakeProcess(timeout=True)
    _stub_popen(monkeypatch, proc)
    with pytest.raises(RuntimeError, match="timed out"):
        obs._download_archive(["evt2"])
    assert proc.killed, "the hung process must be killed"
    assert proc.communicate_calls == 2, "must reap the killed process"


def _refuse_to_run(*args, **kwargs):
    raise AssertionError("the download must not be attempted")


def test_a_completed_download_short_circuits_a_second_call(tmp_path, monkeypatch):
    """A marker written after success is what says "already downloaded"."""
    obs = _obs(tmp_path)
    (obs.workdir / "secondary").mkdir(parents=True)
    obs._archive_download_marker().write_text("done")

    monkeypatch.setattr(observation.subprocess, "Popen", _refuse_to_run)
    monkeypatch.setattr(observation.cda, "get_ocat_local", _refuse_to_run)
    obs._download_archive(["evt2"])


def test_a_partial_download_is_retried_rather_than_trusted(tmp_path, monkeypatch):
    """secondary/ existing used to mean "done", which made a failure sticky.

    download_chandra_obsid creates secondary/ early, and the timeout and error
    paths do not remove it, so an interrupted download left a directory that every
    later attempt read as success -- and the run proceeded with missing files until
    someone deleted it by hand.
    """
    obs = _obs(tmp_path)
    (obs.workdir / "secondary").mkdir(parents=True)
    assert not obs._archive_download_marker().exists()

    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    attempted = []
    _stub_popen(
        monkeypatch,
        FakeProcess(on_communicate=lambda: attempted.append(True)),
    )
    obs._download_archive(["evt2"])

    assert attempted, "a download with no completion marker must be retried"
    assert obs._archive_download_marker().exists()


def test_a_failed_download_leaves_no_completion_marker(tmp_path, monkeypatch):
    """So the next attempt retries instead of inheriting the failure."""
    obs = _obs(tmp_path)
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    _stub_popen(monkeypatch, FakeProcess(returncode=1))

    with pytest.raises(Exception, match="exit code"):
        obs._download_archive(["evt2"])

    assert not obs._archive_download_marker().exists()


def test_asol_is_linked_into_secondary_after_download(tmp_path, monkeypatch):
    """CDA puts asol in primary/; the rest of the pipeline looks in secondary/."""
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    obs = _obs(tmp_path)

    def create_layout():
        (obs.workdir / "primary").mkdir(parents=True, exist_ok=True)
        (obs.workdir / "secondary").mkdir(parents=True, exist_ok=True)
        (obs.workdir / "primary" / "acisf08007_000N001_asol1.fits").write_text("asol")

    _stub_popen(monkeypatch, FakeProcess(stdout=b"ok\n", on_communicate=create_layout))
    obs._download_archive(["evt2"])
    link = obs.workdir / "secondary" / "acisf08007_000N001_asol1.fits"
    assert link.is_symlink(), "asol should be linked into secondary/"
    assert link.resolve().read_text() == "asol"


# --- _fix_archive_bpix -----------------------------------------------------


def _write_gz(path, payload=b"BPIX"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as fh:
        fh.write(payload)


def test_bpix_gz_is_decompressed_into_primary_as_a_real_file(tmp_path):
    """fluximage cannot follow a symlink to the bpix file, so it must be real."""
    obs = _obs(tmp_path)
    _write_gz(obs.workdir / "secondary" / "acisf08007_000N001_bpix1.fits.gz")
    (obs.workdir / "primary").mkdir(parents=True, exist_ok=True)

    obs._fix_archive_bpix()

    dest = obs.workdir / "primary" / "acisf08007_000N001_bpix1.fits"
    assert dest.is_file() and not dest.is_symlink()
    assert dest.read_bytes() == b"BPIX"


def test_bpix_stale_gz_symlink_is_removed(tmp_path):
    """A leftover .gz symlink would make the primary/*_bpix1.fits* glob ambiguous."""
    obs = _obs(tmp_path)
    src = obs.workdir / "secondary" / "acisf08007_000N001_bpix1.fits.gz"
    _write_gz(src)
    primary = obs.workdir / "primary"
    primary.mkdir(parents=True, exist_ok=True)
    stale = primary / "acisf08007_000N001_bpix1.fits.gz"
    stale.symlink_to(src)

    obs._fix_archive_bpix()

    assert not stale.exists() and not stale.is_symlink()
    assert len(list(primary.glob("*_bpix1.fits*"))) == 1, "glob must be unambiguous"


def test_fix_archive_bpix_is_idempotent(tmp_path):
    """Its docstring promises it is safe to call on every run."""
    obs = _obs(tmp_path)
    _write_gz(obs.workdir / "secondary" / "acisf08007_000N001_bpix1.fits.gz")
    (obs.workdir / "primary").mkdir(parents=True, exist_ok=True)

    obs._fix_archive_bpix()
    dest = obs.workdir / "primary" / "acisf08007_000N001_bpix1.fits"
    first = dest.stat().st_mtime_ns
    obs._fix_archive_bpix()

    assert dest.read_bytes() == b"BPIX"
    assert dest.stat().st_mtime_ns == first, "should not rewrite an existing file"


def test_uncompressed_bpix_is_copied_into_primary(tmp_path):
    """An already-decompressed bpix in secondary/ must still reach primary/.

    download_chandra_obsid does not always gzip its output, and the branch that
    handles the plain-.fits case is the one that has never run.
    """
    obs = _obs(tmp_path)
    src = obs.workdir / "secondary" / "acisf08007_000N001_bpix1.fits"
    src.parent.mkdir(parents=True, exist_ok=True)
    src.write_bytes(b"PLAIN BPIX")
    (obs.workdir / "primary").mkdir(parents=True, exist_ok=True)

    obs._fix_archive_bpix()

    dest = obs.workdir / "primary" / "acisf08007_000N001_bpix1.fits"
    assert dest.is_file(), "plain bpix should be placed in primary/"
    assert dest.read_bytes() == b"PLAIN BPIX"


# --- _get_mica_obspar ------------------------------------------------------


def test_mica_obspar_miss_returns_none(tmp_path, monkeypatch):
    """mica returns None for an obsid it has no entry for; that is a miss."""
    obs = _obs(tmp_path)
    from mica.archive import obspar as mica_obspar

    monkeypatch.setattr(mica_obspar, "get_obspar", lambda obsid, **kw: None)
    assert obs._get_mica_obspar() is None


def test_mica_obspar_failure_is_treated_as_a_miss(tmp_path, monkeypatch):
    """No $SKA, or an unsynced archive, must fall back to downloading."""
    obs = _obs(tmp_path)
    from mica.archive import obspar as mica_obspar

    def boom(obsid, **kw):
        raise KeyError("SKA")

    monkeypatch.setattr(mica_obspar, "get_obspar", boom)
    assert obs._get_mica_obspar() is None


def test_mica_obspar_hit_returns_the_dict(tmp_path, monkeypatch):
    obs = _obs(tmp_path)
    from mica.archive import obspar as mica_obspar

    monkeypatch.setattr(
        mica_obspar,
        "get_obspar",
        lambda obsid, **kw: {"instrume": "ACIS", "obsid": obsid},
    )
    assert obs._get_mica_obspar()["instrume"] == "ACIS"


# --- the obspar has exactly two sources, and CDA is not one of them ---------


def test_archive_mode_refuses_an_obspar_download(tmp_path, monkeypatch):
    """CDA publishes no obspar, so asking for one is a mistake, not a download.

    download_chandra_obsid's filetypes are asol, bpix, ... vv, vvref -- there is
    no observation-parameter file among them. _download_archive ignoring ftypes
    meant this request quietly fetched a full evt2 and five other filetypes and
    then reported that "the download produced nothing".
    """
    obs = _obs(tmp_path)
    monkeypatch.setattr(observation.subprocess, "Popen", _refuse_to_run)
    monkeypatch.setattr(observation.cda, "get_ocat_local", _refuse_to_run)

    with pytest.raises(observation.ObsparUnavailable) as excinfo:
        obs._download_archive(["obspar"])

    message = str(excinfo.value)
    assert "mica" in message
    assert "arc5gl" in message


def test_archive_mode_still_downloads_the_filetypes_it_can_supply(
    tmp_path, monkeypatch
):
    """Only the impossible request is refused; the rest is unchanged."""
    past = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d %H:%M:%S")
    _stub_ocat(monkeypatch, past)
    _stub_ciao(monkeypatch)
    seen = {}
    _stub_popen(
        monkeypatch, FakeProcess(on_communicate=lambda: seen.setdefault("ran", True))
    )

    _obs(tmp_path)._download_archive(["evt2"])

    assert seen.get("ran")


def test_get_obspar_names_its_two_real_sources_when_both_are_absent(
    tmp_path, monkeypatch
):
    """The old message blamed the download; the download was never possible."""
    from mica.archive import obspar as mica_obspar

    obs = _obs(tmp_path)
    monkeypatch.setattr(mica_obspar, "get_obspar", lambda obsid, **kw: None)
    monkeypatch.setattr(observation.subprocess, "Popen", _refuse_to_run)
    monkeypatch.setattr(observation.cda, "get_ocat_local", _refuse_to_run)

    with pytest.raises(observation.ObsparUnavailable) as excinfo:
        obs.get_obspar()

    message = str(excinfo.value)
    assert str(obs.obsid) in message
    assert "mica" in message and "arc5gl" in message


# --- cleaned workdirs stay reusable ----------------------------------------


def test_cleanup_marks_make_images_as_not_done(tmp_path):
    """cleanup_downloads deletes images/, so the task that made them is not done.

    Its docstring keeps cache/ deliberately, "so reruns skip completed steps".
    That left make_images' stored result claiming success while its outputs were
    gone, and the columns derived from them -- pileup, acis_streak, grating_arm --
    came back as measured zeros on any later recompute.
    """
    from astromon.task import ReturnValue

    obs = _obs(tmp_path)
    images = obs.workdir / "images"
    images.mkdir(parents=True)
    (images / "something.img").write_bytes(b"x")
    observation.make_images.set_result(
        obs, ReturnValue(return_code=observation.ReturnCode.OK, msg="made them")
    )

    obs.cleanup_downloads()

    assert not images.exists()
    assert observation.make_images.get_result(obs).return_code == (
        observation.ReturnCode.INVALID
    )


def test_cleanup_leaves_the_detection_results_valid(tmp_path):
    """Only make_images is invalidated -- keeping the .src files is the point.

    cleanup_downloads keeps sources/*.src so detection does not have to re-run,
    so invalidating make_images must not cascade into the detection tasks.
    """
    from astromon.task import ReturnValue

    obs = _obs(tmp_path)
    (obs.workdir / "images").mkdir(parents=True)
    for task in (observation.make_images, observation.wavdetect):
        task.set_result(
            obs, ReturnValue(return_code=observation.ReturnCode.OK, msg="done")
        )

    obs.cleanup_downloads()

    assert observation.wavdetect.get_result(obs).return_code == (
        observation.ReturnCode.OK
    )


def test_cleanup_clears_the_download_marker_so_a_retry_redownloads(tmp_path):
    """cleanup_downloads deletes the raw evt2 that a later retry needs back.

    The marker written by ``_download_archive`` says "already downloaded" and
    used to survive cleanup untouched, so a later retry's ``_download_archive``
    call short-circuited on it and skipped re-fetching the raw evt2 -- leaving
    make_images permanently unable to find its input ("Missing input files for
    task make_images: events") no matter how many times the obsid was retried.
    Documented in production as a ~437-obsid cluster that needed
    clear_invalid_make_images.py to unstick; clearing the marker here means
    a retry never falls into that trap in the first place.
    """
    obs = _obs(tmp_path)
    secondary = obs.workdir / "secondary"
    secondary.mkdir(parents=True)
    obs._archive_download_marker().write_text("download_chandra_obsid completed\n")

    primary = obs.workdir / "primary"
    primary.mkdir(parents=True)
    (primary / "acisf8007N004_evt2.fits.gz").write_bytes(b"raw")
    (primary / "8007_evt2_filtered.fits.gz").write_bytes(b"filtered")

    obs.cleanup_downloads()

    assert not (primary / "acisf8007N004_evt2.fits.gz").exists()
    assert not obs._archive_download_marker().exists(), (
        "the marker must be cleared so a later _download_archive() call"
        " re-fetches the raw evt2 instead of trusting a stale 'done' signal"
    )
