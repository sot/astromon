import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table
from cxotime import CxoTime

from astromon import db, utils
from astromon.web import celmon


def _fake_matches(n_years=5):
    """Cross-match rows for two detectors with clearly different offset distributions.

    ACIS-S gets small radial offsets (dr in [0.1, 0.3)), HRC-I gets large ones
    (dr in [1.0, 2.0)), so each detector's own quantiles are easy to tell apart
    from the all-detector aggregate and from each other.
    """
    rng = np.random.default_rng(0)
    n = 40
    now = CxoTime()
    time = now - rng.uniform(0, n_years - 1, 2 * n) * 365.25 * u.day

    acis_dr = rng.uniform(0.1, 0.3, n)
    hrc_dr = rng.uniform(1.0, 2.0, n)
    dr = np.concatenate([acis_dr, hrc_dr])
    angle = rng.uniform(0, 2 * np.pi, 2 * n)

    return Table(
        {
            "obsid": np.arange(2 * n),
            "x_id": np.ones(2 * n, dtype=int),
            "detector": ["ACIS-S"] * n + ["HRC-I"] * n,
            "caldb_version": ["4.10.0"] * (2 * n),
            "time": time,
            "dr": dr,
            "dy": dr * np.cos(angle),
            "dz": dr * np.sin(angle),
        }
    )


@pytest.fixture
def _no_plotting(monkeypatch):
    """create_figures_cal's plots write PNGs; the bug is only about `result`."""
    monkeypatch.setattr(celmon, "plot_offsets_history", lambda *args, **kwargs: None)
    monkeypatch.setattr(celmon, "plot_cdf", lambda *args, **kwargs: None)


def test_create_figures_cal_uses_each_detectors_own_quantiles(
    tmp_path, monkeypatch, _no_plotting
):
    """Each detector's entry in the result dict must reflect its own quantiles.

    create_figures_cal's per-detector loop used to call result.update(...) with
    `quantiles` from *before* cdf_(m) recomputed it for the current detector, so
    the first detector got the all-detector aggregate and every later detector got
    the previous detector's quantiles (one iteration stale). With the fix, each
    detector's quantiles must match calling cdf_() on that detector's own subset,
    and must differ from the aggregate and from the other detector.
    """
    matches = _fake_matches()
    monkeypatch.setattr(db, "get_cross_matches", lambda **kwargs: matches)
    monkeypatch.setattr(
        utils,
        "get_calalign_offsets",
        lambda all_matches, **kwargs: Table(
            {"after_caldb": np.zeros(len(all_matches), dtype=bool)}
        ),
    )

    result = celmon.create_figures_cal(tmp_path)

    acis = matches[matches["detector"] == "ACIS-S"]
    hrc = matches[matches["detector"] == "HRC-I"]
    _, _, acis_quantiles = celmon.cdf_(acis)
    _, _, hrc_quantiles = celmon.cdf_(hrc)
    expected_acis = {f"q{100 * q['q']:.0f}": q["offset"] for q in acis_quantiles}
    expected_hrc = {f"q{100 * q['q']:.0f}": q["offset"] for q in hrc_quantiles}

    assert result["ACIS_S"] == pytest.approx(expected_acis)
    assert result["HRC_I"] == pytest.approx(expected_hrc)
    assert result["ACIS_S"] != pytest.approx(result["HRC_I"])
