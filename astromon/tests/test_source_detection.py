from unittest.mock import patch

import numpy as np
import scipy.optimize
from astropy.table import Table

from astromon.source_detection import _fit, fit_gaussian_2d


def _point_source_events(yag0=0.0, zag0=0.0, sigma=0.5, n=500, seed=42):
    rng = np.random.default_rng(seed)
    yag = rng.normal(yag0, sigma, n)
    zag = rng.normal(zag0, sigma, n)
    return Table({"y_angle": yag, "z_angle": zag})


def _non_converging_result(x0):
    """A scipy.optimize.OptimizeResult as returned by a failed minimize() call."""
    return scipy.optimize.OptimizeResult(
        x=np.asarray(x0, dtype=float),
        success=False,
        message="Non-convergent (forced for test)",
        hess_inv=1e10 * np.eye(len(x0)),
        fun=np.inf,
        nit=0,
    )


class TestFitDoesNotRaiseOnNonConvergence:
    """_fit used to raise RuntimeError when scipy.optimize.minimize did not converge,
    which meant fit_gaussian_2d's own `if not result.success: return fail_value` check
    was dead code -- the exception had already propagated past it. That in turn crashed
    _fit_gaussian_sources's per-source loop (astromon/observation.py), aborting fit
    detection for every other source in the same obsid instead of recording just the
    one failed fit. _fit must return the failed OptimizeResult instead, so the caller's
    existing fail_value branch can do its job.
    """

    def test_fit_returns_failed_result_instead_of_raising(self):
        events = _point_source_events(n=500)
        source = {"y_angle": 0.0, "z_angle": 0.0}

        with patch.object(
            scipy.optimize,
            "minimize",
            side_effect=lambda _fn, x0, **kwargs: _non_converging_result(x0),
        ):
            result = _fit(events, source)

        assert result.success is False

    def test_fit_gaussian_2d_returns_fail_value_on_non_convergence(self):
        events = _point_source_events(n=500)
        source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": 3}

        with patch.object(
            scipy.optimize,
            "minimize",
            side_effect=lambda _fn, x0, **kwargs: _non_converging_result(x0),
        ):
            result = fit_gaussian_2d(events, source)

        assert result["fit_ok"] is False
        assert result["COMPONENT"] == 3
