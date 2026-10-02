"""Tests for astromon.scripts.reprocess_pointing."""

import tempfile
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from astromon import db
from astromon.scripts import reprocess_pointing


def test_recompute_angles_preserves_rows_missing_from_pointing(monkeypatch):
    """An obsid absent from `pointing` must keep its existing angles, not be zeroed.

    get_corrected_pointing can leave some obsids unresolved (logged as
    n_missing). recompute_angles used to pre-zero the whole y_angle/z_angle/
    r_angle arrays and only fill in resolved obsids, so an unresolved obsid's
    real angles were silently overwritten with 0.0 -- indistinguishable from a
    source located exactly at the aimpoint.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "reprocess_pointing_test.h5"
        monkeypatch.setattr(reprocess_pointing, "DBFILE", str(dbfile))

        t = Table(np.zeros(2, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
        t["obsid"] = [1001, 1002]
        t["id"] = [1, 1]
        t["ra"] = [10.0, 20.0]
        t["dec"] = [5.0, -5.0]
        t["y_angle"] = [111.0, 222.0]
        t["z_angle"] = [333.0, 444.0]
        t["r_angle"] = [
            np.sqrt(111.0**2 + 333.0**2),
            np.sqrt(222.0**2 + 444.0**2),
        ]
        db.save("astromon_xray_src", t, str(dbfile))

        # Only obsid 1001 resolves to a corrected pointing; 1002 is missing.
        pointing = {1001: (10.0, 5.0, 0.0)}

        result = reprocess_pointing.recompute_angles("astromon_xray_src", pointing)

        resolved = result[result["obsid"] == 1001][0]
        unresolved = result[result["obsid"] == 1002][0]

        # Resolved obsid was actually recomputed (not required to match the
        # original values, just to have gone through radec_to_yagzag).
        assert np.isfinite(resolved["y_angle"])
        assert np.isfinite(resolved["z_angle"])

        # Unresolved obsid must keep its original angles, not 0.0.
        assert float(unresolved["y_angle"]) == 222.0
        assert float(unresolved["z_angle"]) == 444.0
        assert float(unresolved["r_angle"]) == pytest.approx(
            np.sqrt(222.0**2 + 444.0**2), rel=1e-5
        )
