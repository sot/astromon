"""
Offline regression check: does the locally-built Tycho-2 index (see
astromon/tycho2_local_index.py) agree with HEASARC's independently-sourced Hipparcos
Main Catalog?

This is a data-integrity check on the Tycho-2 index *build*, not on any particular
x-ray cross-match (that is covered separately by test_tycho2_local_index.py's held-out
fixture). It is meant to be re-run whenever the index is rebuilt from the ESO mirror
files, to catch a silent parsing/column error independent of CDS/VizieR -- which is
otherwise unverifiable while CDS Strasbourg/VizieR is down. It intentionally does not
run a live query: the Hipparcos comparison table is a one-time bulk snapshot fetched
from HEASARC's `hipparcos` table (built by HEASARC directly from CDS I/239 in 2000,
independent of VizieR) via scripts/parse_hipparcos.py (not part of this repository),
not re-fetched at test time.

Both this test and test_tycho2_local_index.py are skipped cleanly if their inputs
aren't found locally:

    - set ASTROMON_TYCHO2_INDEX / ASTROMON_HIPPARCOS_INDEX to the built parquet files, or
    - place them at the documented default locations under
      $SKA/data/astromon/tycho2_local_index/.

Tycho-2 main-catalog positions are quoted at each star's mean epoch (J2000 for most
stars), while Hipparcos positions are quoted at the fixed epoch J1991.25 -- a ~9 year
gap. Comparing raw positions without accounting for this (as an earlier ad hoc check of
a handful of *supplement*-file stars happened not to need, since those are already
quoted at J1991.25) gives a false sense of large disagreement, scaling with each star's
proper motion. Positions are propagated to J1991.25 with the Tycho-2 star's own proper
motion (via propagate_positions) before comparing.

Stars that are resolved components of a multi-star system (ccdm_n_comp > 1) are
excluded: Hipparcos reports a different, blended astrometric/photometric solution for
those, not a bug (see test_tycho2_local_index.KNOWN_CONE_EDGE_CASES for a worked
example of the same effect).
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from astromon.tycho2_local_index import load_star_table, propagate_positions

HIPPARCOS_EPOCH = 1991.25  # J1991.25, the Hipparcos Catalogue's fixed epoch

# Thresholds have margin below the full-catalog agreement rates observed when this test
# was written (99.3% of single-component HIP-linked stars within 0.5 arcsec after
# epoch propagation; 95%+ within 0.1 mag in VT). They are meant to catch gross
# corruption (wrong columns, byte-offset drift, broken parsing), not to pin the exact
# long tail of complex proper motion / imperfect ccdm_n_comp cases.
MIN_STARS_COMPARED = 1000
POSITION_TOLERANCE_ARCSEC = 0.5
MIN_FRAC_POSITION_OK = 0.98
MAG_TOLERANCE = 0.3
MIN_FRAC_MAG_OK = 0.90


def _default_path(env_var: str, relative: str) -> Path | None:
    env_path = os.environ.get(env_var)
    if env_path:
        return Path(env_path)
    if "SKA" in os.environ:
        default = Path(os.environ["SKA"]) / "data" / "astromon" / relative
        if default.exists():
            return default
    return None


@pytest.fixture(scope="module")
def tycho2_stars() -> pd.DataFrame:
    path = _default_path(
        "ASTROMON_TYCHO2_INDEX", "tycho2_local_index/tyc2_stars.parquet"
    )
    if path is None or not path.exists():
        pytest.skip(
            "local Tycho-2 index not found; set the ASTROMON_TYCHO2_INDEX environment "
            "variable to a built tyc2_stars.parquet file (see "
            "astromon/tycho2_local_index.py) to run this test"
        )
    return load_star_table(str(path))


@pytest.fixture(scope="module")
def hipparcos_stars() -> pd.DataFrame:
    path = _default_path(
        "ASTROMON_HIPPARCOS_INDEX", "tycho2_local_index/hipparcos_stars.parquet"
    )
    if path is None or not path.exists():
        pytest.skip(
            "local Hipparcos snapshot not found; set the ASTROMON_HIPPARCOS_INDEX "
            "environment variable to a built hipparcos_stars.parquet file (see "
            "scripts/parse_hipparcos.py) to run this test"
        )
    return pd.read_parquet(path)


def test_hip_linked_stars_agree_with_hipparcos(
    tycho2_stars: pd.DataFrame, hipparcos_stars: pd.DataFrame
) -> None:
    tyc = tycho2_stars[tycho2_stars["HIP"].notna()].copy()
    tyc["HIP"] = tyc["HIP"].astype("int64")

    merged = tyc.merge(hipparcos_stars, on="HIP", suffixes=("_tyc", "_hip"))
    merged = merged[merged["ccdm_n_comp"] <= 1]
    assert len(merged) >= MIN_STARS_COMPARED, (
        f"too few single-component HIP-linked stars to check ({len(merged)}); "
        "expected a full Tycho-2 index and Hipparcos snapshot"
    )

    for_propagation = merged.rename(
        columns={"pmRA_tyc": "pmRA", "pmDE_tyc": "pmDE"}
    )
    ra_prop, dec_prop = propagate_positions(for_propagation, HIPPARCOS_EPOCH)
    sep_arcsec = 3600 * np.hypot(
        (ra_prop - merged["ra"].to_numpy()) * np.cos(np.radians(merged["dec"])),
        dec_prop - merged["dec"].to_numpy(),
    )
    frac_position_ok = float(np.mean(sep_arcsec < POSITION_TOLERANCE_ARCSEC))
    assert frac_position_ok >= MIN_FRAC_POSITION_OK, (
        f"only {frac_position_ok:.1%} of {len(merged)} single-component HIP-linked "
        f"stars agree in position to better than {POSITION_TOLERANCE_ARCSEC} arcsec "
        "after propagating to the Hipparcos J1991.25 epoch -- possible local Tycho-2 "
        "index corruption"
    )

    mag_diff = (merged["VT_tyc"] - merged["VT_hip"]).dropna().abs()
    frac_mag_ok = float((mag_diff < MAG_TOLERANCE).mean())
    assert frac_mag_ok >= MIN_FRAC_MAG_OK, (
        f"only {frac_mag_ok:.1%} of {len(mag_diff)} single-component HIP-linked stars "
        f"agree in VT magnitude to better than {MAG_TOLERANCE} mag -- possible local "
        "Tycho-2 index corruption"
    )
