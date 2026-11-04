"""Quaia and GaiaVarStar are no longer queried, stored or matched.

The absolute-astrometry calibration dropped both (Jean Connelly, 2026-11-04): Quaia's
quasars are Gaia DR3 quasar candidates that GaiaQSO, MilliquasGaia or DESIV161 mostly
hold too, and GaiaVarStar's stars are variables.
"""

import pytest

from astromon import cross_match
from astromon.scripts import get_cat_obs_data
from astromon.scripts.maintenance import requery_cat_src

RETIRED = ("Quaia", "GaiaVarStar")


@pytest.mark.parametrize(
    "name",
    [
        "get_quaia",
        "get_quaia_candidates",
        "get_gaia_var_stars",
        "get_gaia_var_stars_by_obsid",
    ],
)
def test_the_catalog_getters_are_gone(name):
    assert not hasattr(cross_match, name)
    assert not hasattr(get_cat_obs_data, name)
    assert not hasattr(requery_cat_src, name)


@pytest.mark.parametrize("catalog", RETIRED)
def test_no_catalog_cache_is_kept(catalog):
    assert catalog not in cross_match.CATALOG_CACHE_PATHS


@pytest.mark.parametrize("catalog", RETIRED)
def test_requery_does_not_know_them(catalog):
    assert catalog not in requery_cat_src.ALL_CATALOGS
    assert catalog not in requery_cat_src.CATALOG_GETTERS
    assert catalog not in requery_cat_src.BATCHED_GETTERS


def test_requery_catalogs_are_the_default_selections_plus_rough_matches():
    assert set(requery_cat_src.ALL_CATALOGS) == set(cross_match.DEFAULT_CATALOGS) | {
        "USNO-B1.0",
        "2MASS",
        "SDSS",
    }
