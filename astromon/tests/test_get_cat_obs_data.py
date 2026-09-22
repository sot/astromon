"""Tests for the per-obsid pipeline entry point in astromon.scripts.get_cat_obs_data."""


def test_split_versions_separates_peak_gaussian():
    """peak_gaussian_detect is run for its .src output, not saved as a method.

    Its fitted position never won on match count; what it carries is the
    peak_offset diagnostic, which get_sources reads off the .src file. So it is
    run but not persisted.
    """
    from astromon.scripts.get_cat_obs_data import _split_versions

    saved, diagnostic = _split_versions(
        ("celldetect", "gaussian_detect", "peak_gaussian_detect")
    )
    assert saved == ["celldetect", "gaussian_detect"]
    assert diagnostic == ["peak_gaussian_detect"]


def test_split_versions_preserves_order_and_handles_absence():
    from astromon.scripts.get_cat_obs_data import _split_versions

    assert _split_versions(("gaussian_detect", "celldetect")) == (
        ["gaussian_detect", "celldetect"],
        [],
    )
    assert _split_versions(("peak_gaussian_detect",)) == ([], ["peak_gaussian_detect"])
    assert _split_versions(()) == ([], [])
