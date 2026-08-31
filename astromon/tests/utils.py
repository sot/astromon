from pathlib import Path

from astropy.table import Table, vstack

from astromon import db

DATA_DIR = Path(__file__).parent / "data"


def create_h5(dbfile="test_data.h5"):
    # utility to create the h5 file from current test data
    dbfile = Path(dbfile)
    names = [
        "astromon_xcorr",
        "astromon_cat_src",
        "astromon_xray_src",
        "astromon_obs",
        "astromon_regions",
    ]
    for name in names:
        # t = db.get_table(name, dbfile)  # fails, file does not exist
        db.save(
            name, Table.read(DATA_DIR / f"{name}.ecsv"), dbfile
        )  # populate the table


def save_test_data(dbfile, datadir=None):
    # utility to create new test data from a dbfile
    names = [
        "astromon_xcorr",
        "astromon_cat_src",
        "astromon_xray_src",
        "astromon_obs",
        "astromon_regions",
    ]
    if datadir is None:
        datadir = DATA_DIR
    for name in names:
        table = db.get_table(name, dbfile=dbfile)
        table.write(datadir / f"{name}.ecsv", overwrite=True)


def minimal_obs(
    obsid: int = 99900,
    ra: float = 83.82,
    dec: float = -5.39,
    target: str = "test",
    category_id: int = 50,
) -> Table:
    """One-row astromon_obs table for use in cross_match tests."""
    return Table(
        {
            "obsid": [obsid],
            "detector": ["ACIS-S"],
            "target": [target],
            "grating": ["NONE"],
            "sim_z": [-190.143],
            "date_obs": ["2020-01-01T00:00:00"],
            "tstart": [0.0],
            "ascdsver": ["10.0"],
            "ra": [ra],
            "dec": [dec],
            "roll": [0.0],
            "category_id": [category_id],
            "version": ["10.0"],
        }
    )


def minimal_xray_src(
    obsid: int = 99900,
    ra: float = 83.82,
    dec: float = -5.39,
    snr: float = 10.0,
    net_counts: float = 500.0,
) -> Table:
    """One celldetect X-ray source for use in cross_match tests."""
    return Table(
        {
            "obsid": [obsid],
            "id": [1],
            "ra": [ra],
            "dec": [dec],
            "net_counts": [net_counts],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "r_angle": [10.0],
            "snr": [snr],
            "near_neighbor_dist": [60.0],
            "pileup": [0.0],
            "acis_streak": [False],
            "caldb_version": ["4.10.0"],
            "detect_method": ["celldetect"],
        }
    )


def minimal_cat_src(
    obsid: int = 99900,
    celldetect_x_id: int = 1,
    catalog: str = "Tycho2",
    name: str = "test-star",
    ra: float = 83.82028,
    dec: float = -5.39,
    mag: float = 8.0,
) -> Table:
    """One catalog source positioned ~1 arcsec from the xray source."""
    return Table(
        {
            "obsid": [obsid],
            "id": [0],
            "celldetect_x_id": [celldetect_x_id],
            "catalog": [catalog],
            "name": [name],
            "ra": [ra],
            "dec": [dec],
            "mag": [mag],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "separation": [1.0],
        }
    )


def two_method_xray_src(obsid=99900, ra=83.82, dec=-5.39):
    """The same source id detected by both methods, as the pipeline stores it."""
    both = vstack([minimal_xray_src(obsid, ra=ra, dec=dec)] * 2)
    both["detect_method"] = ["celldetect", "gaussian_detect"]
    return both
