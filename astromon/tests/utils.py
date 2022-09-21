from pathlib import Path

from astropy.table import Table

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
