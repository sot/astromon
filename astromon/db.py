import logging
import os
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import tables
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table.table import Table
from cxotime import CxoTime
from ska_helpers.retry import tables_open_file

from astromon import observation
from astromon.utils import MissingTableException

__all__ = ["get_table", "get_cross_matches", "save", "add_regions", "remove_regions"]

if "ASTROMON_FILE" in os.environ:
    FILE = Path(os.environ["ASTROMON_FILE"])
elif "NO_DEFAULT_ASTROMON" in os.environ:
    FILE = None
else:
    FILE = Path(os.environ["SKA"]) / "data" / "astromon" / "astromon.h5"

ASTROMON_XCORR_DTYPE = np.dtype(
    [
        ("select_name", "S14"),
        ("obsid", np.int32),
        ("c_id", np.int32),
        ("x_id", np.int32),
        ("dy", np.float32),
        ("dz", np.float32),
        ("dr", np.float32),
    ]
)


ASTROMON_XRAY_SRC_DTYPE = np.dtype(
    [
        ("obsid", np.int32),
        ("id", np.int32),
        ("ra", np.float64),
        ("dec", np.float64),
        ("net_counts", np.float32),
        ("y_angle", np.float32),
        ("z_angle", np.float32),
        ("r_angle", np.float32),
        ("snr", np.float32),
        ("near_neighbor_dist", np.float32),
        ("pileup", np.float32),
        ("acis_streak", np.int32),
        ("caldb_version", "S10"),
    ]
)


ASTROMON_CAT_SRC_DTYPE = np.dtype(
    [
        ("obsid", np.int32),
        ("id", np.int32),
        ("x_id", np.int32),
        ("catalog", "S16"),
        ("name", "S24"),
        ("ra", np.float64),
        ("dec", np.float64),
        ("separation", np.float32),
        ("mag", np.float32),
        ("y_angle", np.float32),
        ("z_angle", np.float32),
    ]
)


ASTROMON_OBS_DTYPE = np.dtype(
    [
        ("obsid", np.int32),
        ("version", np.float32),
        ("detector", "S6"),
        ("target", "S28"),
        ("grating", "S4"),
        ("sim_z", np.float32),
        ("date_obs", "S20"),
        ("tstart", np.float32),
        ("ascdsver", "S32"),
        ("ra", np.float64),
        ("dec", np.float64),
        ("roll", np.float64),
        ("category_id", np.int32),
    ]
)


ASTROMON_REGION_DTYPE = np.dtype(
    [
        ("region_id", np.int32),
        ("ra", np.float64),
        ("dec", np.float64),
        ("radius", np.float32),
        ("obsid", np.int32),
        ("user", "S50"),
        ("comments", "S200"),
    ]
)


ASTROMON_META_DTYPE = np.dtype(
    [
        ("last_region_id", np.int32),
    ]
)


DTYPES = {
    "astromon_xcorr": ASTROMON_XCORR_DTYPE,
    "astromon_xray_src": ASTROMON_XRAY_SRC_DTYPE,
    "astromon_cat_src": ASTROMON_CAT_SRC_DTYPE,
    "astromon_obs": ASTROMON_OBS_DTYPE,
    "astromon_regions": ASTROMON_REGION_DTYPE,
    "astromon_meta": ASTROMON_META_DTYPE,
}


def create_table(table_name):
    """
    Create an empty table using standard dtypes.

    Known table names:

    - astromon_xcorr
    - astromon_xray_src
    - astromon_cat_src
    - astromon_obs
    - astromon_regions
    - astromon_meta

    Parameters
    ----------
    table_name: str
        Name of the table to retrieve, which specifies the dtype.

    Returns
    -------
    :any:`astropy.table.Table`
    """
    return table.Table(names=DTYPES[table_name].names, dtype=DTYPES[table_name])


def get_table(table_name, dbfile=None):
    """
    Get an entire table from the DB file.

    Known table names:

    - astromon_xcorr
    - astromon_xray_src
    - astromon_cat_src
    - astromon_obs
    - astromon_regions
    - astromon_meta

    Parameters
    ----------
    table_name: str
        Name of the table to retrieve. If the requested table is not in the file,
        an exception will be raised.
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`

    Returns
    -------
    :any:`astropy.table.Table`
    """
    # logger = logging.getLogger('astromon')
    # if not Path(dbfile).exists():
    #     raise RuntimeError(f'Astromon DB file does not exist {dbfile}')

    with connect(dbfile) as con:
        try:
            res = con.get_node(f"/{table_name}")[:]
            if table_name in DTYPES:
                res = res.astype(DTYPES[table_name])
            result = table.Table(res)
            result.convert_bytestring_to_unicode()
            set_formats(result)
            if "obsid" in result.colnames:
                result.add_index("obsid")
        except tables.NoSuchNodeError:
            names = sorted(set([n.name for n in con.root] + list(DTYPES.keys())))
            raise MissingTableException(
                f"{table_name} not in file. Available tables: {names}"
            ) from None

    return result


@contextmanager
def connect(dbfile=None, mode="r"):
    """
    Context manager that returns a DB connection (or an HDF5 file).

    Parameters
    ----------
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`

    Returns
    -------
    :any:`tables.File <tables.file.File>` or :any:`sqlite3.Connection`
    """
    if dbfile is None:
        dbfile = FILE

    if isinstance(dbfile, tables.file.File):
        yield dbfile
    else:
        logger = logging.getLogger("astromon")
        if dbfile is None:
            dbfile = FILE
        dbfile = Path(str(dbfile)).absolute()

        if mode == "r+" and not dbfile.exists():
            mode = "w"
        logger.debug(f"{dbfile} open")
        h5 = tables_open_file(dbfile, mode, delay=1, tries=10)
        try:
            yield h5
        except Exception:
            if h5.isopen:
                h5.close()
                logger.debug(f"{dbfile} closed (1)")
            raise
        finally:
            if h5.isopen:
                h5.close()
                logger.debug(f"{dbfile} closed (2)")


def save(table_name, data, dbfile):
    """
    Insert data into a table, deleting previous entries for the same OBSID.

    If the table does not exist, it is created using pre-existing table definitions.

    Parameters
    ----------
    table_name: str
        The name of the table.
    data: :any:`astropy.table.Table`
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`
    """
    with connect(dbfile, mode="r+") as h5:
        # sanity checks: assert that file is open for writing
        assert h5.isopen, f"{h5.filename} is not open"
        assert h5.mode in ["r+", "w"], f"{h5.filename} is not open for writing"
        assert isinstance(data, table.Table), "input to _save_hdf5 must be a table"

        if table_name in DTYPES:
            dtype = DTYPES[table_name]
            names = [n for n in dtype.names if n in data.dtype.names]
            if len(names) == 0:
                raise Exception("Input data has no columns in common with table in DB")
            missing = [name for name in dtype.names if name not in data.dtype.names]
            if missing:
                raise Exception(
                    f"Saving table {table_name} with missing columns: {', '.join(missing)}"
                )
            b = np.zeros(len(data), dtype=dtype)
            b[names] = data[names].as_array().astype(dtype[names])
            data = b
        else:
            data = data.as_array()

        if table_name in h5.root:
            node = h5.get_node(f"/{table_name}")
            if "obsid" in data.dtype.names:
                # remove rows for these obsids
                obsids = np.unique(data["obsid"])
                data_out = node[:].astype(dtype)
                data_out = data_out[~np.in1d(data_out["obsid"], obsids)]
                # append current data
                data = np.concatenate((data_out, data))
            h5.remove_node(node)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            h5.create_table("/", table_name, data)


def remove_regions(regions, dbfile=None):
    """
    Remove exclusion regions (by ID) from the astromon_regions table.

    Raises an exception if the astromon_regions table is not in the dbfile.

    Parameters
    ----------
    regions: list
        list of integer region ID.
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`
    """
    with connect(dbfile, mode="r+") as h5:
        all_regions = get_table("astromon_regions", h5)
        all_regions = all_regions[~np.in1d(all_regions["region_id"], regions)]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if "astromon_regions" in h5.root:
                h5.remove_node("/astromon_regions")
            h5.create_table("/", "astromon_regions", all_regions.as_array())


def add_regions(regions, dbfile=None):
    """
    Add exclusion regions to the astromon_regions table.

    A unique region ID is automatically generated.
    Raises an exception if the astromon_regions table is not in the dbfile.

    Parameters
    ----------
    regions: `astropy.table.Table`-compatible
        This parameter gets converted to an `astropy.table.Table`.
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`
    """
    with connect(dbfile, mode="r+") as h5:
        logger = logging.getLogger("astromon")
        logger.info(f"Adding regions: {regions}")
        all_regions = get_table("astromon_regions", h5)
        if "astromon_meta" in h5.root:
            meta = Table(h5.root.astromon_meta[:], dtype=DTYPES["astromon_meta"])
        else:
            meta = Table(np.zeros(1, dtype=DTYPES["astromon_meta"]))
        rid = meta["last_region_id"][0] + 1
        regions = Table(regions)
        names = [n for n in all_regions.dtype.names if n in regions.dtype.names]
        b = np.zeros(len(regions), dtype=all_regions.dtype)
        b[names] = regions[names].as_array().astype(all_regions.dtype[names])
        b["region_id"] = np.arange(rid, rid + len(b))
        b = Table(b)
        all_regions = table.vstack([all_regions, b])
        meta["last_region_id"][0] = b["region_id"][-1]
        save("astromon_regions", all_regions, h5)
        save("astromon_meta", meta, h5)
        return b


def set_formats(dat):
    """
    Sets format of columns with float dtype to show 2 decimals, except ra/dec/pileup (4 decimals).

    Parameters
    ----------
    dat: `astropy.table.Table`
    """
    fmts = {
        "ra": ".4f",
        "x_ra": ".4f",
        "c_ra": ".4f",
        "dec": ".4f",
        "x_dec": ".4f",
        "c_dec": ".4f",
        "pileup": ".4f",
    }
    for col in dat.itercols():
        if (
            hasattr(col, "name")
            and col.name in dat.colnames
            and hasattr(dat[col.name], "dtype")
            and col.dtype.kind == "f"
        ):
            dat[col.name].info.format = fmts.get(col.name, ".2f")


def get_cross_matches(name="astromon_21", dbfile=None, **kwargs):
    """
    Get a standard cross-match of observations, x-ray sources and catalog counterparts in `dbfile`.

    A *standard* cross-match is a pre-computed cross-match between x-ray sources and catalog
    counterparts. The name of the standard cross-match specifies the algorithm and the set of
    parameters that have been used to cross-match.

    If you want a *non-standard* cross-match, with your own set of parameters, refer to the
    :any:`compute_cross_matches <cross_match.compute_cross_matches>` function.

    This function returns a :any:`Table <astropy:astropy.table.Table>` indexed by OBSID and adds a
    few columns on the fly:

    - time
    - c_loc
    - x_loc

    Parameters
    ----------
    name: str
        name of the standard cross-match.
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`
    kwargs: dict
        All extra arguments are passed directly to :any:`cross_match.filter_matches`. These are
        some of the allowed arguments:

        - snr
        - dr
        - r_angle
        - start
        - stop
        - r_angle_grating
        - near_neighbor_dist
        - sim_z
        - exclude_regions
        - exclude_bad_targets
        - exclude_categories

        Other extra arguments can be passed. In these cases, the argument name determines on which
        column to filter. Each argument is interpreted according to the type of value passed:

        - list values cause a row to be selected if the value at that row and column is in the list.
        - other values cause a row to be selected if the value at that row and column is equal to
          the argument value

    Returns
    -------
    :any:`astropy.table.Table`
    """
    from astromon.cross_match import filter_matches

    matches = get_table("astromon_xcorr", dbfile)
    matches = matches[matches["select_name"] == name]
    astromon_cat_src = get_table("astromon_cat_src", dbfile)
    astromon_xray_src = get_table("astromon_xray_src", dbfile)
    astromon_obs = get_table("astromon_obs", dbfile)

    astromon_obs["category"] = [
        observation.ID_CATEGORY_MAP[cat] for cat in astromon_obs["category_id"]
    ]

    astromon_cat_src.remove_column("x_id")
    astromon_cat_src.rename_columns(
        ["id", "ra", "dec", "y_angle", "z_angle"],
        ["c_id", "c_ra", "c_dec", "c_y_angle", "c_z_angle"],
    )
    astromon_xray_src.rename_columns(
        ["id", "ra", "dec", "y_angle", "z_angle"],
        ["x_id", "x_ra", "x_dec", "x_y_angle", "x_z_angle"],
    )
    matches = table.join(matches, astromon_obs, keys=["obsid"])
    matches = table.join(matches, astromon_cat_src, keys=["obsid", "c_id"])
    matches = table.join(matches, astromon_xray_src, keys=["obsid", "x_id"])

    matches["time"] = CxoTime(matches["date_obs"])
    matches["c_loc"] = SkyCoord(matches["c_ra"] * u.deg, matches["c_dec"] * u.deg)
    matches["x_loc"] = SkyCoord(matches["x_ra"] * u.deg, matches["x_dec"] * u.deg)

    set_formats(matches)
    if kwargs:
        ok = filter_matches(matches, **kwargs)
        matches = matches[ok]
    matches.add_index("obsid")
    return matches
