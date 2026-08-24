import hashlib
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

__all__ = [
    "get_table",
    "get_cross_matches",
    "save",
    "add_regions",
    "update_regions",
    "remove_regions",
    "sync_regions",
    "get_regions",
    "is_in_excluded_region",
]

if "ASTROMON_FILE" in os.environ:
    FILE = Path(os.environ["ASTROMON_FILE"])
elif "NO_DEFAULT_ASTROMON" in os.environ:
    FILE = None
else:
    FILE = Path(os.environ["SKA"]) / "data" / "astromon" / "astromon.h5"

# astromon_xcorr.x_id plus detect_method is the final per-method match key: it
# joins to astromon_xray_src on (obsid, id, detect_method). Unlike
# astromon_cat_src.celldetect_x_id, this one really does identify a pairing.
ASTROMON_XCORR_DTYPE = np.dtype(
    [
        ("select_name", "S14"),
        ("obsid", np.int32),
        ("c_id", np.int32),
        ("x_id", np.int32),
        ("dy", np.float32),
        ("dz", np.float32),
        ("dr", np.float32),
        ("detect_method", "S24"),
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
        # Ratio of fitted Gaussian sigma to PSF ECF=90% radius.
        # Values > 1 indicate emission broader than the PSF.
        ("psfratio", np.float32),
        # Fraction of source counts within 2" vs 10" of the fitted position.
        # Near 1 for point sources; near 0 for extended emission.
        ("concentration_ratio", np.float32),
        # True if source falls inside a grating arm mask region (HEG, MEG, or LETG),
        # excluding the zero-order circle. Used to reject dispersed-spectrum detections.
        ("grating_arm", np.int32),
        # True for the single brightest source (highest SNR) in the observation.
        # Used to allow the calibration target into cross-matching even when it is
        # also flagged acis_streak=True (e.g. grating spectrum detected as a streak).
        ("brightest", np.int32),
        # Source detection method that produced this row, e.g. "celldetect"
        # or "gaussian_detect".
        ("detect_method", "S24"),
        # Angular distance in arcsec from this row's position to the
        # peak-seeded Gaussian fit of the same source. A stability diagnostic:
        # when the peak-seeded fit disagrees with the centroid, the source is
        # extended, confused, or faint. NaN if not computed.
        ("peak_offset", np.float32),
    ]
)


# astromon_cat_src.celldetect_x_id is a celldetect-scoped anchor: the celldetect
# source this catalogue row was queried around, and what `separation` is measured
# from. It is provenance, not a match key -- there is one row per catalogue source
# per observation and no detect_method column, so it cannot express a per-method
# pairing. The pairing lives in astromon_xcorr, which does have detect_method.
# get_cross_matches drops this column before joining for that reason.
ASTROMON_CAT_SRC_DTYPE = np.dtype(
    [
        ("obsid", np.int32),
        ("id", np.int32),
        ("celldetect_x_id", np.int32),
        ("catalog", "S16"),
        ("name", "S32"),
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
        ("region_id_str", "S10"),
        ("last_modified", "S21"),
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


def create_empty_tables(dbfile=None, table_names=None):
    """Create any of the known tables that are missing from `dbfile`, with zero rows.

    Idempotent: tables that already exist are left untouched, data and all. Use this
    to initialize a new database so that incremental writers can pass
    ``expect_existing=True`` to :func:`save` and have a missing table treated as
    corruption rather than as a first write.

    Parameters
    ----------
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`
    table_names: list of str, optional
        Tables to create. Defaults to all of :data:`DTYPES`.

    Returns
    -------
    list of str
        Names of the tables that were actually created.
    """
    if table_names is None:
        table_names = list(DTYPES)
    created = []
    with connect(dbfile, mode="r+") as h5:
        for table_name in table_names:
            if table_name in h5.root:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                h5.create_table("/", table_name, np.zeros(0, dtype=DTYPES[table_name]))
            created.append(table_name)
    return created


def missing_column_fill(dtype: np.dtype, name: str):
    """Value for a column required by *dtype* but absent from the data.

    Float columns get NaN, because 0.0 is a meaningful measurement for the ones
    that can be missing (psfratio, concentration_ratio) and zero-filling makes an
    unmeasured source indistinguishable from a genuinely compact one -- silently
    admitting it into threshold cuts. Everything else gets numpy's zero value:
    there is no integer NaN, and 0/"" is the right conservative default for the
    flag and label columns.

    This is the single fill policy for the package; get_sources and the pipeline
    save path both go through it so a column cannot mean different things
    depending on which door it came in.
    """
    return np.nan if dtype[name].kind == "f" else np.zeros(1, dtype=dtype[name])[0]


def conform_to_dtype(data: table.Table, table_name: str) -> table.Table:
    """Add any columns `table_name`'s schema requires but `data` lacks.

    Fills them per :func:`missing_column_fill`. Returns `data` (modified in place).
    """
    dtype = DTYPES[table_name]
    for name in dtype.names:
        if name in data.colnames:
            continue
        former = RENAMED_COLUMNS.get(name)
        if former is not None and former in data.colnames:
            data[name] = data[former].astype(dtype[name])
            continue
        fill = missing_column_fill(dtype, name)
        data[name] = np.full(len(data), fill, dtype=dtype[name])
    return data


#: Columns that have been renamed, as {new name: former name}. _cast_to_dtype
#: matches fields by name and fills anything absent via missing_column_fill --
#: which returns 0 for integers -- so without this a file written before a rename
#: would read the new column as a silent zero rather than its stored value. Kept
#: permanently: backups and older copies of the database stay readable, and the
#: next rename gets the same protection for free.
RENAMED_COLUMNS = {"celldetect_x_id": "x_id"}


def _cast_to_dtype(arr: np.ndarray, dtype: np.dtype) -> np.ndarray:
    """Cast a structured numpy array to *dtype*.

    Fields missing from `arr` are taken from their former name if they have one
    (see :data:`RENAMED_COLUMNS`), and otherwise filled per
    :func:`missing_column_fill`.
    """
    if arr.dtype == dtype:
        return arr
    result = np.zeros(len(arr), dtype=dtype)
    for name in dtype.names:
        if name in arr.dtype.names:
            result[name] = arr[name]
            continue
        former = RENAMED_COLUMNS.get(name)
        if former is not None and former in arr.dtype.names:
            result[name] = arr[former]
            continue
        result[name] = missing_column_fill(dtype, name)
    return result


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
                res = _cast_to_dtype(res, DTYPES[table_name])
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


def _validated_replace_keys(
    table_name, replace_keys, data_names, ignore_obsid, select_name_key
):
    """Validate save()'s replace_keys against its legacy flags and the data columns."""
    if ignore_obsid:
        raise ValueError("replace_keys cannot be combined with ignore_obsid=True")
    if select_name_key:
        raise ValueError("replace_keys cannot be combined with select_name_key=True")
    missing_keys = [k for k in replace_keys if k not in data_names]
    if missing_keys:
        raise ValueError(
            f"replace_keys columns not in {table_name}: {', '.join(missing_keys)}"
        )
    return tuple(replace_keys)


def save(  # noqa: PLR0912
    table_name,
    data,
    dbfile,
    ignore_obsid=False,
    select_name_key=False,
    replace_keys=None,
    expect_existing=False,
):
    """
    Insert data into a table, deleting previous entries for the same OBSID.

    If the table does not exist, it is created using pre-existing table definitions.

    If `ignore_obsid` is `False`, and `data` has an "obsid" column, this function replaces all rows
    in the table whose OBSID is in `data["obsid"]`. It does not replace the whole table. If `data`
    has all rows of a given obsid removed, those entries are NOT removed from the table. For
    example, the following code has no effect:

        data = db.get_table("astromon_xray_src", dbfile)
        data = data[data['obsid'] != 12345]
        db.save("astromon_xray_src", data, dbfile)

    To remove all entries for a given obsid, set `ignore_obsid=True`. If you do that,
    the entire table is replaced by `data`. The following removes all entries for obsid 12345:

        data = db.get_table("astromon_xray_src", dbfile)
        data = data[data['obsid'] != 12345]
        db.save("astromon_xray_src", data, dbfile, ignore_obsid=True)

    Parameters
    ----------
    table_name: str
        The name of the table.
    data: :any:`astropy.table.Table`
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`
    ignore_obsid: bool
        If True, do not consider obsid to decide which rows to keep in the table.
    select_name_key: bool
        When True and the table has both ``detect_method`` and ``select_name`` columns,
        key on ``(obsid, detect_method, select_name)`` instead of ``(obsid, detect_method)``.
        This lets independent catalog backfills update one select_name without clobbering rows
        from other select_names for the same obsid. Equivalent to
        ``replace_keys=("obsid", "detect_method", "select_name")``.
    replace_keys: tuple of str, optional
        Explicit key columns for the delete-before-append: existing rows whose key-column
        values match any row of `data` are replaced, everything else is kept. For example
        ``replace_keys=("catalog",)`` on ``astromon_cat_src`` rebuilds one catalog without
        touching the others, and ``replace_keys=("obsid", "catalog")`` replaces only the
        (obsid, catalog) pairs present in `data`. Cannot be combined with `ignore_obsid`
        or `select_name_key`.
    expect_existing: bool
        When True, raise :any:`MissingTableException` instead of creating `table_name`
        if it is not already in the file. Incremental writers should set this: `save`
        replaces a table by removing the node and re-creating it, which is not atomic,
        so a process killed in that window leaves the table missing. Without this flag
        the next `save` recreates it holding only its own rows, silently discarding
        every other obsid. Initialize a new database with :func:`create_empty_tables`
        so that a missing table always means corruption rather than a first write.
    """
    with connect(dbfile, mode="r+") as h5:
        if not h5.isopen:
            raise RuntimeError(f"{h5.filename} is not open")
        if h5.mode not in ["r+", "w"]:
            raise RuntimeError(f"{h5.filename} is not open for writing")
        if not isinstance(data, table.Table):
            raise TypeError("input to _save_hdf5 must be a table")

        if table_name in DTYPES:
            dtype = DTYPES[table_name]
            # A column under its former name is present, not missing: resolve it
            # rather than rejecting the write. Fixtures and external dumps predate
            # the rename, and the values are the right ones either way.
            source_of = {
                name: (
                    name
                    if name in data.dtype.names
                    else RENAMED_COLUMNS.get(name, name)
                )
                for name in dtype.names
            }
            names = [n for n in dtype.names if source_of[n] in data.dtype.names]
            if len(names) == 0:
                raise Exception("Input data has no columns in common with table in DB")
            missing = [
                name for name in dtype.names if source_of[name] not in data.dtype.names
            ]
            if missing:
                raise Exception(
                    f"Saving table {table_name} with missing columns: {', '.join(missing)}"
                )
            b = np.zeros(len(data), dtype=dtype)
            for name in names:
                b[name] = data[source_of[name]].data.astype(dtype[name])
            data = b
        else:
            data = data.as_array()

        if replace_keys is not None:
            replace_keys = _validated_replace_keys(
                table_name,
                replace_keys,
                data.dtype.names,
                ignore_obsid,
                select_name_key,
            )

        if expect_existing and table_name not in h5.root:
            raise MissingTableException(
                f"{table_name} is missing from {h5.filename} but expect_existing=True. "
                "The table may have been destroyed by a process killed during a save; "
                "check for data loss before re-running. Use create_empty_tables() to "
                "initialize a new database."
            )

        if table_name in h5.root:
            node = h5.get_node(f"/{table_name}")
            if replace_keys is not None:
                stored = node[:]
                data_out = (
                    _cast_to_dtype(stored, dtype) if table_name in DTYPES else stored
                )
                key_cols = list(replace_keys)
                data_out = data_out[
                    ~np.isin(data_out[key_cols], np.unique(data[key_cols]))
                ]
                data = np.concatenate((data_out, data))
            elif not ignore_obsid and "obsid" in data.dtype.names:
                data_out = _cast_to_dtype(node[:], dtype)
                has_detect = (
                    "detect_method" in data.dtype.names
                    and "detect_method" in data_out.dtype.names
                )
                has_select = (
                    "select_name" in data.dtype.names
                    and "select_name" in data_out.dtype.names
                )
                # Replace rows matching the incoming data on the key columns:
                # (obsid) by default, (obsid, detect_method) when the table has
                # detect_method so different methods coexist, and additionally
                # select_name with select_name_key=True so independent catalog
                # backfills don't clobber each other's rows.
                key_cols = ["obsid"]
                if has_detect:
                    key_cols.append("detect_method")
                    if select_name_key and has_select:
                        key_cols.append("select_name")
                # One vectorized membership test on the structured key subset
                # instead of one full-table scan per unique key.
                data_out = data_out[
                    ~np.isin(data_out[key_cols], np.unique(data[key_cols]))
                ]
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
        sel = ~np.isin(all_regions["region_id"], regions)
        all_regions = all_regions[sel]

        # We just removed some regions. In the process, maybe all rows for a given obsid were
        # removed. In the following, if ignore_obsids is False, these rows would not be removed.
        save("astromon_regions", all_regions, dbfile, ignore_obsid=True)


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

        if "astromon_meta" in h5.root and len(h5.root.astromon_meta) > 0:
            meta = Table(h5.root.astromon_meta[:], dtype=DTYPES["astromon_meta"])
        else:
            meta = Table(np.zeros(1, dtype=DTYPES["astromon_meta"]))

        # Take the floor from whichever of the counter and the data is further
        # along. An absent or zero-row astromon_meta means "the counter is
        # unknown", not "no region_id has ever been allocated" -- reading it as
        # the latter restarts numbering at 1 and collides with existing rows,
        # and because remove_regions() deletes by region_id a single duplicate
        # makes one removal take out two regions. Trusting only the table would
        # be wrong in the other direction: sync_regions can leave the counter
        # above max(region_id) once rows are removed, and those retired ids must
        # not be handed out again.
        max_existing = int(np.max(all_regions["region_id"])) if len(all_regions) else 0
        rid = max(int(meta["last_region_id"][0]), max_existing) + 1
        regions = Table(regions)
        names = [n for n in all_regions.dtype.names if n in regions.dtype.names]
        b = np.zeros(len(regions), dtype=all_regions.dtype)
        b[names] = regions[names].as_array().astype(all_regions.dtype[names])
        b["region_id"] = np.arange(rid, rid + len(b))
        b["last_modified"] = CxoTime.now().date
        # we need to ensure region_id_str are unique. Using the auto-incrementing region_id does
        # not work when we want to synchronize two or more DBs, because adding different regions
        # to each DB could generate the same region_id and cause a collision when syncing.
        # With a 40-bit random number/string, the probability of collision is less than 10^9.
        region_id_str = [
            hashlib.sha256(str(row).encode()).hexdigest()[:10] for row in b
        ]
        while np.any(sel := np.isin(region_id_str, all_regions["region_id_str"])):
            region_id_str[sel] = [
                hashlib.sha256(str(row).encode()).hexdigest()[:10]
                for row in region_id_str[sel]
            ]
        b["region_id_str"] = region_id_str
        b = Table(b)
        all_regions = table.vstack([all_regions, b])
        meta["last_region_id"][0] = b["region_id"][-1]
        save("astromon_regions", all_regions, h5)
        save("astromon_meta", meta, h5)
        return b


def update_regions(regions, dbfile=None):
    """
    Update exclusion regions in the astromon_regions table.

    Things to keep in mind:

        - last_modified is updated automatically.
        - Changes region_id_str are ignored.
        - Entries that do not differ from the values in the table are not updated
          (last_modified is not set).
        - This does not remove regions. Use :any:`remove_regions` for that.

    Raises an exception if:

        - the astromon_regions table is not in the dbfile.
        - any of the region_id in `regions` does not exist in the table.

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
        logger.info(f"Updating regions: {regions}")

        all_regions = get_table("astromon_regions", h5)

        _, from_idx, to_idx = np.intersect1d(
            regions["region_id_str"], all_regions["region_id_str"], return_indices=True
        )

        last_modified = CxoTime.now().date
        for to_i, from_i in zip(to_idx, from_idx, strict=True):
            if not np.array_equal(all_regions[to_i], regions[from_i]):
                all_regions[to_i] = regions[from_i]
                all_regions[to_i]["last_modified"] = last_modified

        save("astromon_regions", all_regions, h5)

        return regions


def get_regions(obsid=None, dbfile=None, radius=5 * u.arcmin):
    """
    Get exclusion regions from the astromon_regions table.

    Parameters
    ----------
    obsid: int
        If provided, only return regions with obsid=0 or centered around the pointing
        of the given obsid, as stored in the astromon_obs table.
    dbfile: :any:`pathlib.Path`
        File where tables are stored.
        The default is `$ASTROMON_FILE` or `$SKA/data/astromon/astromon.h5`
    radius: :any:`astropy.units.Quantity`
        Maximum distance from the nominal pointing of the observation to include a region.
        Only used if obsid is provided. Default is 5 arcmin.
    Returns
    -------
    :any:`astropy.table.Table`
    """
    regions = get_table("astromon_regions", dbfile)
    if obsid is not None:
        astromon_obs = get_table("astromon_obs", dbfile)
        obs_row = astromon_obs[astromon_obs["obsid"] == obsid][0]
        obs_loc = SkyCoord(obs_row["ra"] * u.deg, obs_row["dec"] * u.deg)

        loc = SkyCoord(regions["ra"] * u.deg, regions["dec"] * u.deg)
        regions = regions[
            (obs_loc.separation(loc) < radius) | (regions["obsid"] == obsid)
        ]
    return regions


def _ensure_unique_region_ids(regions, last_region_id):
    """Give every row a distinct ``region_id``.

    Rows keep their id on a first-come basis, so the destination's existing rows
    (which come first in the merge) stay stable and only later duplicates move.
    New ids continue above both `last_region_id` and the largest id in the table,
    so a retired id is never reissued.

    Returns
    -------
    tuple
        ``(regions, last_region_id)`` with the counter advanced to the new
        maximum id.
    """
    if len(regions) == 0:
        return regions, int(last_region_id)

    ids = np.asarray(regions["region_id"], dtype=int)
    next_id = max(int(last_region_id), int(ids.max())) + 1
    seen = set()
    for i, region_id in enumerate(ids):
        if int(region_id) in seen:
            regions["region_id"][i] = next_id
            seen.add(next_id)
            next_id += 1
        else:
            seen.add(int(region_id))
    new_last = max(int(last_region_id), int(np.max(regions["region_id"])))
    return regions, new_last


def sync_regions(from_file, to_file, remove=False):
    """
    Sync exclusion regions from one astromon DB file to another.

    Parameters
    ----------
    from_file: :any:`pathlib.Path`
        Source DB file.
    to_file: :any:`pathlib.Path`
        Destination DB file.
    """
    from_regions = get_table("astromon_regions", from_file)
    to_regions = get_table("astromon_regions", to_file)
    meta = get_table("astromon_meta", to_file)
    # update records that are in both, and last_modified is newer in from_file
    _, from_idx, to_idx = np.intersect1d(
        from_regions["region_id_str"], to_regions["region_id_str"], return_indices=True
    )
    sel = from_regions["last_modified"][from_idx] > to_regions["last_modified"][to_idx]
    up_to_idx = to_idx[sel]
    up_from_idx = from_idx[sel]
    for i1, i2 in zip(up_to_idx, up_from_idx, strict=True):
        to_regions[i1] = from_regions[i2]
    # remove records that are in to_file but not in from_file
    missing = ~np.isin(to_regions["region_id_str"], from_regions["region_id_str"])
    if remove and np.any(missing):
        to_regions = to_regions[~missing]
    # add records that are only in from_file
    new = ~np.isin(from_regions["region_id_str"], to_regions["region_id_str"])
    if np.any(new):
        to_regions = table.vstack([to_regions, from_regions[new]])

    # region_id_str is what identifies a region across files, but region_id is
    # what remove_regions() deletes by, so the merged table has to come out with
    # unique ids. Renumber against the whole merged set: reassigning only the
    # rows that clash with the destination, from a counter that ignores the other
    # incoming ids, lands on ids those incoming rows already hold and turns one
    # duplicate into several. A duplicate already present in either input is
    # repaired here too.
    to_regions, last_region_id = _ensure_unique_region_ids(
        to_regions, meta["last_region_id"][0]
    )
    meta["last_region_id"][0] = last_region_id

    to_regions.sort("region_id")

    save("astromon_regions", to_regions, to_file)
    save("astromon_meta", meta, to_file)


def is_in_excluded_region(position, obsid=None, regions=None, dbfile=None):
    """
    Returns a mask to remove x-ray sources that are close to excluded regions.

    Parameters
    ----------
    position: :any:`astropy.coordinates.SkyCoord`
        Array of :class:`astropy.coordinates.SkyCoord` objects.
    obsid: array-like
        Array of OBSIDs corresponding to each position.
    """
    squeeze = not np.shape(position)
    position = np.atleast_1d(position)

    if obsid is not None:
        obsid = np.atleast_1d(obsid)
        if obsid.shape != position.shape:
            raise ValueError("obsid and position must have the same shape")

    if regions is None:
        regions = get_table("astromon_regions", dbfile)
    ii, jj = np.broadcast_arrays(
        np.arange(len(position))[None, :], np.arange(len(regions))[:, None]
    )
    i, j = ii.flatten(), jj.flatten()
    regions["loc"] = SkyCoord(regions["ra"] * u.deg, regions["dec"] * u.deg)
    in_region = (
        position[i].separation(regions["loc"][j]) < regions["radius"][j] * u.arcsec
    )
    if obsid is not None:
        in_region &= (regions["obsid"][j] <= 0) | (regions["obsid"][j] == obsid[i])
    in_region = in_region.reshape(ii.shape)
    result = np.any(in_region, axis=0)
    if squeeze:
        return result.squeeze()
    return result


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

    # The anchor is celldetect-scoped provenance, not the match key: xcorr already
    # carries (x_id, detect_method), which is the per-method pairing. Dropping it
    # here keeps it from colliding with the xray_src x_id below, and keeps callers
    # from mistaking it for the matched source.
    astromon_cat_src.remove_column("celldetect_x_id")
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
    matches = table.join(
        matches, astromon_xray_src, keys=["obsid", "x_id", "detect_method"]
    )

    matches["time"] = CxoTime(matches["date_obs"])
    matches["c_loc"] = SkyCoord(matches["c_ra"] * u.deg, matches["c_dec"] * u.deg)
    matches["x_loc"] = SkyCoord(matches["x_ra"] * u.deg, matches["x_dec"] * u.deg)

    set_formats(matches)
    if kwargs:
        ok = filter_matches(matches, **kwargs)
        matches = matches[ok]
    matches.add_index("obsid")
    return matches
