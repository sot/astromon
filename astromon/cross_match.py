""" """

import logging
import re
from pathlib import Path

import numpy as np
import requests
from astropy import coordinates as coords
from astropy import io, table
from astropy import units as u
from astroquery.vizier import Vizier
from cxotime import CxoTime
from Ska.DBI import DBI

import astromon
from astromon import db, observation, utils

logger = logging.getLogger("astromon")


SIM_Z = {"ACIS-I": -233.587, "ACIS-S": -190.143, "HRC-I": 126.983, "HRC-S": 250.466}


def _get_vizier(source, ra, dec, time, radius):
    """
    This fetches the vizier url, but it doesn't parse the result.
    """
    time = CxoTime(time)
    url = (
        "http://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
        "-source={source}"
        "&-sort=_r"
        "&-c={ra:.5f}+{dec:.4f}"
        "&-out.add=_RA(J2000,J{year})"
        "&-c.rs={radius}"
    )
    opts = {
        "source": source,
        "ra": float(ra),
        "dec": float(dec),
        "year": time.frac_year,
        "radius": float(radius),
    }
    rc = requests.get(url.format(**opts))
    output = rc.content.decode()
    n_rows = len(
        [
            line
            for line in output.split("\n")
            if not re.match("#", line) and line.strip()
        ]
    )
    if n_rows:
        return io.ascii.read(output, data_start=3)


CROSS_MATCH_DTYPE = np.dtype(
    [
        # ('obsid', 'int'),
        # ('id', 'int'),
        ("catalog", "<U24"),
        ("name", "<U16"),
        ("ra", float),
        ("dec", float),
        ("mag", float),
    ]
)


def get_vizier(
    pos,
    catalog,
    cat_identifier,
    name_cols,
    columns,
    radius=3 * u.arcsec,
    logging_tag="",
    use_astroquery=True,
    raw=False,
):
    if use_astroquery:
        vizier = Vizier(
            columns=[
                f"_RA(J2000,{pos.obstime.frac_year:8.3f})",
                f"_DE(J2000,{pos.obstime.frac_year:8.3f})",
            ],
        )
        vizier_result = vizier.query_region(
            pos, radius=radius, catalog=cat_identifier, cache=False
        )
        vizier_result = list(vizier_result)
    else:
        vizier_result = [
            _get_vizier(
                cat_identifier,
                p.ra / u.deg,
                p.dec / u.deg,
                pos.obstime,
                radius / u.arcsec,
            )
            for p in pos
        ]
        vizier_result = [r for r in vizier_result if r is not None]

    if raw:
        return vizier_result

    if len(vizier_result) == 0:
        logger.debug(f"{logging_tag} {catalog:>24s} has no results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    vizier_result = table.vstack(vizier_result, metadata_conflicts="silent")

    vizier_result["catalog"] = catalog

    vizier_result["name"] = [
        "-".join([str(row[n]) for n in name_cols]) for row in vizier_result
    ]

    result = table.Table(data=np.zeros(len(vizier_result), dtype=CROSS_MATCH_DTYPE))
    for col in CROSS_MATCH_DTYPE.names:
        src_col = columns.get(col, col)
        if src_col in vizier_result.colnames:
            result[col] = vizier_result[src_col]
        else:
            result[col] = table.MaskedColumn(
                dtype=CROSS_MATCH_DTYPE[col],
                length=len(vizier_result),
                mask=np.ones(len(vizier_result)),
            )
    logger.debug(f"{logging_tag} {catalog:>24s} has {len(vizier_result)} results")
    return result


VIZIER_CATALOGS = {
    "Tycho2": {
        "catalog": "Tycho2",
        "cat_identifier": "I/259/tyc2",
        "name_cols": ["TYC1", "TYC2", "TYC3"],
        "columns": {
            "ra": "_RAJ2000_{time.frac_year:.3f}",
            "dec": "_DEJ2000_{time.frac_year:.3f}",
            "mag": "VTmag",
        },
    },
    # https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/323
    "ICRS": {
        "catalog": "ICRS",
        "cat_identifier": "I/323",
        "name_cols": ["ICRF"],
        "columns": {"ra": "_RAJ2000", "dec": "_DEJ2000"},
    },
    # https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/284
    "USNO-B1.0": {
        "catalog": "USNO-B1.0",
        "cat_identifier": "USNO-B1.0",
        "name_cols": ["USNO-B1.0"],
        "columns": {"ra": "_RAJ2000", "dec": "_DEJ2000", "mag": "R1mag"},
    },
    # http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/239
    "HIP": {
        "catalog": "HIP",
        "cat_identifier": "I/239/hip_main",
        "name_cols": ["HIP"],
        "columns": {"ra": "_RAJ2000", "dec": "_DEJ2000", "mag": "Vmag"},
    },
    # http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/311
    # https://heasarc.gsfc.nasa.gov/W3Browse/all/hipnewcat.html
    "HIP2": {
        "catalog": "HIP",
        "cat_identifier": "I/311/hip2",
        "name_cols": ["HIP"],
        "columns": {
            "ra": "_RAJ2000",
            "dec": "_DEJ2000",
            "mag": "Vmag",
        },
    },
    # http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/322
    "UCAC4": {
        "catalog": "UCAC4",
        "cat_identifier": "I/322",
        "name_cols": ["UCAC4"],
        "columns": {"ra": "_RAJ2000", "dec": "_DEJ2000", "mag": "f.mag"},
    },
    "2MASS": {
        "catalog": "2MASS",
        "cat_identifier": "II/246/out",
        "name_cols": ["_2MASS"],
        "columns": {"ra": "_RAJ2000", "dec": "_DEJ2000", "mag": "Kmag"},
    },
    "SDSS": {
        "catalog": "SDSS",
        "cat_identifier": "II/294",
        "name_cols": ["SDSS"],
        "columns": {
            "ra": "_RAJ2000_{time.frac_year:.3f}",
            "dec": "_DEJ2000_{time.frac_year:.3f}",
            "mag": "rmag",
        },
    },
    # https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/345/gaia2
    "Gaia2": {
        "catalog": "Gaia2",
        "cat_identifier": "I/345/gaia2",
        "name_cols": ["Source"],
        "columns": {"ra": "_RA_ICRS", "dec": "_DE_ICRS", "mag": "Gmag"},
    },
}


def _get(
    catalog,
    time,
    pos=None,
    ra=None,
    dec=None,
    radius=3 * u.arcsec,
    logging_tag="",
    raw=False,
):
    if (ra is None or dec is None) and pos is None:
        raise Exception("pos or ra/dec required")
    if pos is None:
        pos = coords.SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs", obstime=time)
    elif ra is not None or dec is not None:
        raise Exception("ra/dec not required if giving pos")
    params = VIZIER_CATALOGS[catalog].copy()
    columns = {}
    for name, col in params["columns"].items():
        columns[name] = col.format(time=time)
    params["columns"] = columns
    return get_vizier(pos, radius=radius, logging_tag=logging_tag, **params, raw=raw)


def rough_match(
    sources,
    time,
    radius=3 * u.arcsec,
    catalogs=("Tycho2", "ICRS", "USNO-B1.0", "2MASS", "SDSS"),
):
    """
    Find sources in a set of :ref:`catalogs <catalog-list>` around the x-ray sources given.

    Find sources in a set of :ref:`catalogs <catalog-list>` around the x-ray sources given
    in `sources`,  within an angular separation of at most `radius`
    (in :any:`astropy units <astropy:astropy-units>`).

    This function uses the 'ra', 'dec', 'id' and 'obsid' columns of the `sources` table.
    The obsid column is optional and is used only to check that the query corresponds to a single
    OBSID.

    .. Note::

        The current implementation of this function creates a temporary table with the
        cartesian product of the x-ray sources and the catalog sources tables
        (length ~ n_x_ray_sources^2). This is intended to be used with few x-ray sources at a time.

        The alternative to this decision would have been to use a custom join function, which can
        lead to wrong results if the sources happen to be closer than `2*radius`. See warning on
        :ref:`Custom Join functions <astropy:astropy-table-join-functions>`. Using a cartesian
        join allows one to get rough matches without a requirement in `near_neighbor_dist`.

    Parameters
    ----------
    sources: :any:`astropy.table.Table`
        Table of x-ray sources as returned by
        :any:`db.get_table("astromon_xray_src") <db.get_table>` or
        :any:`observation.Observation.get_sources`. Required.
    time: :any:`CxoTime <cxotime:index>`-compatible
        Time of the observation. Required to account for proper motion.
    radius: :any:`astropy.units.Quantity`
        Radius around the x-ray sources to look for counterparts.
    catalogs: list
        A list of catalog names.
        Default is ('Tycho2', 'ICRS', 'USNO-B1.0', '2MASS', 'SDSS')
    """
    if len(sources) == 0:
        return []

    if "obsid" in sources.dtype.names:
        if len(np.unique(sources["obsid"])) > 1:
            raise Exception("rough_match only handles one OBSID at a time")
        logging_tag = f"OBSID={sources['obsid'][0]} "
    else:
        logging_tag = ""

    logger.debug(f"{logging_tag} rough_match started")
    pos = coords.SkyCoord(
        ra=sources["ra"], dec=sources["dec"], unit="deg", frame="icrs", obstime=time
    )

    res = [
        _get(pos=pos, time=time, radius=radius, catalog=name, logging_tag=logging_tag)
        for name in catalogs
    ]
    res = table.vstack(list(res), metadata_conflicts="silent")

    if len(sources) and len(res):
        sources["coord_xray"] = coords.SkyCoord(
            sources["ra"], sources["dec"], unit="deg"
        )
        res["coord_cat"] = coords.SkyCoord(res["ra"], res["dec"], unit="deg")
        res = table.join(res, sources[["coord_xray", "id"]], join_type="cartesian")
        sep = res["coord_xray"].separation(res["coord_cat"])
        res["separation"] = sep.to_value(unit="arcsec")
        res.remove_columns(["coord_xray", "coord_cat"])
        res = res[sep < radius]
        res.rename_column("id", "x_id")
    else:
        dtype = _join_dtype(res.dtype, sources[["id"]].dtype, [])
        res = table.Table(dtype=dtype)
        res.rename_column("id", "x_id")
        res["separation"] = table.Column(dtype=np.float32)

    return res


@utils.logging_call_decorator
def do_sql_cross_match(selection_name):
    sql_script = (
        Path(astromon.__file__).parent / "sql" / "x-corr" / f"{selection_name}.sql"
    )
    if not sql_script.exists():
        logging.error(f"File {sql_script} does not exist")
        sql_script = Path(selection_name)
    if not sql_script.exists():
        logging.error(f"File {sql_script} does not exist")
        msg = f"{selection_name} is not a known selection"
        logging.error(msg)
        raise Exception(msg)

    with open(sql_script) as fh:
        sql_query = fh.read()
        dbi = DBI("sqlite", db.FILE, numpy=False)
        return table.Table(dbi.fetchall(sql_query))


def get_bad_target_mask(matches):
    """
    Returns a mask to remove targets that are known to give cross-matches with large residuals.

    This is not used in standard processing and is here for convenience.

    Parameters
    ----------
    matches: :any:`astropy.table.Table`
        Table of astromon cros-matches as returned by :any:`cross_match` and
        :any:`db.get_cross_matches`. Required.
    """
    ok = np.ones(len(matches), dtype=bool)
    bad_targets = [
        "RW Aur",
        "Tau Boo",
        "70 OPH",
        "16 Cyg",
        "M87",
        "Orion",
        "HD 97950",
        "HD4915",
    ]
    bad_targets = [x.replace(" ", "").lower() for x in bad_targets]
    for ii, target in enumerate(matches["target"]):
        trgt = target.replace(" ", "").lower()
        for bad_target in bad_targets:
            if trgt.startswith(bad_target):
                ok[ii] = False
    return ~ok


def get_excluded_regions_mask(matches, regions=None):
    """
    Returns a mask to remove x-ray sources that are close to excluded regions.

    Parameters
    ----------
    matches: :any:`astropy.table.Table`
        Table of astromon cros-matches as returned by :any:`cross_match` and
        :any:`db.get_cross_matches`. Required.
    """
    if regions is None:
        regions = db.get_table("astromon_regions")
    ii, jj = np.broadcast_arrays(
        np.arange(len(matches))[None, :], np.arange(len(regions))[:, None]
    )
    i, j = ii.flatten(), jj.flatten()
    loc = coords.SkyCoord(regions["ra"] * u.deg, regions["dec"] * u.deg)
    in_region = (
        matches["x_loc"][i].separation(loc[j])
        < regions["radius"][j] * u.arcsec
    )
    in_region &= (regions["obsid"][j] <= 0) | (
        regions["obsid"][j] == matches["obsid"][i]
    )
    in_region = in_region.reshape(ii.shape)
    return np.any(in_region, axis=0)


def filter_matches(  # noqa: PLR0912
    matches,
    snr=None,
    dr=None,
    r_angle=None,
    start=None,
    stop=None,
    r_angle_grating=None,
    near_neighbor_dist=None,
    sim_z=None,
    exclude_regions=False,
    exclude_bad_targets=False,
    exclude_categories=(),
    **kwargs,
):
    """
    Return a mask to filter out cross-matches based on some common criteria.

    Parameters
    ----------
    matches: :any:`astropy.table.Table`
        Table of cross-matches as returned by :any:`db.get_cross_matches` or :any:`cross_match`
    snr: float
        Filter matches based on signal-to-noise. Selects matches['snr'] <= snr.
    dr: float
        Filter matches based on angular offset. Selects matches['dr'] <= dr.
    r_angle: float
        Filter matches based on distance to optical axis.
        Selects matches['r_angle'] <= r_angle (apply to non-grating observations only).
    start: :any:`CxoTime <cxotime:index>`
        Filter matches based on time. Selects matches['tstart'] >= start.
    stop: :any:`CxoTime <cxotime:index>`
        Filter matches based on time. Selects matches['tstart'] <= stop.
    r_angle_grating: float
        Filter matches based on distance to optical axis (apply to grating observations only).
        Selects matches['r_angle_grating'] <= r_angle_grating.
    near_neighbor_dist: float
        Filter matches based on distance to closest neighbor.
        Selects matches['near_neighbor_dist'] <= near_neighbor_dist.
    sim_z: float
        Maximum allowed SIM-Z in mm.
    exclude_bad_targets: bool
        Default is False.
    exclude_regions: bool
        Default is False.
    exclude_categories: tuple
        A list of observation categories to exclude. Default is empty.
    kwargs: dict
        The keys in kwargs determine on which columns to filter. The values of kwargs are used
        according to their type:
        - list types cause a row to be selected if the column value is in the list.
        - otherwise rows are selected if the column value is equal to the kwargs value
    """
    ok = np.ones(len(matches), dtype=bool)

    if dr is not None:
        ok &= matches["dr"] <= dr

    if snr is not None:
        ok &= matches["snr"] >= snr

    if r_angle is not None:
        ok &= matches["r_angle"] <= r_angle

    if start is not None:
        ok &= matches["tstart"] >= CxoTime(start).secs

    if stop is not None:
        ok &= matches["tstart"] <= CxoTime(stop).secs

    if r_angle_grating is not None:
        ok &= matches["r_angle_grating"] <= r_angle_grating

    if near_neighbor_dist is not None:
        ok &= matches["near_neighbor_dist"] <= near_neighbor_dist

    if sim_z is not None:
        sim_z_offset = np.zeros(len(matches))
        for det, det_sim_z in SIM_Z.items():
            msk = matches["detector"] == det
            sim_z_offset[msk] = matches["sim_z"][msk] - det_sim_z
        ok &= np.abs(sim_z_offset) <= sim_z

    if exclude_categories:
        ok &= ~np.in1d(matches["category"], exclude_categories)

    for key, val in kwargs.items():
        if isinstance(val, list):
            ok &= np.isin(matches[key], val)
        else:
            ok &= matches[key] == val

    if exclude_bad_targets:
        ok &= ~get_bad_target_mask(matches)

    if exclude_regions:
        ok &= ~get_excluded_regions_mask(matches)

    return ok


def compute_cross_matches(
    name=None,
    astromon_obs=None,
    astromon_xray_src=None,
    astromon_cat_src=None,
    dbfile=None,
    exclude_regions=False,
    exclude_bad_targets=False,
    logging_tag="",
    **kwargs,
):
    """
    Cross-match x-ray sources with catalog counterparts.

    The arguments to this function are relayed to the actual implementation. The default
    algorithm is the `name="simple"` cross-match
    (see the documentation of :any:`simple_cross_match` for details).

    Parameters
    ----------
    name: str
        Can be the name one of the :any:`standard cross-matches <CROSS_MATCHES>`
        or the name of the algorithm to use (only 'simple' is implemented). Default: 'simple'.
    astromon_obs: :any:`astropy.table.Table`
        Table of astromon observations as returned by
        :any:`db.get_table("astromon_obs") <db.get_table>`.  The default is to get it from `dbfile`.
    astromon_xray_src: :any:`astropy.table.Table`
        Table of x-ray sources as returned by
        :any:`db.get_table("astromon_xray_src") <db.get_table>` or
        :any:`observation.Observation.get_sources`.  The default is to get it from `dbfile`.
    astromon_cat_src: :any:`astropy.table.Table`
        Table of catalog counterparts as returned by
        :any:`db.get_table("astromon_cat_src") <db.get_table>`
        or :any:`rough_match`. The default is to get it from `dbfile`.
    dbfile: str
        HDF5 or SQLite file where the tables are. The default is to let :any:`db.get_table` decide.
    exclude_bad_targets: bool
        Default is False.
    exclude_regions: bool
        Default is False.
    logging_tag: str
        A string to prepend to the logging messages.
    kwargs: dict
        Arguments passed to the algorithm implementation. Ignored if `name` is one of the
        :any:`standard cross-matches <CROSS_MATCHES>`. kwargs is implementation dependant,
        but the arguments of the 'simple' algorithm are:

            - **catalogs**: A list of catalog names. The order matters.
              Default is ['ICRS', 'Tycho2']
            - **snr**: Minimum signal-to-noise ratio of x-ray sources to consider.
              Default is 3.
            - **r_angle**: Maximum r_angle (in arcsec).
              Default is 120 arcsec (2 arcmin).
            - **dr**: Maximum separation between x-ray source and catalog counterpart (in arcsec).
              Default is 3 arcsec.
            - **r_angle_grating**: Maximum r_angle in the case of grating observations (in arcsec).
              Default is 24 arcsec (0.4 arcmin).
            - **near_neighbor_dist**: Only consider x-ray sources with no other x-ray source
              within this radial distance (arcsec). Default is 6 arcsec.
            - **start**: Only consider observations after this :any:`CxoTime <cxotime:index>`.
              Default is to consider all.
    """
    if astromon_xray_src is None:
        astromon_xray_src = db.get_table("astromon_xray_src", dbfile)
    if astromon_cat_src is None:
        astromon_cat_src = db.get_table("astromon_cat_src", dbfile)
    if astromon_obs is None:
        astromon_obs = db.get_table("astromon_obs", dbfile)

    if name is None:
        name = "default" if not kwargs else "simple"

    if name in CROSS_MATCHES_ARGS:
        if kwargs:
            logger.warning(
                f"calling astromon.cross_match with {name=}. kwargs are ignored."
            )
        args = CROSS_MATCHES_ARGS[name].copy()
        method = CROSS_MATCH_METHODS[args.pop("method")]
    elif name in CROSS_MATCH_METHODS:
        method = CROSS_MATCH_METHODS[name]
        args = dict(name=name, **kwargs)
    else:
        raise Exception(f'Unknown x-matching name: "{name}"')

    result = method(
        astromon_obs,
        astromon_xray_src,
        astromon_cat_src,
        logging_tag=logging_tag,
        **args,
    )

    result["time"] = CxoTime(result["date_obs"])
    result["c_loc"] = coords.SkyCoord(result["c_ra"], result["c_dec"], unit="deg")
    result["x_loc"] = coords.SkyCoord(result["x_ra"], result["x_dec"], unit="deg")

    result["category"] = [
        observation.ID_CATEGORY_MAP[cat] for cat in result["category_id"]
    ]

    if exclude_bad_targets:
        result = result[~get_bad_target_mask(result)]
    if exclude_regions:
        regions = db.get_table("astromon_regions", dbfile=dbfile)
        result = result[~get_excluded_regions_mask(result, regions=regions)]

    result.add_index("obsid")
    _set_formats(result)

    return result


def _join_dtype(dtype1, dtype2, keys):
    """
    Get the dtype resulting from the join of two tables with the given dtypes.
    """
    dtype1 = {name: dtype1[name] for name in dtype1.names}
    dtype2 = {name: dtype2[name] for name in dtype2.names}
    for key in keys:
        del dtype2[key]
    dtype1.update(dtype2)
    return np.dtype([(name, dtype1[name]) for name in dtype1])


def simple_cross_match(
    astromon_obs,
    astromon_xray_src,
    astromon_cat_src,
    name="",
    catalogs=("ICRS", "Tycho2"),
    snr=3,
    r_angle=120.0,
    dr=3,
    r_angle_grating=24.0,
    near_neighbor_dist=6.0,
    start=None,
    stop=None,
    logging_tag="",
):
    """
    The simplest cross-match of x-ray sources with catalog counterparts.

    This algorithm selects pairs of x-ray sources and catalog counterparts according to these
    criteria:

    - Observations done after `start`.
    - X-ray sources with signal-over-noise ratio > `snr`.
    - X-ray sources at most `r_angle` off-axis (grating observations) or `r_angle_grating`
      arcsec (non-grating observations).
    - Angular separation between X-ray and catalog counterpart less than `dr` arcsec.
    - X-ray sources that are at most `near_neighbor_dist` arcsec from the closest x-ray source.
    - Counterparts from catalogs included in `catalog`.

    The selected pairs are sorted according to catalog and angular separation.
    The order of precedence for the catalogs is the order in the `catalogs` argument.
    The first pair for each x-ray source is selected as the catalog match.

    .. Warning::

        This algorithm does not check whether two or more sources are paired with the same
        counterpart within a radius `dr`. To prevent this from happening, `near_neighbor_dist`
        should be at least twice `dr`.

    .. Note::

        Remember that this cross-match function is applied over a set of *rough matches* resulting
        from applying :any:`rough_match <cross_match.rough_match>`. In consequence, `dr` should be
        smaller than the one used in :any:`rough_match <cross_match.rough_match>`.

    Parameters
    ----------
    astromon_obs: :any:`astropy.table.Table`
        Table of astromon observations as returned by
        :any:`db.get_table("astromon_obs") <db.get_table>`. Required.
    astromon_xray_src: :any:`astropy.table.Table`
        Table of x-ray sources as returned by
        :any:`db.get_table("astromon_xray_src") <db.get_table>` or
        :any:`observation.Observation.get_sources`. Required.
    astromon_cat_src: :any:`astropy.table.Table`
        Table of catalog counterparts as returned by
        :any:`db.get_table("astromon_cat_src") <db.get_table>`
        or :any:`rough_match`. Required.
    name: str
        The name of the algorithm to use or the name of a standard set of arguments.
        Default is 'simple'
    catalogs: list
        A list of catalog names. The order matters. Default is ['ICRS', 'Tycho2']
    snr: float
        Minimum signal-to-noise ratio of x-ray sources to consider.
        Default is 3.
    r_angle: float
        Maximum r_angle (in arcsec).
        Default is 120 arcsec (2 arcmin).
    dr: float
        Maximum separation between x-ray source and catalog counterpart (in arcsec).
        Default is 3 arcsec.
    r_angle_grating: float
        Maximum r_angle in the case of grating observations (in arcsec).
        Default is 24 arcsec (0.4 arcmin).
    near_neighbor_dist: float
        Only consider x-ray sources with no other x-ray source within this radial distance (arcsec).
        Default is 6 arcsec.
    start:  :any:`CxoTime <cxotime:index>`-compatible timestamp
        Only consider observations after this time.
        Default is to consider all.
    stop:  :any:`CxoTime <cxotime:index>`-compatible timestamp
        Only consider observations before this time.
        Default is to consider all.
    logging_tag: str
        A string to prepend to the logging messages.
    """

    astromon_xray_src = astromon_xray_src[
        (astromon_xray_src["snr"] > snr)
        & (astromon_xray_src["near_neighbor_dist"] > near_neighbor_dist)
    ]
    astromon_cat_src = astromon_cat_src[np.in1d(astromon_cat_src["catalog"], catalogs)]

    # I can rename and drop these to match the standard in astromon.db because I just made copies
    astromon_cat_src.rename_columns(
        ["id", "ra", "dec", "y_angle", "z_angle"],
        ["c_id", "c_ra", "c_dec", "c_y_angle", "c_z_angle"],
    )
    astromon_xray_src.rename_columns(
        ["id", "ra", "dec", "y_angle", "z_angle"],
        ["x_id", "x_ra", "x_dec", "x_y_angle", "x_z_angle"],
    )
    if "name" in astromon_xray_src.colnames:
        astromon_xray_src.remove_column("name")

    if start is not None:
        date_obs = CxoTime(astromon_obs["date_obs"].astype(str))
        astromon_obs = astromon_obs[date_obs > start]

    if stop is not None:
        date_obs = CxoTime(astromon_obs["date_obs"].astype(str))
        astromon_obs = astromon_obs[date_obs <= stop]

    if len(astromon_xray_src) == 0 or len(astromon_cat_src) == 0:
        logger.debug(f"{logging_tag} No xray or cat sources")
        dtype = _join_dtype(astromon_obs.dtype, astromon_xray_src.dtype, ["obsid"])
        dtype = _join_dtype(dtype, astromon_cat_src.dtype, ["obsid", "x_id"])
        return table.Table(dtype=dtype)

    matches = table.join(
        astromon_obs,
        astromon_xray_src,
        keys=["obsid"],
    )
    matches = table.join(
        matches, astromon_cat_src, keys=["obsid", "x_id"], table_names=["xray", "cat"]
    )

    matches["dz"] = matches["x_z_angle"] - matches["c_z_angle"]
    matches["dy"] = matches["x_y_angle"] - matches["c_y_angle"]
    matches["dr"] = np.sqrt(matches["dy"] ** 2 + matches["dz"] ** 2)
    matches["cat_order"] = np.full(len(matches), 200)
    for i, k in enumerate(catalogs):
        matches["cat_order"][matches["catalog"] == k] = i

    matches = matches[
        ((matches["grating"] == "NONE") & (matches["r_angle"] < r_angle))
        | ((matches["grating"] != "NONE") & (matches["r_angle"] < r_angle_grating))
    ]
    matches = matches[matches["dr"] < dr]

    if len(matches) == 0:
        return matches

    mg = matches.group_by(["obsid", "x_id"])
    indices = mg.groups.indices
    idxs = []
    for i0, i1 in zip(indices[:-1], indices[1:], strict=True):
        idxs_sort = np.lexsort((mg["dr"][i0:i1], mg["cat_order"][i0:i1]))
        idxs.append(i0 + idxs_sort[0])
    result = mg[idxs]

    result["select_name"] = name

    return result


CROSS_MATCHES_ARGS = {
    "astromon_21": {
        "name": "astromon_21",
        "method": "simple",
        "catalogs": ["ICRS", "Tycho2"],
        "snr": 3,
        "r_angle": 120.0,
        "r_angle_grating": 120.0,
        "near_neighbor_dist": 6.0,
        "dr": 3.0,
    },
    "astromon_22": {
        "name": "astromon_22",
        "method": "simple",
        "catalogs": ["ICRS", "Tycho2"],
        "snr": 3,
        "r_angle": 120.0,
        "r_angle_grating": 24.0,
        "near_neighbor_dist": 6.0,
        "dr": 3.0,
    },
}
"""
*Standard* cross-match arguments.
"""
CROSS_MATCHES_ARGS["default"] = CROSS_MATCHES_ARGS["astromon_21"]

"""
*Standard* cross-matches (the keys of :any:`CROSS_MATCHES_ARGS`).
"""
CROSS_MATCHES = sorted(CROSS_MATCHES_ARGS.keys())

"""
Available cross-matching algorithms to be used in :any:`compute_cross_matches`.
"""
CROSS_MATCH_METHODS = {"simple": simple_cross_match}


def _set_formats(dat):
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
