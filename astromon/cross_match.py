""" """

import logging
import re
import threading
import time as _time
import urllib.request
from pathlib import Path

import numpy as np
import requests
from astropy import coordinates as coords
from astropy import io, table
from astropy import units as u
from astroquery.vizier import Vizier
from cxotime import CxoTime
from Ska.DBI import DBI
from ska_helpers.retry import retry

import astromon
from astromon import db, observation, utils

logger = logging.getLogger("astromon")


class _CallTimeoutError(Exception):
    """Raised by `_run_with_timeout` when the call doesn't finish in time."""


def _run_with_timeout(func, args=(), kwargs=None, timeout=120):
    """Run func in a daemon thread; raise _CallTimeoutError if it exceeds timeout seconds."""
    if kwargs is None:
        kwargs = {}
    outcome = {}

    def _target():
        try:
            outcome["value"] = func(*args, **kwargs)
        except BaseException as e:
            outcome["error"] = e

    thread = threading.Thread(target=_target, daemon=True)
    thread.start()
    thread.join(timeout=timeout)
    if thread.is_alive():
        raise _CallTimeoutError(
            f"{getattr(func, '__name__', func)} did not complete within {timeout}s"
        )
    if "error" in outcome:
        raise outcome["error"]
    return outcome["value"]


SIM_Z = {"ACIS-I": -233.587, "ACIS-S": -190.143, "HRC-I": 126.983, "HRC-S": 250.466}

# Vizier queries occasionally fail with a transient connection reset unrelated to the query
# itself.
retry_on_connection_error = retry(
    exceptions=(requests.exceptions.ConnectionError, requests.exceptions.Timeout),
    tries=3,
    delay=1,
    backoff=2,
)


@retry_on_connection_error
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
        ("name", "<U32"),
        ("ra", float),
        ("dec", float),
        ("mag", float),
    ]
)


@retry_on_connection_error
def _query_vizier_region(vizier, pos, radius, cat_identifier):
    return vizier.query_region(pos, radius=radius, catalog=cat_identifier, cache=False)


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
        vizier_result = _query_vizier_region(vizier, pos, radius, cat_identifier)
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


def _cache_is_stale(cache_path: Path, max_age_days: int | None) -> bool:
    """True if `cache_path` is absent, or exists but is older than `max_age_days`.

    Shared staleness policy for the locally cached catalogs (RFC, ICRF3, Milliquas,
    Quaia): absent means "never downloaded", `max_age_days=None` means "never
    force a refresh".
    """
    return not cache_path.exists() or (
        max_age_days is not None
        and (_time.time() - cache_path.stat().st_mtime) > max_age_days * 86400
    )


_CATALOG_READ_CACHE: dict[Path, tuple[float, table.Table]] = {}


def _read_catalog_cached(cache_path: Path, reader) -> table.Table:
    """Read a cached catalog file, memoized on (path, mtime).

    The per-obsid pipeline calls each catalog getter once per obsid; without this
    the multi-million-row FITS/ECSV/RFC files are re-read and re-parsed on every
    call (measured ~0.6 s per obsid in total). Keying on mtime means a refreshed
    download (e.g. a new RFC release) is picked up on the next call. Callers must
    treat the returned table as read-only.
    """
    mtime = cache_path.stat().st_mtime
    entry = _CATALOG_READ_CACHE.get(cache_path)
    if entry is None or entry[0] != mtime:
        entry = (mtime, reader(cache_path))
        _CATALOG_READ_CACHE[cache_path] = entry
    return entry[1]


def _local_catalog_near(
    catalog: table.Table,
    catalog_name: str,
    pos: coords.SkyCoord,
    radius: u.Quantity,
    logging_tag: str = "",
) -> table.Table:
    """Cone-search a locally cached (name, ra, dec[, mag]) catalog around *pos*.

    Shared by the catalogs that are downloaded once and cached whole (RFC, ICRF3,
    Quaia, GaiaAGN, GaiaQSO) rather than queried live per-position like VizieR or
    Gaia TAP. If the catalog has a ``mag`` column its values are used (masked
    where non-finite); otherwise `mag` is left fully masked.

    A cheap flat-sky bounding cone runs first so that the multi-million-row catalogs
    never have to become a `SkyCoord` in full, then the survivors get an exact
    per-source separation cut. Only sources within `radius` of at least one position
    in `pos` are returned.

    Parameters
    ----------
    catalog
        Full local catalog table with columns 'name', 'ra' (deg), 'dec' (deg),
        and optionally 'mag'.
    catalog_name
        Value to fill in the output 'catalog' column (e.g. 'RFC').
    pos
        Positions to search around (scalar or array SkyCoord).
    radius
        Search radius. Sources farther than this from every position in `pos`
        are excluded.
    logging_tag
        Prefix for log messages.

    Returns
    -------
    table.Table
        Catalog sources with columns matching :data:`CROSS_MATCH_DTYPE`.
    """
    if not pos.isscalar and len(pos) == 0:
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    # Build a bare SkyCoord from the coordinates: this handles a scalar `pos`
    # uniformly and drops frame attributes (e.g. obstime) that would otherwise have
    # to match when comparing against the catalog positions.
    pos_sc = coords.SkyCoord(
        ra=np.atleast_1d(pos.ra.deg), dec=np.atleast_1d(pos.dec.deg), unit="deg"
    )

    # Bounding cone: a cheap flat-sky pre-filter so the multi-million-row catalogs
    # never become a SkyCoord in full. Any input position works as the cone
    # reference, so use the first one -- a coordinate-wise mean would land at RA=180
    # for a field straddling the RA=0/360 meridian and exclude the whole field. The
    # cone radius adapts to how far the positions actually spread (a hardcoded pad
    # fails for crowded/cluster fields spanning many arcminutes, e.g. Perseus).
    center = pos_sc[0]
    max_offset_deg = float(np.max(center.separation(pos_sc).deg))
    # Over-inclusion here is free -- the exact cut below removes the extras -- while
    # under-inclusion is a silent miss, so keep the flat-sky bound deliberately loose.
    search_rad = 1.05 * (radius.to_value(u.deg) + max_offset_deg) + 1e-4
    cos_dec = np.cos(np.radians(center.dec.deg))
    # Wrap the RA difference into (-180, 180] before scaling. Without this, a source
    # just across the RA=0/360 meridian looks ~360 deg away instead of ~0, so every
    # counterpart on the other side of the meridian is silently dropped.
    dra = ((catalog["ra"] - center.ra.deg + 180.0) % 360.0 - 180.0) * cos_dec
    ddec = catalog["dec"] - center.dec.deg
    nearby_mask = np.sqrt(dra**2 + ddec**2) < search_rad
    nearby = catalog[nearby_mask]

    if len(nearby) > 0:
        # Exact per-source cut. The bounding cone is sized to the whole field, so on
        # its own it admits sources arcminutes from any input position, and the
        # callers that use these catalogs directly (get_gaia_agn,
        # get_gaia_qso_candidates, get_quaia_candidates) feed astromon_cat_src with
        # no later radius cut of their own.
        nearby_sc = coords.SkyCoord(nearby["ra"], nearby["dec"], unit="deg")
        _, sep_to_nearest, _ = nearby_sc.match_to_catalog_sky(pos_sc)
        nearby = nearby[sep_to_nearest < radius]

    if len(nearby) == 0:
        logger.debug(f"{logging_tag}{catalog_name} has 0 results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    if "mag" in nearby.colnames:
        mag_vals = np.asarray(nearby["mag"], dtype=np.float32)
        mag = table.MaskedColumn(
            mag_vals, dtype=CROSS_MATCH_DTYPE["mag"], mask=~np.isfinite(mag_vals)
        )
    else:
        mag = table.MaskedColumn(
            dtype=CROSS_MATCH_DTYPE["mag"],
            length=len(nearby),
            mask=np.ones(len(nearby), dtype=bool),
        )
    output = table.Table(
        {
            "catalog": np.full(len(nearby), catalog_name),
            "name": np.asarray(nearby["name"], dtype=CROSS_MATCH_DTYPE["name"]),
            "ra": np.asarray(nearby["ra"], dtype=float),
            "dec": np.asarray(nearby["dec"], dtype=float),
            "mag": mag,
        }
    )
    logger.debug(f"{logging_tag}{catalog_name} has {len(output)} results")
    return output


# Where the locally cached whole-catalog files live (RFC, ICRF3, Milliquas,
# Quaia, GaiaAGN, GaiaQSO). Resolved in utils so that observation.ARCHIVE_DIR,
# which sits in the same tree, cannot disagree with it. Honours
# ASTROMON_DATA_DIR; see utils.
_ASTROMON_DATA_DIR = utils.ASTROMON_DATA_DIR

_RFC_INDEX_URL = "http://astrogeo.org/sol/rfc/"
# The page where astrogeo.org states which release is current. Read this rather
# than listing _RFC_INDEX_URL: a quarterly directory appears there, holding only a
# landing page, before its data files are published.
_RFC_LANDING_URL = "http://astrogeo.org/rfc/"
# Matches only a link *into* a release directory, so that image and PDF paths for
# other releases (/rfc/rfc_2026b_map.pdf) are not mistaken for the announcement.
_RFC_ANNOUNCED_RE = re.compile(r'href="/sol/rfc/(rfc_\d{4}[a-d])', re.IGNORECASE)
_RFC_CACHE_PATH = _ASTROMON_DATA_DIR / "rfc_catalog.txt"
# How often to re-check astrogeo.org for a newer quarterly release. The production
# pipeline (astromon-cross-match) runs once a day, so checking more often than that
# buys nothing, and checking daily means we are never more than a day behind a new
# release.
_RFC_CHECK_INTERVAL_DAYS = 1


def _rfc_release_marker(cache_path: Path) -> Path:
    """Path to the small sidecar file recording which release `cache_path` holds."""
    return cache_path.with_name(cache_path.name + ".release")


@retry_on_connection_error
def _discover_latest_rfc_release() -> str:
    """Return the RFC release that astrogeo.org currently announces.

    astrogeo.org publishes RFC as dated quarterly releases (rfc_2026a, rfc_2026b,
    ...) with no stable "latest" alias, so the current one has to be discovered.

    It has to be read from the announcement page, not inferred from the directory
    listing at `_RFC_INDEX_URL`. A release's directory is created there -- with a
    landing page and no data files -- before the release is published: rfc_2026c
    existed like that while rfc_2026b was still current. Taking the
    lexicographically-largest directory name therefore names a release whose
    catalogue is not there yet, every download from it 404s, and `get_rfc` falls
    back to the cache while logging what looks like a transient network problem.
    Because the fallback never writes the release marker, the check is then due on
    every subsequent call.
    """
    resp = requests.get(_RFC_LANDING_URL, timeout=30)
    resp.raise_for_status()
    releases = sorted(set(_RFC_ANNOUNCED_RE.findall(resp.text)))
    if not releases:
        raise RuntimeError(f"No RFC release announced at {_RFC_LANDING_URL}")
    return releases[-1]


def _parse_rfc(text: str) -> table.Table:
    """Parse RFC fixed-width ASCII into a Table with columns name, ra, dec."""
    names: list[str] = []
    ras: list[float] = []
    decs: list[float] = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        name = line[0:14].strip()
        ra_h = int(line[26:28])
        ra_m = int(line[29:31])
        ra_s = float(line[32:41])
        sign = line[42]
        dec_d = int(line[43:45])
        dec_m = int(line[46:48])
        dec_s = float(line[49:57])
        ra_deg = 15.0 * (ra_h + ra_m / 60.0 + ra_s / 3600.0)
        dec_deg = dec_d + dec_m / 60.0 + dec_s / 3600.0
        if sign == "-":
            dec_deg = -dec_deg
        names.append(name)
        ras.append(ra_deg)
        decs.append(dec_deg)
    return table.Table({"name": names, "ra": np.array(ras), "dec": np.array(decs)})


def get_rfc(
    cache_path: Path | None = None,
    max_age_days: int | None = _RFC_CHECK_INTERVAL_DAYS,
) -> table.Table:
    """Return the RFC catalog, downloading/refreshing from astrogeo.org as needed.

    astrogeo.org has no stable URL for "the current RFC release": each quarter gets
    its own directory (rfc_2026a, rfc_2026b, ...) and old releases are not updated
    in place.  So this discovers the current release name from the astrogeo.org
    directory listing (a small HTML fetch) and only downloads the ~5 MB catalog
    file itself when that release differs from what is already cached — tracked in
    a small sidecar file next to `cache_path`.

    The release check only runs when `cache_path` is absent or older than
    `max_age_days`.  The production pipeline runs once a day, so the default of 1
    day means the check happens at most once per pipeline run, not once per obsid.

    If the release check or download fails and a cached copy already exists, this
    logs a warning and falls back to the (possibly stale) cached copy rather than
    raising, so a transient network problem does not take down the whole
    cross-match run.

    Parameters
    ----------
    cache_path
        Override the default cache location (default: ``$ASTROMON_DATA_DIR/rfc_catalog.txt``).
    max_age_days
        How often to check astrogeo.org for a newer release, in days. Default 1.
        Pass ``None`` to check only when `cache_path` is absent.
    """
    if cache_path is None:
        cache_path = _RFC_CACHE_PATH
    cache_path = Path(cache_path)
    release_marker = _rfc_release_marker(cache_path)

    if _cache_is_stale(cache_path, max_age_days) or not release_marker.exists():
        try:
            latest_release = _discover_latest_rfc_release()
        except Exception as e:
            if cache_path.exists():
                logger.warning(
                    f"RFC release check failed ({e}); using cached copy at {cache_path}"
                )
                return _read_catalog_cached(
                    cache_path, lambda p: _parse_rfc(p.read_text())
                )
            raise

        cached_release = (
            release_marker.read_text().strip() if release_marker.exists() else None
        )
        if cached_release == latest_release and cache_path.exists():
            logger.debug(f"RFC release unchanged ({latest_release}); using cached copy")
            cache_path.touch()  # reset the check-interval clock
        else:
            url = f"{_RFC_INDEX_URL}{latest_release}/{latest_release}_cat.txt"
            logger.info(f"Downloading RFC {latest_release} from {url} → {cache_path}")
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                urllib.request.urlretrieve(url, cache_path)
            except Exception as e:
                if cache_path.exists():
                    logger.warning(
                        f"RFC download failed ({e}); using stale cached copy at {cache_path}"
                    )
                    return _read_catalog_cached(
                        cache_path, lambda p: _parse_rfc(p.read_text())
                    )
                raise
            release_marker.write_text(latest_release)
            logger.info(
                f"  Saved {cache_path.stat().st_size / 1e6:.1f} MB ({latest_release})"
            )
    else:
        logger.debug(f"Loading cached RFC from {cache_path}")
    return _read_catalog_cached(cache_path, lambda p: _parse_rfc(p.read_text()))


def _get_rfc(
    pos: coords.SkyCoord, radius: u.Quantity, logging_tag: str = ""
) -> table.Table:
    """Return RFC sources near *pos* within *radius*, formatted like VizieR results."""
    return _local_catalog_near(get_rfc(), "RFC", pos, radius, logging_tag)


# ---------------------------------------------------------------------------
# ICRF3 — International Celestial Reference Frame 3rd realisation (S/X band)
# Charlot et al. 2020, A&A 644, A159  (VizieR: J/A+A/644/A159)
# 4536 radio sources with VLBI-measured positions.  Cached locally like RFC.
# ---------------------------------------------------------------------------

_ICRF3_VIZIER_ID = "J/A+A/644/A159/table10"
_ICRF3_CACHE_PATH = _ASTROMON_DATA_DIR / "icrf3_catalog.ecsv"
_ICRF3_MAX_AGE_DAYS = 365  # ICRF3 is updated very rarely


def get_icrf3(
    cache_path: Path | None = None,
    max_age_days: int | None = _ICRF3_MAX_AGE_DAYS,
) -> table.Table:
    """Return the ICRF3 S/X-band catalog, downloading/caching from VizieR as needed.

    The catalog (4536 sources) is cached at `cache_path` and refreshed when
    absent or older than `max_age_days`.

    Parameters
    ----------
    cache_path
        Override the default cache location (default: ``$ASTROMON_DATA_DIR/icrf3_catalog.ecsv``).
    max_age_days
        Refresh threshold in days.  Default 365.  Pass ``None`` to download
        only when the file is absent.

    Returns
    -------
    table.Table
        Table with columns ``name`` (str), ``ra`` (deg), ``dec`` (deg).

    Examples
    --------
    >>> t = get_icrf3()
    >>> len(t) > 4000
    True
    """
    if cache_path is None:
        cache_path = _ICRF3_CACHE_PATH
    cache_path = Path(cache_path)
    if _cache_is_stale(cache_path, max_age_days):
        logger.info(f"Downloading ICRF3 from VizieR {_ICRF3_VIZIER_ID} → {cache_path}")
        vizier = Vizier(row_limit=-1, columns=["ICRF", "_RAJ2000", "_DEJ2000"])
        result = vizier.get_catalogs(_ICRF3_VIZIER_ID)
        if not result:
            raise RuntimeError(
                f"VizieR returned no data for ICRF3 catalog {_ICRF3_VIZIER_ID}"
            )
        raw = result[0]
        icrf3 = table.Table(
            {
                "name": np.asarray(raw["ICRF"], dtype="U32"),
                "ra": np.asarray(raw["_RAJ2000"], dtype=float),
                "dec": np.asarray(raw["_DEJ2000"], dtype=float),
            }
        )
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        icrf3.write(cache_path, format="ascii.ecsv", overwrite=True)
        logger.info(f"  Saved {len(icrf3)} ICRF3 sources to {cache_path}")
    else:
        logger.debug(f"Loading cached ICRF3 from {cache_path}")
        icrf3 = _read_catalog_cached(
            cache_path, lambda p: table.Table.read(p, format="ascii.ecsv")
        )
    return icrf3


def _get_icrf3(
    pos: coords.SkyCoord, radius: u.Quantity, logging_tag: str = ""
) -> table.Table:
    """Return ICRF3 sources near *pos* within *radius*, formatted like VizieR results."""
    return _local_catalog_near(get_icrf3(), "ICRF3", pos, radius, logging_tag)


# ---------------------------------------------------------------------------
# ICRF2 — International Celestial Reference Frame 2nd realisation
# Ma et al. 2009, IERS Technical Note 35  (VizieR: I/323)
# 3414 radio sources with VLBI-measured positions. Superseded by ICRF3, but kept
# available (not in any standard CROSS_MATCHES_ARGS hierarchy) as a reference
# baseline for comparing against the pre-ICRF3/RFC catalog set. Cached locally
# like RFC/ICRF3.
# ---------------------------------------------------------------------------

_ICRF2_VIZIER_ID = "I/323"
_ICRF2_CACHE_PATH = _ASTROMON_DATA_DIR / "icrf2_catalog.ecsv"
_ICRF2_MAX_AGE_DAYS = 365  # ICRF2 is superseded and will not be updated again


def get_icrf2(
    cache_path: Path | None = None,
    max_age_days: int | None = _ICRF2_MAX_AGE_DAYS,
) -> table.Table:
    """Return the ICRF2 catalog, downloading/caching from VizieR as needed.

    The catalog (3414 sources) is cached at `cache_path` and refreshed when
    absent or older than `max_age_days`.

    Parameters
    ----------
    cache_path
        Override the default cache location (default: ``$ASTROMON_DATA_DIR/icrf2_catalog.ecsv``).
    max_age_days
        Refresh threshold in days.  Default 365.  Pass ``None`` to download
        only when the file is absent.

    Returns
    -------
    table.Table
        Table with columns ``name`` (str), ``ra`` (deg), ``dec`` (deg).

    Examples
    --------
    >>> t = get_icrf2()
    >>> len(t) > 3000
    True
    """
    if cache_path is None:
        cache_path = _ICRF2_CACHE_PATH
    cache_path = Path(cache_path)
    if _cache_is_stale(cache_path, max_age_days):
        logger.info(f"Downloading ICRF2 from VizieR {_ICRF2_VIZIER_ID} → {cache_path}")
        vizier = Vizier(row_limit=-1, columns=["ICRF", "_RAJ2000", "_DEJ2000"])
        result = vizier.get_catalogs(_ICRF2_VIZIER_ID)
        if not result:
            raise RuntimeError(
                f"VizieR returned no data for ICRF2 catalog {_ICRF2_VIZIER_ID}"
            )
        raw = result[0]
        icrf2 = table.Table(
            {
                "name": np.asarray(raw["ICRF"], dtype="U32"),
                "ra": np.asarray(raw["_RAJ2000"], dtype=float),
                "dec": np.asarray(raw["_DEJ2000"], dtype=float),
            }
        )
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        icrf2.write(cache_path, format="ascii.ecsv", overwrite=True)
        logger.info(f"  Saved {len(icrf2)} ICRF2 sources to {cache_path}")
    else:
        logger.debug(f"Loading cached ICRF2 from {cache_path}")
        icrf2 = _read_catalog_cached(
            cache_path, lambda p: table.Table.read(p, format="ascii.ecsv")
        )
    return icrf2


def _get_icrf2(
    pos: coords.SkyCoord, radius: u.Quantity, logging_tag: str = ""
) -> table.Table:
    """Return ICRF2 sources near *pos* within *radius*, formatted like VizieR results."""
    return _local_catalog_near(get_icrf2(), "ICRF2", pos, radius, logging_tag)


_MILLIQUAS_CACHE_PATH = _ASTROMON_DATA_DIR / "milliquas_catalog.fits"
_MILLIQUAS_MAX_AGE_DAYS = 180  # update twice a year; v8 released 2023-07-30


def get_milliquas(
    cache_path: Path | None = None,
    max_age_days: int | None = _MILLIQUAS_MAX_AGE_DAYS,
) -> table.Table:
    """Return the Milliquas v8 spectroscopic AGN catalog, downloading and caching as needed.

    Filters to spectroscopically-confirmed types (Q, A, B, N) — the same subset
    used by :func:`get_milliquas_gaia`.  Cached at ``cache_path`` as a FITS binary
    table and refreshed when absent or older than ``max_age_days``.

    Parameters
    ----------
    cache_path
        Override the default cache location
        (default: ``$ASTROMON_DATA_DIR/milliquas_catalog.fits``).
    max_age_days
        Refresh threshold in days.  Default 180.  Pass ``None`` to download
        only when the file is absent.

    Returns
    -------
    table.Table
        Table with columns ``name`` (str), ``ra`` (deg), ``dec`` (deg).

    Examples
    --------
    >>> t = get_milliquas()
    >>> len(t) > 800_000
    True
    """
    import gzip  # noqa: PLC0415 — standard library, import deferred to avoid startup cost

    if cache_path is None:
        cache_path = _MILLIQUAS_CACHE_PATH
    cache_path = Path(cache_path)

    if not _cache_is_stale(cache_path, max_age_days):
        logger.debug(f"Loading cached Milliquas from {cache_path}")
        return _read_catalog_cached(
            cache_path, lambda p: table.Table.read(p, format="fits")
        )

    logger.info(
        f"Downloading Milliquas v8 from CDS → {cache_path} (~80 MB compressed, one-time)"
    )
    url = "https://cdsarc.cds.unistra.fr/ftp/VII/294/catalog.dat.gz"
    resp = requests.get(url, timeout=600, stream=True)
    resp.raise_for_status()

    raw_bytes = b"".join(resp.iter_content(chunk_size=65536))
    raw_text = gzip.decompress(raw_bytes).decode("latin-1")

    # Parse fixed-width format per the CDS ReadMe:
    #   cols  1-11 (0-indexed  0-10): RAdeg  F11.7 deg
    #   cols 13-23 (0-indexed 12-22): DEdeg  F11.7 deg
    #   cols 26-50 (0-indexed 25-49): Name   A25
    #   cols 52-55 (0-indexed 51-54): Type   A4
    ra_vals: list[float] = []
    dec_vals: list[float] = []
    name_vals: list[str] = []

    for line in raw_text.splitlines():
        if len(line) < 55:
            continue
        type_str = line[51:55].strip()
        if not type_str or type_str[0] not in "QABN":
            continue
        try:
            ra_vals.append(float(line[0:11]))
            dec_vals.append(float(line[12:23]))
            name_vals.append(line[25:50].strip())
        except ValueError:
            continue

    mq = table.Table(
        {
            "name": np.asarray(name_vals, dtype="U25"),
            "ra": np.asarray(ra_vals, dtype=float),
            "dec": np.asarray(dec_vals, dtype=float),
        }
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    mq.write(cache_path, format="fits", overwrite=True)
    logger.info(f"Saved {len(mq)} Milliquas spectroscopic sources to {cache_path}")
    return mq


# Quaia G<20.0: Gaia DR3 + unWISE quasar catalog (Storey-Fisher et al. 2024, ApJ 964, 69).
# v1.0.0 released 2023-12-18; tied to Gaia DR3, no update expected before Gaia DR4.
_QUAIA_URL = "https://zenodo.org/records/10403370/files/quaia_G20.0.fits"
_QUAIA_CACHE_PATH = _ASTROMON_DATA_DIR / "quaia_catalog.fits"
# None = download only when absent; the catalog is static until Gaia DR4.
_QUAIA_MAX_AGE_DAYS: int | None = None

# Local FITS caches for the two Gaia DR3 AGN/QSO catalog tables.  Both are
# static until Gaia DR4, so no automatic refresh is needed.
_GAIA_AGN_CACHE_PATH = _ASTROMON_DATA_DIR / "gaia_agn_catalog.fits"
_GAIA_QSO_CACHE_PATH = _ASTROMON_DATA_DIR / "gaia_qso_catalog.fits"


def get_quaia(
    cache_path: Path | None = None,
    max_age_days: int | None = _QUAIA_MAX_AGE_DAYS,
) -> table.Table:
    """Return the Quaia G<20.0 Gaia-unWISE quasar catalog, downloading/caching as needed.

    Quaia (Storey-Fisher et al. 2024, ApJ 964, 69) is built from Gaia DR3 quasar
    candidates and unWISE infrared photometry.  The G<20.0 sample (755,850 sources)
    is used: it is a clean subset of the full G<20.5 catalog with tighter magnitude
    constraints and better-characterised redshifts.

    The catalog is cached at ``cache_path`` as a FITS binary table containing only
    the three columns needed for cross-matching (``name``, ``ra``, ``dec``).  It
    is never automatically refreshed because the catalog is tied to Gaia DR3 and
    will not change until Gaia DR4.

    Parameters
    ----------
    cache_path
        Override the default cache location
        (default: ``$ASTROMON_DATA_DIR/quaia_catalog.fits``).
    max_age_days
        Refresh threshold in days.  Default ``None`` (never auto-refresh).

    Returns
    -------
    table.Table
        Table with columns ``name`` (str), ``ra`` (deg), ``dec`` (deg).

    Examples
    --------
    >>> t = get_quaia()
    >>> len(t) > 700_000
    True
    """
    import tempfile  # noqa: PLC0415

    if cache_path is None:
        cache_path = _QUAIA_CACHE_PATH
    cache_path = Path(cache_path)

    if not _cache_is_stale(cache_path, max_age_days):
        logger.debug(f"Loading cached Quaia from {cache_path}")
        return _read_catalog_cached(
            cache_path, lambda p: table.Table.read(p, format="fits")
        )

    logger.info(
        f"Downloading Quaia G<20.0 from Zenodo → {cache_path} (~100 MB, one-time)"
    )
    resp = requests.get(_QUAIA_URL, timeout=600, stream=True)
    resp.raise_for_status()

    with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        for chunk in resp.iter_content(chunk_size=65536):
            tmp.write(chunk)

    try:
        raw = table.Table.read(tmp_path, format="fits")
        # Strip to the three columns needed; name uses the Gaia source_id.
        quaia = table.Table(
            {
                "name": np.asarray(
                    [f"Gaia DR3 {sid}" for sid in raw["source_id"]], dtype="<U30"
                ),
                "ra": np.asarray(raw["ra"], dtype=float),
                "dec": np.asarray(raw["dec"], dtype=float),
            }
        )
    finally:
        tmp_path.unlink(missing_ok=True)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    quaia.write(cache_path, format="fits", overwrite=True)
    logger.info(f"Saved {len(quaia):,} Quaia G<20.0 sources to {cache_path}")
    return quaia


def get_quaia_candidates(
    pos: coords.SkyCoord,
    radius: u.Quantity = 3 * u.arcsec,
    logging_tag: str = "",
) -> table.Table:
    """Return Quaia sources near *pos* within *radius*.

    Loads the locally cached Quaia catalog (downloading it on first call via
    :func:`get_quaia`) and applies a bounding-cone pre-filter before the exact
    per-source separation check.  No network call is made after the initial download.

    Parameters
    ----------
    pos
        Positions to search around (scalar or array SkyCoord).
    radius
        Search radius.  Default 3 arcsec.
    logging_tag
        Prefix for log messages.

    Returns
    -------
    table.Table
        Catalog sources with columns matching :data:`CROSS_MATCH_DTYPE`.
    """
    if not pos.isscalar and len(pos) == 0:
        # Guard before get_quaia(), which triggers a ~100 MB one-time download on
        # first call: an empty position array should never pay that cost.
        return table.Table(dtype=CROSS_MATCH_DTYPE)
    return _local_catalog_near(get_quaia(), "Quaia", pos, radius, logging_tag)


VIZIER_CATALOGS = {
    "Tycho2": {
        "catalog": "Tycho2",
        "cat_identifier": "I/259/tyc2",
        "name_cols": ["TYC1", "TYC2", "TYC3"],
        "columns": {
            "ra": "_RAJ2000/{time.frac_year:.3f}",
            "dec": "_DEJ2000/{time.frac_year:.3f}",
            "mag": "VTmag",
        },
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
            "mag": "Hpmag",
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
        "name_cols": ["2MASS"],
        "columns": {"ra": "_RAJ2000", "dec": "_DEJ2000", "mag": "Kmag"},
    },
    "SDSS": {
        "catalog": "SDSS",
        "cat_identifier": "II/294",
        "name_cols": ["SDSS"],
        "columns": {
            "ra": "_RAJ2000/{time.frac_year:.3f}",
            "dec": "_DEJ2000/{time.frac_year:.3f}",
            "mag": "rmag",
        },
    },
    # https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/345/gaia2
    "Gaia2": {
        "catalog": "Gaia2",
        "cat_identifier": "I/345/gaia2",
        "name_cols": ["Source"],
        "columns": {"ra": "_RAJ2000", "dec": "_DEJ2000", "mag": "Gmag"},
    },
}


@retry(
    exceptions=(
        requests.exceptions.HTTPError,
        requests.exceptions.Timeout,
        requests.exceptions.ConnectionError,
        _CallTimeoutError,
    ),
    tries=3,
    delay=10,
    backoff=2,
    max_delay=120,
    logger=logger,
)
def _execute_gaia_tap_query(query):
    """Execute a single Gaia TAP query with retry and hard timeout."""
    # Imported here rather than at module scope: astroquery instantiates GaiaClass() at
    # import time, which contacts gea.esac.esa.int for status messages with no timeout.
    # At module scope that made "import astromon" (and so pytest collection) hang whenever
    # the ESA server was slow or unresponsive.
    from astroquery.gaia import Gaia  # noqa: PLC0415

    def _query():
        job = Gaia.launch_job(query, verbose=False, upload_resource=None)
        return job.get_results()

    return _run_with_timeout(_query, timeout=120)


def _apply_proper_motion(ra, dec, pmra, pmdec, obs_epoch_yr, catalog_epoch_yr=2016.0):
    """Propagate (ra, dec) from catalog_epoch_yr to obs_epoch_yr using proper motion.

    Parameters
    ----------
    ra, dec : array-like, degrees
    pmra : array-like, mas/yr (mu_alpha * cos(delta) — Gaia convention)
    pmdec : array-like, mas/yr
    obs_epoch_yr : float, target epoch in Julian years
    catalog_epoch_yr : float, catalog reference epoch (Gaia DR3 = J2016.0)

    Returns
    -------
    ra_corr, dec_corr : ndarray, degrees
    """
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    pmra = np.where(np.isfinite(pmra), pmra, 0.0)
    pmdec = np.where(np.isfinite(pmdec), pmdec, 0.0)
    delta_yr = obs_epoch_yr - catalog_epoch_yr
    mas_to_deg = 1.0 / 3_600_000.0
    ra_corr = ra + pmra * delta_yr * mas_to_deg / np.cos(np.radians(dec))
    dec_corr = dec + pmdec * delta_yr * mas_to_deg
    return ra_corr, dec_corr


def _download_gaia_catalog_to_cache(
    adql_query: str,
    cache_path: Path,
    label: str,
    name_prefix: str,
) -> table.Table:
    """Download a full Gaia DR3 table via a single async TAP query and cache it as FITS.

    Stores four columns — name, ra, dec, mag — enough for cross-matching.
    The ``name`` column is built as ``"{name_prefix}-{source_id}"``.

    For tables exceeding the 3M-row Gaia archive cap, use
    :func:`_download_gaia_qso_catalog_chunked` instead.

    Parameters
    ----------
    adql_query : str
        ADQL SELECT returning source_id, ra, dec, phot_g_mean_mag.
    cache_path : Path
        Destination FITS file.
    label : str
        Human-readable name used in log messages.
    name_prefix : str
        Prepended to source_id in the ``name`` column.
    """
    from astroquery.gaia import Gaia  # noqa: PLC0415

    logger.info(
        f"Downloading {label} from Gaia archive → {cache_path} (one-time, may take several minutes)"
    )
    job = Gaia.launch_job_async(adql_query, verbose=False)
    raw = job.get_results()

    mag_col = np.asarray(raw["phot_g_mean_mag"], dtype=np.float32)
    cat = table.Table(
        {
            "name": np.asarray(
                [f"{name_prefix}-{sid}" for sid in raw["source_id"]], dtype="<U32"
            ),
            "ra": np.asarray(raw["ra"], dtype=float),
            "dec": np.asarray(raw["dec"], dtype=float),
            "mag": mag_col,
        }
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cat.write(cache_path, format="fits", overwrite=True)
    logger.info(f"Saved {len(cat):,} {label} sources to {cache_path}")
    return cat


def _download_gaia_qso_catalog_chunked(cache_path: Path) -> table.Table:
    """Download the full gaiadr3.qso_candidates catalog in source_id chunks.

    The Gaia archive hard-caps async results at 3M rows regardless of ROW_LIMIT.
    qso_candidates has ~6.65M rows, so we split into five non-overlapping
    source_id ranges that each stay well under the cap, then concatenate.

    Parameters
    ----------
    cache_path : Path
        Destination FITS file.
    """
    from astroquery.gaia import Gaia  # noqa: PLC0415

    # Each tuple is (lo, hi) for the WHERE clause on g.source_id.
    # Row counts verified against the live archive:
    #   chunk 1: 1,840,645  chunk 2: 720,322  chunk 3: 2,763,765
    #   chunk 4:   461,403  chunk 5: 863,027  total: 6,649,162
    chunks: list[tuple[int | None, int | None]] = [
        (None, 3_500_000_000_000_000_000),
        (3_500_000_000_000_000_000, 4_500_000_000_000_000_000),
        (4_500_000_000_000_000_000, 5_000_000_000_000_000_000),
        (5_000_000_000_000_000_000, 5_500_000_000_000_000_000),
        (5_500_000_000_000_000_000, None),
    ]

    logger.info(
        f"Downloading GaiaQSO (qso_candidates) in {len(chunks)} chunks → {cache_path}"
    )
    raw_chunks = []
    for idx, (lo, hi) in enumerate(chunks):
        if lo is not None and hi is not None:
            filter_clause = f"WHERE g.source_id >= {lo} AND g.source_id < {hi}"
        elif lo is not None:
            filter_clause = f"WHERE g.source_id >= {lo}"
        elif hi is not None:
            filter_clause = f"WHERE g.source_id < {hi}"
        else:
            filter_clause = ""

        query = f"""
            SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag
            FROM gaiadr3.gaia_source g
            INNER JOIN gaiadr3.qso_candidates q ON g.source_id = q.source_id
            {filter_clause}
        """
        logger.info(
            f"  chunk {idx + 1}/{len(chunks)}: {filter_clause or '(no filter)'}"
        )
        job = Gaia.launch_job_async(query, verbose=False)
        chunk = job.get_results()
        logger.info(f"  → {len(chunk):,} rows")
        raw_chunks.append(chunk)

    raw = table.vstack(raw_chunks)
    logger.info(f"  total before dedup: {len(raw):,} rows")

    # Deduplicate in case chunk boundaries produced any overlap.
    _, keep = np.unique(raw["source_id"], return_index=True)
    if len(keep) < len(raw):
        logger.info(f"  deduplicated {len(raw) - len(keep)} duplicate rows")
        raw = raw[sorted(keep)]

    mag_col = np.asarray(raw["phot_g_mean_mag"], dtype=np.float32)
    cat = table.Table(
        {
            "name": np.asarray(
                [f"GaiaQSO-{sid}" for sid in raw["source_id"]], dtype="<U32"
            ),
            "ra": np.asarray(raw["ra"], dtype=float),
            "dec": np.asarray(raw["dec"], dtype=float),
            "mag": mag_col,
        }
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cat.write(cache_path, format="fits", overwrite=True)
    logger.info(f"Saved {len(cat):,} GaiaQSO sources to {cache_path}")
    return cat


def get_gaia_agn_catalog(cache_path: Path | None = None) -> table.Table:
    """Return the full Gaia DR3 AGN cross-ID catalog, downloading and caching as needed.

    Downloads ``gaiadr3.agn_cross_id`` (~1.6 M sources) on first call and caches
    the result at ``cache_path`` as a FITS binary table with columns
    ``name``, ``ra``, ``dec``, ``mag``.  No automatic refresh is performed
    because the catalog is tied to Gaia DR3 and will not change before Gaia DR4.

    Parameters
    ----------
    cache_path
        Override the default cache location
        (default: ``$ASTROMON_DATA_DIR/gaia_agn_catalog.fits``).
    """
    if cache_path is None:
        cache_path = _GAIA_AGN_CACHE_PATH
    cache_path = Path(cache_path)
    if cache_path.exists():
        logger.debug(f"Loading cached GaiaAGN catalog from {cache_path}")
        return _read_catalog_cached(
            cache_path, lambda p: table.Table.read(p, format="fits")
        )

    return _download_gaia_catalog_to_cache(
        adql_query="""
            SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag
            FROM gaiadr3.gaia_source g
            INNER JOIN gaiadr3.agn_cross_id a ON g.source_id = a.source_id
        """,
        cache_path=cache_path,
        label="GaiaAGN (agn_cross_id)",
        name_prefix="GaiaAGN",
    )


def get_gaia_qso_catalog(cache_path: Path | None = None) -> table.Table:
    """Return the full Gaia DR3 QSO candidates catalog, downloading and caching as needed.

    Downloads ``gaiadr3.qso_candidates`` (~6.6 M sources) on first call and caches
    the result at ``cache_path`` as a FITS binary table with columns
    ``name``, ``ra``, ``dec``, ``mag``.  No automatic refresh is performed
    because the catalog is tied to Gaia DR3 and will not change before Gaia DR4.

    Parameters
    ----------
    cache_path
        Override the default cache location
        (default: ``$ASTROMON_DATA_DIR/gaia_qso_catalog.fits``).
    """
    if cache_path is None:
        cache_path = _GAIA_QSO_CACHE_PATH
    cache_path = Path(cache_path)
    if cache_path.exists():
        logger.debug(f"Loading cached GaiaQSO catalog from {cache_path}")
        return _read_catalog_cached(
            cache_path, lambda p: table.Table.read(p, format="fits")
        )

    return _download_gaia_qso_catalog_chunked(cache_path)


def get_gaia_agn(
    pos: coords.SkyCoord,
    radius: u.Quantity = 3 * u.arcsec,
    logging_tag: str = "",
) -> table.Table:
    """Return Gaia DR3 AGN cross-ID sources near *pos* within *radius*.

    Loads the locally cached ``gaiadr3.agn_cross_id`` catalog (~1.6 M sources,
    downloaded on first call via :func:`get_gaia_agn_catalog`) and applies a
    bounding-cone pre-filter before the exact per-source separation check.
    No network call is made after the initial one-time download.

    Parameters
    ----------
    pos : astropy.coordinates.SkyCoord
        Positions to search around (array or scalar).
    radius : astropy.units.Quantity
        Search radius. Default 3 arcsec.
    logging_tag : str
        Prefix for log messages.

    Returns
    -------
    astropy.table.Table
        Catalog sources with columns matching CROSS_MATCH_DTYPE.
    """
    if not pos.isscalar and len(pos) == 0:
        logger.debug(f"{logging_tag}GaiaAGN has 0 results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)
    try:
        catalog = get_gaia_agn_catalog()
    except Exception as exc:
        logger.warning(f"{logging_tag}GaiaAGN catalog unavailable: {exc}")
        return table.Table(dtype=CROSS_MATCH_DTYPE)
    return _local_catalog_near(catalog, "GaiaAGN", pos, radius, logging_tag)


def get_gaia_qso_candidates(
    pos: coords.SkyCoord,
    radius: u.Quantity = 3 * u.arcsec,
    logging_tag: str = "",
) -> table.Table:
    """Return Gaia DR3 QSO candidate sources near *pos* within *radius*.

    Uses ``gaiadr3.qso_candidates`` (~6.6 M sources), which extends the
    confirmed cross-identifications in ``gaiadr3.agn_cross_id`` (~1.6 M, used
    by :func:`get_gaia_agn`) with an additional ~5 M sources classified from
    Gaia astrometry, photometry, and colour.  The extra depth is most useful
    at G > 20, covering faint AGN in extragalactic fields that Quaia (G < 20)
    and agn_cross_id miss.

    Loads the locally cached catalog (downloading it on first call via
    :func:`get_gaia_qso_catalog`) and applies a bounding-cone pre-filter
    before the exact per-source separation check.  No network call is made
    after the initial one-time download.

    Parameters
    ----------
    pos : astropy.coordinates.SkyCoord
        Positions to search around (array or scalar).
    radius : astropy.units.Quantity
        Search radius. Default 3 arcsec.
    logging_tag : str
        Prefix for log messages.

    Returns
    -------
    astropy.table.Table
        Catalog sources with columns matching CROSS_MATCH_DTYPE.
    """
    if not pos.isscalar and len(pos) == 0:
        logger.debug(f"{logging_tag}GaiaQSO has 0 results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)
    try:
        catalog = get_gaia_qso_catalog()
    except Exception as exc:
        logger.warning(f"{logging_tag}GaiaQSO catalog unavailable: {exc}")
        return table.Table(dtype=CROSS_MATCH_DTYPE)
    return _local_catalog_near(catalog, "GaiaQSO", pos, radius, logging_tag)


def get_gaia_var_stars(
    pos: coords.SkyCoord,
    radius: u.Quantity = 3 * u.arcsec,
    max_ruwe: float = 1.4,
    batch_size: int = 50,
    logging_tag: str = "",
) -> table.Table:
    """Query Gaia DR3 rotation-modulation variable stars around given sky positions via TAP.

    Positions are corrected from the Gaia DR3 reference epoch (J2016.0) to the
    observation epoch using each source's proper motion before returning.

    Sources with poor astrometric solutions (RUWE >= max_ruwe) are excluded so
    that unreliable proper motions do not degrade the position correction.

    Parameters
    ----------
    pos : astropy.coordinates.SkyCoord
        Positions to search around (array or scalar). ``pos.obstime`` is used
        as the target epoch for proper motion correction; if not set, no
        correction is applied and a warning is logged.
    radius : astropy.units.Quantity
        Search radius. Default 3 arcsec.
    max_ruwe : float
        Maximum allowed RUWE (renormalised unit weight error). Sources with
        RUWE >= max_ruwe are excluded. Default 1.4 (standard Gaia quality cut).
    batch_size : int
        Number of positions per TAP CONTAINS-OR query.  Default 50.  Increase
        to 200–500 when querying many positions at once (bulk screening) to
        reduce round-trip count; lower if the TAP server times out.
    logging_tag : str
        Prefix for log messages.

    Returns
    -------
    astropy.table.Table
        Catalog sources with columns matching CROSS_MATCH_DTYPE, positions
        propagated to the observation epoch.
    """
    if not pos.isscalar and len(pos) == 0:
        logger.debug(f"{logging_tag} GaiaVarStar has no results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    if pos.obstime is not None:
        obs_epoch_yr = (
            pos.obstime.jyear
            if hasattr(pos.obstime, "jyear")
            else float(pos.obstime.decimalyear)
        )
    else:
        logger.warning(
            f"{logging_tag} GaiaVarStar: pos.obstime not set, skipping proper motion correction"
        )
        obs_epoch_yr = None

    radius_deg = radius.to(u.deg).value
    all_results = []

    pos_ra_deg = np.atleast_1d(pos.ra.deg)
    pos_dec_deg = np.atleast_1d(pos.dec.deg)
    for batch_start in range(0, len(pos_ra_deg), batch_size):
        batch_ra = pos_ra_deg[batch_start : batch_start + batch_size]
        batch_dec = pos_dec_deg[batch_start : batch_start + batch_size]
        where_clauses = [
            f"CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', {ra}, {dec}, {radius_deg}))=1"
            for ra, dec in zip(batch_ra, batch_dec, strict=True)
        ]
        query = f"""
        SELECT g.source_id, g.ra, g.dec, g.pmra, g.pmdec, g.phot_g_mean_mag
        FROM gaiadr3.gaia_source g
        INNER JOIN gaiadr3.vari_rotation_modulation v ON g.source_id = v.source_id
        WHERE g.ruwe < {max_ruwe}
        AND ({" OR ".join(where_clauses)})
        """
        try:
            result = _execute_gaia_tap_query(query)
            if result is not None and len(result) > 0:
                all_results.append(result)
            logger.debug(
                f"{logging_tag} GaiaVarStar batch got {len(result) if result else 0} sources"
            )
        except Exception as e:
            logger.warning(
                f"{logging_tag} GaiaVarStar batch failed after all retries: {e}"
            )
            continue

    if not all_results:
        logger.debug(f"{logging_tag} GaiaVarStar has no results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    combined = table.vstack(all_results)
    _, unique_idx = np.unique(combined["source_id"], return_index=True)
    combined = combined[sorted(unique_idx)]

    ra = np.array(combined["ra"], dtype=float)
    dec = np.array(combined["dec"], dtype=float)
    pmra = np.array(combined["pmra"], dtype=float)
    pmdec = np.array(combined["pmdec"], dtype=float)

    if obs_epoch_yr is not None:
        ra, dec = _apply_proper_motion(ra, dec, pmra, pmdec, obs_epoch_yr)

    combined["catalog"] = "GaiaVarStar"
    combined["name"] = [f"GaiaVarStar-{sid}" for sid in combined["source_id"]]
    # The query always selects phot_g_mean_mag, so the column is always present.
    combined["mag"] = combined["phot_g_mean_mag"].astype(np.float32)

    output = table.Table(data=np.zeros(len(combined), dtype=CROSS_MATCH_DTYPE))
    output["ra"] = ra
    output["dec"] = dec
    for col in ("catalog", "name", "mag"):
        output[col] = combined[col]
    logger.debug(f"{logging_tag} GaiaVarStar has {len(output)} results")
    return output


# Radius used to positionally match Milliquas optical positions to Gaia DR3 sources.
# Milliquas positions are typically accurate to ~0.3 arcsec; 1.5 arcsec gives comfortable margin.
_MILLIQUAS_GAIA_XMATCH_ARCSEC = 1.5
_MILLIQUAS_GAIA_XMATCH_DEG = _MILLIQUAS_GAIA_XMATCH_ARCSEC / 3600.0


def get_milliquas_gaia(  # noqa: PLR0915
    pos, radius=3 * u.arcsec, logging_tag=""
):
    """Query Milliquas v8 (VII/294) + Gaia DR3 for spectroscopic AGN around sky positions.

    Cone-searches Milliquas v8 around each position in `pos`, filters to
    spectroscopically-confirmed types (Q, A, B, N), then requires a Gaia DR3
    positional match within 1.5 arcsec.  Sources without a Gaia DR3 counterpart
    are dropped.  Photometric candidates (K-type) and the ~66,000 radio/X-ray
    association candidates are excluded by the type filter.

    Parameters
    ----------
    pos : astropy.coordinates.SkyCoord
        Positions to search around (array or scalar).
    radius : astropy.units.Quantity
        Search radius. Default 3 arcsec.
    logging_tag : str
        Prefix for log messages.

    Returns
    -------
    astropy.table.Table
        Catalog sources with columns matching CROSS_MATCH_DTYPE.
    """
    if not pos.isscalar and len(pos) == 0:
        logger.debug(f"{logging_tag} MilliquasGaia has no results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    # Step 1: VizieR cone search on Milliquas v8 (VII/294/catalog).
    vizier = Vizier(
        columns=["RAJ2000", "DEJ2000", "Name", "Type", "Rmag"], row_limit=-1
    )
    try:
        mq_result = _query_vizier_region(vizier, pos, radius, "VII/294/catalog")
        mq_tables = list(mq_result)
    except Exception as e:
        logger.warning(f"{logging_tag} Milliquas VizieR query failed: {e}")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    if not mq_tables:
        logger.debug(f"{logging_tag} MilliquasGaia has no results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    mq = table.vstack(mq_tables, metadata_conflicts="silent")

    # Deduplicate by Name so the same Milliquas source only appears once even if it
    # fell within radius of multiple X-ray positions.
    _, unique_idx = np.unique(np.array(mq["Name"], dtype=str), return_index=True)
    mq = mq[sorted(unique_idx)]

    # Step 2: Keep only spectroscopic types (first char Q, A, B, or N).
    # This excludes K-type photometric candidates and all radio/X-ray association entries.
    type_strs = [str(t).strip() for t in mq["Type"]]
    type_mask = np.array([bool(t) and t[0] in "QABN" for t in type_strs])
    mq = mq[type_mask]

    if len(mq) == 0:
        logger.debug(f"{logging_tag} MilliquasGaia has no spectroscopic-type results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    # Step 3: Upgrade positions to Gaia DR3 for sources with a Gaia counterpart.
    # Batch query Gaia DR3 with circles around each filtered Milliquas source,
    # then assign each Gaia result to the nearest Milliquas source within 1.5 arcsec.
    mq_ra = np.array(mq["RAJ2000"], dtype=float)
    mq_dec = np.array(mq["DEJ2000"], dtype=float)
    final_ra = mq_ra.copy()
    final_dec = mq_dec.copy()
    gaia_matched_dist = np.full(len(mq), np.inf)

    batch_size = 50
    for batch_start in range(0, len(mq), batch_size):
        batch_end = batch_start + batch_size
        b_ra = mq_ra[batch_start:batch_end]
        b_dec = mq_dec[batch_start:batch_end]
        circles = " OR ".join(
            f"CONTAINS(POINT('ICRS', g.ra, g.dec),"
            f" CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {_MILLIQUAS_GAIA_XMATCH_DEG:.8f}))=1"
            for ra, dec in zip(b_ra, b_dec, strict=True)
        )
        query = f"""
        SELECT g.source_id, g.ra, g.dec
        FROM gaiadr3.gaia_source g
        WHERE {circles}
        """
        try:
            gaia_batch = _execute_gaia_tap_query(query)
        except Exception as e:
            logger.warning(f"{logging_tag} MilliquasGaia Gaia TAP batch failed: {e}")
            continue

        if gaia_batch is None or len(gaia_batch) == 0:
            continue

        gaia_sc = coords.SkyCoord(
            np.array(gaia_batch["ra"], dtype=float),
            np.array(gaia_batch["dec"], dtype=float),
            unit="deg",
        )
        mq_batch_sc = coords.SkyCoord(b_ra, b_dec, unit="deg")

        # Vectorized nearest-Milliquas match for all Gaia rows at once (each
        # scalar SkyCoord index plus per-row separation() is ~1000x slower).
        nearest_idx, sep2d, _ = coords.match_coordinates_sky(gaia_sc, mq_batch_sc)
        for gi, (nearest_i, sep_arcsec) in enumerate(
            zip(np.atleast_1d(nearest_idx), np.atleast_1d(sep2d.arcsec), strict=True)
        ):
            nearest = int(nearest_i)
            if (
                sep_arcsec < _MILLIQUAS_GAIA_XMATCH_ARCSEC
                and sep_arcsec < gaia_matched_dist[batch_start + nearest]
            ):
                gaia_matched_dist[batch_start + nearest] = sep_arcsec
                final_ra[batch_start + nearest] = float(gaia_batch["ra"][gi])
                final_dec[batch_start + nearest] = float(gaia_batch["dec"][gi])

    # Step 4: Keep only sources with a confirmed Gaia DR3 position.
    gaia_mask = np.isfinite(gaia_matched_dist)
    n_dropped = int((~gaia_mask).sum())
    if n_dropped:
        logger.debug(
            f"{logging_tag} MilliquasGaia: dropping {n_dropped} sources "
            "without Gaia DR3 counterpart"
        )
    mq = mq[gaia_mask]
    final_ra = final_ra[gaia_mask]
    final_dec = final_dec[gaia_mask]

    if len(mq) == 0:
        logger.debug(
            f"{logging_tag} MilliquasGaia has no results with Gaia DR3 positions"
        )
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    rmag = mq["Rmag"]
    rmag_arr = (
        rmag.filled(np.nan).astype(np.float32)
        if hasattr(rmag, "filled")
        else np.array(rmag, dtype=np.float32)
    )

    output = table.Table(data=np.zeros(len(mq), dtype=CROSS_MATCH_DTYPE))
    output["ra"] = final_ra
    output["dec"] = final_dec
    output["catalog"] = "MilliquasGaia"
    output["name"] = [str(n).strip()[:32] for n in mq["Name"]]
    output["mag"] = rmag_arr

    logger.debug(
        f"{logging_tag} MilliquasGaia: {len(output)} results with Gaia DR3 positions"
    )
    return output


def _dedupe_desi_by_delchi2(desi: table.Table) -> table.Table:
    """One row per physical object, keeping the highest-delChi2 observation.

    DESI's V/161 spectroscopic catalog spans multiple survey programs
    (backup/bright/dark), and the same physical target is often observed more
    than once across them -- sometimes under the *same* TargetID (a repeat
    observation in a different program), sometimes under a *different* TargetID
    entirely (the same Legacy Survey DR9 imaging position re-targeted in a later
    survey phase). Neither case is a second real source.

    delChi2 -- redrock's own margin between the best-fit redshift and the
    runner-up -- separates them cleanly: measured duplicate pairs in the live
    database differ by 3x-20x, never a close call, and the higher-delChi2
    observation always has the longer exposure time too, so it is not a noisy
    proxy. Applied twice: first within a TargetID's own repeat observations,
    then across different TargetIDs that land on the same position.

    Must run after the OType/ZWARN filter, not before: a TargetID's better
    observation can have ZWARN=0 while a worse one for the *same* TargetID has
    ZWARN != 0, and deduplicating first risked keeping the bad one by table
    order and then filtering out the TargetID entirely.
    """
    if len(desi) == 0:
        return desi

    def _keep_best(t: table.Table, key: np.ndarray) -> table.Table:
        order = np.argsort(key, kind="stable")
        sorted_key = key[order]
        group_start = np.concatenate(
            ([0], np.flatnonzero(np.diff(sorted_key) != 0) + 1)
        )
        group_end = np.concatenate((group_start[1:], [len(order)]))
        delchi2 = np.asarray(t["delChi2"])[order]
        best_in_group = [
            start + np.argmax(delchi2[start:end])
            for start, end in zip(group_start, group_end, strict=True)
        ]
        return t[sorted(order[best_in_group])]

    desi = _keep_best(desi, np.asarray(desi["TargetID"]))
    desi = _keep_best(desi, _group_by_position(desi["RAICRS"], desi["DEICRS"]))
    return desi


# DESI's own documented positional accuracy (see get_desi_v161_candidates'
# docstring): repeat observations of one object can differ by this much even
# though they share the same underlying imaging source, so anything closer
# than this is the same object, not two coincidentally close ones.
_DESI_POSITION_DUPLICATE_ARCSEC = 0.1


def _group_by_position(
    ra: np.ndarray,
    dec: np.ndarray,
    tolerance_arcsec: float = _DESI_POSITION_DUPLICATE_ARCSEC,
) -> np.ndarray:
    """A group id per row: connected components of "within tolerance of".

    Not round-then-group: DESI's repeat targetings of one object share the
    identical Legacy Survey DR9 imaging position to a fraction of an arcsec, but
    a fraction of an arcsec is exactly small enough to straddle a rounding
    boundary -- 190.544187 and 190.544196 are 0.033" apart, inside a 0.05"
    tolerance, but round to different values at 5 decimal places (...4419 vs
    ...4420) and so land in different buckets under naive rounding. This
    happened for real: 12 of 122 known duplicates survived the first version of
    this fix for exactly that reason.

    Union-find on pairwise separation has no boundary to straddle. Quadratic in
    the number of rows, which is fine here -- one query's candidate list is a
    handful to a few hundred rows, not the whole database.
    """
    ra, dec = np.asarray(ra), np.asarray(dec)
    n = len(ra)
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        dra = (ra - ra[i]) * np.cos(np.radians(dec[i]))
        ddec = dec - dec[i]
        close = np.flatnonzero(np.hypot(dra, ddec) * 3600 < tolerance_arcsec)
        for j in close:
            ri, rj = find(i), find(int(j))
            if ri != rj:
                parent[ri] = rj

    return np.array([find(i) for i in range(n)])


def get_desi_v161_candidates(
    pos: coords.SkyCoord,
    radius: u.Quantity = 3 * u.arcsec,
    logging_tag: str = "",
) -> table.Table:
    """Query DESI EDR spectroscopic catalog (VizieR V/161) for QSOs and galaxies.

    Cone-searches the DESI Early Data Release spectroscopic catalog around each
    position in *pos*, keeping only targets with:

    * ``OType`` in ``{'QSO', 'GALAXY'}`` — spectroscopically confirmed quasars
      or galaxies (excludes STAR-typed objects and blank types).
    * ``ZWARN == 0`` — reliable redshift / spectral-type classification; no
      known quality issues with the fiber or template fit.

    DESI positions come from Legacy Survey DR9 imaging and are accurate to
    ≲ 0.1 arcsec for point-like sources — comparable to Gaia, so no position
    upgrade step is needed.  The DESI EDR footprint (~4100 sq deg in discrete
    survey-validation tiles) is the primary coverage limitation; sources outside
    the footprint return no results.

    Parameters
    ----------
    pos : astropy.coordinates.SkyCoord
        Positions to search around (array or scalar).
    radius : astropy.units.Quantity
        Search radius.  Default 3 arcsec.
    logging_tag : str
        Prefix for log messages.

    Returns
    -------
    astropy.table.Table
        Catalog sources with columns matching :data:`CROSS_MATCH_DTYPE`.

    Notes
    -----
    VizieR catalog V/161 (DESI Collaboration 2023): the DESI Early Data
    Release (EDR) spectroscopic redshift catalog, containing ~3.5 M spectra
    from the commissioning, survey-validation, and early main-survey phases.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> pos = SkyCoord([36.45], [−4.59], unit='deg')
    >>> t = get_desi_v161_candidates(pos, radius=3 * u.arcsec)  # doctest: +SKIP
    """
    if not pos.isscalar and len(pos) == 0:
        logger.debug(f"{logging_tag} DESIV161 has no results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    vizier = Vizier(
        columns=["RAICRS", "DEICRS", "TargetID", "OType", "ZWARN", "delChi2"],
        row_limit=-1,
    )

    try:
        result = _query_vizier_region(vizier, pos, radius, "V/161")
        tables = list(result)
    except Exception as e:
        logger.warning(f"{logging_tag} DESIV161 VizieR query failed: {e}")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    if not tables:
        logger.debug(f"{logging_tag} DESIV161 has no results")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    desi = table.vstack(tables, metadata_conflicts="silent")

    # Filter first: spectroscopic type QSO or GALAXY, reliable redshift. Must
    # happen before deduplication -- see _dedupe_desi_by_delchi2.
    type_arr = np.array([str(t).strip() for t in desi["OType"]])
    zwarn_arr = np.array(desi["ZWARN"], dtype=int)
    keep = np.isin(type_arr, ["QSO", "GALAXY"]) & (zwarn_arr == 0)
    desi = desi[keep]

    if len(desi) == 0:
        logger.debug(f"{logging_tag} DESIV161 has no QSO/GALAXY results with ZWARN=0")
        return table.Table(dtype=CROSS_MATCH_DTYPE)

    # Collapse DESI's own repeat targeting of one physical object -- both within
    # a TargetID's repeat observations and across different TargetIDs at the
    # same position -- to the single best (highest-delChi2) row.
    desi = _dedupe_desi_by_delchi2(desi)

    output = table.Table(data=np.zeros(len(desi), dtype=CROSS_MATCH_DTYPE))
    output["ra"] = np.array(desi["RAICRS"], dtype=float)
    output["dec"] = np.array(desi["DEICRS"], dtype=float)
    output["catalog"] = "DESIV161"
    output["name"] = [f"DESIV161-{tid}" for tid in desi["TargetID"]]
    output["mag"] = table.MaskedColumn(
        dtype=np.float32, length=len(desi), mask=np.ones(len(desi), dtype=bool)
    )

    logger.debug(
        f"{logging_tag} DESIV161: {len(output)} QSO/GALAXY results with ZWARN=0"
    )
    return output


# Locally cached whole catalogs served by the bounding-cone prefilter instead of a
# live VizieR query. rough_match/_get dispatch here first; add new local catalogs
# to this mapping rather than branching in _get.
LOCAL_CATALOG_GETTERS = {
    "RFC": _get_rfc,
    "ICRF3": _get_icrf3,
    "ICRF2": _get_icrf2,
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
    if getter := LOCAL_CATALOG_GETTERS.get(catalog):
        return getter(pos=pos, radius=radius, logging_tag=logging_tag)
    params = VIZIER_CATALOGS[catalog].copy()
    columns = {}
    for name, col in params["columns"].items():
        columns[name] = col.format(time=time)
    params["columns"] = columns
    return get_vizier(pos, radius=radius, logging_tag=logging_tag, **params, raw=raw)


def remap_x_id_to_sources(
    candidates: table.Table, version_sources: table.Table
) -> table.Table:
    """Return a copy of `candidates` with ``x_id`` re-pointed at `version_sources`.

    ``astromon_cat_src`` stores one ``x_id`` per catalog source, but the pipeline
    writes one match per (catalog source, detect_method) pair -- so for every method
    beyond the one the stored ``x_id`` was assigned from, it has to be re-matched to
    that method's source ids first. The pipeline does that per version in memory and
    never persists the result.

    Which means anything recomputing matches from the stored table has to redo this
    rather than trust the stored ``x_id``. Skipping it silently yields one row per
    catalog source instead of one per (source, method) -- half the matches, for two
    methods -- and only the selections carrying a ``detect_method_filter`` come out
    whole, because they want exactly one method anyway.

    Parameters
    ----------
    candidates
        Catalog candidates with ``ra``, ``dec`` and an ``x_id`` column.
    version_sources
        The x-ray sources of a single detect method, with ``ra``, ``dec`` and ``id``.
    """
    remapped = candidates.copy()
    if len(candidates) == 0 or len(version_sources) == 0:
        return remapped
    nearest, _, _ = coords.SkyCoord(
        candidates["ra"], candidates["dec"], unit="deg"
    ).match_to_catalog_sky(
        coords.SkyCoord(version_sources["ra"], version_sources["dec"], unit="deg")
    )
    remapped["x_id"] = np.asarray(version_sources["id"])[nearest]
    return remapped


def rough_match(
    sources,
    time,
    radius=3 * u.arcsec,
    catalogs=("Tycho2", "USNO-B1.0", "2MASS", "SDSS"),
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
        Default is ('Tycho2', 'USNO-B1.0', '2MASS', 'SDSS')
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


_BAD_SOURCES_FILE = Path(__file__).parent / "data" / "bad_sources.ecsv"


def load_bad_sources(bad_sources_file: Path | None = None) -> table.Table:
    """Load the bad-sources exclusion table from an ECSV file.

    The file has three kinds of rows controlled by the ``excl_type`` column:

    ``'target'``
        Exclude all cross-matches from any observation whose target name starts
        with ``target_prefix`` (comparison is space-stripped and case-insensitive).
    ``'cat_source'``
        Exclude any cross-match to a catalog source whose name exactly equals
        ``cat_name``, regardless of obsid.
    ``'obsid_src'``
        Exclude the specific detected X-ray source identified by ``(obsid, x_id)``.

    Parameters
    ----------
    bad_sources_file:
        Path to the ECSV file.  Defaults to ``astromon/data/bad_sources.ecsv``
        bundled with the package.

    Returns
    -------
    astropy.table.Table
    """
    path = Path(bad_sources_file or _BAD_SOURCES_FILE)
    return _read_catalog_cached(
        path, lambda p: table.Table.read(p, format="ascii.ecsv")
    )


def get_bad_source_mask(
    matches: table.Table,
    bad_sources_file: Path | None = None,
) -> np.ndarray:
    """Return a boolean mask marking rows of *matches* that should be excluded.

    Applies all three exclusion types defined in the bad-sources file (see
    :func:`load_bad_sources`).  The returned mask is ``True`` for rows that
    are excluded (same convention as the legacy :func:`get_bad_target_mask`).

    An empty `matches` table returns an empty mask rather than raising.

    Parameters
    ----------
    matches:
        Table of astromon cross-matches as returned by :func:`db.get_cross_matches`.
        May have zero rows.
    bad_sources_file:
        Optional override path for the bad-sources ECSV file.

    Returns
    -------
    numpy.ndarray of bool, shape ``(len(matches),)``
    """
    bad = load_bad_sources(bad_sources_file)
    exclude = np.zeros(len(matches), dtype=bool)

    if len(matches) == 0:
        # The masks below are built from list comprehensions, which yield an empty
        # float64 array for a zero-row table and so cannot be or-ed into a bool one.
        # filter_matches and compute_cross_matches both call this on tables that can
        # legitimately be empty (all rows removed by the dr / r_angle cuts).
        return exclude

    targets_norm = np.array(
        [t.replace(" ", "").lower() for t in np.array(matches["target"])]
    )
    target_prefixes = [
        row["target_prefix"].replace(" ", "").lower()
        for row in bad
        if row["excl_type"] == "target"
    ]
    for prefix in target_prefixes:
        exclude |= np.array([t.startswith(prefix) for t in targets_norm])

    cat_names_excluded = {
        row["cat_name"]
        for row in bad
        if row["excl_type"] == "cat_source" and row["cat_name"]
    }
    if cat_names_excluded and "name" in matches.colnames:
        exclude |= np.array(
            [n in cat_names_excluded for n in np.array(matches["name"])]
        )

    obsid_src_excluded = {
        (int(row["obsid"]), int(row["x_id"]))
        for row in bad
        if row["excl_type"] == "obsid_src"
    }
    if (
        obsid_src_excluded
        and "obsid" in matches.colnames
        and "x_id" in matches.colnames
    ):
        obsids = np.array(matches["obsid"])
        x_ids = np.array(matches["x_id"])
        exclude |= np.array(
            [
                (int(o), int(x)) in obsid_src_excluded
                for o, x in zip(obsids, x_ids, strict=True)
            ]
        )

    return exclude


def get_bad_target_mask(matches: table.Table) -> np.ndarray:
    """Return an exclusion mask based on observation target name.

    .. deprecated::
        Use :func:`get_bad_source_mask` instead, which reads from
        ``astromon/data/bad_sources.ecsv`` and also supports catalog-source
        and obsid-level exclusions.
    """
    return get_bad_source_mask(matches)


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
    in_region = matches["x_loc"][i].separation(loc[j]) < regions["radius"][j] * u.arcsec
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
        ok &= ~np.isin(matches["category"], exclude_categories)

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
              Default is ['RFC', 'Tycho2']
            - **snr**: Minimum signal-to-noise ratio of x-ray sources to consider.
              Default is 3.  NOTE: SNR scales differ by detect_method (gaussian_detect SNR
              is ~13x celldetect SNR for the same source), so this threshold is only
              meaningful per method - the standard gaussian selections use 9.5, the
              quantile-equivalent of celldetect's 3.
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

    detect_method_filter = args.pop("detect_method_filter", None)
    if (
        detect_method_filter is not None
        and "detect_method" in astromon_xray_src.dtype.names
    ):
        astromon_xray_src = astromon_xray_src[
            astromon_xray_src["detect_method"] == detect_method_filter
        ]

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
    catalogs=("RFC", "Tycho2"),
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
    - X-ray sources not flagged as lying on an ACIS readout streak (``acis_streak=False``),
      or sources flagged as the brightest in the observation (``brightest=True``), which
      overrides the streak exclusion for the calibration target.
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
        A list of catalog names. The order matters. Default is ['RFC', 'Tycho2']
    snr: float
        Minimum signal-to-noise ratio of x-ray sources to consider.
        Default is 3.  SNR scales differ by detect_method (gaussian_detect SNR is ~13x
        celldetect SNR for the same source), so choose this per method: the standard
        gaussian selections use 9.5, the quantile-equivalent of celldetect's 3.
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

    arm_flag = (
        astromon_xray_src["grating_arm"].astype(bool)
        if "grating_arm" in astromon_xray_src.dtype.names
        else np.zeros(len(astromon_xray_src), dtype=bool)
    )
    brightest_flag = (
        astromon_xray_src["brightest"].astype(bool)
        if "brightest" in astromon_xray_src.dtype.names
        else np.zeros(len(astromon_xray_src), dtype=bool)
    )
    astromon_xray_src = astromon_xray_src[
        (astromon_xray_src["snr"] > snr)
        & (astromon_xray_src["near_neighbor_dist"] > near_neighbor_dist)
        & (~astromon_xray_src["acis_streak"].astype(bool) | brightest_flag)
        & (~arm_flag)
    ]
    astromon_cat_src = astromon_cat_src[np.isin(astromon_cat_src["catalog"], catalogs)]

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


# Catalog hierarchy for the astromon_24-27 selections. Order matters: for each
# X-ray source the match from the earliest listed catalog wins.
_MATCH_HIERARCHY = [
    "ICRF3",
    "RFC",
    "GaiaAGN",
    "Quaia",
    "DESIV161",
    "MilliquasGaia",
    "GaiaQSO",
    "GaiaVarStar",
    "Tycho2",
]
# astromon_23 predates the hierarchy reorder and keeps its original ordering.
_MATCH_HIERARCHY_23 = [
    "RFC",
    "GaiaAGN",
    "GaiaQSO",
    "Quaia",
    "MilliquasGaia",
    "DESIV161",
    "GaiaVarStar",
    "Tycho2",
]


def _simple_match_args(name, catalogs, **overrides):
    """One CROSS_MATCHES_ARGS entry with the shared 'simple'-method defaults."""
    args = {
        "name": name,
        "method": "simple",
        "catalogs": list(catalogs),
        "snr": 3,
        "r_angle": 120.0,
        "r_angle_grating": 120.0,
        "near_neighbor_dist": utils.NEAR_NEIGHBOR_DIST_ARCSEC,
        "dr": 3.0,
    }
    args.update(overrides)
    return args


# Notes on the non-default values:
#
# detect_method_filter (astromon_24-27): compute_cross_matches discards xray
# sources from other methods before the join, so each select_name carries a
# single clean detect_method.
#
# snr=9.5 (astromon_25/27): SNR floors are method-specific because the SNR
# *scales* are not comparable: for the same source, gaussian_detect SNR is ~13x
# celldetect SNR (median; the ratio is not even constant) - the two tools use
# different noise definitions.  snr=9.5 on gaussian SNR is the same quantile as
# snr=3 on celldetect SNR, so the gaussian selections use 9.5 to match the
# strictness celldetect has always had.  With a shared snr=3 the gaussian
# samples admit faint junk: matches below gaussian SNR 9.5 are ~10x enriched in
# dr > 1.5" (replicates on peak_gaussian_detect, survives alignment-drift
# correction).  See the off-axis-angle-limit notebook in absolute_astrometry
# (2026-08).
#
# r_angle=160 (astromon_26/27): extended off-axis matches using existing
# astromon_cat_src coverage.  The 160-180" bin shows noticeably elevated
# centroid scatter (2.4% dr > 2" vs ~0.6-1.1% at smaller radii), so the limit
# is 160" rather than 180".  Note: astromon_cat_src was populated with a 3"
# rough_match radius, so sources where the PSF spread pushes the true
# counterpart beyond 3" are not represented.
CROSS_MATCHES_ARGS = {
    args["name"]: args
    for args in [
        _simple_match_args("astromon_21", ["RFC", "Tycho2"]),
        _simple_match_args("astromon_22", ["RFC", "Tycho2"], r_angle_grating=24.0),
        # Single-catalog selections, mostly for per-catalog diagnostics.
        _simple_match_args("gaia_agn", ["GaiaAGN"], r_angle_grating=24.0),
        _simple_match_args("gaia_qso", ["GaiaQSO"], r_angle_grating=24.0),
        _simple_match_args("gaia_var_star", ["GaiaVarStar"], r_angle_grating=24.0),
        _simple_match_args("milliquas_gaia", ["MilliquasGaia"], r_angle_grating=24.0),
        _simple_match_args("desi_v161", ["DESIV161"], r_angle_grating=24.0),
        _simple_match_args("quaia", ["Quaia"]),
        _simple_match_args("rfc", ["RFC"]),
        _simple_match_args("tycho2", ["Tycho2"]),
        _simple_match_args("icrf3", ["ICRF3"]),
        # Not part of any combined hierarchy (astromon_23-27): ICRF2 is superseded by
        # ICRF3/RFC and kept only as a reference baseline for the pre-ICRF3/RFC
        # catalog set, so no reprocessing of the hierarchical selections is needed.
        _simple_match_args("icrf2", ["ICRF2"]),
        _simple_match_args("astromon_23", _MATCH_HIERARCHY_23),
        _simple_match_args(
            "astromon_24", _MATCH_HIERARCHY, detect_method_filter="celldetect"
        ),
        _simple_match_args(
            "astromon_25",
            _MATCH_HIERARCHY,
            detect_method_filter="gaussian_detect",
            snr=9.5,
        ),
        _simple_match_args(
            "astromon_26",
            _MATCH_HIERARCHY,
            detect_method_filter="celldetect",
            r_angle=160.0,
            r_angle_grating=160.0,
        ),
        _simple_match_args(
            "astromon_27",
            _MATCH_HIERARCHY,
            detect_method_filter="gaussian_detect",
            snr=9.5,
            r_angle=160.0,
            r_angle_grating=160.0,
        ),
    ]
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
