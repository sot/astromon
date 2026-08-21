Astromon Components
===================

.. _cross_matching_algorithms:

Cross-Matching Algorithms
-------------------------

The process of finding counterparts for detected x-ray sources proceeds in two major steps:

- Rough match: Find all counterparts within a radius around the x-ray sources. These *rough* matches
  constitute a superset of all *reasonable* cross-matches, and are collected within 3 arcsec of each
  x-ray source. They come from two paths, both of which write to ``astromon_cat_src``:

  - :any:`rough_match <astromon.cross_match.rough_match>` encapsulates the VizieR-backed catalogs
    plus the locally cached radio catalogs. The pipeline calls it with RFC, ICRF3, Tycho2,
    USNO-B1.0, 2MASS and SDSS.
  - The AGN, quasar and variable-star catalogs are queried directly by the
    ``astromon.scripts.get_cat_obs_data`` pipeline script, since each needs its own service
    (Gaia TAP) or post-processing (epoch propagation, type and quality cuts) rather than a plain
    VizieR cone search.

  See :ref:`catalog-list` for both groups.
- Cross-match. Pair each x-ray source with at most one counterpart from the rough matches. This is
  done by the :any:`compute_cross_matches  <astromon.cross_match.compute_cross_matches>` function,
  which in turn delegates to specific implementation functions.

The default cross-match set is identified as ``"astromon_21"`` (see
:ref:`pre-computed-queries`). This set is produced using the `"simple"` matching
algorithm implemented in :any:`simple_cross_match
<astromon.cross_match.simple_cross_match>`.

.. _catalog-list:

Available catalogs
^^^^^^^^^^^^^^^^^^

Catalog names below are the values stored in the ``catalog`` column of ``astromon_cat_src``, and
are the names accepted in the ``catalogs`` argument of the :ref:`pre-computed queries
<pre-computed-queries>`.

Rough-match catalogs
""""""""""""""""""""

Available through :any:`rough_match <astromon.cross_match.rough_match>`. Those marked *(pipeline)*
are the ones the standard pipeline actually queries; the rest are available but not currently used
by any cross-match set.

- **RFC** *(pipeline)*. The `Radio Fundamental Catalogue
  <https://astrogeo.smce.nasa.gov/sol/rfc/>`_: ~22,800 compact radio sources with VLBI positions,
  typically accurate to well under a milliarcsecond. Downloaded from astrogeo.org as fixed-width
  ASCII and cached at ``$SKA/data/astromon/rfc_catalog.txt``, refreshed every 30 days.
  `Petrov & Kovalev 2025, ApJS 276, 38 <https://doi.org/10.3847/1538-4365/ad8c36>`_.
- **ICRF3** *(pipeline)*. `Third realization of the International Celestial Reference Frame
  <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/644/A159>`_ (S/X band), 4536 VLBI
  sources defining the celestial reference frame. Fetched from VizieR
  (``J/A+A/644/A159/table10``) and cached at ``$SKA/data/astromon/icrf3_catalog.ecsv``, refreshed
  yearly. `Charlot et al. 2020, A&A 644, A159 <https://doi.org/10.1051/0004-6361/202038368>`_.
  This supersedes ICRF2, which is no longer used.
- `Tycho2 <https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/259/tyc2>`_ *(pipeline)*.
  Reference catalog of 2.5 million stars observed by the Tycho instrument abord the ESA Hipparcos
  satellite. Astrometric accuracy ~25 mas with stars down to ~11.5 mag. More information available
  in the `Guide to the Tycho2 catalog (PDF) <http://www.astro.ku.dk/~cf/CD/docs/guide.pdf>`_.
- `USNO-B1.0 <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/284>`_ *(pipeline)*.
  Monet, D.G. et al. (2003), "The USNO-B Catalog", The Astronomical Journal, vol. 125, no. 2,
  pp. 984-993.
- `2MASS <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/246>`_ *(pipeline)*.
  Two micron all-sky survey: 162,213,354 million point sources from 19,600 square degrees of sky.
- `SDSS <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/294>`_ *(pipeline)*.
  The SDSS Photometric Catalog, Release 7.
- `HIP <http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/239>`_ and
  `HIP2 <http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/311>`_. The original and re-reduced
  Hipparcos catalogs.
- `UCAC4 <http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/322>`_. The fourth USNO CCD Astrograph
  Catalog.
- `Gaia2 <https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/345/gaia2>`_. Gaia DR2 sources.

AGN, quasar and variable-star catalogs
""""""""""""""""""""""""""""""""""""""

Queried directly by the pipeline. These are the catalogs added to improve the supply of
astrometrically reliable extragalactic counterparts, which anchor the absolute astrometry solution
far better than stars do.

- **GaiaAGN**. Gaia DR3 sources with a confirmed AGN cross-identification
  (``gaiadr3.agn_cross_id``, ~2.2 M sources), queried via Gaia TAP.
- **GaiaQSO**. Gaia DR3 quasar candidates (``gaiadr3.qso_candidates``, ~6.6 M sources), queried via
  Gaia TAP. Extends GaiaAGN with ~4.4 M sources classified from Gaia astrometry, photometry and
  colour, adding depth at G > 20 where the confirmed samples run out.
- **Quaia**. The `Gaia-unWISE quasar catalog <https://zenodo.org/records/10403370>`_, G < 20.0
  sample (755,850 sources), built from Gaia DR3 quasar candidates and unWISE infrared photometry.
  Cached at ``$SKA/data/astromon/quaia_catalog.fits``; static until Gaia DR4, so never
  auto-refreshed. `Storey-Fisher et al. 2024, ApJ 964, 69 <https://doi.org/10.3847/1538-4357/ad1328>`_.
- **MilliquasGaia**. `Milliquas v8 <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=VII/294>`_
  restricted to spectroscopically confirmed types (Q, A, B, N) and further required to have a Gaia
  DR3 counterpart within 1.5 arcsec, which supplies the Gaia position. Photometric candidates and
  the radio/X-ray association candidates are excluded.
- **DESIV161**. `DESI Early Data Release
  <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=V/161>`_ spectroscopic catalog, keeping
  ``OType`` of QSO or GALAXY with ``ZWARN == 0``. Positions come from Legacy Survey DR9 imaging and
  are accurate to about 0.1 arcsec. Coverage is limited to the ~4100 sq deg EDR footprint.
- **GaiaVarStar**. Gaia DR3 rotation-modulation variable stars, queried via Gaia TAP and restricted
  to RUWE < 1.4 so that unreliable astrometric solutions do not corrupt the proper-motion
  correction. Positions are propagated from the Gaia DR3 epoch (J2016.0) to the observation epoch.

.. Note::

   Only GaiaVarStar and the Tycho2/SDSS rough matches apply a proper-motion or epoch correction.
   The radio and extragalactic catalogs (RFC, ICRF3, GaiaAGN, GaiaQSO, Quaia, MilliquasGaia,
   DESIV161) describe sources with no measurable proper motion, so their catalog positions are
   used as-is.

Database
--------

The main data product of Astromon is a database that contains the following tables:

* Observations (astromon_obs). Observation info from OBSPARS.
* X-ray sources (astromon_xray_src). Point-sources detected in Chandra observations.
* Catalog sources (astromon_cat_src). Optical/Radio catalog entries in the vicinity of the X-ray
  sources
* Cross-matches (astromon_xcorr). Pairs of X-ray and catalog sources, where there is at most one
  catalog source associated with a given detected X-ray source.
* Excluded regions (astromon_regions). Circular regions within which all x-ray sources are excluded
  from analysis.

Astromon data can be accessed through the :any:`astromon.db` module using convenience functions that
return the tables as an :ref:`Astropy Table <astropy:astropy-table>` ::

    from astromon import db
    matches = db.get_cross_matches()
    observations = db.get_table('astromon_obs')
    xray_src = db.get_table('astromon_xray_src')
    cat_src = db.get_table('astromon_cat_src')
    xcorr = db.get_table('astromon_xcorr')
    regions = db.get_table('astromon_regions')

.. Note::

    One can modify the behavior of :any:`astromon.db` and specify a different database file by
    defining the ASTROMON_FILE environmental variable.

Observations
------------

The access to proprietary information related to observations is encapsulated in the
:any:`astromon.observation` module. The main component is the Observation class, which encapsulates
all calls to arc5gl, CIAO scripts, etc

The following statement instanciates an observation::

    from astromon.observation import Observation
    obs = Observation(8008)

Upon creation, if the files corresponding to the observation are not available locally, it will use
arc5gl to download them. The files are stored locally in a temporary directory unless a working
directory is given to the constructor::

    from astromon.observation import Observation
    obs = Observation(8008, workdir='./astromon/work')

in which case the files will persist after the observation instance is deleted.

The Observation class encapsulates several common operations with CIAO. The following runs a
prescribed sequence of commands that should result in the creation of files with detected sources::

    obs.process()
