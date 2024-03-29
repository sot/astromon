Astromon Components
===================

.. _cross_matching_algorithms:

Cross-Matching Algorithms
-------------------------

The process of finding counterparts for detected x-ray sources proceeds in two major steps:

- Rough match: Find all counterparts within a radius around the x-ray sources. These *rough* matches
  constitute a superset of all *reasonable* cross-matches. This is done by
  :any:`rough_match <astromon.cross_match.rough_match>`, which encapsulates all queries to the
  Tycho2, ICRS, USNO-B1.0, 2MASS, SDSS catalogs and selects matches within 3 arcsec.
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

These are the catalogs available for rough matches:

- `Tycho2 <https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/259/tyc2>`_. Reference catalog of
  2.5 million stars observed by the Tycho instrument abord the ESA Hipparcos satellite.
  Astrometric accuracy ~25 mas with stars down to ~11.5 mag. More information available in the
  `Guide to the Tycho2 catalog (PDF) <http://www.astro.ku.dk/~cf/CD/docs/guide.pdf>`_.
- `ICRF2 <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/323>`_.
  Sources from The Second Realization of the International Celestial Reference Frame by Very
  Long Baseline Interferometry `Ma et al. 2009, IERS Technical Note No. 35 (pdf)
  <http://cdsarc.u-strasbg.fr/ftp/cats/I/323/tn35.pdf>`_
- `USNO-B1.0 <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/284>`_.
  Monet, D.G. et al. (2003), "The USNO-B Catalog", The Astronomical Journal, vol. 125, no. 2,
  pp. 984-993.
- `2MASS <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/246>`_.
  Two micron all-sky survey: 162,213,354 million point sources from 19,600 square degrees of sky.
- `SDSS <https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/294>`_. The SDSS Photometric Catalog,
  Release 7

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
