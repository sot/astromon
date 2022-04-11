Astromon Components
===================

Astromon is a package to monitor the astrometric accuracy of the Chandra X-ray observatory.

Database
--------

The main data product of Astromon is a database that contains the following tables:

* Observations (astromon_obs). Observation info from OBSPARS.
* X-ray sources (astromon_xray_src). Point-sources detected in Chandra observations.
* Catalog sources (astromon_cat_src). Optical/Radio catalog entries in the vicinity of the X-ray
  sources
* Cross-matches (astromon_xcorr). Pairs of X-ray and catalog sources, where there is at most one
  catalog source associated with a given detected X-ray source. 

Astromon data can be accessed through the :any:`astromon.db` module using convenience functions that
return the tables as an :ref:`Astropy Table <astropy:astropy-table>` ::

    from astromon import db
    matches = db.get_cross_matches()
    observations = db.get('astromon_obs')
    xray_src = db.get('astromon_xray_src')
    cat_src = db.get('astromon_cat_src')
    xcorr = db.get('astromon_xcorr')

.. Note::

    One can modify the behavior of :any:`astromon.db` and specify a different database file by
    defining the ASTROMON_FILE environmental variable.

A more elaborate example to do some analysis, correlating several MSIDs with astromon data, and
excluding observations of some specific targets::

    from astropy.table import join
    from astromon import db

    matches = db.get_cross_matches()
    matches.add_index('obsid')

    ok = np.ones(len(matches), dtype=bool)
    excl_targets = ['RW Aur', 'Tau Boo', '70 OPH', '16 Cyg', 'M87', 'Orion', 'HD 97950' 'HD4915']
    excl_targets = [x.replace(' ', '').lower() for x in excl_targets]
    for ii, target in enumerate(matches['target']):
        target = target.replace(' ', '').lower()
        for bad_target in excl_targets:
            if target.startswith(bad_target):
                ok[ii] = False
    matches['ok'] = ok


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


Cross-Matching Algorithms
-------------------------

The first step in the matching process is to find all counterparts within a radius around the
x-ray sources. These *rough* matches constitute a superset of all *reasonable* cross-matches.
This is done by :any:`astromon.cross_match.rough_match`, which encapsulates all queries to various
catalogs:

.. _catalog-list:

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

The actual cross-match, which yields at most one catalog source per x-ray source is performed by
:any:`astromon.cross_match.cross_match`. This in turn delegates to specific
implementation functions that implement algorithms described below.

Simple Matching
^^^^^^^^^^^^^^^

The algorithm currently used is the ``"astromon_21"`` algorithm, which does the following:

    - x-ray sources with ``SNR < 3`` are discarded
    - observations are filtered based on date
    - **x-ray sources with other x-ray sources within a given radius are discarded**
    - observations, x-ray sources and catalog sources are joined on OBSID and x-ray source ID.
    - **matches with source-counterpart separation larger than a given value are discarded**
    - for each x-ray source, the candidate matches are sorted according to catalog precedence and
      angular separation, and the first match is selected. The catalogs considered are: (ICRF2,
      Tycho2).

**Bold statements are the ones I do not see in the source.**