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

Astromon data can be accessed through the astromon.db module using convenience functions that
return the tables as an astropy.table.Table::

    from astromon import db
    matches = db.get_cross_matches()
    observations = db.get('astromon_obs')
    xray_src = db.get('astromon_xray_src')
    cat_src = db.get('astromon_cat_src')
    xcorr = db.get('astromon_xcorr')

NOTE: One can modify the behavior of astromon.db and specify a different database file by defining
the ASTROMON_FILE environmental variable.

A more elaborate example to do some analysis, correlating several MSIDs with astromon data, and
excluding observations of some specific targets::

    import os
    os.environ['ASTROMON_FILE'] = '/Users/javierg/SAO/ska/data/astromon/ASTROMON_table.h5'

    from astropy.table import join
    from astromon import db, telemetry

    matches = db.get_cross_matches()
    matches = join(matches, telemetry.get(matches['obsid'], '2013:001', '2022:001'), keys=['obsid'])
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
astromon.observation module. The main component is the Observation class, which encapsulates all
calls to arc5gl, CIAO scripts, etc

The following statement instanciates an observation::

    from astromon.observation import Observation
    obs = Observation(8008)

Upon creation, if the files corresponding to the observation are not available locally, it will use
arc5gl to download them. The files are stored locally in a temporary directory unless a working
directory is given to the constructor::

    from astromon.observation import Observation
    obs = Observation(8008, workdir='/Users/javierg/SAO/ska/data/astromon')

in which case the files will persist after the observation instance is deleted.

The Observation class encapsulates several common operations with CIAO. The following runs a
prescribed sequence of commands that should result in the creation of files with detected sources::

    obs.process()

Cross-Matching Algorithms
-------------------------

The very first step in the matching process is to find all counterparts within a radius around the
x-ray source. This is what we call a "rough match" and is done by astromon.cross_match.rough_match,
which encapsulates all queries to various standard catalogs.

The actual cross-match, which yields at most one catalog source per x-ray source is performed by
astromon.cross_match.rough_match.cross_match. This in turn delegates to specific implementation
functions.

Simple
^^^^^^

The algorithm currently used is a simple algorithm where:

    - x-ray sources with SNR below a given value are discarded
    - observations are filtered based on date
    - x-ray sources with other x-ray sources within a given radius are discarded
    - observations, x-ray sources and catalog sources are joined on OBSID and x-ray source ID.
    - matches with source-counterpart separation larger than a given value are discarded
    - catalog sources are selected from specific catalogs in a prescribed precedence. The first
      counterpart found is taken to be the match.