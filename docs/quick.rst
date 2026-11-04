Quick Start
===========

.. _pre-computed-queries:

Pre-computed queries
--------------------

A common usage of the astromon data products is to get cross-matches with
:any:`db.get_cross_matches`. It returns a table that results from joining a few underlying
tables. Without a ``name``, it computes the *default selection* from the stored observations,
x-ray sources and catalog candidates, so the result always reflects what the database holds::

    from astromon import get_cross_matches
    matches = get_cross_matches()  # the default selection

With a ``name``, it reads that selection's pre-computed rows from ``astromon_xcorr``. The pipeline
stores one, ``"astromon_21"`` (flight's, read by the quarterly status number). Older databases may
still hold rows of retired selections (``"astromon_22"`` to ``"astromon_27"`` and the
single-catalog sets), which can be read by name but are no longer computed.

The arguments behind both are in :any:`cross_match.CROSS_MATCHES_ARGS
<astromon.cross_match.CROSS_MATCHES_ARGS>`::

    'astromon_21': {          # stored
        'name': 'astromon_21',
        'method': 'simple',
        'catalogs': ['RFC', 'Tycho2'],
        'snr': 3,
        'r_angle': 120.,
        'r_angle_grating': 120.,
        'near_neighbor_dist': 6.,
        'dr': 3.,
    },
    'default': {              # computed on demand
        'name': 'default',
        'method': 'simple',
        'catalogs': ['ICRF3', 'RFC', 'GaiaAGN', 'DESIV161', 'MilliquasGaia', 'GaiaQSO', 'Tycho2'],
        'detect_method_filter': 'celldetect',
        'snr': 3,
        'r_angle': 160.,
        'r_angle_grating': 160.,
        'near_neighbor_dist': 6.,
        'dr': 3.,
    }

Catalogs are listed in priority order: where an x-ray source has counterparts in more than one
catalog, the first matching catalog in the list wins. The default selection is the
absolute-astrometry calibration's.

Here is a more elaborate example to do some analysis, filtering on some columns,
and excluding x-ray sources in known problematic regions::

    from astromon import get_cross_matches
    matches = get_cross_matches(
        exclude_regions=True,
        exclude_bad_targets=True,
        sim_z=4,
        snr=10,
        dr=2,
        catalog=['Tycho2']
    )

.. Warning::
    Supplying additional filter arguments to ``get_cross_matches`` requires attention to
    ensure that the filter arguments are a *strict subset* of the selection's own
    arguments. In most cases you should use a custom cross-match as
    shown in the next section.

Custom cross-matches
--------------------

It is possible to re-compute cross-matches using the
:any:`compute_cross_matches <cross_match.compute_cross_matches>` function. The following uses
default arguments (the arguments to compute the standard cross-match described in
:ref:`cross_matching_algorithms`)::

    from astromon import cross_match
    matches = cross_match.compute_cross_matches()

and this is an example of calling it with some custom arguments (see the notes and warnings for the
:any:`compute_cross_matches <cross_match.compute_cross_matches>` and :any:`simple_cross_match
<cross_match.simple_cross_match>` functions)::

    from astromon import cross_match
    matches = cross_match.compute_cross_matches(
        catalogs=('USNO-B1.0'),
        snr=5,
        r_angle=60.,
        dr=2,
        r_angle_grating=6.,
        near_neighbor_dist=5.,
        start='2022:001',
        exclude_regions=True
    )

Filtering the List of Cross-matches
-----------------------------------

There are convenience methods to help filter x-ray sources within known problematic regions
(:any:`get_excluded_regions_mask <cross_match.get_excluded_regions_mask>`), in known difficult
targets (:any:`get_bad_target_mask <cross_match.get_bad_target_mask>`), or based on common criteria
(:any:`filter_matches <cross_match.filter_matches>`). These functions are used internally by
:any:`get_cross_matches <db.get_cross_matches>`, but they can be used directly::

    from astropy.table import join
    from astromon import cross_match, db

    matches = db.get_cross_matches(name='astromon_21')
    not_excluded = cross_match.get_excluded_regions_mask(matches)
    selected = cross_match.filter_matches(matches, snr=10, dr=2, catalog=['Tycho2'])
    matches = matches[not_excluded & selected]

.. Note::
    Filtering on certain fields might not do what you expect. For example, filtering on `catalog`
    removes sources that might have been matched with your selected catalog if you had re-computed
    the cross-match.


Astromon File Location
----------------------

With all these functions, there are two ways one can specify a custom DB file location::

    from astromon import db
    matches = db.get_cross_matches(dbfile='some_file.h5')

or::

    import os
    os.environ['ASTROMON_FILE'] = 'some_file.h5'
    from astromon import db
    matches = db.get_cross_matches()
