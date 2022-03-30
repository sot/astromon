Quick Start
===========

The most common usage of the astromon data products is the list of cross-matches::

    from astromon import db
    matches = db.get_cross_matches()

This table is the result of joining a few underlying tables. One can easily get those tables too::

    from astromon import db
    observations = db.get('astromon_obs')
    xray_src = db.get('astromon_xray_src')
    cat_src = db.get('astromon_cat_src')
    xcorr = db.get('astromon_xcorr')

With all these functions, there are two ways one can optionally specify the DB file to use::

    from astromon import db
    matches = db.get_cross_matches(dbfile='some_file.h5')

or::

    import os
    os.environ['ASTROMON_FILE'] = 'some_file.h5'
    from astromon import db
    matches = db.get_cross_matches()
