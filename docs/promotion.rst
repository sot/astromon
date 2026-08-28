.. _promotion:

Promoting a Dev Database to Production
=======================================

Astromon development happens against a separate *dev* database file (for example
``/Volumes/Black/data/astromon/astromon.h5``), which is where new catalogs, detection-method
changes, and backfill/cleanup scripts get tried out and validated before they reach the pipeline
that Chandra operations and the public celmon pages (``astromon.web.celmon``) depend on.
*Promotion* is the act of making that dev database (or the changes validated in it) the new
production database at ``$SKA/data/astromon/astromon.h5``.

There is no automated promotion tool. This page describes the manual procedure and the things that
have caused problems before.

How production updates itself day to day
-----------------------------------------

It helps to know what runs against the production file every day, since promotion has to fit
around it without corrupting an in-progress cron cycle. The ``astromon_update`` task
(``task_schedules/task_schedule_update.cfg``) runs on a cron schedule and, each time:

#. Copies ``$SKA/data/astromon/astromon.h5`` to a working copy under
   ``$SKA/data/astromon/tmp``.
#. Runs ``astromon-cross-match`` (``astromon.scripts.get_cat_obs_data``) against that working
   copy, which processes any newly-public obsids and appends their sources, catalog candidates and
   cross-matches.
#. Rsyncs the new per-obsid image/archive products into ``$SKA/data/astromon/archive``.
#. Moves the working copy back over ``$SKA/data/astromon/astromon.h5``.
#. Regenerates the celmon pages with ``astromon-web-pages`` into ``$SKA/data/astromon/web``.

Promotion has to replace ``$SKA/data/astromon/astromon.h5`` outside of that cycle, or the daily
job's copy/move dance can race with it (the daily job overwriting the file mid-swap, or promotion
overwriting the daily job's new obsids). Pause the ``astromon_update`` cron task before starting,
and confirm no cycle is mid-run, before touching the production file.

Pre-promotion checklist
------------------------

These apply to the dev file before it (or its content) is promoted:

#. **Code is already deployed.** Any code the dev database's content depends on -- new catalog
   getters, new ``CROSS_MATCHES_ARGS`` selections (for example ``astromon_25``), changed SNR
   defaults, detection-time exclusions -- must already be released to the production Ska3
   environment. A promoted database with ``select_name="astromon_25"`` rows is useless if the
   deployed ``astromon.cross_match`` does not know what ``astromon_25`` means yet.
#. **Validate the reprocessing.** Use the catalog validation notebook's checks (offset
   distributions, magnitude coverage, celldetect-vs-gaussian_detect comparison) and a celmon dry
   run against the dev file (``astromon-web-pages --db-file <dev-file> --out <scratch-dir>``) to
   confirm the pages look right *before* promoting, not after.
#. **Sync exclusion regions.** Run ``astromon.scripts.maintenance.sync_regions_from_primary`` to
   pull in any region excluded in production since the dev file was branched off it::

       python -m astromon.scripts.maintenance.sync_regions_from_primary \
           --primary $SKA/data/astromon/astromon.h5 \
           --dev /Volumes/Black/data/astromon/astromon.h5 \
           --dry-run

   Drop ``--dry-run`` once the report looks right. Regions are the one piece of state that
   production can gain independently of a dev branch (someone excludes a bad target while dev work
   is ongoing), so this step is not optional -- skipping it silently reverts that exclusion on
   promotion.
#. **Confirm maintenance backfills are complete, not just present.** The dev database accumulates
   its own backlog of one-off backfill/cleanup scripts as code changes land (see
   ``astromon/scripts/maintenance/``). Before promoting, check that every backfill relevant to what
   changed has actually been run against *this* dev file -- a script existing in the repo does not
   mean it has been applied to the file you are about to promote.

Promotion steps
----------------

#. Pause the ``astromon_update`` cron task.
#. Back up the current production file (keep it -- this is the rollback path)::

       cp $SKA/data/astromon/astromon.h5 $SKA/data/astromon/astromon.h5.bak-<date>

#. Copy the validated dev file over production::

       cp /Volumes/Black/data/astromon/astromon.h5 $SKA/data/astromon/astromon.h5

#. Regenerate the celmon pages manually and inspect them before resuming the cron task::

       astromon-web-pages --db-file $SKA/data/astromon/astromon.h5 --out $SKA/data/astromon/web

#. Resume the ``astromon_update`` cron task.

Rollback
--------

If something looks wrong after promotion, restore the backup made above and regenerate the pages
again::

    cp $SKA/data/astromon/astromon.h5.bak-<date> $SKA/data/astromon/astromon.h5
    astromon-web-pages --db-file $SKA/data/astromon/astromon.h5 --out $SKA/data/astromon/web

Any obsids processed by ``astromon_update`` between the promotion and the rollback are lost and
need to be reprocessed once the file is stable again.
