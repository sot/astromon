This directory contains scripts to calculate CALALIGN offsets in regular time bins and to create
new CALALIGN files with the new values. The following commands first estimate the median offset
in 6-month intervals starting after 2013:182, saves the results in `offsets.json`, and creates
CALALIGN files in the  `files` directory.

```
./calc_median.py --start 2013:182 --out offsets.json
./set_calalign_offsets.py --offsets offsets.json
```

NOTE: The offsets from `calc_median.py` are (x-ray-location - cat-source-location), so the CALALIGN
correction should be the negative of the offset.

 The `check` subdirectory contains scripts to run the aspect
pipeline on a list of OBSIDs, while explicitly setting the CALALIGN file to use, and to produce
plots of the offsets.

```
./reproject_obsids.py  # runs the aspect pipeline with the new and old CALALIGN
./combine_obsids.py    # combines the output from the previous step into summary FITS files
./make_plots.py        # Make the plots
```

To run the pipeline, these script use the script in https://github.com/sot/runasp