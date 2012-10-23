"""
Compute the Celestial location radius RMS corresponding to the PRD requirement
of 1.0 arcsec.
"""

import numpy as np
import asciitable

# Read using Tab instead of Rdb because the RDB 2nd-line header is wrong.
dat = asciitable.read('/proj/sot/ska/data/astromon/standard_xcorr/plot.rdb',
                      Reader=asciitable.Tab,
                      data_start=2, guess=False)
ok = dat['status_id'] == ''
dat = dat[ok]

# For all data
print 'All', np.sqrt(np.mean(dat['dr'] ** 2))
# 0.410

# For data since 2010
ok = dat['date_obs'] > '2009'
print 'Since 2010:001', np.sqrt(np.mean(dat['dr'][ok] ** 2))
# 0.517
