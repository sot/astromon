"""
Compute the Celestial location radius RMS corresponding to the PRD requirement
of 1.0 arcsec.
"""

import numpy as np
import asciitable
from Chandra.Time import DateTime

# Read using Tab instead of Rdb because the RDB 2nd-line header is wrong.
dat = asciitable.read(
    '/proj/sot/ska/data/astromon/standard_xcorr/plot.rdb',
    Reader=asciitable.Tab,
    data_start=2,
    guess=False
)
ok = dat['status_id'] == ''
dat = dat[ok]

start = DateTime() - (5 * 365)
stop = DateTime()
ok = ((DateTime(dat['date_obs']).date > start.date)
      & (DateTime(dat['date_obs']).date < stop.date))
print("{} to {}".format(start.date, stop.date))

print("N srcs: {}".format(len(dat[ok])))
print('RMS radius {}'.format(np.sqrt(np.mean(dat[ok]['dr'] ** 2))))
print("90 percentile radius = {} arcsec".format(np.percentile(dat[ok]['dr'], 90)))
print("99 percentile radius = {} arcsec".format(np.percentile(dat[ok]['dr'], 99)))

for detector in ['ACIS-S', 'ACIS-I', 'HRC-S', 'HRC-I']:
    det = dat[ok]['detector'] == detector
    print(
        "90 percentile radius for {} is {} arcsec".format(
            detector,
            np.percentile(dat[ok]['dr'][det], 90)
        )
    )

print("{:.1f} percent outside a 1 arcsec radius".format(
    100.0 * np.count_nonzero(dat[ok]['dr'] > 1.0) / len(dat[ok]['dr']))
)


print("Worst case is {:.1f}".format(np.max(dat[ok]['dr'])))
