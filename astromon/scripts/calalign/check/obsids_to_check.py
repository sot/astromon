from pathlib import Path

ALIGN_DIR = Path('/proj/sot/ska/jgonzalez/aca_cal_align/update_2022-feb/files')

ALIGN = [
    'pcadD2013-01-19alignN0010.fits',
    'pcadD2013-07-02alignN0010.fits',
    'pcadD2014-01-01alignN0010.fits',
    'pcadD2014-07-02alignN0010.fits',
    'pcadD2015-01-01alignN0010.fits',
    'pcadD2015-07-02alignN0010.fits',
    'pcadD2016-01-01alignN0010.fits',
    'pcadD2016-07-02alignN0010.fits',
    'pcadD2017-01-01alignN0010.fits',
    'pcadD2017-07-02alignN0010.fits',
    'pcadD2018-01-01alignN0010.fits',
    'pcadD2018-07-02alignN0010.fits',
    'pcadD2019-01-01alignN0010.fits',
    'pcadD2019-07-02alignN0010.fits',
    'pcadD2020-01-01alignN0010.fits',
    'pcadD2020-07-02alignN0010.fits',
    'pcadD2021-01-01alignN0010.fits',
    'pcadD2021-07-02alignN0010.fits'
]


# detector, obsid, alignment file index
OBS = [
    ('hrc-s', 15479, 0),
    ('hrc-s', 14606, 1),
    ('hrc-s', 16397, 3),
    ('hrc-s', 17387, 4),
    ('hrc-s', 19855, 8),
    ('hrc-s', 19794, 9),
    ('hrc-s', 19444, 10),
    ('hrc-s', 20685, 11),
    ('hrc-s', 21741, 13),
    ('hrc-s', 22793, 14),
    ('hrc-s', 22822, 15),
    ('hrc-s', 23382, 16),
    ('hrc-s', 25067, 17),
    ('hrc-i', 17351, 5),
    ('hrc-i', 18913, 8),
    ('hrc-i', 21058, 10),
    ('hrc-i', 20663, 11),
    ('hrc-i', 21762, 12),
    ('hrc-i', 21701, 13),
    ('hrc-i', 21700, 14),
    ('hrc-i', 22658, 15),
    ('hrc-i', 24525, 16),
    ('hrc-i', 24546, 17),
    # ('acis-s', 14534, 0),  this onegets stuck
    ('acis-s', 14533, 1),
    ('acis-s', 15713, 2),
    ('acis-s', 15715, 3),
    ('acis-s', 15718, 4),
    ('acis-s', 16653, 5),
    ('acis-s', 16654, 6),
    ('acis-s', 17732, 7),
    ('acis-s', 17734, 8),
    ('acis-s', 18999, 9),
    ('acis-s', 20148, 10),
    ('acis-s', 20709, 11),
    ('acis-s', 19616, 12),
    ('acis-s', 21156, 13),
    ('acis-s', 22152, 14),
    ('acis-s', 22340, 15),
    ('acis-s', 22348, 16),
    ('acis-s', 23411, 17),
    ('acis-s', 24322, 17),
    ('acis-i', 13276, 0),
    ('acis-i', 14536, 2),
    ('acis-i', 16138, 3),
    ('acis-i', 16822, 4),
    ('acis-i', 17141, 5),
    ('acis-i', 17133, 6),
    ('acis-i', 17138, 7),
    ('acis-i', 17134, 8),
    ('acis-i', 19281, 9),
    ('acis-i', 20269, 10),
    ('acis-i', 20361, 11),
    ('acis-i', 21138, 12),
    ('acis-i', 21247, 13),
    ('acis-i', 22693, 14),
    ('acis-i', 22316, 15),
    ('acis-i', 23815, 17),
]

# OBS = OBS[:1]
# OBS = [
#    ('hrc-i', 20663, 11)
# ]

# OBS = [
#     ('acis-s', 18999, 9),
#     ('hrc-i', 21058, 10),
#     ('hrc-s', 23382, 16),
#     ('hrc-s', 19855, 8),
# ]
