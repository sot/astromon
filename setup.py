# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import setup
import os

try:
    from testr.setup_helper import cmdclass
except ImportError:
    cmdclass = {}

entry_points = {
    'console_scripts': [
        'astromon-process-obs=astromon.observation:main',
        'astromon-cat-obs-data=astromon.scripts.get_cat_obs_data:main',
        'astromon-web-pages=astromon.web.celmon:main',
    ]
}


data_files = [(
    os.path.join('share', 'astromon'),
    ['task_schedules/task_schedule_update.cfg']
)]


setup(name='astromon',
      author='Javier Gonzalez',
      description='Tools used for absolute astrometry monitoring',
      author_email='javier.gonzalez@cfa.harvard.edu',
      packages=[
          'astromon',
          'astromon.scripts',
          'astromon.scripts.calalign',
          'astromon.scripts.calalign.check',
          'astromon.web'
      ],
      package_data={
          'astromon': ['sql/x-corr/*.sql', 'sql/tables/*.sql'],
          'astromon.web': ['templates/*']
      },
      data_files=data_files,
      license=("New BSD/3-clause BSD License\nCopyright (c) 2019"
               " Smithsonian Astrophysical Observatory\nAll rights reserved."),
      url='https://sot.github.io/astromon',
      entry_points=entry_points,
      use_scm_version=True,
      setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
      zip_safe=False,
      # tests_require=['pytest'],
      cmdclass=cmdclass,
      )
