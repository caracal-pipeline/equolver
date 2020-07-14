#!/usr/bin/env python

import os

try:
    from setuptools import setup
except ImportError as e:
    from distutils.core import setup

requirements = [
#    'numpy>=1.14',
#    'python-casacore',
#    'scipy',
#    'astropy',
]

PACKAGE_NAME = 'equolver'
__version__ = '0.0.0'

setup(name=PACKAGE_NAME,
      version=__version__,
      description="Development Status :: 5 - Production/Stable",
      author="Gyula Jozsa",
      author_email="gigjozsast@gmail.com",
      url="https://github.com/gigjozsa/equolver",
      packages=[PACKAGE_NAME],
      python_requires='>=3.5',
      install_requires=requirements,
      include_package_data=True,
      # package_data - any binary or meta data files should go into MANIFEST.in
      #scripts=["bin/" + j for j in os.listdir("bin")],
      license=["GNU GPL v3"],
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: POSIX :: Linux",
          "Programming Language :: Python",
          "Topic :: Scientific/Engineering :: Astronomy"
      ]
      )
