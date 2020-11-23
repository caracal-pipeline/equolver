#!/usr/bin/env python

import os

try:
    from setuptools import setup
except ImportError as e:
    from distutils.core import setup

requirements = [
    'numpy',
    'scipy',
    'astropy',
    'radio_beam',
    'pyfftw',
    'bokeh',
]

PACKAGE_NAME = 'equolver'
__version__ = '0.0.4'

setup(name=PACKAGE_NAME,
      version=__version__,
      description="Development Status :: 4 - Beta",
      author="Gyula Jozsa",
      author_email="gigjozsast@gmail.com",
      url="https://github.com/gigjozsa/equolver",
      packages=[PACKAGE_NAME],
      python_requires='>=3.6',
      install_requires=requirements,
      include_package_data=True,
      # package_data - any binary or meta data files should go into MANIFEST.in
      scripts=["bin/" + j for j in os.listdir("bin")],
      license=["BSD 3-Clause License"],
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: BSD License",
          "Operating System :: POSIX :: Linux",
          "Programming Language :: Python",
          "Topic :: Scientific/Engineering :: Astronomy"
      ]
      )
