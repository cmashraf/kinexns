from __future__ import absolute_import, division, print_function
from os.path import join as pjoin

#  Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ''
_version_extra = 'dev'
#  _version_extra = ''  #  Uncomment this for full releases

#  Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

#  Description should be a one-liner:
description = "kinexns: Building and optimizing chemical kinetic models"
#  Long description will go up on the pypi page
long_description = """

kinexns
========
kinexns is a tool to build kinetic models with given reaction mechanism
and optimize the model parameters by solving ODEs.

It contains software implementations of an analysis of some simple data, but
more importantly, it contains infrastructure for testing, documentation,
continuous integration and deployment, which can be easily adapted
to use in other projects.


License
=======
``kinexns`` is licensed under the terms of the MIT license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective holders.

"""

NAME = "kinexns"
MAINTAINER = "Chowdhury Ashraf"
MAINTAINER_EMAIL = "cmashraf@uw.edu"
DESCRIPTION = "Building and optimizing chemical kinetic models"
LONG_DESCRIPTION = long_description
URL = "https://github.com/cmashraf/kinexns"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Chowdhury Ashraf"
AUTHOR_EMAIL = "cmashraf@uw.edu"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGE_DATA = {'kinexns': [pjoin('data', '*')]}
REQUIRES = ["numpy"]
