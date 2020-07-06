import os
from setuptools import setup, find_packages
from kinexns.version import *

PACKAGES = find_packages()

# Get version and release info, which is all stored in kmpy/version.py
ver_file = os.path.join('kinexns', 'version.py')
with open(ver_file) as f:
    exec(f.read())

#HERE = path.abspath(path.dirname(__file__))

opts = dict(name="kinexns",
            maintainer="Chowdhury Ashraf",
            maintainer_email="cmashraf@uw.edu",
            description="Building and optimizing kinetic chemical models",
            url="https://github.com/cmashraf/kinexns",
            download_url="https://github.com/cmashraf/kinexns",
            license="MIT",
            classifiers=CLASSIFIERS,
            author="Chowdhury Ashraf",
            author_email="cmashraf@uw.edu",
            platforms="OS independent",
            version=VERSION,
            packages=PACKAGES,
            package_data=PACKAGE_DATA,
            install_requires=['numpy',
                              'pandas',
                              'spotpy',
                              'rdkit=2020.03.1.0=py37h65625ec_1',
                              'matplotlib',
                              'salib',
                              'assimulo',
                              'scipy',
                             ] ,
           )


if __name__ == '__main__':
    setup(**opts)
