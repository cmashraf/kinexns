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
            description="Building and optimizing kinetic models for bio-fuel pyrolysis",
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
            install_requires=REQUIRES,
            requires=REQUIRES)


if __name__ == '__main__':
    setup(install_requires=['numpy', 'pandas', 'spotpy', 'matplotlib'])
