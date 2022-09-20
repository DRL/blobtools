import pip
from setuptools import setup, find_packages

__version__ = '1.1'

# Get the long description from the README file
with open('README.md', 'r') as readme:
    long_description = readme.read()

# get the dependencies and installs
with open('requirements.txt', 'r') as requirements:
    reqs = requirements.read().splitlines()

setup(
    name='blobtools',
    version=__version__,
    description='A modular command-line solution for visualisation, quality control and taxonomic partitioning of genome datasets',
    long_description=long_description,
    url='https://github.com/DRL/blobtools',
    download_url='https://github.com/DRL/blobtools/tarball/' + __version__,
    license='GnuGPL3',
    classifiers=[
      'Development Status :: 4 - Beta',
      'Operating System :: POSIX',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Topic :: Scientific/Engineering :: Visualization',
      'Programming Language :: Python :: 3',
    ],
    keywords='Bioinformatics visualisation genome assembly QC',
    packages=find_packages(exclude=['docs', 'tests*']),
    include_package_data=True,
    author='Dominik R Laetsch',
    entry_points={
        'console_scripts': [
            "blobtools=lib.interface:main",
            ],
        },
    author_email='dominik.laetsch@gmail.com'
)
