import pip
from setuptools import setup, find_packages
from setuptools.command.install import install
from pip._internal.req import parse_requirements
from codecs import open
from os import path

__version__ = '1.1'

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# get the dependencies and installs
install_reqs = parse_requirements(path.join(here, 'requirements.txt'), session=False)
reqs = [str(ir.req) for ir in install_reqs]

class OverrideInstall(install):

    """
    Emulate sequential install of pip install -r requirements.txt
    To fix numpy bug in scipy, scikit in py2
    """

    def run(self):
        for req in reqs:
            pip._internal.main(["install", req])

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
      'Programming Language :: Python :: 2.7',
    ],
    keywords='Bioinformatics visualisation genome assembly QC',
    packages=find_packages(exclude=['docs', 'tests*']),
    include_package_data=True,
    author='Dominik R Laetsch',
    cmdclass={'install': OverrideInstall},
    author_email='dominik.laetsch@gmail.com'
)
