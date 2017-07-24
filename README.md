BlobTools v1.0
===============================
A modular command-line solution for visualisation, quality control and taxonomic partitioning of genome datasets

- Discussions, questions and answers: [BlobTools GoogleGroup](https://groups.google.com/forum/#!forum/blobtools)
- Issues, bug reports and feature requests: [GitHub issues](https://github.com/DRL/blobtools/issues)
- Documentation: [blobtools.readme.io](https://blobtools.readme.io) (under reconstruction)

Dependencies
------------
- UNIX system
- Python 2.7+
- ```pip```

Installation
------------

    $ git clone https://github.com/DRL/blobtools.git
    $ cd blobtools
    $ ./install

The installation script will:
- install Python dependencies through ```pip```
- download and install a copy of [samtools-1.5](http://www.htslib.org/download/) into the folder ```blobtools/samtools/```
- download a copy of [NCBI TaxDump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) and create a nodesDB.txt file in ```blobtools/data/```. BlobTools will use this file for linking TaxIDs to NCBI taxonomies
- create a BlobTools executable (```blobtools```) in the main directory

Usage
-----

    $ ./blobtools -h
