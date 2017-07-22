BlobTools
===============================
A modular command-line solution for visualisation, quality control and taxonomic partitioning of genome datasets

Dependencies
------------
* functional UNIX installation with:
- bash
- wget
- tar
- Python2.7+
- pip

Installation
------------

    $ git clone https://github.com/DRL/blobtools.git
    $ cd blobtools
    $ ./install

The installation script will:
- install python dependencies through PIP
- download and install a copy of [samtools-1.5](http://www.htslib.org/download/) into the folder ```blobtools/samtools/```
- download a copy of [NCBI TaxDump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) and create a nodesDB.txt file which BlobTools will use for linking TaxIDs to taxonomies
- create a BlobTools executable

Usage
-----

    $ ./blobtools -h

Documentation
-------------

[blobtools.readme.io](https://blobtools.readme.io)
