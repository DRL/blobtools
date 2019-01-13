#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: blobtools [<command>] [<args>...] [--help] [--version]

commands:

    create        create a BlobDB
    view          generate tabular view, CONCOCT input or COV files from BlobDB
    plot          generate a BlobPlot from a BlobDB
    covplot       generate a CovPlot from a BlobDB and a COV file

    map2cov       generate a COV file from BAM file
    taxify        generate a BlobTools compatible HITS file (TSV)
    bamfilter     subset paired-end reads from a BAM file
    seqfilter     subset sequences in FASTA file based sequence IDs in list
    nodesdb       create nodesdb based on NCBI Taxdump's names.dmp and nodes.dmp

    -h, --help      show this
    -v, --version   show version number

See 'blobtools <command> --help' for more information on a specific command.

examples:

    # 1. Create a BlobDB
    ./blobtools create -i example/assembly.fna -b example/mapping_1.bam -t example/blast.out -o example/test

    # 2. Generate a tabular view
    ./blobtools view -i example/test.blobDB.json

    # 3. Generate a blobplot
    ./blobtools blobplot -i example/test.blobDB.json

"""

# from __future__ import division
import sys
from subprocess import call
import os
from docopt import docopt

__version__ = "blobtools v1.0"
LIBDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ''))
MAINDIR = os.path.abspath(os.path.join(LIBDIR, '../'))
DATADIR = os.path.abspath(os.path.join(MAINDIR, 'data/'))
# SAMTOOLS = os.path.abspath(os.path.join(MAINDIR, 'samtools/bin/samtools'))
SAMTOOLS = "samtools"

if __name__ == '__main__':
    args = docopt(__doc__,
                  version=__version__,
                  options_first=True)

    argv = [args['<command>']] + args['<args>']
    if args['<command>'] == 'create':
        exit(call(['python', LIBDIR + '/create.py'] + argv))
    elif args['<command>'] == 'view':
        exit(call(['python', LIBDIR + '/view.py'] + argv))
    elif args['<command>'] == 'blobplot' or args['<command>'] == 'plot':
        argv[0] = "blobplot"
        exit(call(['python', LIBDIR + '/blobplot.py'] + argv))
    elif args['<command>'] == 'map2cov':
        exit(call(['python', LIBDIR + '/map2cov.py'] + argv))
    elif args['<command>'] == 'covplot' or args['<command>'] == 'comparecov':
        argv[0] = "covplot"
        exit(call(['python', LIBDIR + '/covplot.py'] + argv))
    elif args['<command>'] == 'seqfilter':
        exit(call(['python', LIBDIR + '/seqfilter.py'] + argv))
    elif args['<command>'] == 'taxify':
        exit(call(['python', LIBDIR + '/taxify.py'] + argv))
    elif args['<command>'] == 'bamfilter':
        exit(call(['python', LIBDIR + '/bamfilter.py'] + argv))
    elif args['<command>'] == 'nodesdb':
        exit(call(['python', LIBDIR + '/nodesdb.py'] + argv))
    else:
        exit(call([os.path.join(MAINDIR, './blobtools'), '-h']))

