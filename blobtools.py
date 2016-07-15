#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: blobtools <command> [<args>...] [--help] [--version]

commands:

    create        create a BlobDB
    view          print BlobDB as a table
    blobplot      plot BlobDB as a blobplot
    covplot       compare BlobDB cov(s) to additional cov file

    bam2cov       generate cov file from bam file
    sumcov        sum coverage from multiple COV files
    seqfilter     filter FASTA sequences based on header in list
    taxify        assign taxids to blast-results based on list

    -h, --help      show this
    -v, --version   show version number

examples:

    1. blobtools create -i assembly.fna -b reads.vs.assembly.bam -t assembly.vs.nt.blast.out -o test
    2. blobtools view -i test.blobDB.json
    3. blobtools blobplot -i test.blobDB.json

"""

from __future__ import division
import sys
from subprocess import call
from os.path import join, dirname
try:
    from docopt import docopt
except ImportError:
    sys.exit("[ERROR]\t: The module docopt is not installed. \n \tPlease run : pip install docopt")


__version__ = "blobtools v0.9.18"
MAINDIR = join(dirname(__file__), '')
DATADIR = join(MAINDIR, 'data/')
SRCDIR = join(MAINDIR, 'src/')
LIBDIR = join(MAINDIR, 'lib/')

if __name__ == '__main__':
    args = docopt(__doc__,
                  version=__version__,
                  options_first=True)

    argv = [args['<command>']] + args['<args>']
    if args['<command>'] == 'create':
        exit(call(['python', SRCDIR + 'create.py'] + argv))
    elif args['<command>'] == 'view':
        exit(call(['python', SRCDIR + 'view.py'] + argv))
    elif args['<command>'] == 'blobplot' or args['<command>'] == 'plot':
        argv[0] = "blobplot"
        exit(call(['python', SRCDIR + 'blobplot.py'] + argv))
    elif args['<command>'] == 'bam2cov':
        exit(call(['python', SRCDIR + 'bam2cov.py'] + argv))
    elif args['<command>'] == 'covplot' or args['<command>'] == 'comparecov':
        argv[0] = "covplot"
        exit(call(['python', SRCDIR + 'covplot.py'] + argv))
    elif args['<command>'] == 'sumcov':
        exit(call(['python', SRCDIR + 'sumcov.py'] + argv))
    elif args['<command>'] == 'seqfilter':
        exit(call(['python', SRCDIR + 'seqfilter.py'] + argv))
    elif args['<command>'] == 'taxify':
        exit(call(['python', SRCDIR + 'taxify.py'] + argv))
    else:
        exit("%r is not a blobtools command. See 'blobtools -h'." % args['<command>'])
