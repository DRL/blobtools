#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: blobtools <command> [<args>...] [--help] [--version]

commands:

    create        create a BlobDB
    view          print BlobDB as a table
    blobplot      plot BlobDB as a blobplot
    covplot       compare BlobDB cov(s) to additional cov file

    map2cov       generate cov file from bam file
    seqfilter     filter FASTA sequences based on header in list
    taxify        assign taxids to blast-results based on list

    -h, --help      show this
    -v, --version   show version number

See 'blobtools <command> --help' for more information on a specific command.

examples:

    # Create a BlobDB
    ./blobtools create -i test.fna -b test.bam -t test.blast.out -o test

    # Generate a tabular view
    ./blobtools view -i test.blobDB.json

    # Generate a blobplot
    ./blobtools blobplot -i test.blobDB.json

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
    elif args['<command>'] == 'map2cov':
        exit(call(['python', SRCDIR + 'map2cov.py'] + argv))
    elif args['<command>'] == 'covplot' or args['<command>'] == 'comparecov':
        argv[0] = "covplot"
        exit(call(['python', SRCDIR + 'covplot.py'] + argv))
    elif args['<command>'] == 'seqfilter':
        exit(call(['python', SRCDIR + 'seqfilter.py'] + argv))
    elif args['<command>'] == 'taxify':
        exit(call(['python', SRCDIR + 'taxify.py'] + argv))
    else:
        exit(call(['./blobtools', '--help']))

