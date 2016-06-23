#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools seqfilter       -i FASTA -l LIST [-v]
                                    [-h|--help]

    Options:
        -h --help                   show this

        -i, --infile FASTA          FASTA file of sequences. Headers are split at whitespaces.
        -l, --list LIST             TXT file containing headers of seuquences to keep
        -v, --invert                Inverts filtering (Sequences whose headers are NOT in the list are written)
"""

from __future__ import division
from docopt import docopt
import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtPlot as BtPlot
from os.path import dirname, isfile


if __name__ == '__main__':
    main_dir = dirname(__file__)
    print main_dir
    #print data_dir
    args = docopt(__doc__)
    fasta_f = args['--infile']
    list_f = args['--list']
    invert = args['--invert']

    # Check input files
    if not isfile(fasta_f):
        BtLog.error('0', fasta_f)
    if not isfile(list_f):
        BtLog.error('0', list_f)

    items = BtIO.parseSet(list_f)
    if not (invert):
        for header, sequence in BtIO.readFasta(fasta_f):
            if header in items:
                print ">%s\n%s\n" % (header, sequence)
    else:
        for header, sequence in BtIO.readFasta(fasta_f):
            if not header in items:
                print ">%s\n%s\n" % (header, sequence)
