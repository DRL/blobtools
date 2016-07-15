#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools seqfilter       -i FASTA -l LIST [-o PREFIX] [-v]
                                    [-h|--help]

    Options:
        -h --help                   show this

        -i, --infile <FASTA>        FASTA file of sequences (Headers are split at whitespaces)
        -l, --list <LIST>           TXT file containing headers of sequences to keep
        -o, --out <PREFIX>          Output prefix
        -v, --invert                Invert filtering (Sequences w/ headers NOT in list)
"""

from __future__ import division
from docopt import docopt
import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtPlot as BtPlot
from os.path import dirname, isfile, basename, splitext



if __name__ == '__main__':
    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    fasta_f = args['--infile']
    list_f = args['--list']
    invert = args['--invert']
    prefix = args['--out']

    out_f = BtIO.getOutFile(fasta_f, prefix, "filtered.fna")

    items = BtIO.parseSet(list_f)
    output = []
    for header, sequence in BtIO.readFasta(fasta_f):
        if header in items:
            if not (invert):
                output.append(">%s\n%s" % (header, sequence))
        else:
            if (invert):
                output.append(">%s\n%s" % (header, sequence))
    with open(out_f, "w") as fh:
        fh.write("".join(output))


