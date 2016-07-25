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

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import lib.BtLog as BtLog
import lib.BtIO as BtIO

def main():
    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    fasta_f = args['--infile']
    list_f = args['--list']
    invert = args['--invert']
    prefix = args['--out']

    output = []
    out_f = BtIO.getOutFile(fasta_f, prefix, "filtered.fna")

    print BtLog.status_d['21'] % list_f
    items = BtIO.parseSet(list_f)

    print BtLog.status_d['22'] % fasta_f
    parsed_items = []
    sequences = 0
    for header, sequence in BtIO.readFasta(fasta_f):
        sequences += 1
        if header in items:
            if not (invert):
                parsed_items.append(header)
                output.append(">%s\n%s\n" % (header, sequence))
        else:
            if (invert):
                parsed_items.add(header)
                output.append(">%s\n%s\n" % (header, sequence))
        BtLog.progress(len(output), 10, len(items), no_limit=True)
    print BtLog.status_d['23'] % ('{:.2%}'.format(parsed_items/sequences), parsed_items, sequences)
    if not len(parsed_items) == len(set(parsed_items)):
        print BtLog.warn_d['8'] % "\n\t".join(list(set([x for x in parsed_items if parsed_items.count(x) > 1])))
    with open(out_f, "w") as fh:
        fh.write("".join(output))

if __name__ == '__main__':
    main()
