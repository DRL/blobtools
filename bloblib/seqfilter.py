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

import bloblib.BtLog as BtLog
import bloblib.BtIO as BtIO

def main():
    args = docopt(__doc__)
    fasta_f = args['--infile']
    list_f = args['--list']
    invert = args['--invert']
    prefix = args['--out']

    output = []
    out_f = BtIO.getOutFile(fasta_f, prefix, "filtered.fna")

    print BtLog.status_d['1'] % ("list", list_f)
    items = BtIO.parseSet(list_f)
    items_count = len(items)
    print BtLog.status_d['22'] % fasta_f
    items_parsed = []
    sequences = 0
    for header, sequence in BtIO.readFasta(fasta_f):
        sequences += 1
        if header in items:
            if not (invert):
                items_parsed.append(header)
                output.append(">%s\n%s\n" % (header, sequence))
        else:
            if (invert):
                items_parsed.append(header)
                output.append(">%s\n%s\n" % (header, sequence))
        BtLog.progress(len(output), 10, items_count, no_limit=True)
    BtLog.progress(items_count, 10, items_count)

    items_parsed_count = len(items_parsed)
    print BtLog.status_d['23'] % ('{:.2%}'.format(items_parsed_count/sequences), "{:,}".format(items_count), "{:,}".format(items_parsed_count), "{:,}".format(sequences))

    items_parsed_count_unique = len(set(items_parsed))
    if not items_parsed_count == items_parsed_count_unique:
        print BtLog.warn_d['8'] % "\n\t\t\t".join(list(set([x for x in items_parsed if items_parsed.count(x) > 1])))

    with open(out_f, "w") as fh:
        print BtLog.status_d['24'] % out_f
        fh.write("".join(output))

if __name__ == '__main__':
    main()
