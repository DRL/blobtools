#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: taxify     -i BLAST -d TAXIDS
                                    [-h|--help]

    Options:
        -h --help                   show this

        -i, --infile BLAST          BLAST file of sequences (has to be in blast format)
        -d, --taxids TAXIDS         TACIX file (should have a flag for diamond/rnacentral)
"""

from __future__ import division
from docopt import docopt
#import lib.BtCore as bt
#import lib.BtLog as BtLog
#import lib.BtIO as BtIO
#import lib.BtPlot as BtPlot
from os.path import dirname, isfile

def parse_taxids(taxid_f):
    taxid_d = {}
    with open(taxid_f) as fh:
        for l in fh:
            line = l.rstrip("\n").split()
            taxid_d[line[0]] = line[3]
    return taxid_d

def print_renamed_blast(blast_f, taxid_d):
    with open(blast_f) as fh:
        for l in fh:
            line = l.rstrip("\n").split()
            print "%s\t%s\t%s" % (line[0], taxid_d[line[3]], "\t".join(line[2:]))

if __name__ == '__main__':
    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    blast_f = args['--infile']
    taxid_f = args['--taxids']

    taxid_d = parse_taxids(taxid_f)
    print_renamed_blast(blast_f, taxid_d)
    # Check input files
    #if not isfile(fasta_f):
    #    BtLog.error('0', fasta_f)
    #if not isfile(list_f):
    #    BtLog.error('0', list_f)

    #items = BtIO.parseSet(list_f)
    #if not (invert):
    #    for header, sequence in BtIO.readFasta(fasta_f):
    #        if header in items:
    #            print ">%s\n%s\n" % (header, sequence)
    #else:
    #    for header, sequence in BtIO.readFasta(fasta_f):
    #        if not header in items:
    #            print ">%s\n%s\n" % (header, sequence)

