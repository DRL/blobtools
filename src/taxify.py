#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: taxify                -i BLAST -t TAXIDS [-o PREFIX] [--diamond] [--blast] [-h|--help]

    Options:
        -h --help                     show this
        -i, --infile <TAXFILE>        Similarity search results of sequences (has to be in BLAST/Diamond format)
        -t, --taxids <TAXIDS>         TAXID file
        -o, --out <PREFIX>            Output prefix
        -e, --evalue <EVALUE>         E-value cutoff (default: 1e-25)
        -b, --bitscore <BITSCORE>
        --rnacentral                  TAXFILE is BLAST and TAXIDS is from rnacentral (SILVA analysis)
        --diamond                     TAXFILE is Diamond output (*.daa)
"""

from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))
import lib.BtLog as BtLog
import lib.BtIO as BtIO

def write_output(blast_f, taxid_d):
    with open(blast_f) as fh:
        for l in fh:
            line = l.rstrip("\n").split()
            print "%s\t%s\t%s" % (line[0], taxid_d[line[3]], "\t".join(line[2:]))

if __name__ == '__main__':
    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    tax_f = args['--infile']
    taxid_f = args['--taxids']
    prefix = args['--out']
    diamond = args['--diamond']
    rnacentral = args['--rnacentral']

    out_f, taxid_d = '', {}
    if (rnacentral):
        out_f = BtIO.getOutFile(tax_f, prefix, "rnacentral.out")
        taxid_d = BtIO.parseDict(taxid_f, 0, 3)
    elif (diamond):
        out_f = BtIO.getOutFile(tax_f, prefix, "diamond.out")
        taxid_d = BtIO.parseDict(taxid_f, 0, 1)
    else:
        BtLog.error('26')

    write_output(tax_f, taxid_d)



