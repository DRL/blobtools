#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: taxify                -i BLAST -t TAXIDS [-o PREFIX] [--diamond] [--blast] [-h|--help]

    Options:
        -h --help                     show this
        -i, --infile <TAXFILE>        Similarity search results of sequences (has to be in BLAST/Diamond format)
        -t, --taxids <TAXIDS>         TAXID file
        -o, --out <PREFIX>            Output prefix
        -e, --evalue <EVALUE>         E-value cutoff (default: 1e-25)
        --rnacentral                  TAXFILE is BLAST and TAXIDS is from rnacentral (SILVA analysis)
        --diamond                     TAXFILE is Diamond output (*.daa)
"""

from __future__ import division
from docopt import docopt
#import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
#import lib.BtPlot as BtPlot
from os.path import dirname, isfile, splitext

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
    out_prefix = args['--out']
    diamond = args['--diamond']
    rnacentral = args['--rnacentral']

    # Check input
    if not isfile(tax_f):
        BtLog.error('0', tax_f)
    if not isfile(taxid_f):
        BtLog.error('0', taxid_f)

    out_f = ''
    if out_prefix:
        out_f = out_prefix
    else:
        out_f = "%s" % (splitext(fasta_f)[0])

    if (rnacentral):
        out_f = "%s.%s" % (out_f, 'rnacentral.out')
        taxid_d = BtIO.parseDict(taxid_f, 0, 3)
    elif (diamond):
        out_f = "%s.%s" % (out_f, 'diamond.out')
        taxid_d = BtIO.parseDict(taxid_f, 0, 1)
    else:
        BtLog.error('26')

    write_output(blast_f, taxid_d)



