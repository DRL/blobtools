#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools taxify          -i BLAST [-t TAXIDS] [--taxid INT]
                                    [-o PREFIX] [--diamond] [--rnacentral] [-h|--help]

    Options:
        -h --help                     show this
        -i, --infile <TAXFILE>        Similarity search results of sequences (BLAST/Diamond TSV format)
        -t, --taxids <TAXIDS>         TAXID file
        --taxid <INT>                 TAXID
        -o, --out <PREFIX>            Output prefix
        --rnacentral                  TAXFILE is BLAST and TAXIDS is from rnacentral (SILVA analysis)
        --diamond                     TAXFILE is Diamond output (*.daa)
"""

from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))
from collections import defaultdict
import lib.BtLog as BtLog
import lib.BtIO as BtIO

def main():
    #print data_dir
    args = docopt(__doc__)
    tax_f = args['--infile']
    taxid_f = args['--taxids']

    prefix = args['--out']
    diamond = args['--diamond']
    rnacentral = args['--rnacentral']

    try:
        taxid = int(args['--taxid'])
    except TypeError:
        BtLog.error('26')

    out_f, taxid_d = '', {}
    if (taxid_f):
        print BtLog.status_d['1'] % ("TAXID file", taxid_f)
        if (rnacentral):
            taxid_d = BtIO.parseDict(taxid_f, 0, 3)
            out_f = BtIO.getOutFile(tax_f, prefix, "rnacentral.out")
        elif (diamond):
            taxid_d = BtIO.parseDict(taxid_f, 0, 1)
            out_f = BtIO.getOutFile(tax_f, prefix, "diamond.out")
        else:
            BtLog.error('26')
    elif (taxid):
        out_f = BtIO.getOutFile(tax_f, prefix, "taxified.out")
        taxid_d = defaultdict(lambda: taxid)
    else:
        BtLog.error('26')


    output = []
    with open(tax_f) as fh:
        for l in fh:
            line = l.rstrip("\n").split()
            output.append("%s\t%s\t%s" % (line[0], taxid_d[line[3]], "\t".join(line[2:])))

    with open(out_f, "w") as fh:
        print BtLog.status_d['24'] % out_f
        fh.write("\n".join(output))

if __name__ == '__main__':
    main()



