#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools taxify          -i BLAST [-d FILE] [-r FILE] [-t INT]
                                    [-o PREFIX] [-h|--help]

    Options:
        -h --help                     show this
        -i, --infile <TAXFILE>        Similarity search results of sequences
                                        (BLAST/Diamond TSV format)
        -d, --diamond <TAXIDS>        Diamond TAXID file
        -r, --rnacentral <TAXIDS>     RNAcentral TAXID file
        -t, --taxid <INT>             TAXID (must be integer)
        -o, --out <PREFIX>            Output prefix
"""

from __future__ import division
from docopt import docopt
from collections import defaultdict

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import bloblib.BtLog as BtLog
import bloblib.BtIO as BtIO

def main():
    args = docopt(__doc__)
    tax_f = args['--infile']
    prefix = args['--out']
    diamond_f = args['--diamond']
    rnacentral_f = args['--rnacentral']
    taxid = args['--taxid']

    out_f, taxid_d = '', {}
    if (taxid):
        try:
            taxid = int(taxid)
        except TypeError:
            BtLog.error('26')
        out_f = BtIO.getOutFile(tax_f, prefix, "tax_%s.out" % taxid)
        taxid_d = defaultdict(lambda: taxid)
    elif (rnacentral_f):
        print BtLog.status_d['1'] % ("TAXID file", rnacentral_f)
        taxid_d = BtIO.parseDict(rnacentral_f, 0, 3)
        out_f = BtIO.getOutFile(tax_f, prefix, "rnacentral.out")
    elif (diamond_f):
        taxid_d = BtIO.parseDict(rnacentral_f, 0, 1)
        out_f = BtIO.getOutFile(diamond_f, prefix, "diamond.out")
    else:
        BtLog.error('26')

    output = []
    with open(tax_f) as fh:
        for l in fh:
            line = l.rstrip("\n").split()
            if (diamond_f):
                output.append("%s\t%s\t%s" % (line[0], taxid_d[line[1]], line[11], "\t".join(line[2:])))
            else:
                output.append("%s\t%s\t%s" % (line[0], taxid_d[line[3]], "\t".join(line[2:])))
    with open(out_f, "w") as fh:
        print BtLog.status_d['24'] % out_f
        fh.write("\n".join(output))

if __name__ == '__main__':
    main()



