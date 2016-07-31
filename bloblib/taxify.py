#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools taxify          -i BLAST [-d FILE] [-r FILE] [-t INT]
                                    [-o PREFIX] [--force] [-h|--help]

    Options:
        -h --help                     show this
        -i, --infile <TAXFILE>        Similarity search results of sequences
                                        (BLAST/Diamond TSV format)
        -d, --diamond <TAXIDS>        Diamond TAXID file
        -r, --rnacentral <TAXIDS>     RNAcentral TAXID file
        -t, --taxid <INT>             TaxID (must be integer)
        --force                       Overwrite existing TaxIDs
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
    force = args['--force']

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
        print BtLog.status_d['1'] % ("TAXID file", diamond_f)
        taxid_d = BtIO.parseDict(diamond_f, 0, 1)
        out_f = BtIO.getOutFile(tax_f, prefix, "diamond.out")
    else:
        BtLog.error('26')

    output = []
    print BtLog.status_d['1'] % ("taxonomy file", tax_f)
    with open(tax_f) as fh:
        for idx, l in enumerate(fh):
            line = l.rstrip("\n").split()
            if diamond_f:
                output.append("%s\t%s\t%s\t%s" % (line[0], taxid_d[line[1]], line[11], "\t".join(line[1:])))
            else:
                if line[1] == 'N/A' or force: # so that it does not overwrite existing taxIDs
                    output.append("%s\t%s\t%s" % (line[0], taxid_d[line[3]], "\t".join(line[2:5])))
                else:
                    print BtLog.warn_d['10'] % (idx+1, line[0], line[1])
                    output.append("%s\t%s\t%s" % (line[0], line[1], "\t".join(line[2:5])))

    if output:
        with open(out_f, "w") as fh:
            print BtLog.status_d['24'] % out_f
            fh.write("\n".join(output))

if __name__ == '__main__':
    main()



