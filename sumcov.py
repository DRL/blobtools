#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools sumcov          -i FASTA -c COV... [-o PREFIX]
                                    [-h|--help]

    Options:
        -h --help                   show this

        -i, --infile <FASTA>          FASTA file of assembly. Headers are split at whitespaces.
        -c, --cov <COV>               COV file(s) to sum up
        -o, --out <PREFIX>            Output prefix
"""

from __future__ import division
from docopt import docopt
import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtPlot as BtPlot
from os.path import dirname, isfile

def parseFasta(infile):
    fasta_order = []
    fasta_dict = {}
    for name in BtIO.readFasta(infile):
        fasta_order.append(name)
        fasta_dict[name] = 0.0
    return fasta_dict, fasta_order

if __name__ == '__main__':
    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    assembly_f = args['--infile']
    cov_fs = args['--cov']
    out_prefix = args['--outprefix']

    if out_prefix:
        out_f = "%s.sum.cov" % (out_prefix)
    else:
        out_f = "%s.sum.cov" % (splitext(fasta_f)[0])

    fasta_dict = {}
    fasta_order = []

    if not isfile(assembly_f):
        BtLog.error('0', assembly_f)
    else:
        fasta_dict, fasta_order = parseFasta(assembly_f)

    for cov_f in cov_fs:
        if not isfile(cov_f):
            BtLog.error('0', cov_f)
        else:
            lib_cov_dict = BtPlot.parseCovFile(cov_f)
            for name in fasta_order:
                fasta_dict[name] = fasta_dict.get(name, 0.0) + lib_cov_dict[name]


    for name in fasta_order:
        print "%s\t%s" % (name, fasta_dict[name])
