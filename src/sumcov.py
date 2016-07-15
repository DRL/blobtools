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
from os.path import dirname, isfile, splitext

if __name__ == '__main__':
    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    assembly_f = args['--infile']
    cov_fs = args['--cov']
    prefix = args['--out']

    out_f = BtIO.getOutFile(assembly_f, prefix, "sum.cov")

    fasta_order = BtIO.parseFastaNameOrder(assembly_f)

    sum_of_reads_mapped = 0
    sum_of_reads_total = 0
    sum_of_base_cov = {}
    sum_of_read_cov = {}

    for cov_f in cov_fs:
        print BtLog.status_d['1'] % ("Cov file", cov_f)
        base_cov_dict, reads_total, reads_mapped, reads_unmapped, read_cov_dict = BtIO.parseCov(cov_f, set(fasta_order))
        sum_of_reads_total += reads_total
        sum_of_reads_mapped += reads_mapped
        for name in fasta_order:
            sum_of_base_cov[name] = sum_of_base_cov.get(name, 0.0) + base_cov_dict[name]
            sum_of_read_cov[name] = sum_of_read_cov.get(name, 0) + read_cov_dict[name]

    source = ",".join(cov_fs)
    BtIO.writeCov(source, fasta_order, sum_of_base_cov, sum_of_reads_total, sum_of_reads_mapped, sum_of_read_cov, out_f)
