#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: blobtools map2cov         -i FASTA [-b BAM...] [-a CAS...]
                                    [-o PREFIX] [-c]
                                    [-h|--help]

    Options:
        -h --help                   show this
        -i, --infile FASTA          FASTA file of assembly. Headers are split at whitespaces.
        -b, --bam <BAM>...          BAM file (requires pysam)
        -a, --cas <CAS>...          CAS file (requires clc_mapping_info in $PATH)
        -o, --output <PREFIX>       Output prefix
        -c, --calculate_cov         Legacy coverage, slower. New default is to estimate coverages 
                                        based on read lengths of first 10K reads.
"""

from __future__ import division
from docopt import docopt

import lib.BtLog as BtLog
import lib.BtCore as BtCore
import lib.interface as interface

def main():
    args = docopt(__doc__)
    fasta_f = args['--infile']
    bam_fs = args['--bam']
    cas_fs = args['--cas']
    prefix = args['--output']
    estimate_cov_flag = True if not args['--calculate_cov'] else False

    # Make covLibs
    cov_libs = [BtCore.CovLibObj('bam' + str(idx), 'bam', lib_f) for idx, lib_f in enumerate(bam_fs)] + \
           [BtCore.CovLibObj('cas' + str(idx), 'cas', lib_f) for idx, lib_f in enumerate(cas_fs)]
    if not (cov_libs):
        BtLog.error('31')
    blobDb = BtCore.BlobDb('cov')
    blobDb.version = interface.__version__
    blobDb.parseFasta(fasta_f, None)
    blobDb.parseCoverage(covLibObjs=cov_libs, estimate_cov=estimate_cov_flag, prefix=prefix)

if __name__ == '__main__':
    main()
