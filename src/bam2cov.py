#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools bam2cov         -i FASTA [-b BAM...] [-a CAS...] [-s SAM...]
                                    [-o PREFIX] [--no_base_cov]
                                    [-h|--help]

    Options:
        -h --help                   show this
        -i, --infile FASTA          FASTA file of assembly. Headers are split at whitespaces.
        -b, --bam <BAM>...          BAM file (requires samtools in $PATH)
        -a, --cas <CAS>...          CAS file (requires clc_mapping_info in $PATH)
        -s, --sam <SAM>...          SAM file
        -o, --output <PREFIX>       Output prefix
        --no_base_cov               only parse read coverage (faster, but ...
                                        can only be used for "blobtools blobplot --noblobs")
"""

from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
#path.append(dirname(dirname(abspath(__file__))))

import blobtools
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtCore as Bt

if __name__ == '__main__':
    args = docopt(__doc__)
    fasta_f = args['--infile']
    bam_fs = args['--bam']
    cas_fs = args['--cas']
    sam_fs = args['--sam']
    prefix = args['--output']
    no_base_cov_flag = args['--no_base_cov']

    # Make covLibs
    cov_libs = [Bt.CovLibObj('bam' + str(idx), 'bam', lib_f) for idx, lib_f in enumerate(bam_fs)] + \
           [Bt.CovLibObj('sam' + str(idx), 'sam', lib_f) for idx, lib_f in enumerate(sam_fs)] + \
           [Bt.CovLibObj('cas' + str(idx), 'cas', lib_f) for idx, lib_f in enumerate(cas_fs)]
    blobDb = Bt.BlobDb('cov')
    blobDb.version = blobtools.__version__
    blobDb.parseFasta(fasta_f, None)
    blobDb.parseCoverage(covLibObjs=cov_libs, no_base_cov=no_base_cov_flag)
    for cov_lib in cov_libs:
        out_f = BtIO.getOutFile(cov_lib.f, prefix, None)
        covView = Bt.ViewObj(name="cov", out_f=out_f, suffix="cov", header="", body=[])
        blobDb.view(viewObjs=[covView], ranks=None, taxrule=None, hits_flag=None, seqs=None, cov_libs=[cov_lib.name])

