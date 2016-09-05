#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools bamfilter  -b FILE [-i FILE] [-e FILE] [-o PREFIX]
                                [--sort] [--keep] [--threads INT] [--gzip]
                                [-h|--help]

    Options:
        -h --help                   show this
        -b, --bam FILE              BAM file (sorted by name)
        -i, --include FILE          List of contigs whose reads are included
                                        (if no FILE is provided ALL mapped reads
                                            are included)
        -e, --exclude FILE          List of contigs whose reads are excluded
                                        (if no FILE is provided ALL mapped reads
                                            are included)
        --sort                      Sort BAM file by name
        --keep                      Keep sorted BAM file (deleted otherwise)
        --threads INT               number of sorting/compression threads
                                    for sorting [default: 2]
        --gzip                      Gzip output
        -o, --out PREFIX            Output prefix
"""

from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import blobtools
import bloblib.BtLog as BtLog
import bloblib.BtIO as BtIO
import bloblib.BtCore as Bt
import bloblib.BtPlot as BtPlot

def main():
    args = docopt(__doc__)
    bam_f = args['--bam']
    include_f = args['--include']
    exclude_f = args['--exclude']
    out_prefix = args['--out']
    gzip = args['--gzip']
    do_sort = args['--sort']
    keep_sorted = args['--keep']
    sort_threads = int(args['--threads'])

    print BtLog.status_d['22'] % bam_f
    out_f = BtIO.getOutFile(bam_f, out_prefix, None)
    if include_f and exclude_f:
        print BtLog.error('43')
    elif include_f:
        sequence_list = BtIO.parseList(include_f)
        BtIO.parseBamForFilter(bam_f, out_f, sequence_list, None, gzip, do_sort, keep_sorted, sort_threads)
    elif exclude_f:
        sequence_list = BtIO.parseList(exclude_f)
        BtIO.parseBamForFilter(bam_f, out_f, None, sequence_list, gzip, do_sort, keep_sorted, sort_threads)
    else:
        BtIO.parseBamForFilter(bam_f, out_f, None, None, gzip, do_sort, keep_sorted, sort_threads)

if __name__ == '__main__':
    main()
