#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools bamfilter  -b FILE [-i FILE] [-e FILE] [-u] [-o PREFIX]
                                [--sort] [--keep] [--threads INT]
                                [-h|--help]

    Options:
        -h --help                   show this
        -b, --bam FILE              BAM file (sorted by name)
        -i, --include FILE          List of contigs whose reads are included
                                    - writes interleaved FASTQs of pairs where at least
                                        one read maps sequences in list
                                        (InUn.fq, InIn.fq, ExIn.fq)
        -e, --exclude FILE          List of contigs whose reads are excluded (outputs reads that do not map to sequences in list)
                                    - writes interleaved FASTQs of pairs where at least
                                        one read does not maps to sequences in list
                                        (InUn.fq, InIn.fq, ExIn.fq)
        -u, --include_unmapped      Include pairs where both reads are unmapped
        --sort                      Sort BAM file by name
        --keep                      Keep sorted BAM file (deleted otherwise)
        --threads INT               number of sorting/compression threads
                                    for sorting [default: 2]
        -o, --out PREFIX            Output prefix
"""

from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import lib.blobtools as blobtools
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtCore as Bt
import lib.BtPlot as BtPlot

def main():
    args = docopt(__doc__)
    bam_f = args['--bam']
    include_f = args['--include']
    exclude_f = args['--exclude']
    out_prefix = args['--out']
    include_unmapped = args['--include_unmapped']
    gzip = None
    do_sort = args['--sort']
    keep_sorted = args['--keep']
    sort_threads = int(args['--threads'])
    out_f = BtIO.getOutFile(bam_f, out_prefix, None)
    if include_f and exclude_f:
        print BtLog.error('43')
    elif include_f:
        sequence_list = BtIO.parseList(include_f)
        BtIO.parseBamForFilter(bam_f, include_unmapped, out_f, sequence_list, None, gzip, do_sort, keep_sorted, sort_threads)
    elif exclude_f:
        sequence_list = BtIO.parseList(exclude_f)
        BtIO.parseBamForFilter(bam_f, include_unmapped, out_f, None, sequence_list, gzip, do_sort, keep_sorted, sort_threads)
    else:
        BtIO.parseBamForFilter(bam_f, include_unmapped, out_f, None, None, gzip, do_sort, keep_sorted, sort_threads)

if __name__ == '__main__':
    main()
