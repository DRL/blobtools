#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: blobtools bamfilter  -b FILE [-i FILE] [-e FILE] [-U] [-n] [-o PREFIX] [-f FORMAT]
                                [-h|--help]

    Options:
        -h --help                   show this
        -b, --bam FILE              BAM file (sorted by name)
        -i, --include FILE          List of contigs whose reads are included
                                    - writes FASTAs of pairs where at least
                                        one read maps sequences in list
                                        (InUn.fq, InIn.fq, ExIn.fq)
        -e, --exclude FILE          List of contigs whose reads are excluded (outputs reads that do not map to sequences in list)
                                    - writes FASTAs of pairs where at least
                                        one read does not maps to sequences in list
                                        (InUn.fq, InIn.fq, ExIn.fq)
        -U, --exclude_unmapped      Include pairs where both reads are unmapped
        -n, --noninterleaved        Use if fw and rev reads should be in separate files  
        -f, --read_format FORMAT    FASTQ = fq, FASTA = fa [default: fa]      
        -o, --out PREFIX            Output prefix
"""

from __future__ import division
from docopt import docopt

import sys

import lib.BtLog as BtLog
import lib.BtIO as BtIO

def main():
    args = docopt(__doc__)
    #print(args)
    bam_f = args['--bam']
    include_f = args['--include']
    exclude_f = args['--exclude']
    out_prefix = args['--out']
    read_format = args['--read_format']
    if not read_format in set(['fq', 'fa']):
        sys.exit("[X] Read format must be fq or fa!")
    noninterleaved = args['--noninterleaved']
    include_unmapped = True
    if args['--exclude_unmapped']:
        include_unmapped = False
    out_f = BtIO.getOutFile(bam_f, out_prefix, None)
    if include_f and exclude_f:
        print(BtLog.error('43'))
    elif include_f:
        sequence_list = BtIO.parseList(include_f)
        BtIO.parseBamForFilter(bam_f, include_unmapped, noninterleaved, out_f, sequence_list, None, read_format)
    elif exclude_f:
        sequence_list = BtIO.parseList(exclude_f)
        BtIO.parseBamForFilter(bam_f, include_unmapped, noninterleaved, out_f, None, sequence_list, read_format)
    else:
        BtIO.parseBamForFilter(bam_f, include_unmapped, noninterleaved, out_f, None, None, read_format)

if __name__ == '__main__':
    main()
