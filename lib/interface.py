#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
blobtools

    Usage:
        ./blobtools [<module>] [<args>...] [-h] [-v]        

    Modules:
        create        create a BlobDB
        view          generate tabular view, CONCOCT input or COV files from BlobDB
        plot          generate a BlobPlot from a BlobDB
        covplot       generate a CovPlot from a BlobDB and a COV file
    
        map2cov       generate a COV file from BAM file
        taxify        generate a BlobTools compatible HITS file (TSV)
        bamfilter     subset paired-end reads from a BAM file
        seqfilter     subset sequences in FASTA file based sequence IDs in list
        nodesdb       create nodesdb based on NCBI Taxdump's names.dmp and nodes.dmp

    Options
        -h, --help      show this
        -v, --version   show version number

    See 'blobtools <command> --help' for more information on a specific command.

    Further documentation is available at https://blobtools.readme.io/

    Examples:

        # 1. Create a BlobDB
        ./blobtools create -i example/assembly.fna -b example/mapping_1.bam -t example/blast.out -o example/test
    
        # 2. Generate a tabular view
        ./blobtools view -i example/test.blobDB.json
    
        # 3. Generate a blobplot
        ./blobtools plot -i example/test.blobDB.json

"""

import sys
from docopt import docopt, DocoptExit
from timeit import default_timer as timer

__version__ = '1.1.1'

def main():
    try:
        start_time = timer()
        try:
            args = docopt(__doc__, version=__version__, options_first=True)
        except DocoptExit:
            print(__doc__)
        else:
            if args['<module>']:
                if args['<module>'] == 'create':
                    import lib.create as create
                    create.main()
                elif args['<module>'] == 'view':
                    import lib.view as view
                    view.main()
                elif args['<module>'] == 'plot':
                    import lib.blobplot as plot
                    plot.main()
                elif args['<module>'] == 'map2cov':
                    import lib.map2cov as map2cov
                    map2cov.main()
                elif args['<module>'] == 'seqfilter':
                    import lib.seqfilter as seqfilter
                    seqfilter.main()
                elif args['<module>'] == 'covplot':
                    import lib.covplot as covplot
                    covplot.main()
                elif args['<module>'] == 'taxify':
                    import lib.taxify as taxify
                    taxify.main()
                elif args['<module>'] == 'bamfilter':
                    import lib.bamfilter as bamfilter
                    bamfilter.main()
                elif args['<module>'] == 'nodesdb':
                    import lib.nodesdb as nodesdb
                    nodesdb.main()
                else:
                    sys.exit("%r is not a blobtools module. See 'blobtools -h'." % args['<module>'])
            else:

                print(__doc__)
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)

