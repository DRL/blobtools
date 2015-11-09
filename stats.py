#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Usage:
    blobtools.py stats  --i <BLOBDB> [--rank <RANK>]
                        [--h|--help] 

    Options:
        --h --help              show this
        --i <BLOBDB>            BlobDB file
        --rank <RANK>           Taxonomic rank used colouring of blobs [default: phylum]
                                (Supported: species, genus, family, order, phylum, superkingdom) 
"""

from __future__ import division
from docopt import docopt
import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtPlot as BtPlot
from os.path import basename, isfile, join, dirname

if __name__ == '__main__':
    RANKS = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']
    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)

    blobdb_f = args['--i']
    rank = args['--rank'] 

    # Does blobdb_f exist ?
    if not isfile(blobdb_f):
        BtLog.error('0', blobdb_f)

    # Are ranks sane ?
    if rank not in RANKS:
        BtLog.error('9', rank)

    # Load BlobDb
    print BtLog.status_d['9'] % blobdb_f
    blobDB = bt.BlobDb('new')
    blobDB.load(blobdb_f)
