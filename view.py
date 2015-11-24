#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools view    -i <BLOBDB> [-r <TAXRULE>] [--rank <TAXRANK>...] [--hits]
                            [--list <LIST>] [--out <OUT>]
                            [--h|--help] 
    
    Options:
        --h --help                  show this
        -i, --input <BLOBDB>        BlobDB file (created with "blobtools forge")
        -o, --out <OUT>             Output file [default: STDOUT]
        -l, --list <LIST>           List of sequence names (comma-separated or file). 
                                    If comma-separated, no whitespaces allowed.
        -r, --taxrule <TAXRULE>     Taxrule used for computing taxonomy (supported: "bestsum", "bestsumorder")
                                    [default: bestsum]
        --rank <TAXRANK>...         Taxonomic rank(s) at which output will be written. 
                                    (supported: 'species', 'genus', 'family', 'order', 
                                    'phylum', 'superkingdom', 'all') [default: phylum]
        -b, --hits                  Displays taxonomic hits from tax files
"""

from __future__ import division
from docopt import docopt
import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
from os.path import basename, isfile, join, dirname


if __name__ == '__main__':
    TAXRULES = ['bestsum', 'bestsumorder']
    RANKS = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom', 'all']

    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    blobdb_f = args['--input']
    out_f = args['--out'] 
    ranks = args['--rank']
    taxrule = args['--taxrule']
    hits_flag = args['--hits']
    seq_list = args['--list']

    # Does blobdb_f exist ?
    if not isfile(blobdb_f):
        BtLog.error('0', blobdb_f)

    # Are ranks sane ?
    for rank in ranks:
        if rank not in RANKS:
            BtLog.error('9', rank)
    if 'all' in ranks:
        ranks = RANKS[0:-1]            

    # Is list a list of sequence names or a file?
    seqs = []
    if (seq_list):
        if isfile(seq_list):
            seqs = BtIO.parseList(seq_list)
        elif "," in seq_list:
            seqs = seq_list.split(",")
        else:
            seqs = [seq_list]

    # Load BlobDb
    blobDB = bt.BlobDb('new')
    blobDB.load(blobdb_f)

    # Is taxrule sane and was it computed?
    if (blobDB.hitLibs) and taxrule not in blobDB.taxrules:
        BtLog.error('11', taxrule, blobDB.taxrules)
    blobDB.view(out_f, ranks, taxrule, hits_flag, seqs)

