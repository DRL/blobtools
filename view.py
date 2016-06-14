#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools view    -i <BLOBDB> [-x <TAXRULE>] [--rank <TAXRANK>...] [--hits]
                            [--list <LIST>] [--out <OUT>] [--concoct] [--notable]
                            [--h|--help]

    Options:
        --h --help                  show this
        -i, --input <BLOBDB>        BlobDB file (created with "blobtools create")
        -o, --out <OUT>             Output prefix
        -l, --list <LIST>           List of sequence names (file).
        -x, --taxrule <TAXRULE>     Taxrule used for computing taxonomy (supported: "bestsum", "bestsumorder")
                                    [default: bestsum]
        -r, --rank <TAXRANK>...     Taxonomic rank(s) at which output will be written.
                                    (supported: 'species', 'genus', 'family', 'order',
                                    'phylum', 'superkingdom', 'all') [default: phylum]
        -b, --hits                  Displays taxonomic hits from tax files that contributed to the taxonomy.
        -c, --concoct               Generate concoct files [default: False]
        -n, --notable               Do not generate table view [default: False]
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
    out_prefix = args['--out']
    ranks = args['--rank']
    taxrule = args['--taxrule']
    hits_flag = args['--hits']
    seq_list_f = args['--list']
    concoct = args['--concoct']
    notable = args['--notable']

    # Does blobdb_f exist ?
    if not isfile(blobdb_f):
        BtLog.error('0', blobdb_f)

    out_f = blobdb_f
    # Was output prefix supplied?
    if out_prefix:
        out_f = "%s.%s" % (out_prefix, basename(out_f))

    # Are ranks sane ?
    for rank in ranks:
        if rank not in RANKS:
            BtLog.error('9', rank)
    if 'all' in ranks:
        temp_ranks = RANKS[0:-1]
        ranks = temp_ranks[::-1]

    # Does seq_list file exist?
    seqs = []
    if (seq_list_f):
        if isfile(seq_list_f):
            seqs = BtIO.parseList(seq_list_f)
        else:
            BtLog.error('0', seq_list_f)

    # Load BlobDb
    blobDB = bt.BlobDb('new')
    blobDB.load(blobdb_f)

    # Is taxrule sane and was it computed?
    if (blobDB.hitLibs) and taxrule not in blobDB.taxrules:
        BtLog.error('11', taxrule, blobDB.taxrules)

    # view(s)
    views = []
    if not (notable):
        tableView = bt.ViewObj(name="table", out_f=out_f, suffix="table.txt", header="", body="")
        views.append(tableView)
    if (concoct):
        concoctTaxView = bt.ViewObj(name="concoct_tax", out_f=out_f, suffix="concoct_taxonomy_info.csv", header="", body=dict())
        views.append(concoctTaxView)
        concoctCovView = bt.ViewObj(name="concoct_cov", out_f=out_f, suffix="concoct_coverage_info.tsv", header="", body="")
        views.append(concoctCovView)
    blobDB.view(views, ranks, taxrule, hits_flag, seqs)
