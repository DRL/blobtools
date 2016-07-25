#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools view    -i <BLOBDB> [-x <TAXRULE>] [--rank <TAXRANK>...] [--hits]
                            [--list <LIST>] [--out <OUT>] [--notable]
                            [--concoct] [--cov] [--experimental]
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
        --concoct                   Generate concoct files [default: False]
        --cov                       Generate cov files [default: False]
        --experimental              Experimental output [default: False]
        -n, --notable               Do not generate table view [default: False]
"""

from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import blobtools
import lib.BtCore as Bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO

if __name__ == '__main__':
    TAXRULES = ['bestsum', 'bestsumorder']
    RANKS = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom', 'all']
    #print data_dir
    args = docopt(__doc__)
    blobdb_f = args['--input']
    prefix = args['--out']
    ranks = args['--rank']
    taxrule = args['--taxrule']
    hits_flag = args['--hits']
    seq_list_f = args['--list']
    concoct = args['--concoct']
    cov = args['--cov']
    notable = args['--notable']
    experimental = args['--experimental']
    # Does blobdb_f exist ?
    if not isfile(blobdb_f):
        BtLog.error('0', blobdb_f)

    out_f = BtIO.getOutFile(blobdb_f, prefix, None)

    # Are ranks sane ?
    if 'all' in ranks:
        temp_ranks = RANKS[0:-1]
        ranks = temp_ranks[::-1]
    else:
        for rank in ranks:
            if rank not in RANKS:
                BtLog.error('9', rank)

    # Does seq_list file exist?
    seqs = []
    if (seq_list_f):
        if isfile(seq_list_f):
            seqs = BtIO.parseList(seq_list_f)
        else:
            BtLog.error('0', seq_list_f)

    # Load BlobDb
    blobDb = Bt.BlobDb('new')
    print BtLog.status_d['9'] % (blobdb_f)
    blobDb.load(blobdb_f)
    blobDb.version = blobtools.__version__

    # Is taxrule sane and was it computed?
    if (blobDb.hitLibs) and taxrule not in blobDb.taxrules:
        BtLog.error('11', taxrule, blobDb.taxrules)

    # view(s)
    viewObjs = []
    print BtLog.status_d['14']
    if not (notable):
        tableView = Bt.ViewObj(name="table", out_f=out_f, suffix="table.txt", body=[])
        viewObjs.append(tableView)
    if (experimental):
        experimentalView = Bt.ExperimentalViewObj(name = "experimental", view_dir=out_f)
        viewObjs.append(experimentalView)
    if (concoct):
        concoctTaxView = Bt.ViewObj(name="concoct_tax", out_f=out_f, suffix="concoct_taxonomy_info.csv", body=dict())
        viewObjs.append(concoctTaxView)
        concoctCovView = Bt.ViewObj(name="concoct_cov", out_f=out_f, suffix="concoct_coverage_info.tsv", body=[])
        viewObjs.append(concoctCovView)
    if (cov):
        for cov_lib_name, covLibDict in blobDb.covLibs.items():
            out_f = BtIO.getOutFile(covLibDict['f'], prefix, None)
            covView = Bt.ViewObj(name="covlib", out_f=out_f, suffix="cov", body=[])
            blobDb.view(viewObjs=[covView], ranks=None, taxrule=None, hits_flag=None, seqs=None, cov_libs=[cov_lib_name], progressbar=True)
    if (viewObjs):
        blobDb.view(viewObjs=viewObjs, ranks=ranks, taxrule=taxrule, hits_flag=hits_flag, seqs=seqs, cov_libs=[], progressbar=True)
    print BtLog.status_d['19']
