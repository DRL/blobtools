#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools plot    -i BLOBDB [-p INT] [-l INT] [-c] [-n] [-s]
                            [-r RANK] [-t TAXRULE] [--label GROUPS...] 
                            [-o PREFIX] [-m] [--sort ORDER] [--hist HIST] [--title]
                            [-h|--help] 

    Options:
        -h --help                   show this
        -i, --infile BLOBDB         BlobDB file
        -p, --plotgroups INT        Number of (taxonomic) groups to plot [default: 7]
        -l, --length INT            Minimum sequence length considered for plotting [default: 100]
        -c, --cindex                Colour blobs by 'c index' [default: False]
        -n, --nohit                 Hide sequences without taxonomic annotation [default: False]
        -s, --noscale               Do not scale sequences by length [default: False]
        -o, --out PREFIX            Output prefix
        -m, --multiplot             Multi-plot. Print plot after addition of each (taxonomic) group 
                                     [default: False]
        --sort <ORDER>              Sort order for plotting [default: span]
                                     span  : plot with decreasing span
                                     count : plot with decreasing count 
        --hist <HIST>               Data for histograms [default: span] 
                                     span  : span-weighted histograms
                                     count : count histograms
        --title                     Add title of BlobDB to plot [default: False]
        -r, --rank RANK             Taxonomic rank used for colouring of blobs [default: phylum]
                                     (Supported: species, genus, family, order, phylum, superkingdom) 
        -t, --taxrule TAXRULE       Taxrule which has been used for computing taxonomy 
                                     (Supported: bestsum, bestsumorder) [default: bestsum]
        --label GROUPS...           Relabel (taxonomic) groups, 
                                     e.g. "Bacteria=Actinobacteria,Proteobacteria"
        --colours COLOURFILE        File containing colours for (taxonomic) groups
        --include GROUPS...         Include only these (taxonomic) groups,
                                     e.g. "Actinobacteria,Proteobacteria"
        --exclude GROUPS..          Exclude these (taxonomic) groups,
                                     e.g. "Actinobacteria,Proteobacteria"
"""

from __future__ import division
from docopt import docopt
import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtPlot as BtPlot
from os.path import dirname, isfile

if __name__ == '__main__':
    TAXRULES = ['bestsum', 'bestsumorder']
    RANKS = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']
    main_dir = dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    blobdb_f = args['--infile']
    rank = args['--rank'] 
    c_index = args['--cindex']
    min_length = int(args['--length'])
    multiplot = args['--multiplot']
    hide_nohits = args['--nohit']
    out_prefix = args['--out']
    max_taxa_plot = int(args['--plotgroups'])
    sort_order = args['--sort']
    taxrule = args['--taxrule']
    hist_type = args['--hist']
    plot_title = args['--title']
    ignore_contig_length = args['--noscale']
    labels = args['--label']

    # Does blobdb_f exist ?
    if not isfile(blobdb_f):
        BtLog.error('0', blobdb_f)

    # Are ranks sane ?
    if rank not in RANKS:
        BtLog.error('9', rank)

    # Are sort_order and hist_type sane?
    if not sort_order in ['span', 'count']:
        BtLog.error('14', sort_order)
    if not hist_type in ['span', 'count']:            
        BtLog.error('15', hist_type)

    # is taxrule provided?
    if taxrule not in TAXRULES:
        BtLog.error('8', taxrule)

    # compute labels if supplied
    label_d = {}
    if (labels):
        try:
            for cluster in labels:
                name, groups = cluster.split("=")
                for group in groups.split(","):
                    if (group):
                        if group in label_d:
                            BtLog.error('16', group)            
                        label_d[group] = name
        except:
            BtLog.error('17', labels)

    # Load BlobDb
    print BtLog.status_d['9'] % blobdb_f
    blobDB = bt.BlobDb('new')
    blobDB.load(blobdb_f)

    title = blobDB.title
    if plot_title:
        plot_title = title

    # Is taxrule sane and was it computed?
    if taxrule not in blobDB.taxrules:
        BtLog.error('11', taxrule, blobDB.taxrules)

    # blobDB.getArrays(rank, c_index, min_length, multiplot, hide_nohits, out_prefix, max_taxa_plot, sort_order, taxrule, hist_type, plot_title)
    data_array, cov_arrays, summary_dict = blobDB.getArrays(rank, min_length, hide_nohits, taxrule, c_index, label_d)
    plot_order = BtPlot.getPlotOrder(summary_dict, sort_order, max_taxa_plot)
    colour_dict = BtPlot.getColourDict(plot_order)
    min_cov, max_cov = BtPlot.getMinMaxCov(cov_arrays)
    if len(cov_arrays) > 1:
        cov_arrays = BtPlot.getSumCov(cov_arrays)
    info = 1
    for cov_lib in cov_arrays:
        cov_array = cov_arrays[cov_lib]
        out_f = "%s.%s.%s" % (title, cov_lib, hist_type)
        if out_prefix:
            out_f = "%s.%s" % (out_prefix, out_f)
        if c_index:
            out_f = "%s.%s" % (out_f, "c_index")
        if labels:
            out_f = "%s.%s" % (out_f, "label_" + "_".join(set([name for name in label_d.values()])))
        out_f = "%s.%s.%s" % (out_f, min_length, taxrule)
        BtPlot.plot(data_array, cov_array, summary_dict, plot_order, colour_dict, min_cov, max_cov, multiplot, hist_type, plot_title, out_f, ignore_contig_length, info)
        info = 0