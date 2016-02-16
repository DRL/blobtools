#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools comparecov  -i BLOBDB -c COV [-p INT] [-l INT] [-n] [-s]
                                [--xlabel XLABEL] [--ylabel YLABEL]
                                [--log] [--xmax FLOAT] [--ymax FLOAT]
                                [-r RANK] [-x TAXRULE] [-o PREFIX] [-m] [--title]
                                [--sort ORDER] [--hist HIST] [--format FORMAT]
                                [-h|--help] 

    Options:
        -h --help                   show this
        -i, --infile BLOBDB         BlobDB file
        -c, --cov COV               COV file used for y-axis
        
        --xlabel XLABEL             Label for x-axis [default: BlobDB_cov]
        --ylabel YLABEL             Label for y-axis [default: CovFile_cov]
        --log                       Plot log-scale axes
        --xmax FLOAT                Maximum values for x-axis [default: 1e10]
        --ymax FLOAT                Maximum values for y-axis [default: 1e10]

        -p, --plotgroups INT        Number of (taxonomic) groups to plot, remaining 
                                     groups are placed in 'other' [default: 7]
        -r, --rank RANK             Taxonomic rank used for colouring of blobs [default: phylum]
        -x, --taxrule TAXRULE       Taxrule which has been used for computing taxonomy 
                                     (Supported: bestsum, bestsumorder) [default: bestsum]
        --sort <ORDER>              Sort order for plotting [default: span]
                                     span  : plot with decreasing span
                                     count : plot with decreasing count 
        --hist <HIST>               Data for histograms [default: span] 
                                     span  : span-weighted histograms
                                     count : count histograms

        --title                     Add title of BlobDB to plot [default: False]
        -l, --length INT            Minimum sequence length considered for plotting [default: 100]
        -n, --nohit                 Hide sequences without taxonomic annotation [default: False]
        -s, --noscale               Do not scale sequences by length [default: False]
        -o, --out PREFIX            Output prefix
        -m, --multiplot             Multi-plot. Print plot after addition of each (taxonomic) group 
                                     [default: False]
        --format FORMAT             Figure format for plot (png, pdf, eps, jpeg, 
                                        ps, svg, svgz, tiff) [default: png]
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
    cov_f = args['--cov']
    x_label = args['--xlabel'] 
    y_label = args['--ylabel']
    scale = args['--log'] 
    x_max = float(args['--xmax'])
    y_max = float(args['--ymax'])
    rank = args['--rank'] 
    min_length = int(args['--length'])
    multiplot = args['--multiplot']
    hide_nohits = args['--nohit']
    out_prefix = args['--out']
    max_group_plot = int(args['--plotgroups'])
    sort_order = args['--sort']
    taxrule = args['--taxrule']
    hist_type = args['--hist']
    plot_title = args['--title']
    ignore_contig_length = args['--noscale']
    #labels = args['--label']
    #colour_f = args['--colours']
    #exclude_groups = args['--exclude']
    format = args['--format'] 
    #no_plot_blobs = args['--noblobs']
    #no_plot_reads = args['--noreads']
    #refcov_f = args['--refcov']
    #catcolour_f = args['--catcolour']

    # Does blobdb_f exist ?
    if not isfile(blobdb_f):
        BtLog.error('0', blobdb_f)

    # Does cov_f exist ?
    if not isfile(cov_f):
        BtLog.error('0', cov_f)
    # parse cov file in dict 
    cov_dict = BtPlot.parseCovFile(cov_f)
    
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
    
    #user_labels = BtPlot.parse_labels(labels)
    #
    #if (exclude_groups):
    #    if "," in exclude_groups:
    #        exclude_groups = exclude_groups.rsplit(",")
    #    else:
    #        exclude_groups = exclude_groups
    #
    #refcov_dict = {}
    #if (refcov_f):
    #    refcov_dict = BtPlot.parseRefCov(refcov_f)
#
    #catcolour_dict = {}
    #if (catcolour_f) and (c_index):
    #    BtLog.error('24')
    #elif (catcolour_f):
    #    catcolour_dict = BtPlot.parseCatColour(catcolour_f)
    #else: 
    #    pass

    # Load BlobDb
    print BtLog.status_d['9'] % blobdb_f
    blobDB = bt.BlobDb('new')
    blobDB.load(blobdb_f)

    # clean cov_dict from coverages below
    #for name in cov_dict:
    #    print name, blobDB.dict_of_blobs[name]

    title = blobDB.title
    if plot_title:
        plot_title = title


    # Is taxrule sane and was it computed?
    if taxrule not in blobDB.taxrules:
        BtLog.error('11', taxrule, blobDB.taxrules)
    
    data_dict, max_cov, cov_libs, cov_libs_total_reads = blobDB.getPlotData(rank, min_length, hide_nohits, taxrule, False, False)
    plotObj = BtPlot.PlotObj(data_dict, cov_libs, cov_libs_total_reads)
    #plotObj.exclude_groups = exclude_groups
    if max_cov < x_max:
        x_max = max_cov
    if max_cov < y_max:
        y_max = max_cov
    
    if (scale):
        scale = 'log'
    else:
        scale = 'linear'

    plotObj.max_cov = max_cov
    plotObj.title = title
    plotObj.format = format
    plotObj.multiplot = multiplot
    plotObj.hist_type = hist_type
    plotObj.ignore_contig_length = ignore_contig_length
    plotObj.max_group_plot = max_group_plot
    plotObj.group_order = BtPlot.getSortedGroups(data_dict, sort_order)
    plotObj.labels.update(plotObj.group_order)
    
    #if (user_labels):
    #    for group, label in user_labels.items():
    #        plotObj.labels.add(label)
    plotObj.group_labels = {group : set() for group in plotObj.group_order}
    plotObj.relabel_and_colour(None, {})
    plotObj.compute_stats()

    info_flag = 1

    for cov_lib in plotObj.cov_libs:
        if (plotObj.title):
            plotObj.title = "%s.%s.%s" % (title, taxrule, cov_lib)

        out_f = "%s.%s.%s.p%s.%s.%s" % (title, hist_type, rank, max_group_plot, scale, cov_lib)
        if out_prefix:
            out_f = "%s.%s" % (out_prefix, out_f)
        #if catcolour_dict:
        #    out_f = "%s.%s" % (out_f, "catcolour")
        if ignore_contig_length:
            out_f = "%s.%s" % (out_f, "noscale")
        #if c_index:
        #    out_f = "%s.%s" % (out_f, "c_index")
        #if exclude_groups:
        #    out_f = "%s.%s" % (out_f, "exclude" + "_".join(exclude_groups))
        #if labels:
        #    out_f = "%s.%s" % (out_f, "label_" + "_".join(set([name for name in user_labels.values()])))
        out_f = "%s.%s.%s" % (out_f, min_length, taxrule)
        plotObj.out_f = out_f
        
        plotObj.plotScatterCov(cov_lib, cov_dict, info_flag, x_label, y_label, scale, x_max, y_max)
        info_flag = 0
    plotObj.write_stats()
