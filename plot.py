#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools plot    -i BLOBDB [-p INT] [-l INT] [-c] [-n] [-s]
                            [-r RANK] [-x TAXRULE] [--label GROUPS...] 
                            [-o PREFIX] [-m] [--sort ORDER] [--hist HIST] [--title]
                            [--colours FILE] [--include FILE] [--exclude FILE]
                            [--format FORMAT] [--noblobs] [--noreads] 
                            [--refcov FILE] [--catcolour FILE]
                            [-h|--help] 

    Options:
        -h --help                   show this
        -i, --infile BLOBDB         BlobDB file
        -p, --plotgroups INT        Number of (taxonomic) groups to plot, remaining 
                                     groups are placed in 'other' [default: 7]
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
        -x, --taxrule TAXRULE       Taxrule which has been used for computing taxonomy 
                                     (Supported: bestsum, bestsumorder) [default: bestsum]
        --label GROUPS...           Relabel (taxonomic) groups (not 'all' or 'other'), 
                                     e.g. "Bacteria=Actinobacteria,Proteobacteria"
        --colours COLOURFILE        File containing colours for (taxonomic) groups
        --exclude GROUPS..          Place these (taxonomic) groups in 'other',
                                     e.g. "Actinobacteria,Proteobacteria"
        --format FORMAT             Figure format for plot (png, pdf, eps, jpeg, 
                                        ps, svg, svgz, tiff) [default: png]
        --noblobs                   Omit blobplot [default: False]
        --noreads                   Omit plot of reads mapping [default: False]
        --refcov FILE               File containing number of "total" and "mapped" reads 
                                     per coverage file. (e.g.: bam0,900,100). If provided, info
                                     will be used in read coverage plot(s). 
        --catcolour FILE            Colour plot based on categories from FILE 
                                     (format : "seq\tcategory"). 
                                    
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
    max_group_plot = int(args['--plotgroups'])
    sort_order = args['--sort']
    taxrule = args['--taxrule']
    hist_type = args['--hist']
    plot_title = args['--title']
    ignore_contig_length = args['--noscale']
    labels = args['--label']
    colour_f = args['--colours']
    exclude_groups = args['--exclude']
    format = args['--format'] 
    no_plot_blobs = args['--noblobs']
    no_plot_reads = args['--noreads']
    refcov_f = args['--refcov']
    catcolour_f = args['--catcolour']

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
    
    user_labels = BtPlot.parse_labels(labels)
    
    if (exclude_groups):
        if "," in exclude_groups:
            exclude_groups = exclude_groups.rsplit(",")
        else:
            exclude_groups = exclude_groups
    
    refcov_dict = {}
    if (refcov_f):
        refcov_dict = BtPlot.parseRefCov(refcov_f)

    catcolour_dict = {}
    if (catcolour_f) and (c_index):
        BtLog.error('24')
    elif (catcolour_f):
        catcolour_dict = BtPlot.parseCatColour(catcolour_f)
    else: 
        pass

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
    
    data_dict, max_cov, cov_libs, cov_libs_total_reads = blobDB.getPlotData(rank, min_length, hide_nohits, taxrule, c_index, catcolour_dict)
    plotObj = BtPlot.PlotObj(data_dict, cov_libs, cov_libs_total_reads)
    plotObj.exclude_groups = exclude_groups
    plotObj.format = format
    plotObj.max_cov = max_cov
    plotObj.title = title
    plotObj.multiplot = multiplot
    plotObj.hist_type = hist_type
    plotObj.ignore_contig_length = ignore_contig_length
    plotObj.max_group_plot = max_group_plot
    plotObj.group_order = BtPlot.getSortedGroups(data_dict, sort_order)
    plotObj.labels.update(plotObj.group_order)
    
    if (user_labels):
        for group, label in user_labels.items():
            plotObj.labels.add(label)
    plotObj.group_labels = {group : set() for group in plotObj.group_order}
    plotObj.relabel_and_colour(colour_f, user_labels)
    plotObj.compute_stats()

    info_flag = 1

    for cov_lib in plotObj.cov_libs:
        if (plotObj.title):
            plotObj.title = "%s.%s" % (plotObj.title, cov_lib)

        out_f = "%s.%s.%s.p%s.%s" % (title, hist_type, rank, max_group_plot, cov_lib)
        if out_prefix:
            out_f = "%s.%s" % (out_prefix, out_f)
        if catcolour_dict:
            out_f = "%s.%s" % (out_f, "catcolour")
        if ignore_contig_length:
            out_f = "%s.%s" % (out_f, "noscale")
        if c_index:
            out_f = "%s.%s" % (out_f, "c_index")
        if exclude_groups:
            out_f = "%s.%s" % (out_f, "exclude" + "_".join(exclude_groups))
        if labels:
            out_f = "%s.%s" % (out_f, "label_" + "_".join(set([name for name in user_labels.values()])))
        out_f = "%s.%s.%s" % (out_f, min_length, taxrule)
        plotObj.out_f = out_f
        if not (no_plot_blobs):
            plotObj.plotBlobs(cov_lib, info_flag)
            info_flag = 0
    plotObj.write_stats()

    if not (no_plot_reads):
        plotObj.plotReadCov(refcov_dict)
