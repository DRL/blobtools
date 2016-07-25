#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools covplot  -i BLOBDB -c COV [--max FLOAT]
                                [--xlabel XLABEL] [--ylabel YLABEL]
<<<<<<< HEAD
                                [--lib COVLIB] [-o PREFIX] [-m]
                                [-p INT] [-l INT] [--cindex] [-n] [-s]
                                [-r RANK] [-x TAXRULE] [--label GROUPS...]
                                [--sort ORDER] [--hist HIST] [--notitle]
                                [--colours FILE] [--include FILE] [--exclude FILE]
                                [--refcov FILE] [--catcolour FILE]
                                [--format FORMAT] [--noblobs] [--noreads] [--legend]
                                [--cumulative]
=======
                                [--max FLOAT]
                                [-r RANK] [-x TAXRULE] [-o PREFIX] [-m] [--title]
                                [--sort ORDER] [--hist HIST] [--format FORMAT]
>>>>>>> FETCH_HEAD
                                [-h|--help]

    Options:
        -h --help                   show this
        -i, --infile BLOBDB         BlobDB file
        -c, --cov COV               COV file to be used in y-axis

        --xlabel XLABEL             Label for x-axis [default: BlobDB_cov]
        --ylabel YLABEL             Label for y-axis [default: CovFile_cov]
        --max FLOAT                 Maximum values for x/y-axis [default: 1e10]

        --lib COVLIB                Plot only certain covlib(s). Separated by ","
        --notitle                   Do not add filename as title to plot
        -p, --plotgroups INT        Number of (taxonomic) groups to plot, remaining
                                     groups are placed in 'other' [default: 7]
        -l, --length INT            Minimum sequence length considered for plotting [default: 100]
        --cindex                    Colour blobs by 'c index' [default: False]
        -n, --nohit                 Hide sequences without taxonomic annotation [default: False]
        -s, --noscale               Do not scale sequences by length [default: False]
        --legend                    Plot legend of blobplot in separate figure
        -m, --multiplot             Multi-plot. Print blobplot for each (taxonomic) group separately
        --cumulative                Print plot after addition of each (taxonomic) group
        --sort <ORDER>              Sort order for plotting [default: span]
                                     span  : plot with decreasing span
                                     count : plot with decreasing count
        --hist <HIST>               Data for histograms [default: span]
                                     span  : span-weighted histograms
                                     count : count histograms
        -r, --rank <RANK>           Taxonomic rank used for colouring of blobs [default: phylum]
                                     (Supported: species, genus, family, order, phylum, superkingdom)
        -x, --taxrule <TAXRULE>     Taxrule which has been used for computing taxonomy
                                     (Supported: bestsum, bestsumorder) [default: bestsum]
        --format FORMAT             Figure format for plot (png, pdf, eps, jpeg,
                                        ps, svg, svgz, tiff) [default: png]
        --noblobs                   Omit blobplot [default: False]
        --noreads                   Omit plot of reads mapping [default: False]
        -o, --out PREFIX            Output prefix
        --label GROUPS...           Relabel (taxonomic) groups, can be used several times.
                                     e.g. "A=Actinobacteria,Proteobacteria"
        --colours COLOURFILE        File containing colours for (taxonomic) groups
        --exclude GROUPS..          Place these (taxonomic) groups in 'other',
                                     e.g. "Actinobacteria,Proteobacteria"
        --refcov <FILE>               File containing number of "total" and "mapped" reads
                                     per coverage file. (e.g.: bam0,900,100). If provided, info
                                     will be used in read coverage plot(s).
        --catcolour <FILE>            Colour plot based on categories from FILE
                                     (format : "seq\tcategory").
"""

from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import blobtools
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtCore as Bt
import lib.BtPlot as BtPlot

if __name__ == '__main__':
    main_dir = dirname(__file__)

    args = docopt(__doc__)
    args = BtPlot.check_input(args)
    blobdb_f = args['--infile']
<<<<<<< HEAD
    cov_f = args['--cov']
=======
>>>>>>> FETCH_HEAD
    rank = args['--rank']
    c_index = args['--cindex']
    min_length = int(args['--length'])
    max_group_plot = int(args['--plotgroups'])
    hide_nohits = args['--nohit']
    taxrule = args['--taxrule']
    c_index = args['--cindex']
    exclude_groups = args['--exclude']
    labels = args['--label']
    colour_f = args['--colours']
    refcov_f = args['--refcov']
    catcolour_f = args['--catcolour']

    multiplot = args['--multiplot']
    out_prefix = args['--out']
    sort_order = args['--sort']
    hist_type = args['--hist']
    no_title = args['--notitle']
    ignore_contig_length = args['--noscale']
<<<<<<< HEAD
    format = args['--format']
    no_plot_blobs = args['--noblobs']
    no_plot_reads = args['--noreads']
=======
    labels = args['--label']
    colour_f = args['--colours']
    exclude_groups = args['--exclude']
    format = args['--format']
    no_plot_blobs = args['--noblobs']
    no_plot_reads = args['--noreads']
    refcov_f = args['--refcov']
    catcolour_f = args['--catcolour']
>>>>>>> FETCH_HEAD
    legend_flag = args['--legend']
    cumulative_flag = args['--cumulative']
    cov_lib_selection = args['--lib']

    x_label = args['--xlabel']
    y_label = args['--ylabel']
    axis_max = float(args['--max'])
<<<<<<< HEAD
=======

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
>>>>>>> FETCH_HEAD

    exclude_groups = BtIO.parseCmdlist(exclude_groups)
    refcov_dict = BtIO.parseReferenceCov(refcov_f)
    user_labels = BtIO.parseCmdLabels(labels)
    catcolour_dict = BtIO.parseCatColour(catcolour_f)
    colour_dict = BtIO.parseColours(colour_f)

    # Load BlobDb
    print BtLog.status_d['9'] % blobdb_f
<<<<<<< HEAD
    blobDb = Bt.BlobDb('blobplot')
    blobDb.version = blobtools.__version__
    blobDb.load(blobdb_f)

    # Generate plot data
    print BtLog.status_d['18']
    data_dict, max_cov, cov_lib_dict = blobDb.getPlotData(rank, min_length, hide_nohits, taxrule, False, False)
    plotObj = BtPlot.PlotObj(data_dict, cov_lib_dict, cov_lib_selection, 'covplot')
    plotObj.cov_y_dict, reads_total, reads_mapped, reads_unmapped, read_cov_dict = BtIO.parseCov(cov_f, set(blobDb.dict_of_blobs))
    plotObj.exclude_groups = exclude_groups
=======
    blobDB = Bt.BlobDb('new')
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

    data_dict, max_cov, cov_lib_dict = blobDB.getPlotData(rank, min_length, hide_nohits, taxrule, False, False)
    plotObj = BtPlot.PlotObj(data_dict, cov_lib_dict, cov_lib_selection)
    plotObj.cov_y_dict, reads_total, reads_mapped, reads_unmapped, read_cov_dict = BtIO.parseCov(covLib.f, set(blobDB.dict_of_blobs))
    plotObj.max_cov = max_cov
    plotObj.title = title
>>>>>>> FETCH_HEAD
    plotObj.format = format
    plotObj.max_cov = max_cov
    plotObj.no_title = no_title
    plotObj.multiplot = multiplot
    plotObj.hist_type = hist_type
    plotObj.ignore_contig_length = ignore_contig_length
    plotObj.max_group_plot = max_group_plot
    plotObj.legend_flag = legend_flag
    plotObj.cumulative_flag = cumulative_flag
    # order by which to plot (should know about user label)
    plotObj.group_order = BtPlot.getSortedGroups(data_dict, sort_order)
    # labels for each level of stats
    plotObj.labels.update(plotObj.group_order)
    # plotObj.group_labels is dict that contains labels for each group : all/other/user_label
    if (user_labels):
        for group, label in user_labels.items():
            plotObj.labels.add(label)
    plotObj.group_labels = {group : set() for group in plotObj.group_order}
    plotObj.relabel_and_colour(colour_dict, user_labels)
    plotObj.compute_stats()
    plotObj.refcov_dict = refcov_dict
    # Plotting
    info_flag = 1

    out_f = ''
    for cov_lib in plotObj.cov_libs:
        out_f = "%s.%s.%s.p%s.%s.%s" % (blobDb.title, taxrule, rank, max_group_plot, hist_type, min_length)
        if catcolour_dict:
            out_f = "%s.%s" % (out_f, "catcolour")
        if ignore_contig_length:
            out_f = "%s.%s" % (out_f, "noscale")
        if c_index:
            out_f = "%s.%s" % (out_f, "c_index")
        if exclude_groups:
            out_f = "%s.%s" % (out_f, "exclude" + "_".join(exclude_groups))
        if labels:
            out_f = "%s.%s" % (out_f, "userlabel_" + "_".join(set([name for name in user_labels.values()])))
        out_f = "%s.%s" % (out_f, "covplot")
        out_f = BtIO.getOutFile(out_f, out_prefix, None)
        if not (no_plot_blobs):
            plotObj.plotScatter(cov_lib, info_flag, out_f)
            info_flag = 0
    #plotObj.write_stats(out_f)
