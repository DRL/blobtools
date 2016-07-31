#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools blobplot  -i BLOBDB
                                [-p INT] [-l INT] [--cindex] [-n] [-s]
                                [-r RANK] [-x TAXRULE] [--label GROUPS...]
                                [--lib COVLIB] [-o PREFIX] [-m]
                                [--sort ORDER] [--hist HIST] [--notitle] [--filelabel]
                                [--colours FILE] [--exclude FILE]
                                [--refcov FILE] [--catcolour FILE]
                                [--format FORMAT] [--noblobs] [--noreads] [--legend]
                                [--cumulative]
                                [-h|--help]

    Options:
        -h --help                   show this
        -i, --infile BLOBDB         BlobDB file (created with "blobtools create")
        --lib COVLIB                Plot only certain covlib(s). Separated by ","
        --notitle                   Do not add filename as title to plot
        --filelabel                 Label axis based on filenames
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
                                     (Supported: species, genus, family, order,
                                        phylum, superkingdom)
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
        --exclude GROUPS            Exclude these (taxonomic) groups (also works for 'other')
                                     e.g. "Actinobacteria,Proteobacteria,other"
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
import bloblib.BtLog as BtLog
import bloblib.BtIO as BtIO
import bloblib.BtCore as BtCore
import bloblib.BtPlot as BtPlot

def main():
    args = docopt(__doc__)
    args = BtPlot.check_input(args)

    blobdb_f = args['--infile']
    rank = args['--rank']
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
    format = args['--format']
    no_plot_blobs = args['--noblobs']
    no_plot_reads = args['--noreads']
    legend_flag = args['--legend']
    cumulative_flag = args['--cumulative']
    cov_lib_selection = args['--lib']

    filelabel = args['--filelabel']

    exclude_groups = BtIO.parseCmdlist(exclude_groups)
    refcov_dict = BtIO.parseReferenceCov(refcov_f)
    user_labels = BtIO.parseCmdLabels(labels)
    catcolour_dict = BtIO.parseCatColour(catcolour_f)
    colour_dict = BtIO.parseColours(colour_f)

    # Load BlobDb
    print BtLog.status_d['9'] % blobdb_f
    blobDb = BtCore.BlobDb('blobplot')
    blobDb.version = blobtools.__version__
    blobDb.load(blobdb_f)

    # Generate plot data
    print BtLog.status_d['18']
    data_dict, max_cov, cov_lib_dict = blobDb.getPlotData(rank, min_length, hide_nohits, taxrule, c_index, catcolour_dict)
    plotObj = BtPlot.PlotObj(data_dict, cov_lib_dict, cov_lib_selection, 'blobplot')
    plotObj.exclude_groups = exclude_groups
    plotObj.version = blobDb.version
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
        plotObj.ylabel = "Coverage"
        plotObj.xlabel = "GC proportion"
        if (filelabel):
            plotObj.ylabel = basename(cov_lib_dict[cov_lib]['f'])
        out_f = "%s.%s.%s.p%s.%s.%s" % (blobDb.title, taxrule, rank, max_group_plot, hist_type, min_length)
        if catcolour_dict:
            out_f = "%s.%s" % (out_f, "catcolour")
        if ignore_contig_length:
            out_f = "%s.%s" % (out_f, "noscale")
        if c_index:
            out_f = "%s.%s" % (out_f, "c_index")
        if exclude_groups:
            out_f = "%s.%s" % (out_f, "exclude_" + "_".join(exclude_groups))
        if labels:
            out_f = "%s.%s" % (out_f, "userlabel_" + "_".join(set([name for name in user_labels.values()])))
        out_f = "%s.%s" % (out_f, "blobplot")
        if (plotObj.cumulative_flag):
            out_f = "%s.%s" % (out_f, "cumulative")
        if (plotObj.multiplot):
            out_f = "%s.%s" % (out_f, "multiplot")
        out_f = BtIO.getOutFile(out_f, out_prefix, None)
        if not (no_plot_blobs):
            plotObj.plotScatter(cov_lib, info_flag, out_f)
            info_flag = 0
        if not (no_plot_reads):
            plotObj.plotBar(cov_lib, out_f)
    plotObj.write_stats(out_f)

if __name__ == '__main__':
    main()
