#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: blobtools covplot  -i BLOBDB -c COV [--max FLOAT]
                                [--xlabel XLABEL] [--ylabel YLABEL]
                                [--lib COVLIB] [-o PREFIX] [-m]
                                [-p INT] [-l INT] [--cindex] [-n] [-s]
                                [-r RANK] [-x TAXRULE] [--label GROUPS...]
                                [--sort ORDER] [--sort_first LABELS]
                                [--hist HIST] [--notitle]
                                [--colours FILE] [--exclude FILE]
                                [--refcov FILE] [--catcolour FILE]
                                [--format FORMAT] [--noblobs] [--noreads] [--legend]
                                [--cumulative]
                                [-h|--help]

    Options:
        -h --help                   show this
        -i, --infile BLOBDB         BlobDB file
        -c, --cov COV               COV file to be used in y-axis

        --xlabel XLABEL             Label for x-axis
        --ylabel YLABEL             Label for y-axis
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
        --sort_first <L1,L2,...>    Labels that should always be plotted first, regardless of sort order
                                     ("no-hit,other,undef" is often a useful setting)
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
                                     (format : "seq,category").
"""

from __future__ import division
from docopt import docopt

from os.path import basename

import lib.interface as interface
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtCore as Bt
import lib.BtPlot as BtPlot

def main():
    args = docopt(__doc__)
    args = BtPlot.check_input(args)
    blobdb_f = args['--infile']
    cov_f = args['--cov']
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
    sort_order = args['--sort']
    sort_first = args['--sort_first']
    multiplot = args['--multiplot']
    out_prefix = args['--out']
    sort_order = args['--sort']
    hist_type = args['--hist']
    no_title = args['--notitle']
    ignore_contig_length = args['--noscale']
    format_plot = args['--format']
    no_plot_blobs = args['--noblobs']
    #no_plot_reads = args['--noreads']
    legend_flag = args['--legend']
    cumulative_flag = args['--cumulative']
    cov_lib_selection = args['--lib']

    xlabel = args['--xlabel']
    ylabel = args['--ylabel']
    axis_max = float(args['--max'])

    exclude_groups = BtIO.parseCmdlist(exclude_groups)
    refcov_dict = BtIO.parseReferenceCov(refcov_f)
    user_labels = BtIO.parseCmdLabels(labels)
    catcolour_dict = BtIO.parseCatColour(catcolour_f)
    colour_dict = BtIO.parseColours(colour_f)

    # Load BlobDb
    print(BtLog.status_d['9'] % blobdb_f)
    blobDb = Bt.BlobDb('blobplot')
    blobDb.version = interface.__version__
    blobDb.load(blobdb_f)

    # Generate plot data
    print(BtLog.status_d['1'] % ('cov_y_axis', cov_f))
    cov_y_dict, reads_total, reads_mapped, reads_unmapped, read_cov_dict = BtIO.parseCov(cov_f, set(blobDb.dict_of_blobs))
    print(BtLog.status_d['18'])
    data_dict, min_cov, max_cov, cov_lib_dict = blobDb.getPlotData(rank, min_length, hide_nohits, taxrule, c_index, catcolour_dict)
    plotObj = BtPlot.PlotObj(data_dict, cov_lib_dict, cov_lib_selection, 'covplot', sort_first)
    # set lowest coverage to 0.01
    for contig in cov_y_dict:
        if cov_y_dict[contig] < 0.1:
            cov_y_dict[contig] = 0.1
    plotObj.cov_y_dict = cov_y_dict
    plotObj.exclude_groups = exclude_groups
    plotObj.version = blobDb.version
    plotObj.format = format_plot
    plotObj.max_cov = axis_max
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
        plotObj.xlabel = basename(cov_lib_dict[cov_lib]['f'])
        plotObj.ylabel = cov_f
        if (ylabel):
            plotObj.ylabel = ylabel
        if (xlabel):
            plotObj.xlabel = xlabel
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
        out_f = "%s.%s" % (out_f, "covplot")
        if (plotObj.cumulative_flag):
            out_f = "%s.%s" % (out_f, "cumulative")
        if (plotObj.multiplot):
            out_f = "%s.%s" % (out_f, "multiplot")
        out_f = BtIO.getOutFile(out_f, out_prefix, None)
        if not (no_plot_blobs):
            plotObj.plotScatter(cov_lib, info_flag, out_f)
            info_flag = 0
    plotObj.write_stats(out_f)

if __name__ == '__main__':
    main()
