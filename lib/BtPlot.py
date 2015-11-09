#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File        : BtPlot.py
Version     : 0.1
Author      : Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs        : ?
To do       : ?
"""
from __future__ import division
from numpy import amax, amin, arange, logspace, where
import numpy as np
import math
import lib.BtLog as BtLog
import matplotlib as mat
from matplotlib import cm
from matplotlib.ticker import NullFormatter
from matplotlib.lines import Line2D
from matplotlib.colors import rgb2hex
mat.use('agg')
import matplotlib.pyplot as plt
from itertools import izip

mat.rcParams.update({'font.size': 30})
mat.rcParams['xtick.major.pad'] = '8'
mat.rcParams['ytick.major.pad'] = '8'
mat.rcParams['lines.antialiased'] = True

FONTSIZE = 24
COLOURMAP = "Set2" # "Set1"
BLACK, GREY, BGGREY, WHITE = '#262626', '#d3d3d3', '#F0F0F5', '#ffffff'
FIGFORMAT = 'png'
nullfmt = NullFormatter()

def n50(list_of_lengths):
    total_span = 0
    sorted_list_of_lengths=sorted(list_of_lengths, reverse=True)
    for contig_length in sorted_list_of_lengths:
        total_span += contig_length
    teoN50 = total_span/2.0
    running_sum = 0
    N50 = 0
    for contig_length in sorted_list_of_lengths:
        running_sum += contig_length
        if teoN50 <= running_sum:
            N50 = contig_length
            break
    return N50

def getPlotOrder(summary_dict, sort_order, max_taxa_plot):
    """ Returns list of tax of size max_taxa_plot sorted by span or count. """
    plot_order = []
    if sort_order == 'span':
        plot_order = sorted(summary_dict, key = lambda x : summary_dict[x]['span_visible'] if summary_dict[x]['span_visible'] > 0 else None, reverse=True)
    elif sort_order == 'count':
        plot_order = sorted(summary_dict, key = lambda x : summary_dict[x]['count_visible'] if summary_dict[x]['count_visible'] > 0 else None, reverse=True)
    else:
        pass
    if len(plot_order) > max_taxa_plot:
        plot_order = plot_order[0:max_taxa_plot]
    return plot_order

def getColourDict(plot_order):
    """ Returns colour dict, Not annotated is always grey. """
    cmap = cm.get_cmap(name=COLOURMAP)
    n_tax = 0
    tax_with_colour = []
    if 'no-hit' in plot_order or 'None' in plot_order:
        n_tax = len(plot_order) - 1
        tax_with_colour = [tax for tax in plot_order if not tax == "no-hit" and not tax == 'None']
    else:
        n_tax = len(plot_order)
        tax_with_colour = [tax for tax in plot_order]
    breaks = [0.0 + x*(1.0-0.0)/n_tax for x in range(n_tax)]
    colour_d = {tax: rgb2hex(cmap(b)) for b, tax in izip(breaks, tax_with_colour)}
    if 'no-hit' in plot_order:
        colour_d['no-hit'] = GREY
    elif 'None' in plot_order:
        colour_d['None'] = GREY
    else:
        pass    
    return colour_d

def getMinMaxCov(cov_arrays):
    max_cov, min_cov = 100.00, 100.00
    for cov_lib in cov_arrays:
        lib_max_cov, lib_min_cov = amax(cov_arrays[cov_lib].astype(float)), amin(cov_arrays[cov_lib].astype(float))
        if lib_max_cov > max_cov:
            max_cov = lib_max_cov
        if lib_min_cov < min_cov:
            min_cov = lib_min_cov
    return min_cov, max_cov

def getSumCov(cov_arrays):
    arrays = [cov_array for cov_array in cov_arrays.values()]
    sum_array = np.sum(arrays, axis=0)
    cov_arrays['sum'] = sum_array
    return cov_arrays

def getSummaryTable(summary_dict, plot_order):
    table_data = []
    for elem in sorted(summary_dict, key = lambda x : summary_dict[x]['span_visible'], reverse=True):
        count_total = 0
        count_visible_perc = ''
        span_total = 0
        span_visible_perc = ''
        if elem in plot_order:
            count_total = summary_dict[elem]['count_total']
            count_visible_perc = '{0:.1%}'.format(summary_dict[elem]['count_visible']/count_total)
            span_total = summary_dict[elem]['span_total']
            span_visible_perc = '{0:.1%}'.format(summary_dict[elem]['span_visible']/span_total)
        else:
            count_total = summary_dict[elem]['count_total']
            count_visible_perc = '{0:.1%}'.format(0.0)
            span_total = summary_dict[elem]['span_total']
            span_visible_perc = '{0:.1%}'.format(0.0)
        table_data.append({
                    'name' : elem, 
                    'count_total' : "{:,}".format(count_total), 
                    'count_visible_perc' : count_visible_perc, 
                    'span_total' : "{:,}".format(span_total), 
                    'span_visible_perc' : span_visible_perc})
    return table_data

def writeTableData(table_data, out_f):
    out_f = "%s.txt" % out_f
    with open(out_f,'w') as fh:
        fh.write("{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:<10}\n".format("COUNT", "visible (%)", "SPAN (nt)", "visible (%)", "group"))
        for d in table_data:
            fh.write("{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:<10}\n".format(d['count_total'], d['count_visible_perc'], d['span_total'], d['span_visible_perc'], d['name'] ))

def set_canvas():
    left, width = 0.1, 0.60
    bottom, height = 0.1, 0.60
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    rect_legend = [left_h, bottom_h, 0.2, 0.2]
    return rect_scatter, rect_histx, rect_histy, rect_legend

def set_format_scatterplot(axScatter, max_cov):
    axScatter.set_xlabel("GC proportion", fontsize=35)
    axScatter.set_ylabel("Coverage", fontsize=35)
    axScatter.grid(True, which="major", lw=2., color=WHITE, linestyle='-') 
    axScatter.set_axisbelow(True)
    axScatter.set_xlim( (0, 1) )
    axScatter.set_ylim( (0.01, max_cov+1000) ) # This sets the max-Coverage so that all libraries + sum are at the same scale
    axScatter.xaxis.labelpad = 20
    axScatter.xaxis.labelpad = 20
    return axScatter

def set_format_hist_x(axHistx, axScatter):
    axHistx.set_xlim( axScatter.get_xlim() )
    axHistx.grid(True, which="major", lw=2., color= WHITE, linestyle='-')
    axHistx.xaxis.set_major_formatter(nullfmt) # no labels since redundant
    axHistx.set_axisbelow(True)
    axHistx.yaxis.labelpad = 20
    return axHistx

def set_format_hist_y(axHisty, axScatter):
    axHisty.set_yscale('log')
    axHisty.yaxis.set_major_formatter(nullfmt) # no labels since redundant
    axHisty.set_ylim( axScatter.get_ylim() )
    axHisty.grid(True, which="major", lw=2., color= WHITE, linestyle='-')
    axHisty.set_axisbelow(True)
    axHisty.xaxis.labelpad = 20
    return axHisty

def plot_ref_legend(axScatter):
    s = 15
    # markersize in scatter is in "points^2", markersize in Line2D is in "points" ... that's why we need math.sqrt()
    ref_1 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(1000/15),  markerfacecolor=GREY))
    ref_2 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(5000/15), markerfacecolor=GREY))
    ref_3 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(10000/15), markerfacecolor=GREY))
    axScatter.legend([ref_1,ref_2,ref_3], ["1,000nt", "5,000nt", "10,000nt"], numpoints=1, loc = 4, fontsize=FONTSIZE)

def plot(data_array, cov_array, summary_dict, plot_order, colour_dict, min_cov, max_cov, multiplot, hist_type, plot_title, out_f, ignore_contig_length, info_flag):
    rect_scatter, rect_histx, rect_histy, rect_legend = set_canvas()
    # Setting up plots and axes
    plt.figure(1, figsize=(35,35), dpi=400)
    axScatter = plt.axes(rect_scatter, axisbg=BGGREY, yscale = 'log')
    axScatter = set_format_scatterplot(axScatter, max_cov)
    axHistx = plt.axes(rect_histx, axisbg=BGGREY)
    axHistx = set_format_hist_x(axHistx, axScatter)
    axHisty = plt.axes(rect_histy, axisbg=BGGREY)
    axHisty = set_format_hist_y(axHisty, axScatter)
    if hist_type == "span":
        axHistx.set_ylabel("Span (kb)")
        axHisty.set_xlabel("Span (kb)", rotation='horizontal')
    elif hist_type == "count":
        axHistx.set_ylabel("Count")
        axHisty.set_xlabel("Count", rotation='horizontal')    
    else: 
        pass
    axScatter.yaxis.get_major_ticks()[0].label1.set_visible(False)
    axScatter.yaxis.get_major_ticks()[1].label1.set_visible(False)
    for xtick in axHisty.get_xticklabels(): # rotate text for ticks in cov histogram 
        xtick.set_rotation(270)
    axLegend = plt.axes(rect_legend, axisbg=WHITE)
    axLegend.xaxis.set_major_locator(plt.NullLocator())
    axLegend.xaxis.set_major_formatter(nullfmt)
    axLegend.yaxis.set_major_locator(plt.NullLocator())
    axLegend.yaxis.set_major_formatter(nullfmt)
    # Setting title
    if (plot_title):
        plt.suptitle(plot_title, fontsize=35, verticalalignment='top')
    # Setting bins for histograms
    top_bins = arange(0, 1, 0.01)
    right_bins = logspace(-2, (int(math.log(max_cov)) + 1), 200, base=10.0)
    # empty handles for big legend
    legend_handles = []
    legend_labels = []
    # counter necessary for multiplot so that PNGs are in order when sorted by name
    i = 0
    # initiate variables for plotting
    s, lw, alpha, color = 0, 0, 0, ''
    for elem in plot_order:
        if (summary_dict[elem]['count_visible']):
            i += 1
            # get indices for those rows in data where the phylum == tax
            idx = where(data_array[:,3].astype(str) == elem, True, False)
            # create np_arrays for length, gc and cov for all contigs in phylum 
            elem_length_array = data_array[idx][:,1].astype(int)
            elem_gc_array = data_array[idx][:,2].astype(float)
            elem_cov_array = cov_array[idx].astype(float)
            # calculate values for legend
            elem_span_in_mb = round(summary_dict[elem]['span_visible']/1000000, 2)
            elem_number_of_seqs = summary_dict[elem]['count_visible']
            elem_n50 = n50(elem_length_array)
            # set variables for plotting
            blob_size_array = []
            s, lw, alpha, colour = 15, 0.5, 1, colour_dict[elem]
            if (ignore_contig_length):
                if not elem == "no-hit":
                    s = 65
                blob_size_array = [s for length in elem_length_array]
            else:
                blob_size_array = [length/s for length in elem_length_array]
            if elem == "no-hit":
                alpha = 0.5
            weights_array = elem_length_array/1000
            # generate label for legend
            fmt_seqs = "{:,}".format(elem_number_of_seqs)
            fmt_span = "{:,}".format(elem_span_in_mb)
            fmt_n50 = "{:,}".format(elem_n50)
            label = "%s (%s;%sMB;%snt)" % (elem, fmt_seqs, fmt_span, fmt_n50)
            if (info_flag):
                print BtLog.info_d['0'] % (elem, fmt_seqs, fmt_span, fmt_n50)
            legend_handles.append(Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=24, markerfacecolor=colour))
            legend_labels.append(label)
            if (hist_type == "span"):
                axHistx.hist(elem_gc_array, weights=weights_array, color = colour, bins = top_bins, histtype='step', lw = 3)
                axHisty.hist(elem_cov_array, weights=weights_array, color = colour, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
            else:           
                axHistx.hist(elem_gc_array, color = colour, bins = top_bins, histtype='step', lw = 3)
                axHisty.hist(elem_cov_array , color = colour, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
            axScatter.scatter(elem_gc_array, elem_cov_array, color = colour, s = blob_size_array, lw = lw, alpha=alpha, edgecolor=BLACK, label=label)
            axLegend.axis('off')

            if (multiplot): 
                axLegend.legend(legend_handles, legend_labels, loc=6, numpoints=1, fontsize=FONTSIZE, frameon=True)
                plot_ref_legend(axScatter)
                m_out_f = "%s.%s.%s.png" % (out_f, i, elem)
                print BtLog.status_d['8'] % m_out_f
                plt.savefig(m_out_f, format="png")
    if not (ignore_contig_length):
        plot_ref_legend(axScatter)
    axLegend.legend(legend_handles, legend_labels, numpoints=1, fontsize=FONTSIZE, frameon=True, loc=6 )
    table_data = getSummaryTable(summary_dict, plot_order)
    # writing data table for plot
    writeTableData(table_data, out_f)
    out_f = "%s.png" % out_f
    print BtLog.status_d['8'] % out_f
    plt.savefig(out_f, format="png") 
    plt.close()  
    
    

if __name__ == "__main__": 
    pass