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
from numpy import array, arange, logspace, mean, std
import math
import lib.BtLog as BtLog
import lib.BtIO as BtIO
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
BLACK, GREY, BGGREY, WHITE = unicode('#262626'), unicode('#d3d3d3'), unicode('#F0F0F5'), unicode('#ffffff')
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

def getSortedGroups(filter_dict, sort_order):
    """ Returns list of sorted groups based on span or count. """
    sorted_groups = []
    if sort_order == 'span':
        sorted_groups = sorted(filter_dict, key = lambda x : filter_dict[x]['span_visible'] if filter_dict[x]['span_visible'] > 0 else None, reverse=True)
    elif sort_order == 'count':
        sorted_groups = sorted(filter_dict, key = lambda x : filter_dict[x]['count_visible'] if filter_dict[x]['count_visible'] > 0 else None, reverse=True)
    else:
        pass
    return sorted_groups

def generateColourDict(colour_groups):
    cmap = cm.get_cmap(name=COLOURMAP)
    n_tax = 0
    tax_with_colour = []
    if 'no-hit' in colour_groups or 'None' in colour_groups:
        n_tax = len(colour_groups) - 1
    else:
        n_tax = len(colour_groups)
    breaks = [0.0 + x*(1.0-0.0)/n_tax for x in range(n_tax)]
    colour_d = {group: rgb2hex(cmap(b)) for b, group in izip(breaks, colour_groups)}
    if 'no-hit' in colour_groups:
        colour_d['no-hit'] = GREY
    if 'None' in colour_groups: 
        colour_d['None'] = GREY
    colour_d['other'] = WHITE
    return colour_d

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

def parse_labels(labels):
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

class PlotObj():
    def __init__(self, data_dict, filter_dict, cov_libs):
        self.labels = {}
        self.cov_libs = cov_libs
        self.data_dict = data_dict
        self.read_cov = {}
        self.exclude_groups = []
        self.count = {group : filter_dict[group]['count'] for group in filter_dict.keys()}
        self.count_visible = {group : filter_dict[group]['count_visible'] for group in filter_dict.keys()}
        self.count_hidden = {group : filter_dict[group]['count_hidden'] for group in filter_dict.keys()}
        self.span = {group : filter_dict[group]['span'] for group in filter_dict.keys()}
        self.span_visible = {group : filter_dict[group]['span_visible'] for group in filter_dict.keys()}
        self.span_hidden = {group : filter_dict[group]['span_hidden'] for group in filter_dict.keys()}
        self.colours = {}
        self.group_order = []
        self.plot_order = []
        self.min_cov = 0.01
        self.max_cov = 0.0
        self.out_f = ''
        self.title = ''
        self.n50 = {}
        self.gc_mean = {}
        self.gc_std = {}
        self.cov_mean = {}
        self.cov_std = {}
        self.max_group_plot = 0

    def stats(self):
        stats_data = []
        for idx, group in enumerate(self.group_order):
            table_data.append({
                        'name' : group, 
                        'label' : self.labels[group],
                        'count_total' : "{:,}".format(self.count[group]), 
                        'count_visible_perc' : '{0:.1%}'.format(self.count_visible[group]/self.count[group]), 
                        'span_total' : "{:,}".format(self.span[group]), 
                        'span_visible_perc' : '{0:.1%}'.format(self.span_visible[group]/self.span[group]),
                        'colour' : self.colours[group],
                        'n50' : self.n50[group],
                        'gc_mean' : self.gc_mean[group],
                        'gc_std' : self.gc_std[group],
                        'cov_mean' : {cov_lib : cov_mean for cov_lib, cov_mean in self.cov_mean[group].items()},
                        'cov_std' : {cov_lib : cov_std for cov_lib, cov_std in self.cov_std[group].items()},
                        'read_cov' : {cov_lib : read_cov for cov_lib, cov_std in self.read_cov[group].items() if not read_cov == 0 else pass}
                        })
        
        out_f = "%s.plot.txt" % self.out_f
        with open(out_f,'w') as fh:
            fh.write("{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:<10}\t{:<10}\t{:<10}\n".format("COUNT", "visible (%)", "SPAN (nt)", "visible (%)", "group", "label", "colour"))
            for d in table_data:
                fh.write("{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:<10}\t{:<10}\t{:<10}\n".format(d['count_total'], d['count_visible_perc'], d['span_total'], d['span_visible_perc'], d['name'], d['label'], d['colour'] ))

    def compute_stats(self):
        other_gc = []
        other_length = []
        other_covs = {cov_lib : [] for cov_lib in self.cov_libs}
        for group, label in self.labels.items():
            group_gc = self.data_dict[group]['gc']
            self.gc_mean[group] = mean(array(group_gc))
            self.gc_std[group] = std(array(group_gc))
            group_length = self.data_dict[group]['length']
            self.n50[group] = n50(group_length)
            for cov_lib, covs in self.data_dict[group]['covs'].items():
                if not group in self.cov_mean:
                    self.cov_mean[group] = {}
                self.cov_mean[group][cov_lib] = mean(array(covs))
                if not group in self.cov_std:
                    self.cov_std[group] = {}
                self.cov_std[group][cov_lib] = std(array(covs))
            if label == 'other':
                other_gc += group_gc
                other_length += group_length 
                for cov_lib, covs in self.data_dict[group]['covs'].items():
                    other_covs[cov_lib] += covs
                self.count['other'] = self.count.get('other', 0) + self.count[group]
                self.count_visible['other'] = self.count_visible.get('other', 0) + self.count_visible[group]
                self.count_hidden['other'] = self.count_hidden.get('other', 0) + self.count_hidden[group]
                self.span['other'] = self.span.get('other', 0) + self.span[group]
                self.span_visible['other'] = self.span_visible.get('other', 0) + self.span_visible[group]
                self.span_hidden['other'] = self.span_hidden.get('other', 0) + self.span_hidden[group]
                self.data_dict['other'] = {'length' : [], 'gc' : [], 'covs' : {covLib : [] for covLib in self.cov_libs} }
                self.data_dict['other']['length'] = other_length
                self.data_dict['other']['gc'] = other_gc
                self.gc_mean['other'] = mean(array(other_gc))
                self.gc_std['other'] = std(array(other_gc))
                self.n50['other'] = n50(other_length)
                for cov_lib, covs in other_covs.items():
                    self.data_dict['other']['covs'][cov_lib] = covs
                    if not 'other' in self.cov_mean:
                        self.cov_mean['other'] = {}
                    self.cov_mean['other'][cov_lib] = mean(array(covs))
                    if not 'other' in self.cov_std:
                        self.cov_std['other'] = {}
                    self.cov_std['other'][cov_lib] = std(array(covs))
    
    def relabel_and_colour(self, colour_f, label_d):
        if (colour_f):
            colour_dict = BtIO.parseColourDict(colour_f)
        else:
            colour_groups = self.group_order[0:self.max_group_plot]
            colour_dict = generateColourDict(colour_groups)
            
        for idx, group in enumerate(self.group_order):
            if (label_d):
                if group in label_d:
                    self.labels[group] = label_d[group] 
                    self.colours[group] = colour_dict[group]
            elif (self.exclude_groups):
                if group in self.exclude_groups:
                    self.labels[group] = 'other'
                    self.colours[group] = colour_dict['other']     
            elif group in colour_dict:
                self.labels[group] = group
                self.colours[group] = colour_dict[group] 
                self.plot_order.append(group)
            elif idx > self.max_group_plot:
                self.labels[group] = 'other'
                self.colours[group] = colour_dict['other'] 
            else:
                self.labels[group] = 'other'
                self.colours[group] = colour_dict['other'] 
        self.colours['other'] = colour_dict['other'] 
        self.plot_order.append('other')

    def writePlotSummaryTable(self):
        table_data = []
        for group in self.group_order:
            table_data.append({
                        'name' : group, 
                        'label' : self.labels[group],
                        'count_total' : "{:,}".format(self.count[group]), 
                        'count_visible_perc' : '{0:.1%}'.format(self.count_visible[group]/self.count[group]), 
                        'span_total' : "{:,}".format(self.span[group]), 
                        'span_visible_perc' : '{0:.1%}'.format(self.span_visible[group]/self.span[group]),
                        'colour' : self.colours[group]
                        })
        
        out_f = "%s.plot.txt" % self.out_f
        with open(out_f,'w') as fh:
            fh.write("{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:<10}\t{:<10}\t{:<10}\n".format("COUNT", "visible (%)", "SPAN (nt)", "visible (%)", "group", "label", "colour"))
            for d in table_data:
                fh.write("{:>10}\t{:>10}\t{:>10}\t{:>10}\t{:<10}\t{:<10}\t{:<10}\n".format(d['count_total'], d['count_visible_perc'], d['span_total'], d['span_visible_perc'], d['name'], d['label'], d['colour'] ))

    def plotBlobs(self, cov_lib, info_flag):
        rect_scatter, rect_histx, rect_histy, rect_legend = set_canvas()
        # Setting up plots and axes
        plt.figure(1, figsize=(35,35), dpi=400)
        axScatter = plt.axes(rect_scatter, axisbg=BGGREY, yscale = 'log')
        axScatter = set_format_scatterplot(axScatter, self.max_cov)
        axHistx = plt.axes(rect_histx, axisbg=BGGREY)
        axHistx = set_format_hist_x(axHistx, axScatter)
        axHisty = plt.axes(rect_histy, axisbg=BGGREY)
        axHisty = set_format_hist_y(axHisty, axScatter)
        if self.hist_type == "span":
            axHistx.set_ylabel("Span (kb)")
            axHisty.set_xlabel("Span (kb)", rotation='horizontal')
        else:
            axHistx.set_ylabel("Count")
            axHisty.set_xlabel("Count", rotation='horizontal')    
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
        if (self.title):
            plt.suptitle(self.title, fontsize=35, verticalalignment='top')
        # Setting bins for histograms
        top_bins = arange(0, 1, 0.01)
        right_bins = logspace(-2, (int(math.log(self.max_cov)) + 1), 200, base=10.0)
        # empty handles for big legend
        legend_handles = []
        legend_labels = []
        # counter necessary for multiplot so that PNGs are in order when sorted by name
        i = 0
        for group in self.plot_order:
            i += 1
            group_length_array = array(self.data_dict[group]['length'])
            group_gc_array = array(self.data_dict[group]['gc'])
            group_cov_array = array(self.data_dict[group]['covs'][cov_lib])
            # calculate values for legend
            group_span_in_mb = round(self.span_visible[group]/1000000, 2)
            group_number_of_seqs = self.count_visible[group]
            group_n50 = self.n50[group]
            blob_size_array = []
            s, lw, alpha, colour = 15, 0.5, 1, self.colours[group]
            if (self.ignore_contig_length):
                if not group == "no-hit":
                    s = 65
                blob_size_array = [s for length in group_length_array]
            else:
                blob_size_array = [length/s for length in group_length_array]
            if group == "no-hit":
                alpha = 0.5
            weights_array = group_length_array/1000
            # generate label for legend
            fmt_seqs = "{:,}".format(group_number_of_seqs)
            fmt_span = "{:,}".format(group_span_in_mb)
            fmt_n50 = "{:,}".format(group_n50)
            label = "%s (%s;%sMB;%snt)" % (group, fmt_seqs, fmt_span, fmt_n50)
            if (info_flag):
                print BtLog.info_d['0'] % (group, fmt_seqs, fmt_span, fmt_n50)
            legend_handles.append(Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=24, markerfacecolor=colour))
            legend_labels.append(label)
            if (self.hist_type == "span"):
                axHistx.hist(group_gc_array, weights=weights_array, color = colour, bins = top_bins, histtype='step', lw = 3)
                axHisty.hist(group_cov_array, weights=weights_array, color = colour, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
            else:           
                axHistx.hist(group_gc_array, color = colour, bins = top_bins, histtype='step', lw = 3)
                axHisty.hist(group_cov_array , color = colour, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
            axScatter.scatter(group_gc_array, group_cov_array, color = colour, s = blob_size_array, lw = lw, alpha=alpha, edgecolor=BLACK, label=label)
            axLegend.axis('off')
            if (self.multiplot): 
                axLegend.legend(legend_handles, legend_labels, loc=6, numpoints=1, fontsize=FONTSIZE, frameon=True)
                plot_ref_legend(axScatter)
                m_out_f = "%s.%s.%s.%s" % (self.out_f, i, group, self.format)
                print BtLog.status_d['8'] % m_out_f
                plt.savefig(m_out_f, format=self.format)
        if not (self.ignore_contig_length):
            plot_ref_legend(axScatter)
        axLegend.legend(legend_handles, legend_labels, numpoints=1, fontsize=FONTSIZE, frameon=True, loc=6 )
        #table_data = getSummaryTable(summary_dict, plot_order)
        # writing data table for plot
        #writeTableData(table_data, out_f)
        self.out_f = "%s.%s" % (self.out_f, self.format)
        print BtLog.status_d['8'] % self.out_f
        plt.savefig(self.out_f, format=self.format) 
        plt.close()  
        
        

    

if __name__ == "__main__": 
    pass