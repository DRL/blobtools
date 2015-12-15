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
COLOURMAP = "Spectral" # "Set1", "Paired", "Set2", "Spectral"
BLACK, GREY, BGGREY, WHITE, DGREY = unicode('#262626'), unicode('#d3d3d3'), unicode('#F0F0F5'), unicode('#ffffff'), unicode('#4d4d4d')
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

def parseRefCov(refcov_f):
    refcov_dict = {}
    with open(refcov_f) as fh:
        for l in fh:
            try:
                cov_lib, reads_total_ref, reads_mapped_ref = l.split(",")
                refcov_dict[cov_lib] = {
                                        'reads_total' : int(reads_total_ref), 
                                        'reads_mapped' : int(reads_mapped_ref)
                                       }
            except:
                BtLog.error('21', refcov_f)
    return refcov_dict

def parseCatColour(catcolour_f):
    catcolour_dict = {}
    with open(catcolour_f) as fh:
        for l in fh:
            try:
                seq_name, category = l.rstrip("\n").split(",")
                catcolour_dict[seq_name] = category
            except:
                BtLog.error('23', catcolour_f)
    return catcolour_dict

def parseCovFile(cov_f):
    cov_dict = {}
    with open(cov_f) as fh:
        for l in fh:
            try:
                seq_name, cov = l.rstrip("\n").split("\t")
                if float(cov) < 0.02:
                    cov_dict[seq_name] = 0.02
                else:
                    cov_dict[seq_name] = float(cov)

            except:
                BtLog.error('25', cov_f)
    return cov_dict


def getSortedGroups(data_dict, sort_order):
    """ Returns list of sorted groups based on span or count. """
    sorted_groups = []
    if sort_order == 'span':
        sorted_groups = sorted(data_dict, key = lambda x : data_dict[x]['span_visible'] if data_dict[x]['span_visible'] > 0 else None, reverse=True)
    elif sort_order == 'count':
        sorted_groups = sorted(data_dict, key = lambda x : data_dict[x]['count_visible'] if data_dict[x]['count_visible'] > 0 else None, reverse=True)
    else:
        pass
    return sorted_groups

def generateColourDict(groups):
    cmap = cm.get_cmap(name=COLOURMAP)
    colour_groups = [group for group in groups if not group == 'no-hit' or not group == 'None']
    n_tax = len(colour_groups)
    breaks = [0.0 + x*(1.0-0.0)/n_tax for x in range(n_tax)]
    colour_d = {group: rgb2hex(cmap(b)) for b, group in izip(breaks, colour_groups)}
    if 'no-hit' in groups:
        colour_d['no-hit'] = GREY
    if 'None' in groups: 
        colour_d['None'] = GREY
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

def set_format_covplot(axScatter, x_max, y_max, x_label, y_label, scale, cov_lib):
    axScatter.set_xlabel(x_label + " " + cov_lib, fontsize=35)
    axScatter.set_ylabel(y_label, fontsize=35)
    axScatter.grid(True, which="major", lw=2., color=WHITE, linestyle='-') 
    axScatter.set_axisbelow(True)
    if scale == 'log':
        axScatter.set_xlim( (0.01, x_max+100) )
        axScatter.set_ylim( (0.01, y_max+100) ) 
    else:
        axScatter.set_xlim( (-10, x_max + (x_max*0.1) ) )
        axScatter.set_ylim( (-10, y_max + (y_max*0.1) ) ) 
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

def set_format_hist_x_cov(axHistx, axScatter, scale):
    axHistx.set_xlim(axScatter.get_xlim())
    axHistx.set_xscale(scale)
    axHistx.grid(True, which="major", lw=2., color= WHITE, linestyle='-')
    axHistx.xaxis.set_major_formatter(nullfmt) # no labels since redundant
    axHistx.set_axisbelow(True)
    axHistx.yaxis.labelpad = 20
    return axHistx

def set_format_hist_y_cov(axHisty, axScatter, scale):
    axHisty.set_ylim(axScatter.get_ylim())
    axHisty.set_yscale(scale)
    axHisty.grid(True, which="major", lw=2., color= WHITE, linestyle='-')
    axHisty.yaxis.set_major_formatter(nullfmt) # no labels since redundant
    axHisty.set_axisbelow(True)
    axHisty.xaxis.labelpad = 20
    return axHisty

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
    label_d = {}
    name, groups = '', ''
    if (labels):
        try:
            for label in labels:
                name, groups = str(label).split("=")
                if "," in groups:
                    for group in groups.split(","):
                        label_d[group] = name
                else:
                    label_d[groups] = name
        except:
            BtLog.error('17', labels)
    return label_d

class PlotReadObj():
    def __init__(self):
        self.values = []
        self.labels = []
        self.colours = []

class PlotObj():
    def __init__(self, data_dict, cov_libs, cov_libs_total_reads_d):
        self.labels = {'all'}
        self.group_labels = {}
        self.cov_libs = cov_libs
        self.data_dict = data_dict
        self.cov_libs_total_reads_dict = cov_libs_total_reads_d
        self.stats = {}
        self.exclude_groups = []
        self.colours = {}
        self.group_order = []
        self.plot_order = []
        self.min_cov = 0.01
        self.max_cov = 0.0
        self.out_f = ''
        self.title = ''
        self.max_group_plot = 0
        self.format = ''

    def get_stats_for_group(self, group):
        stats = { 'name' : group, 
                  'count_total' : "{:,}".format(self.stats[group]['count']), 
                  'count_visible' : "{:,}".format(self.stats[group]['count_visible']), 
                  'count_visible_perc' : '{0:.1%}'.format(self.stats[group]['count_visible']/self.stats[group]['count']) if self.stats[group]['count'] > 0 else '{0:.1%}'.format(0.0), 
                  'span_visible' : "{:,}".format(self.stats[group]['span_visible']),
                  'span_total' : "{:,}".format(self.stats[group]['span']), 
                  'span_visible_perc' : '{0:.1%}'.format(self.stats[group]['span_visible']/self.stats[group]['span']) if self.stats[group]['span'] > 0 else '{0:.1%}'.format(0.0),
                  'colour' : str(self.colours[group] if group in self.colours else None),
                  'n50' : "{:,}".format(self.stats[group]['n50']),
                  'gc_mean' : "{0:.2}".format(self.stats[group]['gc_mean']),
                  'gc_std' : "{0:.2}".format(self.stats[group]['gc_std']),
                  'cov_mean' : {cov_lib : "{0:0.1f}".format(cov_mean) for cov_lib, cov_mean in self.stats[group]['cov_mean'].items()},
                  'cov_std' : {cov_lib : "{0:0.1f}".format(cov_std) for cov_lib, cov_std in self.stats[group]['cov_std'].items()},
                  'reads_mapped' : {cov_lib : "{:,}".format(reads_mapped) for cov_lib, reads_mapped in self.stats[group]['reads_mapped'].items()},
                  'reads_mapped_perc' : {cov_lib : '{0:.1%}'.format(reads_mapped_perc) for cov_lib, reads_mapped_perc in self.stats[group]['reads_mapped_perc'].items()}
                }
        return stats

    def write_stats(self):
        stats = []
        stats.append(self.get_stats_for_group('all'))
        for group in self.plot_order: # group/label/other that has been plotted
            stats.append(self.get_stats_for_group(group))
            if not group in self.group_labels: # it is either a label or "other"
                label = group
                for g, labels in self.group_labels.items():
                    if label in labels:
                        stats.append(self.get_stats_for_group(g))
        
        out_f = "%s.stats.txt" % self.out_f
        with open(out_f, 'w') as fh:
            for cov_lib in sorted(self.cov_libs):
                fh.write("# %s - %s\n" % (self.out_f, cov_lib))
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('Group', 'colour', 'count', 'visible (%)', 'span', 'visible(%)', 'n50', 'GC', 'GC (std)', 'cov_mean', 'cov_std', 'read map', 'read map (%)'))        
                for stat in stats:
                    fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (\
                            stat['name'], stat['colour'], stat['count_visible'], stat['count_visible_perc'], stat['span_visible'], \
                            stat['span_visible_perc'], stat['n50'], stat['gc_mean'], stat['gc_std'], stat['cov_mean'][cov_lib], \
                            stat['cov_std'][cov_lib], stat['reads_mapped'][cov_lib], stat['reads_mapped_perc'][cov_lib]))


    def compute_stats(self):
        stats = {}
        for label in self.labels:
            stats[label] = {
                            'name' : [],
                            'gc' : [], 
                            'length': [], 
                            'covs' : {cov_lib : [] for cov_lib in self.cov_libs},
                            'cov_mean' : {cov_lib : 0.0 for cov_lib in self.cov_libs},
                            'cov_std' : {cov_lib : 0.0 for cov_lib in self.cov_libs},
                            'reads_mapped' : {cov_lib : 0 for cov_lib in self.cov_libs},
                            'reads_mapped_perc' : {cov_lib: 0.0 for cov_lib in self.cov_libs},
                            'n50' : 0,
                            'gc_mean' : 0.0,
                            'gc_std' : 0.0,
                            'groups' : set(),
                            'count' : 0,
                            'span' : 0,
                            'count_visible' : 0,
                            'span_visible' : 0,
                            'count_hidden' : 0,
                            'span_hidden' : 0
                            }

        # gather data
        for group, labels in self.group_labels.items():
            for label in labels:    
                stats[label]['name'] = self.data_dict[group]['name']
                stats[label]['groups'].add(group)
                stats[label]['gc'] = stats[label]['gc'] + self.data_dict[group]['gc']
                stats[label]['length'] = stats[label]['length'] + self.data_dict[group]['length']
                stats[label]['count'] += self.data_dict[group]['count']
                stats[label]['span'] += self.data_dict[group]['span']
                stats[label]['count_visible'] += self.data_dict[group]['count_visible']
                stats[label]['count_hidden'] += self.data_dict[group]['count_hidden']
                stats[label]['span_visible'] += self.data_dict[group]['span_visible']
                stats[label]['span_hidden'] += self.data_dict[group]['span_hidden']
                for cov_lib in self.cov_libs:
                    stats[label]['covs'][cov_lib] = stats[label]['covs'][cov_lib] + self.data_dict[group]['covs'][cov_lib]
                    stats[label]['reads_mapped'][cov_lib] += self.data_dict[group]['reads_mapped'][cov_lib]
  
        for label in stats:
            stats[label]['gc_mean'] = mean(array(stats[label]['gc'])) if stats[label]['count_visible'] > 0 else 0.0
            stats[label]['gc_std'] = std(array(stats[label]['gc'])) if stats[label]['count_visible'] > 0 else 0.0
            stats[label]['n50'] = n50(stats[label]['length']) if stats[label]['count_visible'] > 0 else 0.0
            for cov_lib in self.cov_libs:
                stats[label]['cov_mean'][cov_lib] = mean(array(stats[label]['covs'][cov_lib])) if stats[label]['count_visible'] > 0 else 0.0
                stats[label]['cov_std'][cov_lib] = std(array(stats[label]['covs'][cov_lib])) if stats[label]['count_visible'] > 0 else 0.0
                if (self.cov_libs_total_reads_dict[cov_lib]):
                    stats[label]['reads_mapped_perc'][cov_lib] = stats[label]['reads_mapped'][cov_lib]/self.cov_libs_total_reads_dict[cov_lib]
        self.stats = stats
    
    def relabel_and_colour(self, colour_f, user_labels):
        if (colour_f):
            colour_dict = BtIO.parseColourDict(colour_f)
        else:
            groups = self.group_order[0:self.max_group_plot]
            colour_groups = [group if not (group in user_labels) else user_labels[group] for group in groups]
            colour_dict = generateColourDict(colour_groups)
        for idx, group in enumerate(self.group_order):
            if (self.exclude_groups):
                if group in self.exclude_groups:
                    self.group_labels[group].add('other')
                    self.colours[group] = WHITE     
            elif group in user_labels:
                label = user_labels[group]
                self.group_labels[group].add(label)
                self.group_labels[group].add(group)
                self.colours[label] = colour_dict[label]
                if label not in self.plot_order:
                    self.plot_order.append(label)
            elif group in colour_dict:    
                self.group_labels[group].add(group)
                self.colours[group] = colour_dict[group] 
                self.plot_order.append(group)
            elif idx > self.max_group_plot:
                self.group_labels[group].add('other')
                self.group_labels[group].add(group)
                self.colours['other'] = WHITE
                self.labels.add('other')
            else:
                self.group_labels[group].add('other')
                self.group_labels[group].add(group)
                self.colours['other'] = WHITE
                self.labels.add('other')
            self.group_labels[group].add('all')
        if 'other' in self.labels:
            self.plot_order.append('other')

    def plotScatterCov(self, cov_lib, cov_dict, info_flag, x_label, y_label, scale, x_max, y_max):
        rect_scatter, rect_histx, rect_histy, rect_legend = set_canvas()
        # Setting up plots and axes
        plt.figure(1, figsize=(35,35), dpi=400)
        axScatter = plt.axes(rect_scatter, axisbg=BGGREY, yscale = scale, xscale = scale)
        axScatter = set_format_covplot(axScatter, x_max, y_max, x_label, y_label, scale, cov_lib)
        axHistx = plt.axes(rect_histx, axisbg=BGGREY)
        axHistx = set_format_hist_x_cov(axHistx, axScatter, scale)
        axHisty = plt.axes(rect_histy, axisbg=BGGREY)
        axHisty = set_format_hist_y_cov(axHisty, axScatter, scale)
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
        top_bins = logspace(-2, (int(math.log(self.max_cov)) + 1), 200, base=10.0)
        right_bins = logspace(-2, (int(math.log(self.max_cov)) + 1), 200, base=10.0)
        # empty handles for big legend
        legend_handles = []
        legend_labels = []
        # counter necessary for multiplot so that PNGs are in order when sorted by name
        i = 0
        for group in self.plot_order:
            i += 1
            group_length_array = array(self.stats[group]['length'])
            group_cov_y_array = array([cov_dict[name] for name in self.stats[group]['name']])
            group_cov_x_array = array(self.stats[group]['covs'][cov_lib])
            # calculate values for legend
            if len(group_length_array) > 0:
                group_span_in_mb = round(self.stats[group]['span_visible']/1000000, 2)
                group_number_of_seqs = self.stats[group]['count_visible']
                group_n50 = self.stats[group]['n50']
                blob_size_array = []
                s, lw, alpha, colour = 15, 0.5, 0.5, self.colours[group]
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
                    axHistx.hist(group_cov_x_array, weights=weights_array, color = colour, bins = top_bins, histtype='step', lw = 3)
                    axHisty.hist(group_cov_y_array, weights=weights_array, color = colour, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
                else:           
                    axHistx.hist(group_cov_x_array, color = colour, bins = top_bins, histtype='step', lw = 3)
                    axHisty.hist(group_cov_y_array, color = colour, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
                axScatter.scatter(group_cov_x_array, group_cov_y_array, color = colour, s = blob_size_array, lw = lw, alpha=alpha, edgecolor=BLACK, label=label)
                axLegend.axis('off')
                if (self.multiplot): 
                    axLegend.legend(legend_handles, legend_labels, loc=6, numpoints=1, fontsize=FONTSIZE, frameon=True)
                    plot_ref_legend(axScatter)
                    m_out_f = "%s.%s.%s.compare_cov.%s" % (self.out_f, i, group.replace("/", "_").replace(" ", "_"), self.format)
                    print BtLog.status_d['8'] % m_out_f
                    plt.savefig(m_out_f, format=self.format)
        if not (self.ignore_contig_length):
            plot_ref_legend(axScatter)
        axLegend.legend(legend_handles, legend_labels, numpoints=1, fontsize=FONTSIZE, frameon=True, loc=6 )
        out_f = "%s.%s.compare_cov.%s" % (self.out_f, cov_lib, self.format)
        print BtLog.status_d['8'] % out_f
        plt.savefig(out_f, format=self.format) 
        plt.close()  

    def plotReadCov(self, refcov_dict):
        mat.rcParams.update({'font.size': 24})

        main_columns = 2
        if (refcov_dict):
            main_columns += 2
        group_columns = len(self.plot_order)
        
        #for group in self.stats:
        #    print group, self.stats[group]['reads_mapped']
        for cov_lib in self.cov_libs:
            print cov_lib
            plot_data = {}
            if not self.cov_libs_total_reads_dict[cov_lib] == 0:
                main_plot = PlotReadObj()
                group_plot = PlotReadObj()
                # unmapped (assembly)
                reads_total = self.cov_libs_total_reads_dict[cov_lib]
                reads_unmapped = reads_total - self.stats['all']['reads_mapped'][cov_lib]
                if cov_lib in refcov_dict:
                    reads_total_ref = refcov_dict[cov_lib]['reads_total']
                    reads_mapped_ref = refcov_dict[cov_lib]['reads_mapped'] 
                    reads_unmapped_ref = reads_total_ref - reads_mapped_ref
                    main_plot.labels.append('Unmapped (ref)')
                    main_plot.values.append(reads_unmapped_ref/reads_total_ref)
                    main_plot.colours.append(DGREY)
                    main_plot.labels.append('Mapped (ref)')
                    main_plot.values.append(reads_mapped_ref/reads_total_ref)
                    main_plot.colours.append(DGREY)

                main_plot.labels.append('Unmapped (assembly)')
                main_plot.values.append(reads_unmapped/reads_total)
                main_plot.colours.append(DGREY)
                # mapped (assembly)
                main_plot.labels.append('Mapped (assembly)')
                main_plot.values.append(self.stats['all']['reads_mapped_perc'][cov_lib])
                main_plot.colours.append(DGREY)
                # mapped plotted groups
                for group in self.plot_order:
                    group_plot.labels.append(group)
                    group_plot.values.append(self.stats[group]['reads_mapped_perc'][cov_lib])
                    group_plot.colours.append(self.colours[group])
                
                plot_data[cov_lib] = {'main' : main_plot, 'group' : group_plot}
                x_pos_main = arange(main_columns)
                x_pos_group = arange(len(self.plot_order))
        
                fig = plt.figure(1, figsize=(30, 10), dpi=200)  
                gs = mat.gridspec.GridSpec(1, 2, width_ratios=[main_columns, len(self.plot_order)]) 
                ax_main = plt.subplot(gs[0])
                ax_main.set_axis_bgcolor(BGGREY)
                ax_group = plt.subplot(gs[1])
                ax_group.set_axis_bgcolor(BGGREY)
                rect_group = ax_group.bar(x_pos_group, plot_data[cov_lib]['group'].values, width = 0.5, tick_label=plot_data[cov_lib]['group'].labels, align='center', color = plot_data[cov_lib]['group'].colours)
                for rect_g in rect_group:
                    height_g = float(rect_g.get_height())
                    ax_group.text(rect_g.get_x() + rect_g.get_width()/2., 0.005 + height_g, '{:.1f}%'.format(height_g*100), ha='center', va='bottom')
                rect_main = ax_main.bar(x_pos_main, plot_data[cov_lib]['main'].values, width = 0.5, tick_label=plot_data[cov_lib]['main'].labels, align='center', color = plot_data[cov_lib]['main'].colours)
                for rect_m in rect_main:
                    height_m = float(rect_m.get_height())
                    ax_main.text(rect_m.get_x() + rect_m.get_width()/2., 0.005 + height_m, '{:.1f}%'.format(height_m*100), ha='center', va='bottom')
                ax_main.set_ylim(0, 1.1)
                ax_group.set_ylim(0, 1.1)
                ax_main.set_yticklabels(['{:.0f}%'.format(x*100) for x in ax_main.get_yticks()])
                ax_main.set_xticklabels(plot_data[cov_lib]['main'].labels, rotation=45)
                ax_group.set_yticklabels(['{:.0f}%'.format(x*100) for x in ax_group.get_yticks()])
                ax_group.set_xticklabels(plot_data[cov_lib]['group'].labels, rotation=45)
                ax_main.grid(True,  axis='y', which="major", lw=2., color=WHITE, linestyle='--') 
                ax_group.grid(True,  axis='y', which="major", lw=2., color=WHITE, linestyle='--')
                out_f = "%s.%s.read_cov.%s" % (self.out_f, cov_lib, self.format)
                print BtLog.status_d['8'] % out_f
                plt.tight_layout()
                plt.savefig(out_f, format=self.format)
                plt.close()            
                
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
            group_length_array = array(self.stats[group]['length'])
            group_gc_array = array(self.stats[group]['gc'])
            group_cov_array = array(self.stats[group]['covs'][cov_lib])
            if len(group_length_array) > 0:
                # calculate values for legend
                group_span_in_mb = round(self.stats[group]['span_visible']/1000000, 2)
                group_number_of_seqs = self.stats[group]['count_visible']
                group_n50 = self.stats[group]['n50']
                blob_size_array = []
                s, lw, alpha, colour = 15, 0.5, 0.5, self.colours[group]
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
                    m_out_f = "%s.%s.%s.blobs.%s" % (self.out_f, i, group.replace("/", "_").replace(" ", "_"), self.format)
                    print BtLog.status_d['8'] % m_out_f
                    plt.savefig(m_out_f, format=self.format)
        if not (self.ignore_contig_length):
            plot_ref_legend(axScatter)
        axLegend.legend(legend_handles, legend_labels, numpoints=1, fontsize=FONTSIZE, frameon=True, loc=6 )
        out_f = "%s.%s.blobs.%s" % (self.out_f, cov_lib, self.format)
        print BtLog.status_d['8'] % out_f
        plt.savefig(out_f, format=self.format) 
        plt.close()  
        
        

    

if __name__ == "__main__": 
    pass