#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File        : BtPlot.py
Author      : Dominik R. Laetsch, dominik.laetsch at gmail dot com
"""

from __future__ import division
from numpy import array, arange, logspace, mean, std
import math
import bloblib.BtLog as BtLog
import bloblib.BtIO as BtIO
import bloblib.BtTax as BtTax
import matplotlib as mat
from matplotlib import cm
from matplotlib.ticker import NullFormatter, MultipleLocator, AutoMinorLocator
from matplotlib.lines import Line2D
from matplotlib.colors import rgb2hex
mat.use('agg')
import matplotlib.pyplot as plt
from itertools import izip

mat.rcParams.update({'font.size': 36})
mat.rcParams['xtick.major.pad'] = '8'
mat.rcParams['ytick.major.pad'] = '8'
mat.rcParams['lines.antialiased'] = True

LEGEND_FONTSIZE = 24
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

def set_format_scatterplot(axScatter, **kwargs):
    min_x, max_x = None, None
    min_y, max_y = None, None
    if kwargs['plot'] == 'blobplot':
        min_x, max_x = 0, 1
        major_xticks = MultipleLocator(0.2)
        minor_xticks = AutoMinorLocator(20)
        min_y, max_y = 0.005, kwargs['max_cov']+1000
        axScatter.set_yscale('log')
        axScatter.set_xscale('linear')
        axScatter.xaxis.set_major_locator(major_xticks)
        axScatter.xaxis.set_minor_locator(minor_xticks)
    elif kwargs['plot'] == 'covplot':
        min_x, max_x = 0.005, kwargs['max_cov']+1000
        min_y, max_y = 0.005, kwargs['max_cov']+1000
        axScatter.set_yscale('log')
        axScatter.set_xscale('log')
    else:
        BtLog.error('34' % kwargs['plot'])
    axScatter.set_xlim( (min_x, max_x) )
    axScatter.set_ylim( (min_y, max_y) ) # This sets the max-Coverage so that all libraries + sum are at the same scale
    axScatter.grid(True, which="major", lw=2., color=WHITE, linestyle='-')
    axScatter.set_axisbelow(True)
    axScatter.xaxis.labelpad = 20
    axScatter.yaxis.labelpad = 20
    axScatter.yaxis.get_major_ticks()[0].label1.set_visible(False)
    axScatter.tick_params(axis='both', which='both', direction='out')
    return axScatter

def set_format_hist_x(axHistx, axScatter):
    axHistx.set_xlim(axScatter.get_xlim())
    axHistx.set_xscale(axScatter.get_xscale())
    axHistx.grid(True, which="major", lw=2., color= WHITE, linestyle='-')
    axHistx.xaxis.set_major_locator(axScatter.xaxis.get_major_locator()) # no labels since redundant
    axHistx.xaxis.set_minor_locator(axScatter.xaxis.get_minor_locator())
    axHistx.xaxis.set_major_formatter(nullfmt) # no labels since redundant
    axHistx.set_axisbelow(True)
    axHistx.yaxis.labelpad = 20
    axHistx.tick_params(axis='both', which='both', direction='out')
    return axHistx

def set_format_hist_y(axHisty, axScatter):
    axHisty.set_ylim(axScatter.get_ylim())
    axHisty.set_yscale(axScatter.get_yscale())
    axHisty.grid(True, which="major", lw=2., color= WHITE, linestyle='-')
    axHisty.yaxis.set_major_formatter(nullfmt) # no labels since redundant
    axHisty.set_axisbelow(True)
    axHisty.xaxis.labelpad = 20
    axHisty.tick_params(axis='both', which='both', direction='out')
    return axHisty

def get_ref_label(max_length, max_marker_size, fraction):
    length = int(math.ceil(fraction * max_length / 100.0)) * 100
    string = "%snt" % "{:,}".format(length)
    markersize = length/max_length * max_marker_size
    return length, string, markersize

def plot_ref_legend(axScatter, max_length, max_marker_size, ignore_contig_length):
    if not (ignore_contig_length):
        ref1_length, ref1_string, ref1_markersize = get_ref_label(max_length, max_marker_size, 0.05)
        ref2_length, ref2_string, ref2_markersize = get_ref_label(max_length, max_marker_size, 0.1)
        ref3_length, ref3_string, ref3_markersize = get_ref_label(max_length, max_marker_size, 0.25)
        # markersize in scatter is in "points^2", markersize in Line2D is in "points" ... that's why we need math.sqrt()
        ref_1 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(ref1_markersize),  markerfacecolor=GREY))
        ref_2 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(ref2_markersize), markerfacecolor=GREY))
        ref_3 = (Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=math.sqrt(ref3_markersize), markerfacecolor=GREY))
        axScatter.legend([ref_1,ref_2,ref_3], [ref1_string, ref2_string, ref3_string], numpoints=1, ncol = 3, loc = 8, fontsize=LEGEND_FONTSIZE, borderpad=1.2, labelspacing=1.8, handlelength=1, handletextpad=1)

def plot_legend(fig, axLegend, out_f, legend_flag, format, cumulative_flag):
    if (legend_flag):
        extent = axLegend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        legend_out_f = '%s.%s.%s' % (out_f, "legend", format)
        print BtLog.status_d['8'] % legend_out_f
        fig.savefig('%s' % legend_out_f, bbox_inches=extent, format=format)
        fig.delaxes(axLegend)
    return fig

def check_input(args):
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
    no_title = args['--notitle']
    ignore_contig_length = args['--noscale']
    labels = args['--label']
    colour_f = args['--colours']
    exclude_groups = args['--exclude']
    format = args['--format']
    no_plot_blobs = args['--noblobs']
    no_plot_reads = args['--noreads']
    refcov_f = args['--refcov']
    catcolour_f = args['--catcolour']
    legend_flag = args['--legend']
    cumulative_flag = args['--cumulative']
    cov_lib_selection = args['--lib']

    if 'blobplot' in args or 'covplot' in args:
        # Are ranks sane ?
        if rank not in BtTax.RANKS:
            BtLog.error('9', rank)
        # is taxrule provided?
        if taxrule not in BtTax.TAXRULES:
            BtLog.error('8', taxrule)
        # Are sort_order and hist_type sane?
        if not sort_order in ['span', 'count']:
            BtLog.error('14', sort_order)
        if not hist_type in ['span', 'count']:
            BtLog.error('15', hist_type)
        if (catcolour_f) and (c_index):
            BtLog.error('24')
        if (cumulative_flag) and (multiplot):
            BtLog.error('32')
    return args

class PlotObj():
    def __init__(self, data_dict, cov_lib_dict, cov_lib_selection, plot_type):
        self.labels = {'all'}
        self.plot = plot_type # type of plot
        self.group_labels = {}
        self.cov_lib_dict = cov_lib_dict
        self.cov_libs = self.subselect_cov_libs(cov_lib_dict, cov_lib_selection)
        self.cov_libs_total_reads_dict = self.get_cov_libs_total_reads_dict(cov_lib_dict)
        self.cov_libs_mapped_reads_dict = self.get_cov_libs_mapped_reads_dict(cov_lib_dict)
        self.data_dict = data_dict
        self.stats = {}
        self.exclude_groups = []
        self.version = None
        self.colours = {}
        self.group_order = []
        self.plot_order = []
        self.min_cov = 0.01
        self.max_cov = 0.0
        self.out_f = ''
        self.no_title = ''
        self.max_group_plot = 0
        self.format = ''
        self.legend_flag = ''
        self.cumulative_flag = ''

        self.cov_y_dict = {}
        self.xlabel = None
        self.ylabel = None

        self.refcov_dict = {}

    def subselect_cov_libs(self, cov_lib_dict, cov_lib_selection):
        selected_cov_libs = []
        cov_lib_selection_error = 0
        if (cov_lib_selection):
            if cov_lib_selection == 'covsum':
                selected_cov_libs.append('covsum')
            elif "," in cov_lib_selection:
                selected_cov_libs = cov_lib_selection.split(",")
                if not set(selected_cov_libs).issubset(set(cov_lib_dict.keys())):
                    cov_lib_selection_error = 1
            else:
                selected_cov_libs.append(cov_lib_selection)
                if not cov_lib_selection in cov_lib_dict:
                    cov_lib_selection_error = 1
        else:
            selected_cov_libs = cov_lib_dict.keys()
        if cov_lib_selection_error:
            covlib_string = []
            for covlib in cov_lib_dict:
                cov_lib_f = cov_lib_dict[covlib]['f']
                if not cov_lib_f:
                    cov_lib_f = "sum of coverages from all covlibs"
                covlib_string.append("\t\t%s : %s" % (covlib, cov_lib_f))
            BtLog.error('33', "\n".join(covlib_string))
        return selected_cov_libs

    def get_cov_libs_total_reads_dict(self, cov_lib_dict):
        return { x : cov_lib_dict[x]['reads_total'] for x in self.cov_libs}

    def get_cov_libs_mapped_reads_dict(self, cov_lib_dict):
        return { x : cov_lib_dict[x]['reads_mapped'] for x in self.cov_libs}

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

    def write_stats(self, out_f):
        stats = []
        stats.append(self.get_stats_for_group('all'))
        for group in self.plot_order: # group/label/other that has been plotted
            stats.append(self.get_stats_for_group(group))
            if not group in self.group_labels: # it is either a label or "other"
                label = group
                for g, labels in self.group_labels.items():
                    if label in labels:
                        stats.append(self.get_stats_for_group(g))
        output = []
        output.append('## %s' % self.version)
        for cov_lib, cov_lib_dict in self.cov_lib_dict.items():
           if cov_lib in self.cov_libs:
                output.append("## %s=%s" % (cov_lib, cov_lib_dict['f']))
        fields = ['name', 'colour', 'count_visible', 'count_visible_perc', 'span_visible','span_visible_perc', 'n50', 'gc_mean', 'gc_std']
        header = [field for field in fields]
        for cov_lib in sorted(self.cov_libs):
            header.append('%s_mean' % cov_lib)
            header.append('%s_std' % cov_lib)
            header.append('%s_read_map' % cov_lib)
            header.append('%s_read_map_p' % cov_lib)
        output.append('# %s' % "\t".join(header))
        for stat in stats:
            line = []
            for field in fields:
                line.append(stat[field])
            for cov_lib in sorted(self.cov_libs):
                line.append(stat['cov_mean'][cov_lib])
                line.append(stat['cov_std'][cov_lib])
                line.append(stat['reads_mapped'][cov_lib])
                line.append(stat['reads_mapped_perc'][cov_lib])
            output.append("%s" % "\t".join(line))
        out_f = "%s.stats.txt" % out_f
        with open(out_f, 'w') as fh:
            print BtLog.status_d['24'] % ("%s" % out_f)
            fh.write("\n".join(output))

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

        for group, labels in self.group_labels.items():
            for label in labels:
                stats[label]['name'] = stats[label]['name'] + self.data_dict[group]['name']
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
            stats[label]['gc_mean'] = mean(array(stats[label]['gc'])) if stats[label]['count_visible'] > 0.0 else 0.0
            stats[label]['gc_std'] = std(array(stats[label]['gc'])) if stats[label]['count_visible'] > 0.0 else 0.0
            stats[label]['n50'] = n50(stats[label]['length']) if stats[label]['count_visible'] > 0.0 else 0.0
            for cov_lib in self.cov_libs:
                stats[label]['cov_mean'][cov_lib] = mean(array(stats[label]['covs'][cov_lib])) if stats[label]['count_visible'] > 0.0 else 0.0
                stats[label]['cov_std'][cov_lib] = std(array(stats[label]['covs'][cov_lib])) if stats[label]['count_visible'] > 0.0 else 0.0
                if self.cov_libs_total_reads_dict[cov_lib]:
                    stats[label]['reads_mapped_perc'][cov_lib] = stats[label]['reads_mapped'][cov_lib]/self.cov_libs_total_reads_dict[cov_lib]
        self.stats = stats

    def relabel_and_colour(self, colour_dict, user_labels):

        if (colour_dict):
            groups = self.group_order[0:self.max_group_plot]
            groups_not_in_colour_dict = set(groups) - set(colour_dict.keys())
            for group in groups_not_in_colour_dict:
                colour_dict[group] = WHITE
        else:
            groups = self.group_order[0:self.max_group_plot]
            colour_groups = [group if not (group in user_labels) else user_labels[group] for group in groups]
            colour_dict = generateColourDict(colour_groups)
        for idx, group in enumerate(self.group_order):
            if group in self.exclude_groups:
                pass
            elif group in user_labels:
                label = user_labels[group]
                self.group_labels[group].add(label)
                self.group_labels[group].add(group)
                self.colours[label] = colour_dict[group]
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


    def setupPlot(self, plot):
        if plot == 'blobplot' or plot == 'covplot':
            rect_scatter, rect_histx, rect_histy, rect_legend = set_canvas()
            # Setting up plots and axes
            fig = plt.figure(1, figsize=(35,35), dpi=400)
            axScatter = plt.axes(rect_scatter, axisbg=BGGREY)
            axScatter = set_format_scatterplot(axScatter, max_cov=self.max_cov, plot=plot)
            axScatter.set_xlabel(self.xlabel)
            axScatter.set_ylabel(self.ylabel)
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
            for xtick in axHisty.get_xticklabels(): # rotate text for ticks in cov histogram
                xtick.set_rotation(270)
            axLegend = plt.axes(rect_legend, axisbg=WHITE)
            axLegend.xaxis.set_major_locator(plt.NullLocator())
            axLegend.xaxis.set_major_formatter(nullfmt)
            axLegend.yaxis.set_major_locator(plt.NullLocator())
            axLegend.yaxis.set_major_formatter(nullfmt)
            top_bins, right_bins = None, None
            if plot == 'blobplot':
                top_bins = arange(0, 1, 0.01)
                right_bins = logspace(-2, (int(math.log(self.max_cov)) + 1), 200, base=10.0)
            elif plot == 'covplot':
                top_bins = logspace(-2, (int(math.log(self.max_cov)) + 1), 200, base=10.0)
                right_bins = logspace(-2, (int(math.log(self.max_cov)) + 1), 200, base=10.0)
            else:
                pass
            return fig, axScatter, axHistx, axHisty, axLegend, top_bins, right_bins
        elif plot == 'readcov':
            main_columns = 2
            if (self.refcov_dict):
                main_columns += 2
            group_columns = len(self.plot_order)
            fig = plt.figure(1, figsize=(30, 10), dpi=200)
            gs = mat.gridspec.GridSpec(1, 2, width_ratios=[main_columns, group_columns])

            ax_main = plt.subplot(gs[0])
            ax_main.set_axis_bgcolor(BGGREY)
            ax_main.set_ylim(0, 1.1)
            ax_main.set_yticklabels(['{:.0f}%'.format(x*100) for x in ax_main.get_yticks()])
            ax_main.grid(True,  axis='y', which="major", lw=2., color=WHITE, linestyle='--')

            ax_group = plt.subplot(gs[1])
            ax_group.set_axis_bgcolor(BGGREY)
            ax_group.set_ylim(0, 1.1)
            ax_group.set_yticklabels(['{:.0f}%'.format(x*100) for x in ax_group.get_yticks()])
            ax_group.grid(True,  axis='y', which="major", lw=2., color=WHITE, linestyle='--')

            x_pos_main = arange(main_columns)
            x_pos_group = arange(len(self.plot_order))
            return fig, ax_main, ax_group, x_pos_main, x_pos_group
        else:
            return None

    def plotBar(self, cov_lib, out_f):
        fig, ax_main, ax_group, x_pos_main, x_pos_group = self.setupPlot('readcov')
        ax_main_data = {'labels' : [], 'values' : [], 'colours' : [] }
        ax_group_data = {'labels' : [], 'values' : [], 'colours' : [] }
        reads_total = self.cov_libs_total_reads_dict[cov_lib]
        reads_mapped = self.stats['all']['reads_mapped'][cov_lib]
        reads_unmapped = reads_total - self.stats['all']['reads_mapped'][cov_lib]
        ax_main_data['labels'].append('Unmapped (assembly)')
        ax_main_data['values'].append(reads_unmapped/reads_total)
        ax_main_data['colours'].append(DGREY)
        ax_main_data['labels'].append('Mapped (assembly)')
        ax_main_data['values'].append(reads_mapped/reads_total)
        ax_main_data['colours'].append(DGREY)
        if (self.refcov_dict):
            if cov_lib in self.refcov_dict:
                reads_total_ref = self.refcov_dict[cov_lib]['reads_total']
                reads_mapped_ref = self.refcov_dict[cov_lib]['reads_mapped']
                reads_unmapped_ref = reads_total_ref - reads_mapped_ref
                ax_main_data['labels'].append('Unmapped (ref)')
                ax_main_data['values'].append(reads_unmapped_ref/reads_total_ref)
                ax_main_data['colours'].append(DGREY)
                ax_main_data['labels'].append('Mapped (ref)')
                ax_main_data['values'].append(reads_mapped_ref/reads_total_ref)
                ax_main_data['colours'].append(DGREY)
            else:
                BtLog.error('40', cov_lib)

        # mapped plotted groups
        for group in self.plot_order:
           ax_group_data['labels'].append(group)
           ax_group_data['values'].append(self.stats[group]['reads_mapped_perc'][cov_lib])
           ax_group_data['colours'].append(self.colours[group])
        rect_group = ax_group.bar(x_pos_group, ax_group_data['values'], width = 0.5, tick_label=ax_group_data['labels'], align='center', color = ax_group_data['colours'])
        for rect_g in rect_group:
            height_g = float(rect_g.get_height())
            ax_group.text(rect_g.get_x() + rect_g.get_width()/2., 0.005 + height_g, '{:.2f}%'.format(height_g*100), ha='center', va='bottom', fontsize=LEGEND_FONTSIZE)
        rect_main = ax_main.bar(x_pos_main, ax_main_data['values'], width = 0.5, tick_label=ax_main_data['labels'], align='center', color = ax_main_data['colours'])
        for rect_m in rect_main:
            height_m = float(rect_m.get_height())
            ax_main.text(rect_m.get_x() + rect_m.get_width()/2., 0.005 + height_m, '{:.2f}%'.format(height_m*100), ha='center', va='bottom', fontsize=LEGEND_FONTSIZE)

        ax_main.set_xticklabels(ax_main_data['labels'], rotation=45, ha='center', fontsize=LEGEND_FONTSIZE)
        ax_group.set_xticklabels(ax_group_data['labels'], rotation=45, ha='center', fontsize=LEGEND_FONTSIZE)
        #figsuptitle = fig.suptitle(out_f, verticalalignment='top')
        out_f = "%s.read_cov.%s" % (out_f, cov_lib)
        print BtLog.status_d['8'] % "%s.%s" % (out_f, self.format)
        fig.tight_layout()
        #fig.savefig("%s.%s" % (out_f, self.format), format=self.format,  bbox_extra_artists=(figsuptitle,))
        fig.savefig("%s.%s" % (out_f, self.format), format=self.format)
        plt.close(fig)

    def plotScatter(self, cov_lib, info_flag, out_f):

        fig, axScatter, axHistx, axHisty, axLegend, top_bins, right_bins = self.setupPlot(self.plot)
        # empty handles for big legend
        legend_handles = []
        legend_labels = []
        # marker size scaled by biggest blob (size in points^2)
        max_length = max(array(self.stats['all']['length'])) # length of biggest blob
        max_marker_size = 12500 # marker size for biggest blob, i.e. area of 12500^2 pixel
        for idx, group in enumerate(self.plot_order):
            idx += 1
            lw, alpha = 0.5, 0.8
            if group == 'no-hit':
                alpha = 0.5
            group_length_array = array(self.stats[group]['length'])
            if len(group_length_array) > 0 and group not in self.exclude_groups:
                colour = self.colours[group]
                group_x_array = ''
                group_y_array = ''
                if self.plot == 'blobplot':
                    group_x_array = array(self.stats[group]['gc'])
                    group_y_array = array(self.stats[group]['covs'][cov_lib])
                elif self.plot == 'covplot':
                    group_x_array = array(self.stats[group]['covs'][cov_lib])
                    group_y_array = array([self.cov_y_dict.get(name, 0.02) for name in self.stats[group]['name']])
                else:
                    BtLog.error('34', self.plot)
                marker_size_array = []
                if (self.ignore_contig_length): # no scaling
                    if group == "no-hit":
                        s = 20
                    else:
                        s = 100
                    marker_size_array = [s for length in group_length_array]
                else: # scaling by max_length
                    marker_size_array = [(length/max_length)*max_marker_size for length in group_length_array]
                # generate label for legend
                group_span_in_mb = round(self.stats[group]['span_visible']/1000000, 2)
                group_number_of_seqs = self.stats[group]['count_visible']
                group_n50 = self.stats[group]['n50']
                fmt_seqs = "{:,}".format(group_number_of_seqs)
                fmt_span = "{:,}".format(group_span_in_mb)
                fmt_n50 = "{:,}".format(group_n50)
                label = "%s (%s;%sMB;%snt)" % (group, fmt_seqs, fmt_span, fmt_n50)
                if (info_flag):
                    print BtLog.info_d['0'] % (group, fmt_seqs, fmt_span, fmt_n50)
                legend_handles.append(Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=24, markerfacecolor=colour))
                legend_labels.append(label)

                weights_array = None
                if self.hist_type == "span":
                    weights_array = group_length_array/1000

                axHistx.hist(group_x_array, weights=weights_array, color = colour, bins = top_bins, histtype='step', lw = 3)
                axHisty.hist(group_y_array, weights=weights_array, color = colour, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
                axScatter.scatter(group_x_array, group_y_array, color = colour, s = marker_size_array, lw = lw, alpha=alpha, edgecolor=BLACK, label=label)
                axLegend.axis('off')
                if (self.multiplot):
                    fig_m, axScatter_m, axHistx_m, axHisty_m, axLegend_m, top_bins, right_bins = self.setupPlot(self.plot)
                    legend_handles_m = []
                    legend_labels_m = []
                    legend_handles_m.append(Line2D([0], [0], linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=24, markerfacecolor=colour))
                    legend_labels_m.append(label)
                    axHistx_m.hist(group_x_array, weights=weights_array, color = colour, bins = top_bins, histtype='step', lw = 3)
                    axHisty_m.hist(group_y_array, weights=weights_array, color = colour, bins = right_bins, histtype='step', orientation='horizontal', lw = 3)
                    axScatter_m.scatter(group_x_array, group_y_array, color = colour, s = marker_size_array, lw = lw, alpha=alpha, edgecolor=BLACK, label=label)
                    axLegend_m.axis('off')
                    axLegend_m.legend(legend_handles_m, legend_labels_m, loc=6, numpoints=1, fontsize=LEGEND_FONTSIZE, frameon=True)
                    plot_ref_legend(axScatter_m, max_length, max_marker_size, self.ignore_contig_length)
                    m_out_f = "%s.%s.%s.%s" % (out_f, cov_lib, idx, group.replace("/", "_").replace(" ", "_"))
                    fig_m = plot_legend(fig_m, axLegend_m, m_out_f, self.legend_flag, self.format, self.cumulative_flag)
                    print BtLog.status_d['8'] % "%s.%s" % (m_out_f, self.format)
                    fig_m.savefig("%s.%s" % (m_out_f, self.format), format=self.format)
                    plt.close(fig_m)
                elif (self.cumulative_flag):
                    axLegend.legend(legend_handles, legend_labels, loc=6, numpoints=1, fontsize=LEGEND_FONTSIZE, frameon=True)
                    plot_ref_legend(axScatter, max_length, max_marker_size, self.ignore_contig_length)
                    m_out_f = "%s.%s.%s.%s" % (out_f, cov_lib, idx, group.replace("/", "_").replace(" ", "_"))
                    fig.add_axes(axLegend)
                    fig = plot_legend(fig, axLegend, m_out_f, self.legend_flag, self.format, self.cumulative_flag)
                    if not (self.no_title):
                        fig.suptitle(out_f, fontsize=35, verticalalignment='top')
                    print BtLog.status_d['8'] % "%s.%s" % (m_out_f, self.format)
                    fig.savefig("%s.%s" % (m_out_f, self.format), format=self.format)
                else:
                    pass
        plot_ref_legend(axScatter, max_length, max_marker_size, self.ignore_contig_length)
        axLegend.legend(legend_handles, legend_labels, numpoints=1, fontsize=LEGEND_FONTSIZE, frameon=True, loc=6 )
        out_f = "%s.%s" % (out_f, cov_lib)
        fig.add_axes(axLegend)
        fig = plot_legend(fig, axLegend, out_f, self.legend_flag, self.format, self.cumulative_flag)
        if not (self.no_title):
            fig.suptitle(out_f, fontsize=35, verticalalignment='top')
        print BtLog.status_d['8'] % "%s.%s" % (out_f, self.format)
        fig.savefig("%s.%s" % (out_f, self.format), format=self.format)
        plt.close(fig)

if __name__ == "__main__":
    pass
