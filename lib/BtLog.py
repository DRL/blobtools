#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File        : BtLog.py
Version     : 0.1
Author      : Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs        : ?
To do       : ?
"""
from __future__ import division
import sys
from os.path import basename

def error(message, *argv):
    if argv is None:
        print error_d[message]
    else:
        print error_d[message] % (argv)
    exit(1)

def warn(message, *argv):
    if (argv):
        print warn_d[message] % (argv)
    else:
        print warn_d[message]

def progress(iteration, steps, max_value):
    if int(iteration) % int(steps) == 0:
        sys.stdout.write('\r')
        print "[PROGRESS]\t: %d%%" % (float(int(iteration)/int(max_value)*100)),
        sys.stdout.flush()
    if int(iteration == max_value):
        sys.stdout.write('\r')
        print "[PROGRESS]\t: %d%%" % (100)

error_d = {
    '0' : '[ERROR:0]\t: File %s does not exist.',
    '1' : '[ERROR:1]\t: Please provide coverage information.',
    '2' : '[ERROR:2]\t: Assembly type is not valid (%s).',
    '3' : '[ERROR:3]\t: nodesDB ("--nodes" and "--names" or "--db") not specified.',
    '4' : '[ERROR:4]\t: BlobDB.parseFasta() - no sequences found. Check FASTA file.',
    '5' : '[ERROR:5]\t: Sequence header %s is not unique.',
    '6' : '[ERROR:6]\t: BlobDB.readBam() - sequence header %s in %s was not in FASTA.',
    '7' : '[ERROR:7]\t: Please add "samtools" to you PATH variable.',
    '8' : '[ERROR:8]\t: Unsupported taxrule "%s".',
    '9' : '[ERROR:9]\t: Unsupported taxonomic rank "%s".',
    '10' : '[ERROR:10]\t: Unsupported output format "%s".',
    '11' : '[ERROR:11]\t: Taxrule "%s" was not computed for this BlobDb. Available taxrule(s) : %s.',
    '12' : '[ERROR:12]\t: Please provide an output file.',
    '13' : '[ERROR:13]\t: %s does not appear to be a comma-separated list or a file.',
    '14' : '[ERROR:14]\t: Unsupported sort order for plotting : %s. Must be either "span" or "count".',
    '15' : '[ERROR:15]\t: Unsupported histogram type for plotting : %s. Must be either "span" or "count".',
    '16' : '[ERROR:16]\t: Cluster string "%s" seems to be in the wrong format".',
    '17' : '[ERROR:17]\t: Group "%s" was specified in multiple clusters".',
    '18' : '[ERROR:18]\t: Please provide a tax file in BLAST format.',
    '19' : '[ERROR:19]\t: Sequence %s in file %s is not part of the assembly.',
    '20' : '[ERROR:20]\t: Please add "clc_mapping_info" to you PATH variable.'

}

warn_d = {
    '0' : '[WARN:0]\t: No tax files specified.' 
}
status_d = {
    '1' : '[STATUS]\t: Parsing %s - %s',
    '2' : '... Done',
    '3' : '[STATUS]\t: Creating nodesDB from %s and %s',
    '4' : '[STATUS]\t: Retrieving %s nodesDB from %s',
    '5' : '[STATUS]\t: Store nodesDB in location %s',
    '6' : '[STATUS]\t: Computing taxonomy using taxrule(s) %s',
    '7' : '[STATUS]\t: Writing BlobDB to file %s',
    '8' : '[STATUS]\t: Plotting %s',
    '9' : '[STATUS]\t: Reading BlobDb %s'

}

info_d = {
    '0' : '[INFO]\t: %s : sequences = %s, span = %s MB, N50 = %s nt'
    }

if __name__ == "__main__": 
    pass