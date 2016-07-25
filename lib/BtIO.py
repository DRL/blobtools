#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File        : BtIO.py
Version     : 0.1
Author      : Dominik R. Laetsch, dominik.laetsch at gmail dot com
Bugs        : ?
To do       : ?
"""

from __future__ import division
import re
import subprocess
from os.path import basename, isfile, abspath, splitext, join, isdir
import shutil
import os
import sys
import lib.BtLog as BtLog
from collections import deque


def create_dir(directory="", overwrite=True):
    if (directory):
        if not isdir(directory):
            os.makedirs(directory)
        else:
            if (overwrite):
                shutil.rmtree(directory)           #removes all the subdirectories!
                os.makedirs(directory)
        return directory
    else:
        return None

def parseList(infile):
    if not isfile(infile):
         BtLog.error('0', infile)
    with open(infile) as fh:
        items = []
        for l in fh:
            items.append(l.rstrip("\n"))
    return items

def parseReferenceCov(infile):
    refcov_dict = {}
    if (infile):
        if not isfile(infile):
            BtLog.error('0', infile)
        with open(infile) as fh:
            for l in fh:
                try:
                    cov_lib, reads_total_ref, reads_mapped_ref = l.split(",")
                    refcov_dict[cov_lib] = {'reads_total' : int(reads_total_ref),
                                            'reads_mapped' : int(reads_mapped_ref)}
                except:
                    BtLog.error('21', infile)
    return refcov_dict

def parseCmdlist(temp):
    _list = []
    if (temp):
        if "," in temp:
            _list = temp.split(",")
        else:
            _list.append(temp)
    return _list

def parseCmdLabels(labels):
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

def parseCatColour(infile):
    catcolour_dict = {}
    if (infile):
        if not isfile(infile):
            BtLog.error('0', infile)
        with open(infile) as fh:
            for l in fh:
                try:
                    seq_name, category = l.rstrip("\n").split(",")
                    catcolour_dict[seq_name] = category
                except:
                    BtLog.error('23', infile)
    return catcolour_dict

def parseDict(infile, key, value):
    items = {}
    if (infile):
        if not isfile(infile):
            BtLog.error('0', infile)
        with open(infile) as fh:
            items = {}
            k_idx = int(key)
            v_idx = int(value)
            for l in fh:
                temp = l.rstrip("\n").split()
                items[temp[k_idx]] = temp[v_idx]
    return items

def parseColours(infile):
    items = {}
    if (infile):
        if not isfile(infile):
            BtLog.error('0', infile)
        with open(infile) as fh:
            for l in fh:
                temp = l.rstrip("\n").split(",")
                items[temp[0]] = temp[1]
    return items

def parseSet(infile):
    if not isfile(infile):
         BtLog.error('0', infile)
    with open(infile) as fh:
        items = set()
        for l in fh:
            items.add(l.rstrip("\n").lstrip(">"))
    return items

def parseFastaNameOrder(infile):
    fasta_order = []
    for name, seq in readFasta(infile):
        fasta_order.append(name)
    return fasta_order

def readFasta(infile):
    if not isfile(infile):
         BtLog.error('0', infile)
    with open(infile) as fh:
        header, seqs = '', []
        for l in fh:
            if l[0] == '>':
                if (header):
                    yield header, ''.join(seqs)
                header, seqs = l[1:-1].split()[0], [] # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, ''.join(seqs)

def runCmd(command):
    cmd = command.split() # sanitation
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True,
                         bufsize=-1) # buffersize of system
    return iter(p.stdout.readline, b'')

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def checkBam(infile):
    print BtLog.status_d['10']
    if not isfile(infile):
         BtLog.error('0', infile)
    if not (which('samtools')):
        BtLog.error('7')
    reads_mapped_re = re.compile(r"(\d+)\s\+\s\d+\smapped")
    reads_total_re = re.compile(r"(\d+)\s\+\s\d+\sin total")
    reads_total, reads_mapped = 0, 0
    output = ''
    command = "samtools flagstat " + infile
    for line in runCmd(command):
        output += line
    reads_mapped = int(reads_mapped_re.search(output).group(1))
    reads_total = int(reads_total_re.search(output).group(1))
    # check whether there are reads in BAM
    if not (reads_total) or not (reads_mapped):
          BtLog.error('29' % infile)
    print BtLog.status_d['11'] % ('{:,}'.format(reads_mapped), '{:,}'.format(reads_total), '{0:.1%}'.format(reads_mapped/reads_total))
    return reads_total, reads_mapped

def parseSam(infile, set_of_blobs, no_base_cov_flag):
    if not isfile(infile):
         BtLog.error('0', infile)
    base_cov_dict = {blob : [] for blob in set_of_blobs}
    read_cov_dict = {blob : 0 for blob in set_of_blobs}
    cigar_match_re = re.compile(r"(\d+)M|X|=") # only gets digits before M,X,='s
    reads_total = 0
    reads_mapped = 0
    if not (no_base_cov_flag):
        with open(infile) as fh:
            for line in fh:
                if line.startswith("@"):
                    pass
                else:
                    reads_total += 1
                    match = line.split()
                    if not match[2] == '*':
                        reads_mapped += 1
                        try:
                            base_cov_dict[match[2]].append(sum([int(matching) for matching in cigar_match_re.findall(match[5])]))
                            read_cov_dict[match[2]] += 1
                        except:
                            print BtLog.warn_d['2'] % (match[2])
    else:
        with open(infile) as fh:
            for line in fh:
                if line.startswith("@"):
                    pass
                else:
                    reads_total += 1
                    match = line.split()
                    if not match[2] == '*':
                        reads_mapped += 1
                        try:
                            read_cov_dict[match[2]] += 1
                        except:
                            print BtLog.warn_d['2'] % (match[2])
    base_cov_dict = {seq_name: sum(base_covs) for seq_name, base_covs in base_cov_dict.items()}
    return base_cov_dict, reads_total, reads_mapped, read_cov_dict

def parseBam(infile, set_of_blobs, no_base_cov_flag):
    '''
    checkBam returns reads_total and reads_mapped
    base_cov_dict is list of coverages for each contigs, since list appending should be faster

    '''
    if not isfile(infile):
         BtLog.error('0', infile)
    reads_total, reads_mapped = checkBam(infile)
    progress_unit = int(reads_mapped/1000)
    base_cov_dict = {blob : [] for blob in set_of_blobs}
    #base_cov_dict = {blob : 0 for blob in set_of_blobs}
    read_cov_dict = {blob : 0 for blob in set_of_blobs}
    cigar_match_re = re.compile(r"(\d+)M|X|=") # only gets digits before M,X,='s
    # execute samtools to get only mapped reads (no optial duplicates, no 2nd-ary alignment)
    command = "samtools view -F 1028 -F 4 -F 256 " + infile
    seen_reads = 0
    #import time
    #start = time.time()
    if not (no_base_cov_flag):
        for line in runCmd(command):
            seen_reads += 1
            match = line.split()
            try:
                base_cov_dict[match[2]].append(sum([int(matching) for matching in cigar_match_re.findall(match[5])]))
                #base_cov_dict[match[2]] += sum([int(matching) for matching in cigar_match_re.findall(match[5])])
                read_cov_dict[match[2]] += 1
            except:
                print BtLog.warn_d['2'] % (match[2])
            BtLog.progress(seen_reads, progress_unit, reads_mapped)
    else:
        for line in runCmd(command):
            seen_reads += 1
            match = line.split()
            try:
                read_cov_dict[match[2]] += 1
            except:
                print BtLog.warn_d['2'] % (match[2])
            BtLog.progress(seen_reads, progress_unit, reads_mapped)
    if not int(reads_mapped) == int(seen_reads):
        print BtLog.warn_d['3'] % (reads_mapped, seen_reads)
    base_cov_dict = {seq_name: sum(base_covs) for seq_name, base_covs in base_cov_dict.items()}
    #end = time.time()
    #print (end-start)
    return base_cov_dict, reads_total, reads_mapped, read_cov_dict

def parseCovFromHeader(fasta_type, header):
    '''
    Returns the coverage from the header of a FASTA
    sequence depending on the assembly type
    '''
    ASSEMBLY_TYPES = [None, 'spades', 'velvet', 'platanus']
    if not fasta_type in ASSEMBLY_TYPES:
        BtLog.error('2', ",".join(ASSEMBLY_TYPES[1:]))
    if fasta_type == 'spades':
        spades_match_re = re.compile(r"_cov_(\d+\.*\d*)")
        cov = re.findall(r"_cov_(\d+\.*\d*)", header)
        return float(spades_match_re.findall(header)[0])
    elif fasta_type == 'velvet':
        return float(header.split("_")[-1])
    #elif fasta_type == 'abyss' or fasta_type == 'soap':
    #    temp = header.split(" ")
    #    return float(temp[2]/(temp[1]+1-75))
    elif fasta_type == 'platanus':
        temp = header.rstrip("\n").split("_")
        if len(temp) >= 3:
            return float(temp[2].replace("cov", "")) # scaffold/scaffoldBubble/contig
        else:
            return float(temp[1].replace("cov", "")) # gapClosed
    else:
        pass

def parseCov(infile, set_of_blobs):
    if not isfile(infile):
         BtLog.error('0', infile)
    old_cov_line_re = re.compile(r"^(\S+)\t(\d+\.*\d*)")
    base_cov_dict = {}

    cov_line_re = re.compile(r"^(\S+)\t(\d+\.*\d*)\t(\d+\.*\d*)")
    reads_total = 0
    reads_mapped = 0
    reads_unmapped = 0
    read_cov_dict = {}

    seqs_parsed = 0
    progress_unit = 1
    old_format = 1
    with open(infile) as fh:
        for line in fh:
            if line.startswith("#"):
                old_format = 0
            if old_format == 0:
                if line.startswith('#'):
                    if line.startswith("## Total Reads"):
                        reads_total = int(line.split(" = ")[1])
                    elif line.startswith("## Mapped Reads"):
                        reads_mapped = int(line.split(" = ")[1])
                    elif line.startswith("## Unmapped Reads"):
                        reads_unmapped = int(line.split(" = ")[1])
                    else:
                        pass
                else:
                    match = cov_line_re.search(line)
                    if match:
                        seqs_parsed += 1
                        name, read_cov, base_cov = match.group(1), int(match.group(2)), float(match.group(3))
                        if name not in set_of_blobs:
                            print BtLog.warn_d['2'] % (name, infile)
                        else:
                            read_cov_dict[name] = read_cov
                            base_cov_dict[name] = base_cov
            else:
                match = old_cov_line_re.search(line)
                if match:
                    seqs_parsed += 1
                    name, base_cov = match.group(1), float(match.group(2))
                    if name not in set_of_blobs:
                        print BtLog.warn_d['2'] % (name, infile)
                    else:
                        base_cov_dict[name] = base_cov
            BtLog.progress(seqs_parsed, progress_unit, len(set_of_blobs))
        #BtLog.progress(len(set_of_blobs), progress_unit, len(set_of_blobs))
    return base_cov_dict, reads_total, reads_mapped, reads_unmapped, read_cov_dict

def checkCas(infile):
    print BtLog.status_d['12']
    if not isfile(infile):
         BtLog.error('0', infile)
    if not (which('clc_mapping_info')):
        BtLog.error('20')
    seqs_total_re = re.compile(r"\s+Contigs\s+(\d+)")
    reads_total_re = re.compile(r"\s+Reads\s+(\d+)")
    reads_mapping_re = re.compile(r"\s+Mapped reads\s+(\d+)\s+(\d+.\d+)\s+\%")
    seqs_total, reads_total, reads_mapping, mapping_rate = 0, 0, 0, 0.0
    output = ''
    command = "clc_mapping_info -s " + infile
    for line in runCmd(command):
        output += line
    seqs_total = int(seqs_total_re.search(output).group(1))
    reads_mapped = int(reads_mapping_re.search(output).group(1))
    reads_total = int(reads_total_re.search(output).group(1))
    print BtLog.status_d['11'] % ('{:,}'.format(reads_mapped), '{:,}'.format(reads_total), '{0:.1%}'.format(reads_mapped/reads_total))
    return seqs_total, reads_total, reads_mapped

def parseCas(infile, order_of_blobs):
    if not isfile(infile):
         BtLog.error('0', infile)
    seqs_total, reads_total, reads_mapped = checkCas(infile)
    progress_unit = int(len(order_of_blobs)/100)
    cas_line_re = re.compile(r"\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+.\d{2})\s+(\d+)\s+(\d+.\d{2})")
    command = "clc_mapping_info -n " + infile
    cov_dict = {}
    read_cov_dict = {}
    seqs_parsed = 0
    if (runCmd(command)):
        for line in runCmd(command):
            cas_line_match = cas_line_re.search(line)
            if cas_line_match:
                idx = int(cas_line_match.group(1)) - 1 # -1 because index of contig list starts with zero
                try:
                    name = order_of_blobs[idx]
                    reads = int(cas_line_match.group(3))
                    cov = float(cas_line_match.group(6))
                    cov_dict[name] = cov
                    read_cov_dict[name] = reads
                    seqs_parsed += 1
                except:
                    pass
                BtLog.progress(seqs_parsed, progress_unit, seqs_total)
    return cov_dict, reads_total, reads_mapped, read_cov_dict

def readTax(infile, set_of_blobs):
    '''
    If more fields need to be parsed:
        - change hit_line_re
        - catch matches in variables
        - add as key-value pairs to hitDict
    '''
    if not isfile(infile):
         BtLog.error('0', infile)
    hit_line_re = re.compile(r"^(\S+)\s+(\d+)[\;?\d+]*\s+(\d+\.*\d*)") # TEST TEST , if not split it afterwards
    with open(infile) as fh:
        for line in fh:
            match = hit_line_re.search(line)
            if match:
                hitDict = {
                    'name' : match.group(1),
                    'taxId' : match.group(2), # string because if int, conversion is a nightmare ...
                    'score' : float(match.group(3))
                    }
                if hitDict['name'] not in set_of_blobs:
                    BtLog.error('19', hitDict['name'], infile)
                if hitDict['taxId'] == 'N/A':
                    BtLog.error('22', infile)
                yield hitDict

def getOutFile(base_file, prefix, suffix):
    EXTENSIONS = ['.fasta', '.fa', '.fna', '.txt', '.cov', '.out', '.json']
    out_f, extension = splitext(basename(base_file))
    if not extension in EXTENSIONS:
        out_f = '%s%s' % (out_f, extension)
    if (prefix):
        if prefix.endswith("/"):
            out_f = "%s" % (join(prefix, out_f))
        else:
            out_f = "%s.%s" % (prefix, out_f)
    if (suffix):
        out_f = "%s.%s" % (out_f, suffix)
    return out_f

def parseNodesDB(**kwargs):
    '''
    Parsing names.dmp and nodes.dmp into the 'nodes_db' dict of dicts that
    gets JSON'ed into blobtools/data/nodes_db.json if this file
    does not exist. Nodes_db.json is used if neither "--names" and "--nodes"
    nor "--db" is specified.
    '''
    nodesDB = {}
    names_f = kwargs['names']
    nodes_f = kwargs['nodes']
    nodesDB_f = kwargs['nodesDB']
    nodesDB_default = kwargs['nodesDBdefault']

    if (nodes_f and names_f):
        if not isfile(names_f):
            BtLog.error('0', names_f)
        if not isfile(nodes_f):
            BtLog.error('0', nodes_f)
        print BtLog.status_d['3'] % (nodes_f, names_f)
        try:
            nodesDB = readNamesNodes(names_f, nodes_f)
        except:
            BtLog.error('3', nodes_f, names_f)
    elif (nodesDB_f):
        if not isfile(nodesDB_f):
            BtLog.error('0', nodesDB_f)
        print BtLog.status_d['4'] % (nodesDB_f)
        try:
            nodesDB = readNodesDB(nodesDB_f)
        except:
            BtLog.error('27', nodesDB_f)
    elif (nodesDB_default):
        if not isfile(nodesDB_default):
            BtLog.error('28')
        print BtLog.status_d['4'] % (nodesDB_default)
        try:
            nodesDB = readNodesDB(nodesDB_default)
        except:
            BtLog.error('27', nodesDB_default)
        nodesDB_f = nodesDB_default

    # Write nodesDB if not available
    if not isfile(nodesDB_default):
        writeNodesDB(nodesDB, nodesDB_default)

    return nodesDB, nodesDB_f

def readNamesNodes(names_f, nodes_f):
    nodesDB = {}
    nodes_count = 0
    with open(nodes_f) as fh:
        for line in fh:
            nodes_col = line.split("\t")
            node = {}
            node_id = nodes_col[0]
            node['parent'] = nodes_col[2]
            node['rank'] = nodes_col[4]
            nodesDB[node_id] = node
            nodes_count += 1
    with open(names_f) as fh:
        for line in fh:
            names_col = line.split("\t")
            if names_col[6] == "scientific name":
               nodesDB[names_col[0]]['name'] = names_col[2]
    nodesDB['nodes_count'] = nodes_count
    return nodesDB

def readNodesDB(nodesDB_f):
    nodesDB = {}
    nodesDB_count = 0
    nodes_count = 0
    with open(nodesDB_f) as fh:
        for line in fh:
            if line.startswith("#"):
                nodesDB_count = int(line.lstrip("# nodes_count = ").rstrip("\n"))
            nodes_count += 1
            node, rank, name, parent = line.rstrip("\n").split("\t")
            nodesDB[node] = {'rank' : rank, 'name' : name, 'parent' : parent}
            if (nodesDB_count):
                BtLog.progress(i, 1000, nodesDB_count)
    nodesDB['nodes_count'] = nodes_count
    return nodesDB

def writeNodesDB(nodesDB, nodesDB_f):
    print BtLog.status_d['5'] % nodesDB_f
    nodes_count = nodesDB['nodes_count']
    i = 0
    with open(nodesDB_f, 'w') as fh:
        fh.write("# nodes_count = %s\n" % nodes_count)
        for node in nodesDB:
            if not node == "nodes_count":
                i += 1
                BtLog.progress(i, 1000, nodes_count)
                fh.write("%s\t%s\t%s\t%s\n" % (node, nodesDB[node]['rank'], nodesDB[node]['name'], nodesDB[node]['parent']))

def byteify(input):
    '''
    http://stackoverflow.com/a/13105359
    '''
    if isinstance(input, dict):
        return {byteify(key):byteify(value) for key,value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input



def writeJsonGzip(obj, outfile):
    import json
    import gzip
    with gzip.open(outfile, 'wb') as fh:
        json.dump(obj, fh)

def writeJson(obj, outfile, indent=0, separators=(',', ': ')):
    import json
    with open(outfile, 'w') as fh:
        if (indent):
            json.dump(obj, fh, indent=indent, separators=separators)
        else:
            json.dump(obj, fh)
        #json.dump(obj, fh, indent=4, separators=(',', ': ')) #

def parseJsonGzip(infile):
    import json
    import gzip
    with gzip.open(infile, 'rb') as fh:
        obj = json.loads(fh.read().decode("ascii"))
    return byteify(obj)

def parseJson(infile):
    '''http://artem.krylysov.com/blog/2015/09/29/benchmark-python-json-libraries/'''
    if not isfile(infile):
         BtLog.error('0', infile)
    import time
    start = time.time()
    json_parser = ''
    with open(infile, 'r') as fh:
        print BtLog.status_d['15']
        json_string = fh.read()
    try:
        import ujson as json # fastest
        json_parser = 'ujson'
        print BtLog.status_d['16'] % json_parser
    except ImportError:
        try:
            import simplejson as json # fast
            json_parser = 'simplejson'
        except ImportError:
            import json # default
            json_parser = 'json'
        print BtLog.status_d['17'] % json_parser
    try:
        obj = json.loads(json_string.decode("ascii"))
    except ValueError:
        BtLog.error('37', infile, "BlobDB")
    data = byteify(obj)
    print BtLog.status_d['20'] % (time.time() - start)
    return data

if __name__ == "__main__":
    pass
