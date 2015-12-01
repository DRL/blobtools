#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools bam2cov         -i FASTA -b BAM [-h|--help] 
    
    Options:
        -h --help                   show this
        -i, --infile FASTA          FASTA file of assembly. Headers are split at whitespaces.  
        -b, --bam <BAM>             BAM file (requires samtools in $PATH)
"""

from __future__ import division
import lib.BtLog as BtLog
from docopt import docopt
import re
import subprocess
import os

class Fasta():
    def __init__(self, name, seq):
        self.name = name 
        self.length = len(seq)
        self.n_count = seq.count('N')
        self.agct_count = self.length - self.n_count
        self.cov = 0.0

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

def runCmd(command):
    cmd = command.split() # sanitation
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')

def readFasta(infile):
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

def parseFasta(infile):
    fasta_dict = {}
    for name, seq in readFasta(infile):
        fasta = Fasta(name, seq)
        fasta_dict[fasta.name] = fasta
    return fasta_dict

def checkBam(infile):
    print BtLog.status_d['10']
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
    print BtLog.status_d['11'] % ('{:,}'.format(reads_mapped), '{:,}'.format(reads_total), '{0:.1%}'.format(reads_mapped/reads_total))
    return reads_total, reads_mapped

def readBam(infile, fasta_headers):
    reads_total, reads_mapped = checkBam(infile)
    progress_unit = int(int(reads_total)/1000)
    base_cov_dict = {}
    cigar_match_re = re.compile(r"(\d+)M") # only gets digits before M's
    # execute samtools to get only mapped reads
    command = "samtools view -F 4 " + infile
    # only one counter since only yields mapped reads
    parsed_reads = 0 
    for line in runCmd(command):
        match = line.split("\t")
        if match >= 11:
            seq_name = match[2]
            base_cov = sum([int(matching) for matching in cigar_match_re.findall(match[5])])
            if (base_cov):
                parsed_reads += 1
                if seq_name not in fasta_headers:
                    print BtLog.warn_d['2'] % (seq_name, infile)
                else:
                    base_cov_dict[seq_name] = base_cov_dict.get(seq_name, 0) + base_cov 
        BtLog.progress(parsed_reads, progress_unit, reads_total)
    BtLog.progress(reads_total, progress_unit, reads_total)

    if not int(reads_mapped) == int(parsed_reads):
        print warn_d['3'] % (reads_mapped, parsed_reads)
    return base_cov_dict, reads_total, parsed_reads

def parseBam(bam_f, fasta_dict):
    base_cov_dict, reads_total, reads_mapped = readBam(bam_f, set(fasta_dict.keys()))
    if reads_total == 0:
        print BtLog.warn_d['4'] % bam_f
    for name, base_cov in base_cov_dict.items():
        fasta_dict[name].cov = base_cov / fasta_dict[name].agct_count
    return fasta_dict

def writeCov(fasta_dict, out_f):
    with open(out_f, 'w') as fh:
        for name, fasta_obj in fasta_dict.items():
            fh.write("%s\t%s\n" % (name, fasta_obj.cov))

if __name__ == '__main__':
    args = docopt(__doc__)
    
    fasta_f = args['--infile']
    bam_f = args['--bam']
    out_f = os.path.basename(bam_f) + ".cov"

    fasta_dict = parseFasta(fasta_f)
    fasta_dict = parseBam(bam_f, fasta_dict)
    writeCov(fasta_dict, out_f)

