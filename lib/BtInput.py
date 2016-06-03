#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File        : BtInput.py
Version     : 0.1
Author      : Dominik R. Laetsch, dominik.laetsch at gmail dot com
Bugs        : ?
To do       : ?
"""

from __future__ import division
import re
import subprocess
from os.path import basename, isfile, abspath
import os
import lib.BtLog as BtLog
import lib.BtCore as bt

def validate_input_create(main_dir, args):
    '''
    Accepts:
        - main_dir
        - docopt args
    Returns:
        - title
        - fasta_f
        - fasta_type
        - cov_libs
        - hit_libs
        - nodesDB_f
        - taxrules
        - out_f
    '''
    ASSEMBLY_TYPES = [None, 'spades', 'soap', 'abyss', 'velvet', 'platanus']

    fasta_f = args['--infile']
    fasta_type = args['--type']
    sam_fs = args['--sam']
    bam_fs = args['--bam']
    cov_fs = args['--cov']
    cas_fs = args['--cas']
    hit_fs = args['--taxfile']
    out_f = args['--out']
    if (out_f):
        out_f = "%s.%s" % (os.path.basename(out_f), "BlobDB.json")
    else:
        out_f = "%s" % ("BlobDB.json")
    nodesDB_f = args['--db']
    names_f = args['--names']
    nodes_f = args['--nodes']
    taxrules = args['--taxrule']
    title = args['--title'] if (args['--title']) else out_f

    # Do files exist ?
    files = [x for x in list([fasta_f] + sam_fs + bam_fs + cov_fs + cas_fs + [names_f] + [nodes_f] + hit_fs) if x is not None]
    for f in files:
        if not os.path.isfile(f):
            BtLog.error('0', f)

    # Is taxonomy provided?
    if nodesDB_f == "data/nodesDB.txt":
        nodesDB_f = os.path.join(main_dir, nodesDB_f)
    if not os.path.isfile(nodesDB_f) and not ((names_f) and (nodes_f)):
        BtLog.error('3')
    if not (hit_fs):
        BtLog.error('18')
    # can FASTA parser deal with assemblies
    if not fasta_type in ASSEMBLY_TYPES:
        BtLog.error('2', ",".join(ASSEMBLY_TYPES[1:]))
    # Is coverage provided?
    if not (fasta_type) and not bam_fs and not sam_fs and not cov_fs and not cas_fs:
        BtLog.error('1')
    cov_libs = [bt.CovLibObj('bam' + str(idx), 'bam', lib_f) for idx, lib_f in enumerate(bam_fs)] + \
               [bt.CovLibObj('sam' + str(idx), 'sam', lib_f) for idx, lib_f in enumerate(sam_fs)] + \
               [bt.CovLibObj('cas' + str(idx), 'cas', lib_f) for idx, lib_f in enumerate(cas_fs)] + \
               [bt.CovLibObj('cov' + str(idx), 'cov', lib_f) for idx, lib_f in enumerate(cov_fs)]

    hit_libs = [bt.hitLibObj('tax' + str(idx), 'tax', lib_f) for idx, lib_f in enumerate(hit_fs)]

    return title, fasta_f, fasta_type, cov_libs, hit_libs, taxrules, nodesDB_f, nodes_f, names_f, out_f

if __name__ == "__main__":
    pass
