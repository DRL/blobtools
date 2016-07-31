#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools create     -i FASTA [-y FASTATYPE] [-o PREFIX] [--title TITLE]
                              [-b BAM...] [-s SAM...] [-a CAS...] [-c COV...]
                              [--nodes <NODES>] [--names <NAMES>] [--db <NODESDB>]
                              [-t HITS...] [-x TAXRULE...] [-m INT] [--tax_collision_random]
                              [-h|--help]

    Options:
        -h --help                       show this
        -i, --infile FASTA              FASTA file of assembly. Headers are split at whitespaces.
        -y, --type FASTATYPE            Assembly program used to create FASTA. If specified,
                                        coverage will be parsed from FASTA header.
                                        (Parsing supported for 'spades', 'velvet', 'platanus')
        -t, --hitsfile HITS...          Hits file in format (qseqid\\ttaxid\\tbitscore)
                                        (e.g. BLAST output "--outfmt '6 qseqid staxids bitscore'")
        -x, --taxrule <TAXRULE>...      Taxrule determines how taxonomy of blobs
                                        is computed [default: bestsum]
                                        "bestsum"       : sum bitscore across all
                                                          hits for each taxonomic rank
                                        "bestsumorder"  : sum bitscore across all
                                                          hits for each taxonomic rank.
                                                  - If first <TAX> file supplies hits, bestsum is calculated.
                                                  - If no hit is found, the next <TAX> file is used.
        -m, --min_diff <FLOAT>          Minimal score difference between highest scoring
                                        taxonomies (otherwise "unresolved") [default: 0.0]
        --tax_collision_random          Random allocation of taxonomy if highest scoring
                                        taxonomies have equal scores (otherwise "unresolved")
        --nodes <NODES>                 NCBI nodes.dmp file. Not required if '--db'
        --names <NAMES>                 NCBI names.dmp file. Not required if '--db'
        --db <NODESDB>                  NodesDB file (default: $BLOBTOOLS/data/nodesDB.txt).
        -b, --bam <BAM>...              BAM file(s) (requires samtools in $PATH)
        -s, --sam <SAM>...              SAM file(s)
        -a, --cas <CAS>...              CAS file(s) (requires clc_mapping_info in $PATH)
        -c, --cov <COV>...              COV file(s)
        -o, --out <PREFIX>              BlobDB output prefix
        --title TITLE                   Title of BlobDB [default: output prefix)
"""

from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import blobtools
import bloblib.BtCore as BtCore
import bloblib.BtLog as BtLog
import bloblib.BtIO as BtIO

def main():

    #main_dir = dirname(__file__)
    args = docopt(__doc__)
    fasta_f = args['--infile']
    fasta_type = args['--type']
    sam_fs = args['--sam']
    bam_fs = args['--bam']
    cov_fs = args['--cov']
    cas_fs = args['--cas']
    hit_fs = args['--hitsfile']
    prefix = args['--out']
    nodesDB_f = args['--db']
    names_f = args['--names']
    nodes_f = args['--nodes']
    taxrules = args['--taxrule']
    min_bitscore_diff = float(args['--min_diff'])
    tax_collision_random = args['--tax_collision_random']
    title = args['--title']

    # outfile
    out_f = BtIO.getOutFile("blobDB", prefix, "json")
    if not (title):
        title = out_f

    # coverage
    if not (fasta_type) and not bam_fs and not sam_fs and not cov_fs and not cas_fs:
        BtLog.error('1')
    cov_libs = [BtCore.CovLibObj('bam' + str(idx), 'bam', lib_f) for idx, lib_f in enumerate(bam_fs)] + \
           [BtCore.CovLibObj('sam' + str(idx), 'sam', lib_f) for idx, lib_f in enumerate(sam_fs)] + \
           [BtCore.CovLibObj('cas' + str(idx), 'cas', lib_f) for idx, lib_f in enumerate(cas_fs)] + \
           [BtCore.CovLibObj('cov' + str(idx), 'cov', lib_f) for idx, lib_f in enumerate(cov_fs)]

    # taxonomy
    hit_libs = [BtCore.HitLibObj('tax' + str(idx), 'tax', lib_f) for idx, lib_f in enumerate(hit_fs)]

    # Create BlobDB object
    blobDb = BtCore.BlobDb(title)
    blobDb.version = blobtools.__version__
    # Parse FASTA
    blobDb.parseFasta(fasta_f, fasta_type)

    # Parse nodesDB OR names.dmp, nodes.dmp
    nodesDB_default = join(blobtools.DATADIR, "nodesDB.txt")
    nodesDB, nodesDB_f = BtIO.parseNodesDB(nodes=nodes_f, names=names_f, nodesDB=nodesDB_f, nodesDBdefault=nodesDB_default)
    blobDb.nodesDB_f = nodesDB_f

    # Parse similarity hits
    if (hit_libs):
        blobDb.parseHits(hit_libs)
        blobDb.computeTaxonomy(taxrules, nodesDB, min_bitscore_diff, tax_collision_random)
    else:
        print BtLog.warn_d['0']

    # Parse coverage
    blobDb.parseCoverage(covLibObjs=cov_libs, no_base_cov=None)

    # Generating BlobDB and writing to file
    print BtLog.status_d['7'] % out_f
    BtIO.writeJson(blobDb.dump(), out_f)

if __name__ == '__main__':
    main()
