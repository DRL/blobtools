#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools create     -i FASTA [-y FASTATYPE] [-o OUTFILE] [--title TITLE]
                              [-b BAM...] [-s SAM...] [-a CAS...] [-c COV...]
                              [--nodes <NODES>] [--names <NAMES>] [--db <NODESDB>]
                              [-t TAX...] [-x TAXRULE...]
                              [-h|--help]

    Options:
        -h --help                   show this
        -i, --infile FASTA          FASTA file of assembly. Headers are split at whitespaces.
        -y, --type FASTATYPE        Assembly program used to create FASTA. If specified,
                                    coverage will be parsed from FASTA header.
                                    (Parsing supported for 'spades', 'soap', 'velvet', 'abyss', 'platanus')
        -t, --taxfile TAX...        Taxonomy file in format (qseqid\\ttaxid\\tbitscore)
                                    (e.g. BLAST output "--outfmt '6 qseqid staxids bitscore'")
        -x, --taxrule <TAXRULE>...  Taxrule determines how taxonomy of blobs is computed [default: bestsum]
                                    "bestsum"       : sum bitscore across all hits for each taxonomic rank
                                    "bestsumorder"  : sum bitscore across all hits for each taxonomic rank.
                                                  - If first <TAX> file supplies hits, bestsum is calculated.
                                                  - If no hit is found, the next <TAX> file is used.
        --nodes <NODES>             NCBI nodes.dmp file. Not required if '--db'
        --names <NAMES>             NCBI names.dmp file. Not required if '--db'
        --db <NODESDB>              NodesDB file [default: data/nodesDB.txt].
        -b, --bam <BAM>...          BAM file(s) (requires samtools in $PATH)
        -s, --sam <SAM>...          SAM file(s)
        -a, --cas <CAS>...          CAS file(s) (requires clc_mapping_info in $PATH)
        -c, --cov <COV>...          TAB separated. (seqID\\tcoverage)
        -o, --out <OUT>             BlobDB output prefix
        --title TITLE               Title of BlobDB [default: output prefix)
"""

from __future__ import division
from docopt import docopt
import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtInput as BtInput
import os.path


if __name__ == '__main__':

    main_dir = os.path.dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    #print args

    title, fasta_f, fasta_type, cov_libs, hit_libs, taxrules, nodesDB_f, nodes_f, names_f, out_f = BtInput.validate_input_create(main_dir, args)

    # Create BlobDB object
    blobDb = bt.BlobDb(title)

    # Parse FASTA
    blobDb.parseFasta(fasta_f, fasta_type)

    # Parse Tax
    blobDb.parseHits(hit_libs)

    # Parse nodesDB
    nodesDB, nodesDB_f = BtIO.getNodesDB(nodes=nodes_f, names=names_f, nodesDB=nodesDB_f)
    blobDb.nodesDB_f = nodesDB_f

    if not os.path.isfile(nodesDB_f):
        print BtLog.status_d['5'] % nodesDB_f
        BtIO.writeNodesDB(nodesDB, nodesDB_f)

    # Computing taxonomy based on taxrules
    print BtLog.status_d['6'] % ",".join(taxrules)
    blobDb.computeTaxonomy(taxrules, nodesDB)

    # Parse coverage
    blobDb.parseCovs(cov_libs)

    # Generating BlobDB and writing to file
    print BtLog.status_d['7'] % out_f
    BtIO.writeJson(blobDb.dump(), out_f)
