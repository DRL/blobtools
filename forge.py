#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools forge     --i <FASTA> [--type <ASSEMBLY>] [--out <OUT>] [--title <TITLE>]
                              [--bam <BAM>...] [--sam <SAM>...] [--cas <CAS>...] [--cov <COV>...]  
                              [--nodes <NODES>] [--names <NAMES>] [--db <NODESDB>] 
                              [--tax <TAX>...] [--taxrule <TAXRULE>...]
                              [--h|--help] 
    
    Options:
        --h --help              show this
        --i <FASTA>             FASTA file of assembly 
        --type ASSEMBLY         Assembly program used to create FASTA. If specified, 
                                coverage will be parsed from FASTA header. 
                                (Parsing supported for 'spades', 'soap', 'velvet', 'abyss')
        --tax <TAX>...          Taxonomy file in format (qseqid\\ttaxid\\tbitscore) 
                                (e.g. BLAST output "--outfmt '6 std'")
        --taxrule <TAXRULE>...  Taxrule determines how taxonomy of blobs is computed [default: bestsum]
                                "bestsum"       : sum bitscore across all hits for each taxonomic rank
                                "bestsumorder"  : sum bitscore across all hits for each taxonomic rank. 
                                                  - If first <TAX> file supplies hits these are used. 
                                                  - If no hit is found, the next <TAX> file is used.                                 
        --nodes <NODES>         NCBI nodes.dmp file. Not required if '--db'
        --names <NAMES>         NCBI names.dmp file. Not required if '--db' 
        --db <NODESDB>          NodesDB file [default: data/nodesDB.txt]. 
        --bam <BAM>...          BAM file (requires samtools in $PATH)
        --sam <SAM>...          SAM file
        --cas <CAS>...          CAS file (requires clc_mapping_info in $PATH)
        --cov <COV>...          TAB separated. (seqID\\tcoverage)
        --out <OUT>             BlobDB output file [default: blobDb.json]
        --title TITLE           Title of BlobDB [default: FASTA)  
"""

from __future__ import division
from docopt import docopt
import lib.BtCore as bt
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import os.path


if __name__ == '__main__':
    ASSEMBLY_TYPES = [None, 'spades', 'soap', 'abyss', 'velvet']
    main_dir = os.path.dirname(__file__)
    #print data_dir
    args = docopt(__doc__)
    #print args

    fasta_f = args['--i']
    fasta_type = args['--type']
    
    sam_fs = args['--sam']
    bam_fs = args['--bam']
    cov_fs = args['--cov']
    cas_fs = args['--cas']
    hit_fs = args['--tax']

    out_f = args['--out']
    nodesDB_f = args['--db']
    names_f = args['--names']
    nodes_f = args['--nodes']
    taxrules = args['--taxrule']
    title = args['--title'] if (args['--title']) else os.path.basename(".".join(fasta_f.split('.')[0:-1]))


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

    # Create BlobDB object              
    blobDb = bt.BlobDb(title)

    # Parse FASTA
    blobDb.parseFasta(fasta_f, fasta_type)
    # Parse coverage
    blobDb.parseCovs(cov_libs)

    # Parse Tax
    hitLibs = [bt.hitLibObj('tax' + str(idx), 'tax', lib_f) for idx, lib_f in enumerate(hit_fs)]
    blobDb.parseHits(hitLibs)
    
    # Parse nodesDB
    nodesDB, nodesDB_f = BtIO.getNodesDB(nodes=nodes_f, names=names_f, nodesDB=nodesDB_f)
    blobDb.nodesDB_f = nodesDB_f
        
    if not os.path.isfile(nodesDB_f):
        print BtLog.status_d['5'] % nodesDB_f
        BtIO.writeNodesDB(nodesDB, nodesDB_f)

    # Computing taxonomy based on taxrules
    print BtLog.status_d['6'] % ",".join(taxrules)
    blobDb.computeTaxonomy(taxrules, nodesDB)

    # Writing BlobDB to file
    print BtLog.status_d['7'] % out_f
    BtIO.writeJson(blobDb.dump(), out_f)