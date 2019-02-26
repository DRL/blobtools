#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: blobtools taxify          -f FILE [-a INT] [-b INT] [-c INT]
                                    [-m FILE] [-s INT] [-t INT]
                                    [-i FILE] [-x INT] [-v FLOAT]
                                    [-o PREFIX] [-h|--help]

    Options:
        -h --help                           show this

    Options for similarity search input
        -f, --hit_file <FILE>               BLAST/Diamond similarity search result (TSV format).
                                                Defaults assume "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
        -a, --hit_column_qseqid <INT>       Zero-based column of qseqid in similarity search result [default: 0]
                                                Change if different format than (-outfmt '6')
        -b, --hit_column_sseqid <INT>       Zero-based column of sseqid in similarity search result [default: 1]
                                                Change if different format than (-outfmt '6')
        -c, --hit_column_score <INT>        Zero-based column of (bit)score in similarity search result [default: 11]
                                                Change if different format than (-outfmt '6')
    Options for TaxID mapping file
        -m, --taxid_mapping_file <FILE>     TaxID mapping file (contains seqid and taxid)
        -s, --map_col_sseqid <INT>          Zero-based column of sseqid in TaxID mapping file (it will search for sseqid in this column)
        -t, --map_col_taxid <INT>           Zero-based Column of taxid in TaxID mapping file (it will extract for taxid from this column)

    Options for custom input
        -i, --custom <FILE>                 File containing list of sequence IDs
        -x, --custom_taxid <INT>            TaxID to assign to all sequence IDs in list
        -v, --custom_score <FLOAT>          Score to assign to all sequence IDs in list

    General
        -o, --out <PREFIX>                  Output prefix
"""

from __future__ import division
from docopt import docopt
from collections import defaultdict

import lib.BtLog as BtLog
import lib.BtIO as BtIO

def main():
    args = docopt(__doc__)
    out_f, hit_f, map_f, taxid_d = None, None, None, {}
    hit_f = args['--hit_file']
    hit_col_qseqid = args['--hit_column_qseqid']
    hit_col_sseqid = args['--hit_column_sseqid']
    hit_col_score = args['--hit_column_score']
    map_f = args['--taxid_mapping_file']
    map_col_sseqid = args['--map_col_sseqid']
    map_col_taxid = args['--map_col_taxid']
    #custom_f = args['--custom']
    custom_taxid = args['--custom_taxid']
    #custom_score = args['--custom_score']
    prefix = args['--out']

    try:
        hit_col_qseqid = int(hit_col_qseqid)
        hit_col_sseqid = int(hit_col_sseqid)
        hit_col_score = int(hit_col_score)
    except ValueError:
        BtLog.error('41' % ("--hit_column_qseqid, --hit_column_sseqid and --hit_column_score"))

    if custom_taxid:
        try:
            custom_taxid = int(custom_taxid)
        except TypeError:
            BtLog.error('26')
        out_f = BtIO.getOutFile(hit_f, prefix, "taxID_%s.out" % custom_taxid)
        taxid_d = defaultdict(lambda: custom_taxid)
    elif map_f:
        if map_col_sseqid and map_col_taxid:
            try:
                map_col_sseqid = int(map_col_sseqid)
                map_col_taxid = int(map_col_taxid)
            except ValueError:
                BtLog.error('44')
            print(BtLog.status_d['1'] % ("Mapping file", map_f))
            taxid_d = BtIO.parseDict(map_f, map_col_sseqid, map_col_taxid)
            out_f = BtIO.getOutFile(hit_f, prefix, "taxified.out")
        else:
            BtLog.error('44')
    else:
        BtLog.error('41')

    output = []
    print(BtLog.status_d['1'] % ("similarity search result", hit_f))
    with open(hit_f) as fh:
        for idx, line in enumerate(fh):
            col = line.rstrip("\n").split()
            qseqid = col[hit_col_qseqid]
            sseqid = col[hit_col_sseqid]
            score = col[hit_col_score]
            tax_id = None
            if custom_taxid:
                tax_id = taxid_d[sseqid]
            else:
                if sseqid not in taxid_d:
                    BtLog.warn_d['12'] % (sseqid, map_f)
                tax_id = taxid_d.get(sseqid, "N/A")
            output.append("%s\t%s\t%s\t%s" % (qseqid, tax_id, score, sseqid))
    if output:
        with open(out_f, "w") as fh:
            print(BtLog.status_d['24'] % out_f)
            fh.write("\n".join(output) + "\n")

if __name__ == '__main__':
    main()
