#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools taxify          [-b FILE] [-d FILE]
                                    [-u FILE] [-r FILE] [-s FILE] [-t INT]
                                    [-o PREFIX] [--force] [-h|--help]

    Options:
        -h --help                   show this

        -b, --blast <FILE>          BLAST similarity search result (BLAST TSV format)
                                      (-outfmt '6 qseqid staxids bitscore std')
        -d, --diamond <FILE>        Diamond similarity search result (Diamond TSV format)
                                      (--outfmt tab)

        -u, --uniref <FILE>         UnirefID-to-TaxID mapping file
        -r, --rnacentral <FILE>     RNAcentralID-to-TAXID mapping file
        -s, --swissprot <FILE>      UniprotID-to-TaxID mapping file
        -t, --taxid <INT>           TaxID (must be integer)

        --force                     Overwrite existing TaxIDs
        -o, --out <PREFIX>          Output prefix
"""

from __future__ import division
from docopt import docopt
from collections import defaultdict

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import bloblib.BtLog as BtLog
import bloblib.BtIO as BtIO

def main():
    args = docopt(__doc__)
    blast_f = args['--blast']
    diamond_f = args['--diamond']

    uniref_f = args['--uniref']
    rnacentral_f = args['--rnacentral']
    swissprot_f = args['--swissprot']
    taxid = args['--taxid']

    force = args['--force']
    prefix = args['--out']

    out_f, hit_f, map_f, taxid_d = None, None, None, {}

    # Check if blast_f OR diamond_f is speciefied
    if not (bool(blast_f) + bool(diamond_f) == 1):
        BtLog.error('26')
    elif blast_f:
        hit_f = blast_f
    elif diamond_f:
        hit_f = diamond_f
    else:
        pass

    # Check if taxID or Mapping file is supplied
    if (taxid):
        try:
            taxid = int(taxid)
        except TypeError:
            BtLog.error('26')
        out_f = BtIO.getOutFile(hit_f, prefix, "tax_%s.out" % taxid)
        taxid_d = defaultdict(lambda: taxid)
    elif (bool(uniref_f) + bool(rnacentral_f) + bool(swissprot_f) == 1):
        if uniref_f:
            print BtLog.status_d['1'] % ("ID-to-taxID Mapping file", uniref_f)
            taxid_d = BtIO.parseDict(uniref_f, 0, 1)
            out_f = BtIO.getOutFile(hit_f, prefix, "uniref.out")
            map_f = uniref_f
        elif rnacentral_f:
            print BtLog.status_d['1'] % ("ID-to-taxID Mapping file", rnacentral_f)
            taxid_d = BtIO.parseDict(rnacentral_f, 0, 3)
            out_f = BtIO.getOutFile(hit_f, prefix, "rnacentral.out")
            map_f = rnacentral_f
        elif swissprot_f:
            print BtLog.status_d['1'] % ("ID-to-taxID Mapping file", swissprot_f)
            taxid_d = BtIO.parseDict(swissprot_f, 0, 1)
            out_f = BtIO.getOutFile(hit_f, prefix, "swissprot.out")
            map_f = swissprot_f
        else:
            pass
    else:
        BtLog.error('41')

    output = []
    print BtLog.status_d['1'] % ("hits file", hit_f)

    with open(hit_f) as fh:
        for idx, l in enumerate(fh):
            query_id, bitscore, tax_id, subject_id, rest = None, None, None, None, None
            line = l.rstrip("\n").split()
            query_id = line[0]
            if blast_f:
                bitscore = line[2]
                tax_id = line[1]
                subject_id = line[3]
                rest = "\t".join(line[2:])
            elif diamond_f:
                bitscore = line[11]
                subject_id = line[1]
                rest = "\t".join(line[1:])
            if swissprot_f:
                subject_id = subject_id.split("|")[1]
            if blast_f and not tax_id == "N/A" and not force: # so that it does not overwrite existing taxIDs
                print BtLog.warn_d['10'] % (idx+1, line[0], line[1])
                output.append("%s\t%s\t%s\t%s" % (query_id, tax_id, bitscore, rest))
            else:
                try:
                    tax_id = taxid_d[subject_id]
                except KeyError:
                    print BtLog.warn_d['11'] % (subject_id, map_f)
                    tax_id = "N/A"
                output.append("%s\t%s\t%s\t%s" % (query_id, tax_id, bitscore, rest))

    if output:
        with open(out_f, "w") as fh:
            print BtLog.status_d['24'] % out_f
            fh.write("\n".join(output))
if __name__ == '__main__':
    main()



