#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File        : BtCore.py
Author      : Dominik R. Laetsch, dominik.laetsch at gmail dot com
"""

from __future__ import division
import bloblib.BtLog as BtLog
import bloblib.BtIO as BtIO
import bloblib.BtTax as BtTax
from os.path import abspath, isfile, basename, isdir, join
from os import getcwd, mkdir
import json

class BlobDb():
    '''
    class BlobDB holds all information parsed from files
    '''
    def __init__(self, title):
        self.title = title
        self.assembly_f = ''
        self.dict_of_blobs = {}
        self.order_of_blobs = [] # ordereddict
        self.set_of_taxIds = set()
        self.lineages = {}
        self.length = 0
        self.seqs = 0
        self.n_count = 0
        self.covLibs = {}
        self.hitLibs = {}
        self.nodesDB_f = ''
        self.taxrules = []
        self.version = ''
        self.view_dir = ''

    def view(self, **kwargs):
        # arguments
        viewObjs = kwargs['viewObjs']
        ranks = kwargs['ranks']
        taxrule = kwargs['taxrule']
        hits_flag = kwargs['hits_flag']
        seqs = kwargs['seqs']
        cov_libs = kwargs['cov_libs']
        progress_bar = kwargs['progressbar']
        # Default sequences if no subset
        if not (seqs):
            seqs = self.order_of_blobs
        # Default cov_libs if no subset
        cov_lib_names = cov_libs
        if not (cov_libs):
            cov_lib_names = [covLib for covLib in self.covLibs]
        tax_lib_names = [taxLib for taxLib in sorted(self.hitLibs)]
        lineages = self.lineages
        # setup
        for viewObj in viewObjs:
            if viewObj.name == 'table':
                viewObj.header = self.getTableHeader(taxrule, ranks, hits_flag, cov_lib_names)
            if viewObj.name == 'concoct_cov':
                viewObj.header = self.getConcoctCovHeader(cov_lib_names)
            if viewObj.name == 'covlib':
                viewObj.header = self.getCovHeader(cov_lib_names)
            if viewObj.name == 'experimental':
                viewObj.covs = {cov_lib : [] for cov_lib in cov_lib_names}
                viewObj.covs["covsum"] = []
                for taxrule in self.taxrules:
                    viewObj.tax[taxrule] = {rank : [] for rank in BtTax.RANKS}
        # bodies
        for i, seq in enumerate(seqs):
            if (progress_bar):
                BtLog.progress(i, 1000, len(seqs))
            blob = self.dict_of_blobs[seq]
            for viewObj in viewObjs:
                if viewObj.name == 'table':
                    viewObj.body.append(self.getTableLine(blob, taxrule, ranks, hits_flag, cov_lib_names, tax_lib_names, lineages))
                if viewObj.name == 'concoct_cov':
                    viewObj.body.append(self.getConcoctCovLine(blob, cov_lib_names))
                if viewObj.name == 'experimental':
                    viewObj.names.append(blob['name'])
                    viewObj.gc.append(blob['gc'])
                    viewObj.length.append(blob['length'])
                    cov_sum = 0.0
                    for cov_lib in blob['covs']:
                        viewObj.covs[cov_lib].append(blob['covs'][cov_lib])
                        cov_sum += blob['covs'][cov_lib]
                    viewObj.covs['covsum'].append(cov_sum)
                    for taxrule in blob['taxonomy']:
                        for rank in blob['taxonomy'][taxrule]:
                            viewObj.tax[taxrule][rank].append(blob['taxonomy'][taxrule][rank]['tax'])
                if viewObj.name == 'concoct_tax':
                    for rank in ranks:
                        if not rank in viewObj.body:
                            viewObj.body[rank] = []
                        viewObj.body[rank].append(self.getConcoctTaxLine(blob, rank, taxrule))
                if viewObj.name == 'covlib':
                    viewObj.body.append(self.getCovLine(blob, cov_lib_names))
        if (progress_bar):
            BtLog.progress(len(seqs), 1000, len(seqs))
        for viewObj in viewObjs:
            viewObj.output()

    def getCovHeader(self, cov_lib_names):
        cov_lib_name = cov_lib_names[0]
        header = '## %s\n' % (self.version)
        if isinstance(self.covLibs[cov_lib_name], CovLibObj):
            # CovLibObjs
            header += "## Total Reads = %s\n" % (self.covLibs[cov_lib_name].reads_total)
            header += "## Mapped Reads = %s\n" % (self.covLibs[cov_lib_name].reads_mapped)
            header += "## Unmapped Reads = %s\n" % (self.covLibs[cov_lib_name].reads_total - self.covLibs[cov_lib_name].reads_mapped)
            header += "## Source(s) : %s\n" % (self.covLibs[cov_lib_name].f)
        else:
            # CovLibObjs turned dicts
            header += "## Total Reads = %s\n" % (self.covLibs[cov_lib_name]['reads_total'])
            header += "## Mapped Reads = %s\n" % (self.covLibs[cov_lib_name]['reads_mapped'])
            header += "## Unmapped Reads = %s\n" % (self.covLibs[cov_lib_name]['reads_total'] - self.covLibs[cov_lib_name]['reads_mapped'])
            header += "## Source(s) : %s\n" % (self.covLibs[cov_lib_name]['f'])
        header += "# %s\t%s\t%s\n" % ("contig_id", "read_cov", "base_cov")
        return header

    def getCovLine(self, blob, cov_lib_names):
        if isinstance(blob, BlObj):
            # BlobObj
            return "%s\t%s\t%s\n" % (blob.name, blob.read_cov[cov_lib_names[0]], blob.covs[cov_lib_names[0]])
        else:
            # BlobObj turned dict
            return "%s\t%s\t%s\n" % (blob['name'], blob['read_cov'][cov_lib_names[0]], blob['covs'][cov_lib_names[0]])

    def getConcoctCovHeader(self, cov_lib_names):
        return "contig\t%s\n" % "\t".join(cov_lib_names)

    def getConcoctTaxLine(self, blob, rank, taxrule):
        if taxrule in blob['taxonomy']:
            return "%s,%s\n" % (blob['name'], blob['taxonomy'][taxrule][rank]['tax'])

    def getConcoctCovLine(self, blob, cov_lib_names):
        return "%s\t%s\n" % (blob['name'], "\t".join(map(str, [ blob['covs'][covLib] for covLib in cov_lib_names])))

    def getTableHeader(self, taxrule, ranks, hits_flag, cov_lib_names):
        sep = "\t"
        header = ''
        header += '## %s\n' % (self.version)
        header += "## assembly\t: %s\n" % self.assembly_f
        header += "%s\n" % "\n".join("## coverage\t: " + covLib['name'] + " - " + covLib["f"] for covLib in self.covLibs.values())
        if (self.hitLibs):
            header += "%s\n" % "\n".join("## taxonomy\t: " + hitLib['name'] + " - " + hitLib["f"] for hitLib in self.hitLibs.values())
        else:
            header += "## taxonomy\t: no taxonomy information found\n"
        header += "## nodesDB\t: %s\n" % self.nodesDB_f
        header += "## taxrule\t: %s\n" % taxrule
        header += "##\n"
        header += "# %s" % sep.join(map(str, ["name", "length", "GC", "N"]))
        col = 4
        header += "%s%s" % (sep, sep.join([cov_lib_name for cov_lib_name in cov_lib_names]))
        col += len(cov_lib_names)
        if (len(cov_lib_names)) > 1:
            col += 1
            header += "%s%s" % (sep, "cov_sum")
        for rank in ranks:
            header += "%s%s" % (sep, sep.join([rank + ".t." + str(col + 1), rank + ".s." + str(col + 2), rank + ".c." + str(col + 3)]))
            col += 3
            if hits_flag:
                header += "%s%s" % (sep, rank + ".hits." + str(col + 1))
                col += 1
        return header

    def getTableLine(self, blob, taxrule, ranks, hits_flag, cov_lib_names, tax_lib_names, lineages):
        sep = "\t"
        line = ''
        line += "\n%s" % sep.join(map(str, [ blob['name'], blob['length'], blob['gc'], blob['n_count']  ]))
        line += sep
        line += "%s" % sep.join(map(str, [ blob['covs'][covLib] for covLib in cov_lib_names]))
        if len(cov_lib_names) > 1:
            line += "%s%s" % (sep, sum([ blob['covs'][covLib] for covLib in cov_lib_names]))
        for rank in ranks:
            line += sep
            tax, score, c_index = 'N/A', 'N/A', 'N/A'
            if taxrule in blob['taxonomy']:
                tax = blob['taxonomy'][taxrule][rank]['tax']
                score = blob['taxonomy'][taxrule][rank]['score']
                c_index = blob['taxonomy'][taxrule][rank]['c_index']
            line += "%s" % sep.join(map(str, [tax, score, c_index ]))
            if (hits_flag) and (taxrule in blob['taxonomy']):
                line += sep
                for tax_lib_name in tax_lib_names:
                    if tax_lib_name in blob['hits']:
                        line += "%s=" % tax_lib_name
                        sum_dict = {}
                        for hit in blob['hits'][tax_lib_name]:
                            tax_rank = lineages[hit['taxId']][rank]
                            sum_dict[tax_rank] = sum_dict.get(tax_rank, 0.0) + hit['score']
                        line += "%s" % "|".join([":".join(map(str, [tax_rank, sum_dict[tax_rank]])) for tax_rank in sorted(sum_dict, key=sum_dict.get, reverse=True)])
                    else:
                        line += "%s=no-hit:0.0" % (tax_lib_name)
                    line += ";"
        return line

    def dump(self):
        dump = {'title' : self.title,
                'assembly_f' : self.assembly_f,
                'lineages' : self.lineages,
                'order_of_blobs' : self.order_of_blobs,
                'dict_of_blobs' : {name : blObj.__dict__ for name, blObj in self.dict_of_blobs.items()},
                'length' : self.length,
                'seqs' : self.seqs,
                'n_count' : self.n_count,
                'nodesDB_f' : self.nodesDB_f,
                'covLibs' : {name : covLibObj.__dict__ for name, covLibObj in self.covLibs.items()},
                'hitLibs' : {name : hitLibObj.__dict__ for name, hitLibObj in self.hitLibs.items()},
                'taxrules' : self.taxrules,
                'version' : self.version
                }
        return dump

    def load(self, BlobDb_f):
        blobDict = BtIO.parseJson(BlobDb_f)
        for k, v in blobDict.items():
            setattr(self, k, v)
        self.set_of_taxIds = blobDict['lineages'].keys()

    def getPlotData(self, rank, min_length, hide_nohits, taxrule, c_index, catcolour_dict):
        data_dict = {}
        read_cov_dict = {}
        max_cov = 0.0

        cov_lib_dict = self.covLibs
        cov_lib_names_l = self.covLibs.keys() # does not include cov_sum
        if len(cov_lib_names_l) > 1:
            # more than one cov_lib, cov_sum_lib has to be created
            cov_lib_dict['covsum'] = CovLibObj('covsum', 'covsum', 'Sum of cov in %s' % basename(self.title)).__dict__ # ugly
            cov_lib_dict['covsum']['reads_total'] = sum([self.covLibs[x]['reads_total'] for x in self.covLibs])
            cov_lib_dict['covsum']['reads_mapped'] = sum([self.covLibs[x]['reads_mapped'] for x in self.covLibs])
            cov_lib_dict['covsum']['cov_sum'] = sum([self.covLibs[x]['cov_sum'] for x in self.covLibs])
            cov_lib_dict['covsum']['mean_cov'] = cov_lib_dict['covsum']['cov_sum']/self.seqs
        for blob in self.dict_of_blobs.values():
            name, gc, length, group = blob['name'], blob['gc'], blob['length'], ''
            if (catcolour_dict): # annotation with categories specified in catcolour
                group = str(catcolour_dict[name])
            elif (c_index): # annotation with c_index instead of taxonomic group
                if taxrule not in self.taxrules:
                    BtLog.error('11', taxrule, self.taxrules)
                else:
                    group = str(blob['taxonomy'][taxrule][rank]['c_index'])
            else: # annotation with taxonomic group
                if not (taxrule) or taxrule not in self.taxrules:
                    BtLog.warn_d['9'] % (taxrule, self.taxrules)
                if taxrule in blob['taxonomy']:
                    group = str(blob['taxonomy'][taxrule][rank]['tax'])
            if not group in data_dict:
                data_dict[group] = {
                                    'name' : [],
                                    'length' : [],
                                    'gc' : [],
                                    'covs' : {covLib : [] for covLib in cov_lib_dict.keys()},           # includes cov_sum if it exists
                                    'reads_mapped' : {covLib : 0 for covLib in cov_lib_dict.keys()},    # includes cov_sum if it exists
                                    'count' : 0,
                                    'count_hidden' : 0,
                                    'count_visible' : 0,
                                    'span': 0,
                                    'span_hidden' : 0,
                                    'span_visible' : 0,
                                    }
            data_dict[group]['count'] = data_dict[group].get('count', 0) + 1
            data_dict[group]['span'] = data_dict[group].get('span', 0) + int(length)
            if ((hide_nohits) and group == 'no-hit') or length < min_length: # hidden
                data_dict[group]['count_hidden'] = data_dict[group].get('count_hidden', 0) + 1
                data_dict[group]['span_hidden'] = data_dict[group].get('span_hidden', 0) + int(length)
            else: # visible
                data_dict[group]['count_visible'] = data_dict[group].get('count_visible', 0) + 1
                data_dict[group]['span_visible'] = data_dict[group].get('span_visible', 0) + int(length)
                data_dict[group]['name'].append(name)
                data_dict[group]['length'].append(length)
                data_dict[group]['gc'].append(gc)
                cov_sum = 0.0
                reads_mapped_sum = 0
                for cov_lib in sorted(cov_lib_names_l):
                    cov = float(blob['covs'][cov_lib])
                    if cov < 0.02:
                        cov = 0.02
                    # increase max_cov
                    if cov > max_cov:
                        max_cov = cov
                    # add cov of blob to group
                    data_dict[group]['covs'][cov_lib].append(cov)
                    cov_sum += cov
                    # add readcov
                    if cov_lib in blob['read_cov']:
                        reads_mapped = blob['read_cov'][cov_lib]
                        data_dict[group]['reads_mapped'][cov_lib] += reads_mapped
                        reads_mapped_sum += reads_mapped
                if len(cov_lib_names_l) > 1:
                    if cov_sum < 0.02 :
                        cov_sum = 0.02
                    data_dict[group]['covs']['covsum'].append(cov_sum)
                    if cov_sum > max_cov:
                        max_cov = cov_sum
                    if (reads_mapped_sum):
                        data_dict[group]['reads_mapped']['covsum'] += reads_mapped_sum

        return data_dict, max_cov, cov_lib_dict

    def addCovLib(self, covLib):
        self.covLibs[covLib.name] = covLib
        for blObj in self.dict_of_blobs.values():
            blObj.addCov(covLib.name, 0.0)

    def parseFasta(self, fasta_f, fasta_type):
        print BtLog.status_d['1'] % ('FASTA', fasta_f)
        self.assembly_f = abspath(fasta_f)
        if (fasta_type):
            # Set up CovLibObj for coverage in assembly header
            self.covLibs[fasta_type] = CovLibObj(fasta_type, fasta_type, fasta_f)

        for name, seq in BtIO.readFasta(fasta_f):
            blObj = BlObj(name, seq)
            if not blObj.name in self.dict_of_blobs:
                self.seqs += 1
                self.length += blObj.length
                self.n_count += blObj.n_count

                if (fasta_type):
                    cov = BtIO.parseCovFromHeader(fasta_type, blObj.name)
                    self.covLibs[fasta_type].cov_sum += cov
                    blObj.addCov(fasta_type, cov)

                self.order_of_blobs.append(blObj.name)
                self.dict_of_blobs[blObj.name] = blObj
            else:
                BtLog.error('5', blObj.name)

        if self.seqs == 0 or self.length == 0:
            BtLog.error('1')

    def parseCoverage(self, **kwargs):
        # arguments
        covLibObjs = kwargs['covLibObjs']
        no_base_cov = kwargs['no_base_cov']

        for covLib in covLibObjs:
            self.addCovLib(covLib)
            print BtLog.status_d['1'] % (covLib.name, covLib.f)
            if covLib.fmt == 'bam' or covLib.fmt == 'sam':
                base_cov_dict = {}
                if covLib.fmt == 'bam':
                    base_cov_dict, covLib.reads_total, covLib.reads_mapped, read_cov_dict = BtIO.parseBam(covLib.f, set(self.dict_of_blobs), no_base_cov)
                else:
                    base_cov_dict, covLib.reads_total, covLib.reads_mapped, read_cov_dict = BtIO.parseSam(covLib.f, set(self.dict_of_blobs), no_base_cov)

                if covLib.reads_total == 0:
                    print BtLog.warn_d['4'] % covLib.f

                for name, base_cov in base_cov_dict.items():
                    cov = base_cov / self.dict_of_blobs[name].agct_count
                    covLib.cov_sum += cov
                    self.dict_of_blobs[name].addCov(covLib.name, cov)
                    self.dict_of_blobs[name].addReadCov(covLib.name, read_cov_dict[name])
                # Create COV file for future use
                out_f = BtIO.getOutFile(covLib.f, None, None)
                covView = ViewObj(name="covlib", out_f=out_f, suffix="cov", header="", body=[])
                self.view(viewObjs=[covView], ranks=None, taxrule=None, hits_flag=None, seqs=None, cov_libs=[covLib.name], progressbar=False)

            elif covLib.fmt == 'cas':
                cov_dict, covLib.reads_total, covLib.reads_mapped, read_cov_dict = BtIO.parseCas(covLib.f, self.order_of_blobs)
                if covLib.reads_total == 0:
                    print BtLog.warn_d['4'] % covLib.f
                for name, cov in cov_dict.items():
                    covLib.cov_sum += cov
                    self.dict_of_blobs[name].addCov(covLib.name, cov)
                    self.dict_of_blobs[name].addReadCov(covLib.name, read_cov_dict[name])
                out_f = BtIO.getOutFile(covLib.f, None, None)
                covView = ViewObj(name="covlib", out_f=out_f, suffix="cov", header="", body=[])
                self.view(viewObjs=[covView], ranks=None, taxrule=None, hits_flag=None, seqs=None, cov_libs=[covLib.name], progressbar=False)

            elif covLib.fmt == 'cov':
                base_cov_dict, covLib.reads_total, covLib.reads_mapped, covLib.reads_unmapped, read_cov_dict = BtIO.parseCov(covLib.f, set(self.dict_of_blobs))
                #cov_dict = BtIO.readCov(covLib.f, set(self.dict_of_blobs))
                if not len(base_cov_dict) == self.seqs:
                    print BtLog.warn_d['4'] % covLib.f
                for name, cov in base_cov_dict.items():
                    covLib.cov_sum += cov
                    self.dict_of_blobs[name].addCov(covLib.name, cov)
                    if name in read_cov_dict:
                        self.dict_of_blobs[name].addReadCov(covLib.name, read_cov_dict[name])
            else:
                pass
            covLib.mean_cov = covLib.cov_sum/self.seqs
            if covLib.cov_sum == 0.0:
                print BtLog.warn_d['6'] % (covLib.name)
            self.covLibs[covLib.name] = covLib


    def parseHits(self, hitLibs):
        for hitLib in hitLibs:
            self.hitLibs[hitLib.name] = hitLib
            print BtLog.status_d['1'] % (hitLib.name, hitLib.f)
            # only accepts format 'seqID\ttaxID\tscore'
            for hitDict in BtIO.readTax(hitLib.f, set(self.dict_of_blobs)):
                if ";" in hitDict['taxId']:
                    hitDict['taxId'] = hitDict['taxId'].split(";")[0]
                    print BtLog.warn_d['5'] % (hitDict['name'], hitLib)
                self.set_of_taxIds.add(hitDict['taxId'])
                self.dict_of_blobs[hitDict['name']].addHits(hitLib.name, hitDict)

    def computeTaxonomy(self, taxrules, nodesDB, min_bitscore_diff, tax_collision_random):
        print BtLog.status_d['6'] % ",".join(taxrules)
        tree_lists = BtTax.getTreeList(self.set_of_taxIds, nodesDB)
        self.lineages = BtTax.getLineages(tree_lists, nodesDB)
        self.taxrules = taxrules
        i = 0
        for blObj in self.dict_of_blobs.values():
            i += 1
            BtLog.progress(i, 100, self.seqs)
            for taxrule in taxrules:
                if (blObj.hits):
                    blObj.taxonomy[taxrule] = BtTax.taxRule(taxrule, blObj.hits, self.lineages, min_bitscore_diff, tax_collision_random)
                else:
                    blObj.taxonomy[taxrule] = BtTax.noHit()
        self.set_of_taxIds = set()

    def getBlobs(self):
        for blObj in [self.dict_of_blobs[key] for key in self.order_of_blobs]:
            yield blObj

class BlObj():
    def __init__(self, name, seq):
        self.name = name
        self.length = len(seq)
        self.n_count = seq.count('N')
        self.agct_count = self.length - self.n_count
        self.gc = round(self.calculateGC(seq), 4)
        self.covs = {}
        self.read_cov = {}
        self.hits = {}
        self.taxonomy = {}

    def calculateGC(self, seq):
        return float((seq.count('G') + seq.count('C') ) / self.agct_count \
                     if self.agct_count > 0 else 0.0)

    def addCov(self, lib_name, cov):
        self.covs[lib_name] = float("{0:.3f}".format(cov)) # changed to three decimal digits

    def addReadCov(self, lib_name, read_cov):
        self.read_cov[lib_name] = read_cov

    def addHits(self, hitLibName, hitDict):
        if not hitLibName in self.hits:
            self.hits[hitLibName] = []
        self.hits[hitLibName].append(hitDict)

class CovLibObj():
    def __init__(self, name, fmt, f):
        self.name = name
        self.fmt = fmt
        self.f = abspath(f) if isfile(f) else f # pass file/string/''
        self.cov_sum = 0.0
        self.reads_total = 0
        self.reads_mapped = 0
        self.reads_unmapped = 0
        self.mean_cov = 0.0

class HitLibObj():
    def __init__(self, name, fmt, f):
        self.name = name
        self.fmt = fmt
        self.f = abspath(f) if isfile(f) else f # pass file/string/''


class ViewObj():
    def __init__(self, name='', out_f='', suffix='', header='', body=''):
        self.name = name
        self.out_f = out_f
        self.suffix = suffix
        self.header = header
        self.body = body

    def output(self):
        if isinstance(self.body, dict):
            for category in self.body:
                out_f = "%s.%s.%s" % (self.out_f, category, self.suffix)
                print BtLog.status_d['13'] % (out_f)
                with open(out_f, "w") as fh:
                    fh.write(self.header + "".join(self.body[category]))
        elif isinstance(self.body, list):
            out_f = "%s.%s" % (self.out_f, self.suffix)
            print BtLog.status_d['13'] % (out_f)
            with open(out_f, "w") as fh:
                fh.write(self.header + "".join(self.body))
        else:
            sys.exit("[ERROR] - 001")

class newBlobDb():
    def __init__(self, name='blobdb'):
        # meta
        self.title = title
        self.path = path
        self.version = version
        self.blob_count = 0
        self.tax_hit_count = 0
        self.sources = {}
        self.files = {}
        self.cov_lib = []
        self.tax_lib = []
        self.tax_rule = []
        self.reads_total = {}
        self.reads_mapped = {}
        self.ranks = []

        self.blobs_list = []
        self.blobs_dict = {}

        self.blob_id = []
        self.length = []
        self.gc = []
        self.n_count = []
        self.agct_count = []

        self.cov_base = {}
        self.cov_read = {}
        self.tax = {}
        self.tax_hit = {}
        self.tax_id = {}

        self.meta = {}

    def _set_meta(self):
        self.meta = {
            "title" : self.title,
            "version" : self.version,
            "count" : self.count,
            "sources" : self.sources,
            "files" : self.files,
            "cov_lib" : self.cov_lib,
            "tax_lib" : self.tax_lib,
            "tax_rule" : self.tax_rule,
            "reads_total" : self.reads_total,
            "reads_mapped" : self.reads_mapped,
            "tax_hit_count" : self.tax_hit_count
        }

    def _parse_meta(self, meta_f):
        pass


    def _set_files(self):
        self.files = {
            # primary
            'meta' : "%s" % join(self.path, "meta"),
            'blob_id' : "%s" % join(self.path, "blob_id"),
            'length' : "%s" % join(self.path, "length"),
            'gc' : "%s" % join(self.path, "gc"),
            'n_count' : "%s" % join(self.path, "n_count"),
            'agct_count' : "%s" % join(self.path, "agct_count"),
            'tax_id' : "%s" % join(self.path, tax_id),
            # secondary
            'cov_base' : {cov_lib : "%s" % join(self.path, "cov_base") for cov_lib in self.cov_lib},
            'cov_read' : {cov_lib : "%s" % join(self.path, "cov_read") for cov_lib in self.cov_lib},
            'tax' : {tax_rule : "%s" % join(self.path, "tax") for tax_rule in self.tax_rule},
            'tax_hit' : {tax_lib : "%s" % join(self.path, "tax_hit") for tax_lib in self.tax_lib}
            }

    def _write_output(self):
        primary = ['meta', 'blob_id', 'length', 'gc', 'n_count', 'agct_count', 'tax_id']
        secondary = ['cov_base', 'cov_read', 'tax', 'tax_hit']
        directory = BtIO.create_dir(self.path)
        out_fs = []
        if (directory):
            out_fs = []
            strings = []
            for key in self.files:
                if key in primary:
                    out_f.append(self.files[key])
                    data.append(getattr(self, key))
                elif key in secondary:
                    for key2 in self.files[key]:
                        out_f.append(self.files[key][key2])
                        data.append(getattr(self, key)[key2])
                else:
                    pass
            with tarfile.open(out_f, "a:gz") as tar:
                for out_f, string in zip(out_fs, strings):
                    with open(out_f, 'w') as fh:
                        json.dump(string, fh, indent=1, separators=(',', ' : '))
                    tar.add(out_f)

    def dump(self):
        self._set_files()
        self._set_meta()
        self._write_output()

    def load(self, blobdb):
        try:
            meta_f = blobdb.getmember('meta')
        except KeyError:
            pass

    def yield_blob(self, fields):
        Blob = collections.namedtuple('Blob', fields)
        for idx, blob_id in enumerate(self.blob_id):
            data = [x for x in getattr(self, fields)]
            blob = Blob()

    def get_idxs_by_ids(self, ids):
        # list of blob_ids
        if (ids) and isinstance(ids, list):
            return [int(self.blob_dict[id_]) for _id_ in ids]

    def get_data(self, *key, **kwargs):
        idxs = []
        # sort out which idxs
        if (kwargs['ids']): # delimited by list of ids
            idxs = self.get_idxs_by_ids(kwargs['ids'])
        elif (kwargs['idxs']): # delimited by list of idxs
            idxs = kwargs['idxs']
        else: # all
            idxs = range(self.count)

        # return data based on idxs
        if args[0] == "name":
            return [self.name[idx] for idx in idxs]
        elif args[0] == "length":
            return [self.length[idx] for idx in idxs]
        elif args[0] == "gc":
            return [self.gc[idx] for idx in idxs]
        elif args[0] == "n_count":
            return [self.n_count[idx] for idx in idxs]
        elif args[0] == "cov_base":
            if args[1] in self.cov_lib:
                return [self.cov_base[idx] for idx in idxs]
        elif args[0] == "cov_read":
            if args[1] in self.cov_lib:
                return [self.cov_read[idx] for idx in idxs]
        elif args[0] == "tax":
            if args[1] in self.tax_rule:
                if args[2] in self.ranks:
                    return [self.tax[idx] for idx in idxs]
        elif args[0] == "tax_hit":
            if args[1] in self.tax_lib:
                return [self.tax_hit[idx] for idx in idxs]
        else:
            return None

    # Parse
    def parse_meta(self, meta_f):
        meta = BtIO.read_meta(meta_f)
        for key, value in meta.items():
            setattr(key, value)
        self.name = meta['name']
        self.covlib = meta['covlib']
        self.taxrule = meta['taxrule']
        self.taxlib = meta['taxlib']
        self.files = meta['files']
        self.blobs_count = meta['count']
        self.ranks = meta['ranks']

    def parse_data(self, key, exp_count):
        data = BtIO.read_json_list(self.files[key])
        if not len(data) == exp_count:
            # error
            pass
        else:
            setattr(self, key, data)

    def output(self):
        # meta
        meta = self.get_meta()
        meta_f = join(self.view_dir, "meta.json")
        BtIO.writeJson(meta, meta_f, indent=2)
        # gc
        gc_f = join(self.view_dir, "gc.json")
        print BtLog.status_d['13'] % (gc_f)
        BtIO.writeJson(self.gc, gc_f, indent=1)
        # length
        length_f = join(self.view_dir, "length.json")
        print BtLog.status_d['13'] % (length_f)
        BtIO.writeJson(self.length, length_f, indent=1)
        # names
        names_f = join(self.view_dir, "names.json")
        print BtLog.status_d['13'] % (names_f)
        BtIO.writeJson(self.names, names_f, indent=1)
        # cov
        cov_d = join(self.view_dir, "covs")
        BtIO.create_dir(directory=cov_d)
        for cov_lib, cov in self.covs.items():
            cov_f = join(cov_d, "%s.json" % cov_lib)
            print BtLog.status_d['13'] % (cov_f)
            BtIO.writeJson(cov, cov_f, indent=1)
        # tax
        taxrule_d = join(self.view_dir, "taxrule")
        BtIO.create_dir(directory=taxrule_d)
        for taxrule in self.tax:
            tax_d = join(taxrule_d, taxrule)
            BtIO.create_dir(directory=tax_d)
            for rank in self.tax[taxrule]:
                tax = self.tax[taxrule][rank]
                rank_f = join(tax_d, "%s.json" % rank)
                BtIO.writeJson(tax, rank_f, indent=1)

class ExperimentalViewObj():
    def __init__(self, name='experimental', view_dir=''):
        self.name = name
        self.view_dir = view_dir
        self.length = []
        self.gc = []
        self.names = []
        self.tax = {}
        self.covs = {}
        self.meta = {}
        BtIO.create_dir(self.view_dir)


    def get_meta(self):
        meta = {
              "name": self.view_dir,
              "datatypes":{
                "identifiers": {"name" : "identifiers", "type" : "identifier"},
                "len": {"name" : "length","type" : "integer", "min" : min(self.length), "max" : max(self.length)},
                "gc": {"name" : "gc","type" : "percentage", "min" : min(self.gc), "max" : max(self.gc)},
                "covs": {"name" : "covs","type" : "float", "min" : 0.02, "max" : max(self.covs["covsum"])},
                "taxrule":{"name" : "taxrule","type" : "category"},
              },
              "hierarchy":{
                "covs": [cov for cov in self.covs],
                "taxrule": [taxrule for taxrule in self.tax],
                "besthit": [rank for rank in BtTax.RANKS],
                "bestsum": [rank for rank in BtTax.RANKS],
                #"superkingdom":["bitscore","cindex","evalue"],
                #"phylum":["bitscore","cindex","evalue"]
              },
              "rows" : len(self.names),
              "mainplot":{
                "type":"hexbin",
                "x":{
                  "datatype":["gc"],
                  "bounds":[0,100],
                  "scale":"linear"
                },
                "y":{
                  "datatype": [cov_lib for cov_lib in self.covs],
                  "bounds":[0.001,1000000],
                  "scale":"log10"
                },
                "series":{
                  "datatype":["taxrule","besthit","superkingdom"],
                  "scale":"sqrt"
                }
              }
            }
        for cov_lib in self.covs:
            meta['datatypes'][cov_lib] = {"name" : cov_lib, "type" : "float"}
        for taxrule in self.tax:
            meta['datatypes'][taxrule] = {"name": taxrule, "type":"category"}
        for rank in BtTax.RANKS:
            meta['datatypes'][rank] = {"name": rank, "type":"category", "levels" : 7}
        return meta

    def output(self):
        # meta
        meta = self.get_meta()
        meta_f = join(self.view_dir, "meta.json")
        BtIO.writeJson(meta, meta_f, indent=2)
        # gc
        gc_f = join(self.view_dir, "gc.json")
        print BtLog.status_d['13'] % (gc_f)
        BtIO.writeJson(self.gc, gc_f, indent=1)
        # length
        length_f = join(self.view_dir, "length.json")
        print BtLog.status_d['13'] % (length_f)
        BtIO.writeJson(self.length, length_f, indent=1)
        # names
        names_f = join(self.view_dir, "names.json")
        print BtLog.status_d['13'] % (names_f)
        BtIO.writeJson(self.names, names_f, indent=1)
        # cov
        cov_d = join(self.view_dir, "covs")
        BtIO.create_dir(directory=cov_d)
        for cov_lib, cov in self.covs.items():
            cov_f = join(cov_d, "%s.json" % cov_lib)
            print BtLog.status_d['13'] % (cov_f)
            BtIO.writeJson(cov, cov_f, indent=1)
        # tax
        taxrule_d = join(self.view_dir, "taxrule")
        BtIO.create_dir(directory=taxrule_d)
        for taxrule in self.tax:
            tax_d = join(taxrule_d, taxrule)
            BtIO.create_dir(directory=tax_d)
            for rank in self.tax[taxrule]:
                tax = self.tax[taxrule][rank]
                rank_f = join(tax_d, "%s.json" % rank)
                BtIO.writeJson(tax, rank_f, indent=1)


if __name__ == '__main__':
    pass
