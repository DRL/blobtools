#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import lib.BtLog as BtLog
import lib.BtIO as BtIO
import lib.BtTax as BtTax
from os.path import abspath

class BlobDb():
    '''
    class BlobDB holds all information parsed from files
    '''
    def __init__(self, title):
        self.title = title
        self.assembly_f = ''
        self.dict_of_blobs = {}  
        self.order_of_blobs = []
        self.set_of_taxIds = set()
        self.lineages = {}
        self.length = 0
        self.seqs = 0
        self.n_count = 0
        self.covLibs = {} 
        self.hitLibs = {}
        self.nodesDB_f = ''
        self.taxrules = []

    def view(self, out_f, ranks, taxrule, hits_flag, seqs):
        sep = "\t"
        header = ''
        body = ''
        header += "##\n"
        header += "## assembly : %s\n" % self.assembly_f
        header += "%s\n" % "\n".join("## " + covLib['name'] + " : " + covLib["f"] for covLib in self.covLibs.values())
        header += "%s\n" % "\n".join("## " + hitLib['name'] + " : " + hitLib["f"] for hitLib in self.hitLibs.values())
        header += "## nodesDB : %s\n" % self.nodesDB_f
        header += "##\n"    
        header += "# %s" % sep.join(map(str, [ "name", "length", "GC", "N"  ])) 
        header += "%s%s" % (sep, sep.join([cov_lib_name for cov_lib_name in self.covLibs]))
        for rank in ranks:
            header += "%s%s" % (sep, sep.join([rank + ".t", rank + ".s", rank + ".c"]))
            if hits_flag:
                header += "%s%s" % (sep, rank + ".hits")
        if (seqs):
            for seq in seqs:
                blob = self.dict_of_blobs[seq]
                body += self.getViewLine(blob, taxrule, ranks, sep, hits_flag)
        else:
            for name in self.order_of_blobs:
                blob = self.dict_of_blobs[name]
                body += self.getViewLine(blob, taxrule, ranks, sep, hits_flag)

        if out_f == "STDOUT":
            print header + body
        else:
            with open(out_f, 'w') as fh:
                fh.write(header + body)
    
    def getViewLine(self, blob, taxrule, ranks, sep, hits_flag):
        line = ''
        line += "\n%s" % sep.join(map(str, [ blob['name'], blob['length'], blob['gc'], blob['n_count']  ])) 
        line += sep
        line += "%s" % sep.join(map(str, [ blob['covs'][covLib] for covLib in self.covLibs]))
        for rank in ranks:
            line += sep
            tax = blob['taxonomy'][taxrule][rank]['tax']
            score = blob['taxonomy'][taxrule][rank]['score']
            c_index = blob['taxonomy'][taxrule][rank]['c_index']
            line += "%s" % sep.join(map(str, [tax, score, c_index ]))
            line += sep
            if hits_flag:
                for tax_lib_name in sorted(self.hitLibs):
                    if tax_lib_name in blob['hits']:
                        line += "%s=" % tax_lib_name
                        sum_dict = {}
                        for hit in blob['hits'][tax_lib_name]:
                            tax_rank = self.lineages[hit['taxId']][rank] 
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
                'taxrules' : self.taxrules
                }
        return dump

    def load(self, BlobDb_f):
        blobDict = BtIO.readJson(BlobDb_f)
        self.title = blobDict['title']
        self.assembly_f = blobDict['assembly_f']
        self.nodesDB_f = blobDict['nodesDB_f']
        self.lineages = blobDict['lineages']
        self.set_of_taxIds = blobDict['lineages'].keys()
        self.order_of_blobs = blobDict['order_of_blobs']
        self.dict_of_blobs = blobDict['dict_of_blobs'] # this will probably not work
        self.length = int(blobDict['length'])
        self.seqs = int(blobDict['seqs'])
        self.n_count = int(blobDict['n_count'])
        self.covLibs = blobDict['covLibs']
        self.hitLibs = blobDict['hitLibs']
        self.taxrules = blobDict['taxrules']

    def getArrays(self, rank, min_length, hide_nohits, taxrule, c_index, label_d):
        from numpy import array
        summary_dict = {}
        data_list = []
        cov_dict = {covLib : [] for covLib in self.covLibs}
        for blob in self.dict_of_blobs.values():
            name = blob['name']
            gc = blob['gc']
            length = blob['length']
            tax = ''
            if (c_index):
                tax = str(blob['taxonomy'][taxrule][rank]['c_index'])
            else:
                tax = blob['taxonomy'][taxrule][rank]['tax']
                if label_d and tax in label_d:
                    tax = label_d[tax] 
            if not tax in summary_dict:
                summary_dict[tax] = {'count_total' : 0,
                                     'count_hidden' : 0,
                                     'count_visible' : 0,
                                     'span_total': 0, 
                                     'span_hidden' : 0, 
                                     'span_visible' : 0}
            if ((hide_nohits) and tax == 'no-hit') or length < min_length:
                summary_dict[tax]['count_hidden'] = summary_dict[tax].get('count_hidden', 0) + 1
                summary_dict[tax]['span_hidden'] = summary_dict[tax].get('span_hidden', 0) + length
            else:
                data_list.append([(name), (length), (gc), (tax)])
                for covLib in self.covLibs:
                    cov = float(blob['covs'][covLib])
                    if cov < 0.1:
                        cov = 0.1
                    cov_dict[covLib].append(cov)
                summary_dict[tax]['count_visible'] = summary_dict[tax].get('count_visible', 0) + 1
                summary_dict[tax]['span_visible'] = summary_dict[tax].get('span_visible', 0) + int(length)
            summary_dict[tax]['count_total'] = summary_dict[tax].get('count_total', 0) + 1
            summary_dict[tax]['span_total'] = summary_dict[tax].get('span_total', 0) + int(length)
        data_array = array(data_list)
        cov_arrays = {covLib: array(cov) for covLib, cov in cov_dict.items()}
        return data_array, cov_arrays, summary_dict

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
            
    def parseCovs(self, covLibObjs):
        for covLib in covLibObjs:
            self.addCovLib(covLib)
            print BtLog.status_d['1'] % (covLib.name, covLib.f)
            if covLib.fmt == 'bam' or covLib.fmt == 'sam':
                base_cov_dict = {}
                if covLib.fmt == 'bam':
                    base_cov_dict, covLib.total_reads, covLib.mapped_reads, covLib.read_cov_dict = BtIO.readBam(covLib.f, set(self.dict_of_blobs))
                else:
                    base_cov_dict, covLib.total_reads, covLib.mapped_reads, covLib.read_cov_dict = BtIO.readSam(covLib.f, set(self.dict_of_blobs))    
                for name, base_cov in base_cov_dict.items():
                    cov = base_cov / self.dict_of_blobs[name].agct_count
                    covLib.cov_sum += cov
                    self.dict_of_blobs[name].addCov(covLib.name, cov)
            elif covLib.fmt == 'cas':
                for name, cov in BtIO.readCas(covLib.f, self.order_of_blobs):
                    covLib.cov_sum += cov
                    self.dict_of_blobs[name].addCov(covLib.name, cov)
            elif covLib.fmt == 'cov':
                for name, cov in BtIO.readCov(covLib.f, set(self.dict_of_blobs)):
                    covLib.cov_sum += cov
                    self.dict_of_blobs[name].addCov(covLib.name, cov)
            else:
                pass        
            covLib.mean_cov = covLib.cov_sum/self.seqs
            self.covLibs[covLib.name] = covLib

    def parseHits(self, hitLibs):
        for hitLib in hitLibs:
            self.hitLibs[hitLib.name] = hitLib
            print BtLog.status_d['1'] % (hitLib.name, hitLib.f)
            # only accepts format 'seqID\ttaxID\tscore'
            for hitDict in BtIO.readTax(hitLib.f, set(self.dict_of_blobs)):
                self.set_of_taxIds.add(hitDict['taxId'])
                self.dict_of_blobs[hitDict['name']].addHits(hitLib.name, hitDict)
    
    def computeTaxonomy(self, taxrules, nodesDB):
        tree_lists = BtTax.getTreeList(self.set_of_taxIds, nodesDB)
        self.lineages = BtTax.getLineages(tree_lists, nodesDB)
        self.taxrules = taxrules
        i = 0
        for blObj in self.dict_of_blobs.values():
            i += 1
            BtLog.progress(i, 100, self.seqs)
            for taxrule in taxrules:
                if (blObj.hits):
                    blObj.taxonomy[taxrule] = BtTax.taxRule(taxrule, blObj.hits, self.lineages)
                else:
                    blObj.taxonomy[taxrule] = BtTax.noHit()

    def counts(self):
        count_dict = {
            'seqs'     : self.seqs,
            'length'   : self.length,
            'Ns'       : self.n_count,
            'AvgCov'   : {lib : round(covlibObj.cov_sum/self.seqs, 2) for lib, covlibObj in self.covLibs.items()},
            'GC'       : round(sum([blObj.gc for blObj in self.dict_of_blobs.values()])/self.seqs, 2),
            'MappedReads' : {lib : (covlibObj.mapped_reads) for lib, covlibObj in self.covLibs.items()},
            'TotalReads' : {lib : (covlibObj.total_reads) for lib, covlibObj in self.covLibs.items()}
        }
        print count_dict
    
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
        self.hits = {}
        self.taxonomy = {}
        
    def calculateGC(self, seq):
        return float((seq.count('G') + seq.count('C') ) / self.agct_count \
                     if self.agct_count > 0 else 0.0)

    def addCov(self, name, cov):
        self.covs[name] = cov

    def addHits(self, hitLibName, hitDict):
        if not hitLibName in self.hits:
            self.hits[hitLibName] = []
        self.hits[hitLibName].append(hitDict)

class CovLibObj():
    def __init__(self, name, fmt, f):
        self.name = name
        self.fmt = fmt
        self.f = abspath(f)
        self.cov_sum = 0
        self.total_reads = 0
        self.mapped_reads = 0
        self.read_cov_dict = {} 
        self.mean_cov = 0.0

class hitLibObj():
    def __init__(self, name, fmt, f):
        self.name = name
        self.fmt = fmt
        self.f = abspath(f)

if __name__ == '__main__':
    pass