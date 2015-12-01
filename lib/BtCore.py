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
        self.dict_of_blobs = blobDict['dict_of_blobs'] 
        self.length = int(blobDict['length'])
        self.seqs = int(blobDict['seqs'])
        self.n_count = int(blobDict['n_count'])
        self.covLibs = blobDict['covLibs']
        self.hitLibs = blobDict['hitLibs']
        self.taxrules = blobDict['taxrules']

    def getPlotData(self, rank, min_length, hide_nohits, taxrule, c_index, catcolour_dict):
        data_dict = {}
        read_cov_dict = {}
        max_cov = 0.0
        cov_libs = self.covLibs.keys()
        cov_libs_reads_total = {cov_lib : data['reads_total'] for cov_lib, data in self.covLibs.items()}

        for blob in self.dict_of_blobs.values():
            name, gc, length, group = blob['name'], blob['gc'], blob['length'], ''
            
            if (catcolour_dict): # annotation with categories specified in catcolour
                group = str(catcolour_dict[name])
            elif (c_index): # annotation with c_index instead of taxonomic group
                group = str(blob['taxonomy'][taxrule][rank]['c_index'])
            else: # annotation with taxonomic group
                group = str(blob['taxonomy'][taxrule][rank]['tax'])
            
            if not group in data_dict: 
                data_dict[group] = {
                                    'name' : [], 
                                    'length' : [], 
                                    'gc' : [], 
                                    'covs' : {covLib : [] for covLib in cov_libs}, 
                                    'reads_mapped' : {covLib : 0 for covLib in cov_libs},
                                    'count' : 0,
                                    'count_hidden' : 0,
                                    'count_visible' : 0,
                                    'span': 0, 
                                    'span_hidden' : 0, 
                                    'span_visible' : 0,
                                    }
                if len(cov_libs) > 1:
                    data_dict[group]['covs']['cov_sum'] = []
                    data_dict[group]['reads_mapped']['cov_sum'] = 0

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
            for cov_lib in sorted(cov_libs):
                cov = float(blob['covs'][cov_lib]) 
                cov_sum += cov
                cov = cov if cov > 0.02 else 0.02
                if cov > max_cov:
                    max_cov = cov
                data_dict[group]['covs'][cov_lib].append(cov)
                if cov_lib in blob['read_cov']:
                    reads_mapped = blob['read_cov'][cov_lib]
                    reads_mapped_sum += reads_mapped
                    data_dict[group]['reads_mapped'][cov_lib] += reads_mapped  
            
            if len(cov_libs) > 1:
                cov_sum = cov_sum if cov_sum > 0.02 else 0.02
                data_dict[group]['covs']['cov_sum'].append(cov_sum)
                if cov > max_cov:
                    max_cov = cov
                if (reads_mapped_sum):
                    data_dict[group]['reads_mapped']['cov_sum'] += reads_mapped_sum

            data_dict[group]['count'] = data_dict[group].get('count', 0) + 1
            data_dict[group]['span'] = data_dict[group].get('span', 0) + int(length)

        if len(cov_libs) > 1:
            cov_libs.append('cov_sum')
            for cov_lib, data in self.covLibs.items():
                cov_libs_reads_total['cov_sum'] = cov_libs_reads_total.get('cov_sum', 0) + data['reads_total'] 

        return data_dict, max_cov, cov_libs, cov_libs_reads_total

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
                    base_cov_dict, covLib.reads_total, covLib.reads_mapped, read_cov_dict = BtIO.readBam(covLib.f, set(self.dict_of_blobs))
                else:
                    base_cov_dict, covLib.reads_total, covLib.reads_mapped, read_cov_dict = BtIO.readSam(covLib.f, set(self.dict_of_blobs))    
                if covLib.reads_total == 0:
                    print BtLog.warn_d['4'] % covLib.f
                for name, base_cov in base_cov_dict.items():
                    cov = base_cov / self.dict_of_blobs[name].agct_count
                    covLib.cov_sum += cov
                    self.dict_of_blobs[name].addCov(covLib.name, cov)
                    self.dict_of_blobs[name].read_cov = {covLib.name : read_cov_dict[name]}
            elif covLib.fmt == 'cas':
                cov_dict, covLib.reads_total, covLib.reads_mapped, read_cov_dict = BtIO.readCas(covLib.f, self.order_of_blobs)
                if covLib.reads_total == 0:
                    print BtLog.warn_d['4'] % covLib.f
                for name, cov in cov_dict.items():
                    covLib.cov_sum += cov
                    self.dict_of_blobs[name].addCov(covLib.name, cov)
                    self.dict_of_blobs[name].read_cov = {covLib.name : read_cov_dict[name]}
            elif covLib.fmt == 'cov':
                cov_dict = BtIO.readCov(covLib.f, set(self.dict_of_blobs))
                if not len(cov_dict) == self.seqs:
                    print BtLog.warn_d['4'] % covLib.f
                for name, cov in cov_dict.items():
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
                if ";" in hitDict['taxId']:
                    hitDict['taxId'] = hitDict['taxId'].split(";")[0]
                    print BtLog.warn['5'] % (hitDict['name'], hitLib)
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
            'MappedReads' : {lib : (covlibObj.reads_mapped) for lib, covlibObj in self.covLibs.items()},
            'TotalReads' : {lib : (covlibObj.reads_total) for lib, covlibObj in self.covLibs.items()}
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
        self.read_cov = {}
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
        self.reads_total = 0
        self.reads_mapped = 0
        self.mean_cov = 0.0

class hitLibObj():
    def __init__(self, name, fmt, f):
        self.name = name
        self.fmt = fmt
        self.f = abspath(f)

if __name__ == '__main__':
    pass