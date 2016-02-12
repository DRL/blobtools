#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File        : BtTax.py
Version     : 0.1
Author      : Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs        : ?
To do       : ?
"""
from __future__ import division
RANKS = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']

def noHit():
	return {rank : {'tax' : 'no-hit', 'score' : 0.0, 'c_index' : None} for rank in RANKS}

def getTreeList(taxIds, nodesDB):
    known_tree_lists = {}
    for taxId in taxIds: 
    	if not taxId in known_tree_lists:
    		tree_list = []
    		nextTaxId = [taxId]
    		while nextTaxId:
        		thisTaxId = nextTaxId.pop(0)
        		if (not thisTaxId == '1') and (thisTaxId in nodesDB):
        			parent = nodesDB[thisTaxId]['parent']
        			nextTaxId.append(parent)
        			tree_list.append(thisTaxId)
        		else:
    				tree_list.append('1')
    		known_tree_lists[taxId] = tree_list
    return known_tree_lists

def getLineages(tree_lists, nodesDB):
	lineage = {}
	for tree_list_id, tree_list in tree_lists.items():
		lineage[tree_list_id] = {rank : 'undef' for rank in RANKS}
		for taxId in tree_list:
			node = nodesDB[taxId]
			if node['rank'] in RANKS:
				lineage[tree_list_id][node['rank']] = node['name']
		# traverse ranks again so that undef is "higher_def_rank" + "-" + undef
		def_rank = ''
		for rank in reversed(list(RANKS)):
			if not lineage[tree_list_id][rank] == 'undef':
				def_rank = lineage[tree_list_id][rank]
			else:
				if (def_rank):
					lineage[tree_list_id][rank] = def_rank + "-" + lineage[tree_list_id][rank]
	return lineage

def taxRule(taxrule, hits, lineages):
	taxonomy = { rank : {'tax' : '', 'score' : 0.0, 'c_index' : 0 } for rank in RANKS } 
	tempTax = { rank : {} for rank in RANKS } 
	taxDict = getTaxDict(hits, lineages) # here libs are separated
	if taxrule == 'bestsum':
		for lib in sorted(taxDict):
			for rank in RANKS:
				for tax, score in sorted(taxDict[lib][rank].items()): 
					tempTax[rank][tax] = tempTax[rank].get(tax, 0.0) + score
		for rank in tempTax:
			for tax, score in sorted(tempTax[rank].items(), key=lambda x: x[1], reverse=True):
				if not (taxonomy[rank]['tax']):
					taxonomy[rank]['tax'] = tax
					taxonomy[rank]['score'] = score
				else:
					taxonomy[rank]['c_index'] += 1 
	elif taxrule == 'bestsumorder':
		for lib in sorted(taxDict):
			hit_lib = 0
			for rank in RANKS:
				for tax, score in sorted(taxDict[lib][rank].items(), key=lambda x: x[1], reverse=True):
					if not (taxonomy[rank]['tax']):
						taxonomy[rank]['tax'] = tax
						taxonomy[rank]['score'] = score
						hit_lib = 1	
					else:
						if (hit_lib):
							taxonomy[rank]['c_index'] += 1  
	else:
		pass
	return taxonomy 

def getTaxDict(hits, lineages):
	taxDict = {}
	for lib, hits in hits.items():
		taxDict[lib] = {}
		for hit in hits:
			taxId = hit['taxId']
			score = hit['score']
			for rank in RANKS:
				name = lineages[taxId][rank]
				if not rank in taxDict[lib]:
					taxDict[lib][rank] = {name : 0.0}
				taxDict[lib][rank][name] = taxDict[lib][rank].get(name, 0.0) + score  
	return taxDict

if __name__ == "__main__": 
    pass