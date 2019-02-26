#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File        : BtLog.py
Author      : Dominik R. Laetsch, dominik.laetsch at gmail dot com
"""

from __future__ import division
import sys

def error(message, *argv):
    if argv is None:
        sys.exit(error_d[message])
    else:
        sys.exit(error_d[message] % (argv))
    #exit(1)  # change to exit with the actual ERROR number (different than 0)

error_d = {
    '0': '[ERROR:0]\t: File %s does not exist.',
    '1': '[ERROR:1]\t: Please provide coverage information.',
    '2': '[ERROR:2]\t: Assembly type is not valid (%s).',
    '3': '[ERROR:3]\t: names.dmp/nodes.dmp ("--names", "--nodes") could not be read. %s, %s',
    '4': '[ERROR:4]\t: BlobDB.parseFasta() - no sequences found. Check FASTA file.',
    '5': '[ERROR:5]\t: Sequence header %s is not unique.',
    '6': '[ERROR:6]\t: BlobDB.readBam() - sequence header %s in %s was not in FASTA.',
    '7': '[ERROR:7]\t: Please add "samtools" to you PATH variable.',
    '8': '[ERROR:8]\t: Unsupported taxrule "%s".',
    '9': '[ERROR:9]\t: Unsupported taxonomic rank "%s".',
    '10': '[ERROR:10]\t: Unsupported output format "%s".',
    '11': '[ERROR:11]\t: Taxrule "%s" was not computed for this BlobDb. Available taxrule(s) : %s.',
    '12': '[ERROR:12]\t: Please provide an output file.',
    '13': '[ERROR:13]\t: %s does not appear to be a comma-separated list or a file.',
    '14': '[ERROR:14]\t: Unsupported sort order for plotting : %s. Must be either "span" or "count".',
    '15': '[ERROR:15]\t: Unsupported histogram type for plotting : %s. Must be either "span" or "count".',
    '16': '[ERROR:16]\t: Group "%s" was specified in multiple clusters.',
    '17': '[ERROR:17]\t: Label could not be parsed from "%s".',
    '18': '[ERROR:18]\t: Please provide a tax file in BLAST format.',
    '19': '[ERROR:19]\t: Sequence %s in file %s is not part of the assembly.',
    '20': '[ERROR:20]\t: Please add "clc_mapping_info" to your PATH variable.',
    '21': '[ERROR:21]\t: Refcov file %s does not seem to have the right format.',
    '23': '[ERROR:23]\t: Catcolour file %s does not seem to have the right format.',
    '24': '[ERROR:24]\t: Catcolour file incompatible with c-index colouring.',
    '25': '[ERROR:25]\t: COV file %s does not seem to have the right format.',
    '26': '[ERROR:26]\t: TaxID must be integer.',
    '27': '[ERROR:27]\t: nodesDB ("--db") %s could not be read.',
    '28': '[ERROR:28]\t: Please specify "--names" and "--nodes", or "--db"',
    '29': '[ERROR:29]\t: No mapping reads found in %s',
    '30': '[ERROR:30]\t: The module docopt is not installed. Please install it to run blobtools\n\tpip install docopt',
    '31': '[ERROR:31]\t: Please specify a read mapping file (BAM/SAM/CAS)',
    '32': '[ERROR:32]\t: Choose either --cumulative or --multiplot',
    '33': '[ERROR:33] : CovLib(s) not found. The available covlibs are: \n%s',
    '34': '[ERROR:34] : Invalid plot type : %s',
    '35': '[ERROR:35] : Directory %s could not be created',
    '36': '[ERROR:36] : View %s could not be created',
    '37': '[ERROR:37] : %s does not seem to be a valid %s file',
    '38': '[ERROR:38] : %s is not an integer',
    '39': '[ERROR:39] : Please specify a taxid file (mapping subjects to taxids)',
    '40': '[ERROR:40] : CovLib \'%s\' not specified in refcov file',
    '41': '[ERROR:41] : Please specify either a mapping file or a taxID.',
    '42': '[ERROR:42] : SubjectID %s not found in mapping file %s.',
    '43': '[ERROR:43] : %s could not be found.',
    '44': '[ERROR:44] : Please specify integers for --map_col_sseqid and --map_col_taxid.',
    '45': '[ERROR:45] : Both --min_score and --min_diff must be numbers.',
    '46': '[ERROR:46] : Score in %s must be a float, not \'%s\'.',
    '47': '[ERROR:47] : Cannot create new "--db" file from "--names", "--nodes", "--db" file exists. %s'
}

warn_d = {
    '0': '[-] No tax files specified.',
    '1': '[-] %s not in colour file %s ...',
    '2': '[-] %s is not part of the assembly',
    '3': '\n[-] Based on samtools flagstat: expected %s reads, %s reads were parsed',
    '4': '[-] No coverage data found in %s',
    '5': '[-] Hit for sequence %s in tax file %s has multiple taxIds, only first one is used.',
    '6': '[-] Sum of coverage in cov lib %s is 0.0. Please ignore this warning if "--no_base_cov" was specified.',
    '7': '[-] No taxonomy information found.',
    '8': '[-] Duplicated sequences found :\n\t\t\t%s',
    '9': '[-] Taxrule "%s" was not computed for this BlobDb. Available taxrule(s) : %s. Will proceed without taxonomic annotation ...',
    '10': '[-] Line %s: sequence "%s" already has TaxID "%s". Skipped. (use --force to overwrite)',
    '11': '\n[-] The BAM file appears to be truncated.',
    '12': '[-] sseqid %s not found in ID-to-taxID mapping file %s.',
    '13': '[-] Sequence %s in file %s is not part of the assembly.'
}
status_d = {
    '0': '[+] Nothing to be done. %s',
    '1': '[+] Parsing %s - %s',
    '2': '[+] Done',
    '3': '[+] Creating nodesDB from %s and %s',
    '4': '[+] names.dmp/nodes.dmp not specified. Retrieving nodesDB from %s',
    '5': '[+] Store nodesDB in default location %s',
    '6': '[+] Computing taxonomy using taxrule(s) %s',
    '7': '[+] Generating BlobDB and writing to file %s',
    '8': '[+] Plotting %s',
    '9': '[+] Reading BlobDB %s',
    '10': '[+] \tChecking with \'samtools flagstat\'',
    '11': '[+] \tMapping reads = %s, total reads = %s (mapping rate = %s)',
    '12': '[+] \tChecking with \'clc_mapping_info\'',
    '13': '[+] \tWriting %s',
    '14': '[+] Preparing view(s) ...',
    '15': '[+] \tLoading BlobDB into memory ...',
    '16': '[+] \tDeserialising BlobDB (using \'%s\' module) (this may take a while) ...',
    '17': '[+] \tDeserialising BlobDB (using \'%s\' module) (this may take a while, consider installing the \'ujson\' module) ...',
    '18': '[+] Extracting data for plots ...',
    '19': '[+] Writing output ...',
    '20': '[+] \tFinished in %ss',
    '22': '[+] Filtering %s ...',
    '23': '[+] Filtered %s (pairs=%s) ...',
    '24': '[+] Writing %s',
    '25': '[+] Gzip\'ing %s',
    '26': '[+] Reading %s',
    '27': '[+] Creating nodesDB %s from %s and %s',
    '28': '[+] Store nodesDB in %s',
}

info_d = {
    '0': '[I]\t%s : sequences = %s, span = %s MB, N50 = %s nt'
    }

if __name__ == "__main__":
    pass
