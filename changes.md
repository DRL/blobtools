# blobtools v0.9.19
- numerous bugfixes and refactoring
- version numbers now printed in all output files
- blobtools create
    - BAM/SAM parsing:
        - BAM and SAM parsing now fast(er)
        - Support for X/= in CIGAR string (previously only M's)
    - BAM/SAM/CAS parsing:
        - when parsing mapping files, blobtools automatically writes COV files
    - blobDB JSON parser tries first "ujson" module, then "simplejson" module, before defaulting to "json" module.
        "ujson" is fastest, so "pip install ujson" is encouraged
    - taxonomy calculation
        - if two highest scoring taxonomies for a given contig have EQUAL scores, taxonomy gets set to "unresolved"
            - during plotting, the group "unresolved" gets treated as any other group
            - new flag "--tax_collision_random" allows legacy behaviour (random selection of taxonomy)
        - new flag "--min_diff FLOAT" sets minimal score difference between two highest scoring
            taxonomies of a contig in order for the contig to not be called "unresolved"
- blobtools view
    - refactored and is now fast(er)
    - supported views:
        - table: legacy output, taxonomy columns are now numbered (1-based, for cut'ing)
        - concoct: outputs taxonomy and coverage information necessary to run concoct (https://github.com/BinPro/CONCOCT/)
        - cov: outputs coverage files (for using them in scattercov)
        - experimental: something new we are trying,
- blobtools map2cov (previously bam2cov):
    - supports BAM, SAM and CAS files (uses parsers in BtIO)


- blobtools blobplot
    - changed behaviour of "--multiplot", now plots each group separately, then a final plot all together.
    - new flag "--cumulative". Previous "--multiplot" behaviour: incremental addition of groups
    - new flag "--legend". Plots legend in a separate figure, useful for slides.
    - new flag "--lib". Allows selecting particular covlib(s) to plot.
    - fixed bug causing different plots to share filenames, causing overwriting
    - changed scaling of blobs in plots:
        - previously, blobs had an area of (length/1000*65) pixels^2. But if an assembly had a few very big contigs, these covered the smaller ones.
        - now the scaling is done in the following way:
            - the biggest contig gets plotted as a blob of area 12500 pixels^2
            - all other contigs get scaled accordingly
        - reference scale was changed as well:
            - scale shows area of 0.05, 0.1 and 0.25 of longest contig
    - behaviour of colour_f changed, so that now all groups that are not named in colour_f are colour in white
    - added x-axis minor ticklabels, all ticklabels facing outwards
    - [--exclude] removes group from plotting and stats
    - [--filelabel] : labels coverage axis based on coverage file

- blobtools covplot
    - new flag "--lib". Allows selecting particular covlib(s) to plot.
    - both x/y-axes have the same limits (max-cov + 1000)

- blobtools sumcov
    - removed

- blobtools taxify
    - allows setting taxID to taxfiles (sequence similarity search results) using ID-mapping-file
    - rDNA sequences using RNAcentral
        - ID mapping
            - wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/5.0/id_mapping/id_mapping.tsv.gz
            - gunzip id_mapping.tsv.gz
            - Subset SILVA sequence IDs in RNAcentral (id_mapping.tsv)
                grep SILVA id_mapping.tsv | sort | uniq > id_mapping.SILVA.tsv
            - Create list of sequence IDs:
                cut -f1 id_mapping.SILVA.tsv > id_mapping.SILVA.names.txt
        - Sequences
            - wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/5.0/sequences/rnacentral_active.fasta.gz
            - convert rRNA to rDNA
                perl -ne 'chomp; if(m/^>/){@temp=split(" "); print $temp[0]."\n";} else {$_ =~ tr/U/T/; print $_."\n"}' rnacentral_active.fasta > rnacentral_active.rDNA.fasta
            - Subset SILVA sequences in RNAcentral (rnacentral_active.rDNA.fasta)
                blobtools seqfilter -i rnacentral_active.rDNA.fasta -l id_mapping.SILVA.names.txt
    - Taxify blast result against SILVA sequences in RNAcentral
    -
To do:
- check that create makes cov files when bam
- move all *.py to bloblib
- change imports
- make maxcov a filter in getPlotData so that stats for covplots reflect hidden blobs
- write all messages to logfile (including execution parameters)
-

bwa mem -t 52 -M -R '@RG\tID:nHa.MiSeq.400\tSM:nHa.1\tPL:illumina\tLB:nHa.MiSeq.400\tPU:unit1' \
assembly/nHa.MiSeq.400.clc.se.fna \
/exports/colossus/howardula/howardula_aoronymphium/analysis/1_adapter_trimming/nHa.MiSeq.400.trim.1.fastq.gz \
/exports/colossus/howardula/howardula_aoronymphium/analysis/1_adapter_trimming/nHa.MiSeq.400.trim.2.fastq.gz \
| samtools view -bS - > nHa.MiSeq.400.vs.nHa.MiSeq.400.clc.se.bam ; \
bwa mem -t 52 -M -R '@RG\tID:nDp.HiSeq.200.B\tSM:nDp.1\tPL:illumina\tLB:nDp.HiSeq.200\tPU:unit1' \
assembly/nHa.MiSeq.400.clc.se.fna \
/exports/colossus/howardula/howardula_aoronymphium/analysis/1_adapter_trimming/nDp.HiSeq.200.B.trim.1.fastq.gz \
/exports/colossus/howardula/howardula_aoronymphium/analysis/1_adapter_trimming/nDp.HiSeq.200.B.trim.2.fastq.gz \
| samtools view -bS - > nDp.HiSeq.200.B.vs.nHa.MiSeq.400.clc.se.bam ; exit

Milestones:
- deal with problematic sequences in databases (nr)
    - by type (transposons)
    - find non-informative sequences and ignore
    - take contigs that have "weird colour"
    - extract sequence, blast,
    - overenthusiasticly predicted proteins
    - ribosomal promotor associated proteins (ORF in rRNA)
    - aging genes in plants (ORF in rRNA)
- new blastdb format (will speed everythin up)
- automated testing



Format

cov
    - covs
taxrule
    - name
        - rank
length

