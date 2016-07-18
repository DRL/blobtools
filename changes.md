# blobtools v0.9.19
- numerous bugfixes and refactoring
- version numbers now printed in all output files
- blobtools create
    - BAM/SAM parsing:
        - BAM and SAM parsing now fast(er)
        - Support for X/= in CIGAR string (previously only M's)
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
- blobtools map2cov (previously bam2cov):
    - supports BAM, SAM and CAS files (uses parsers in BtIO)
- blobtools blobplot
    - changed behaviour of "--multiplot", now plots each group separately, then a final plot all together.
    - new flag "--cumulative". Previous "--multiplot" behaviour: incremental addition of groups
    - new flag "--legend". Plots legend in a separate figure, useful for slides.
    - fixed bug causing different plots to share filenames, causing overwriting
    - changed scaling of blobs in plots:
        - previously, blobs had an area of (length/1000*65) pixels^2. But if an assembly had a few very big contigs, these covered the smaller ones.
        - now the scaling is done in the following way:
            - the biggest contig gets plotted as a blob of area 12500 pixels^2
            - all other contigs get scaled accordingly
        - reference scale was changed as well:
            - scale shows area of 0.05, 0.1 and 0.25 of longest contig

