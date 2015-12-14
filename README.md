> Once we were blobs in the sea...

> -Terry Pratchett, A Hat Full of Sky

# blobtools
Application for the visualisation of (draft) genome assemblies and general assembly QC using TAGC (Taxon-annotated Gc-Coverage) plots [Kumar et al. 2012](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3843372/pdf/fgene-04-00237.pdf).

## Requirements

```
- Python 2.7+
- Matplotlib 1.5 
- Docopt 
- NCBI Taxonomy (names.dmp and nodes.dmp) ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
- Virtualenv (recommended), see http://docs.python-guide.org/en/latest/dev/virtualenvs/
```

## Installation 
- Recommended 
```
# install virtualenv 
pip install virtualenv

# clone blobtools into folder
git clone https://github.com/DRL/blobtools.git

# create virtual environment for blobtools  
cd blobtools/
virtualenv blob_env

# activate virtual environment
source blob_env/bin/activate

# install matplotlib
(blob_env) $ pip install matplotlib

# install docopt
(blob_env) $ pip install docopt 

# run
(blob_env) $ ./blobtools -h
```
- Basic
```
# clone blobtools into folder
$ git clone https://github.com/DRL/blobtools.git

# install matplotlib
$ pip install matplotlib

# install docopt
$ pip install docopt

# run
$ ./blobtools -h
```

## Doc 
### blobtools 
- main executable
```
usage: blobtools <command> [<args>...] [--help]

commands:
  create    create a BlobDB
  view      print BlobDB
  plot      plot BlobDB as a blobplot

-h --help    show this
```

### blobtools create 
- create a BlobDb JSON file
```
usage: blobtools create     -i FASTA [-y FASTATYPE] [-o OUTFILE] [--title TITLE]
                              [-b BAM...] [-s SAM...] [-a CAS...] [-c COV...]  
                              [--nodes <NODES>] [--names <NAMES>] [--db <NODESDB>] 
                              [-t TAX...] [-r TAXRULE...]
                              [-h|--help] 
    
    Options:
        -h --help                   show this
        -i, --infile FASTA          FASTA file of assembly. Headers are split at whitespaces.  
        -y, --type FASTATYPE        Assembly program used to create FASTA. If specified, 
                                    coverage will be parsed from FASTA header. 
                                    (Parsing supported for 'spades', 'soap', 'velvet', 'abyss')
        -t, --taxfile TAX...        Taxonomy file in format (qseqid\ttaxid\tbitscore) 
                                    (e.g. BLAST output "--outfmt '6 qseqid staxids bitscore'")
        -x, --taxrule <TAXRULE>...  Taxrule determines how taxonomy of blobs is computed [default: bestsum]
                                    "bestsum"       : sum bitscore across all hits for each taxonomic rank
                                    "bestsumorder"  : sum bitscore across all hits for each taxonomic rank. 
                                                  - If first <TAX> file supplies hits, bestsum is calculated. 
                                                  - If no hit is found, the next <TAX> file is used.                                 
        --nodes <NODES>             NCBI nodes.dmp file. Not required if '--db'
        --names <NAMES>             NCBI names.dmp file. Not required if '--db' 
        --db <NODESDB>              NodesDB file [default: data/nodesDB.txt]. 
        -b, --bam <BAM>...          BAM file (requires samtools in $PATH)
        -s, --sam <SAM>...          SAM file
        -a, --cas <CAS>...          CAS file (requires clc_mapping_info in $PATH)
        -c, --cov <COV>...          TAB separated. (seqID\tcoverage)
        -o, --out <OUT>             BlobDB output prefix 
        --title TITLE               Title of BlobDB [default: FASTA) 
```

### blobtools view 
- generate table output from a blobDB file
```
usage: blobtools view    -i <BLOBDB> [-r <TAXRULE>] [--rank <TAXRANK>...] [--hits]
                            [--list <LIST>] [--out <OUT>]
                            [--h|--help] 
    
    Options:
        --h --help                  show this
        -i, --input <BLOBDB>        BlobDB file (created with "blobtools forge")
        -o, --out <OUT>             Output file [default: STDOUT]
        -l, --list <LIST>           List of sequence names (comma-separated or file). 
                                    If comma-separated, no whitespaces allowed.
        -r, --taxrule <TAXRULE>     Taxrule used for computing taxonomy (supported: "bestsum", "bestsumorder")
                                    [default: bestsum]
        --rank <TAXRANK>...         Taxonomic rank(s) at which output will be written. 
                                    (supported: 'species', 'genus', 'family', 'order', 
                                    'phylum', 'superkingdom', 'all') [default: phylum]
        -b, --hits                  Displays taxonomic hits from tax files
```

### blobtools plot 
- generate a blobplot from a blobDB file
```
usage: blobtools plot    -i BLOBDB [-p INT] [-l INT] [-c] [-n] [-s]
                            [-r RANK] [-x TAXRULE] [--label GROUPS...] 
                            [-o PREFIX] [-m] [--sort ORDER] [--hist HIST] [--title]
                            [--colours FILE] [--include FILE] [--exclude FILE]
                            [--format FORMAT] [--noblobs] [--noreads] [--refcov FILE]
                            [-h|--help] 

    Options:
        -h --help                   show this
        -i, --infile BLOBDB         BlobDB file
        -p, --plotgroups INT        Number of (taxonomic) groups to plot, remaining 
                                     groups are placed in 'other' [default: 7]
        -l, --length INT            Minimum sequence length considered for plotting [default: 100]
        -c, --cindex                Colour blobs by 'c index' [default: False]
        -n, --nohit                 Hide sequences without taxonomic annotation [default: False]
        -s, --noscale               Do not scale sequences by length [default: False]
        -o, --out PREFIX            Output prefix
        -m, --multiplot             Multi-plot. Print plot after addition of each (taxonomic) group 
                                     [default: False]
        --sort <ORDER>              Sort order for plotting [default: span]
                                     span  : plot with decreasing span
                                     count : plot with decreasing count 
        --hist <HIST>               Data for histograms [default: span] 
                                     span  : span-weighted histograms
                                     count : count histograms
        --title                     Add title of BlobDB to plot [default: False]
        -r, --rank RANK             Taxonomic rank used for colouring of blobs [default: phylum]
                                     (Supported: species, genus, family, order, phylum, superkingdom) 
        -x, --taxrule TAXRULE       Taxrule which has been used for computing taxonomy 
                                     (Supported: bestsum, bestsumorder) [default: bestsum]
        --label GROUPS...           Relabel (taxonomic) groups (not 'all' or 'other'), 
                                     e.g. "Bacteria=Actinobacteria,Proteobacteria"
        --colours COLOURFILE        File containing colours for (taxonomic) groups
        --exclude GROUPS..          Place these (taxonomic) groups in 'other',
                                     e.g. "Actinobacteria,Proteobacteria"
        --format FORMAT             Figure format for plot (png, pdf, eps, jpeg, 
                                        ps, svg, svgz, tiff) [default: png]
        --noblobs                   Omit blobplot [default: False]
        --noreads                   Omit plot of reads mapping [default: False]
        --refcov FILE               File containing number of "total" and "mapped" reads 
                                    per coverage file. (e.g.: bam0,900,100). If provided, info
                                    will be used in read coverage plot(s). 
```
## Additional features

### blobtools bam2cov
- extract base-coverage for each contig from BAM file
```
usage: blobtools bam2cov         -i FASTA -b BAM [-h|--help] 
    
    Options:
        -h --help                   show this
        -i, --infile FASTA          FASTA file of assembly. Headers are split at whitespaces.  
        -b, --bam <BAM>             BAM file (requires samtools in $PATH)
```
