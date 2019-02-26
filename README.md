BlobTools v1.1
===============================
A modular command-line solution for visualisation, quality control and taxonomic partitioning of genome datasets

- Discussions, questions and answers: [BlobTools GoogleGroup](https://groups.google.com/forum/#!forum/blobtools)
- Issues, bug reports and feature requests: [GitHub issues](https://github.com/DRL/blobtools/issues)
- Documentation: [blobtools.readme.io](https://blobtools.readme.io)
- Citation: [Laetsch DR and Blaxter ML, 2017](https://f1000research.com/articles/6-1287/v1)

![](https://github.com/DRL/blobtools/blob/master/example/blobplot.png)

Obtaining BlobTools
------------
- **Option A**: Download latest [release](https://github.com/DRL/blobtools/releases/latest)
- **Option B**: Clone repository
  ```
  git clone https://github.com/DRL/blobtools.git
  ```

Entering directory
------------
  ```
  cd blobtools
  ```

Install dependencies
------------
- **Option A**: Create [Conda](https://conda.io/en/latest/miniconda.html) environment

  ```
  conda create -n blobtools
  conda activate blobtools
  conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
  conda install -c bioconda pysam --update-deps
  ```
  *Tip*: Check if samtools exists by executing the command 'samtools' in the commandline. If samtools complains about dependencies, simply run the pysam install twice.

- **Option B**: Install dependencies via PIP
  ```
  python setup.py install --user
  ```

Download NCBI taxdump and create nodesdb
------------
  ```
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
  tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
  ./blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp
  ```

Create blobplot
------------
  ```
  ./blobtools create -i example/assembly.fna -b example/mapping_1.sorted.bam -t example/blast.out -o example/test && \
  ./blobtools view -i example/test.blobDB.json && \
  ./blobtools plot -i example/test.blobDB.json
  ```
Usage
-----
```
    ./blobtools --help
```

Docker
------

A docker container can be build using the following command:
```
     docker build -t drl/blobtools .
```
This docker image can be run with sample data as follows:
```
     docker run -v $PWD/example:/example/  -t  drl/blobtools ./blobtools create -i /example/assembly.fna -b /example/mapping_1.sorted.bam -t /example/blast.out -o /example/test
```
