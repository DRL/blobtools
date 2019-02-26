FROM continuumio/miniconda3
MAINTAINER Nick Waters <nickp60@gmail.com>
RUN conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
RUN conda install -c bioconda pysam --update-deps
RUN git clone https://github.com/DRL/blobtools.git
WORKDIR blobtools

RUN ./blobtools -h
# RUN ./blobtools create -i example/assembly.fna -b example/mapping_1.bam -t example/blast.out -o example/test
RUN wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
RUN tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
RUN ./blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp
