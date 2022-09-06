FROM python:3.7.7

LABEL maintainer="Farzaneh Khajouei <fkhajoue@broadinstitute.org>" \
  software="sctools  v.1.0.0" \
  description="A collection of tools for single cell data. Splitting fastq files based on cellbarcodes and other tools to compute metrics on single cell data using barcodes and UMIs."


RUN apt-get update && apt-get upgrade -y && apt-get install -y patch libhdf5-dev vim apt-utils
RUN mkdir /sctools/

COPY . /sctools 

ARG htslib_version="1.13"

RUN cd /sctools/fastqpreprocessing &&\
    wget https://github.com/khajoue2/libStatGen/archive/refs/tags/v1.0.15.broad.tar.gz &&\
    wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 &&\
    tar -zxvf v1.0.15.broad.tar.gz &&\
    tar -jxvf htslib-${htslib_version}.tar.bz2 &&\
    mv libStatGen-1.0.15.broad libStatGen 

RUN cd /sctools/fastqpreprocessing &&\
    wget http://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz &&\
    tar -xvf gzstream.tgz 

RUN cd /sctools/fastqpreprocessing &&\
    make -C libStatGen 

RUN cd /sctools/fastqpreprocessing && make -C htslib-${htslib_version}/ && make -C gzstream

RUN cd /sctools/fastqpreprocessing && mkdir bin obj && make install

RUN cp /sctools/fastqpreprocessing/bin/* /usr/local/bin/

WORKDIR usr/local/bin/sctools

