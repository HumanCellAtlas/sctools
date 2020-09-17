FROM python:3.7.7

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="python 3.7.7" \
  description="python 3.7.7 with pysam, sctools, requests, and a basic science stack"

COPY requirements.txt .

RUN apt-get update && apt-get install -y patch && apt-get install -y libhdf5-dev

RUN pip3 install -r requirements.txt

RUN mkdir /sctools/

COPY .  /sctools 

RUN  pip3 install /sctools

ARG libStatGen_version="v1.0.14"

RUN wget https://github.com/HumanCellAtlas/sctools/archive/kmk-fastqprocessing.zip 

RUN unzip kmk-fastqprocessing.zip && \
    cd sctools-kmk-fastqprocessing/fastqpreprocessing &&\
    wget https://github.com/statgen/libStatGen/archive/${libStatGen_version}.tar.gz &&\
    tar -zxvf ${libStatGen_version}.tar.gz &&\
    patch libStatGen/fastq/FastQFile.cpp patches/FastQFile.cpp.patch &&\
    patch libStatGen/Makefile patches/Makefile.patch &&\
    make -C libStatGen &&\
    mkdir src/obj &&\
    make -C src/ 

RUN cp sctools-kmk-fastqprocessing/fastqpreprocessing/src/fastqprocess /usr/local/bin/

WORKDIR usr/local/bin/sctools


