FROM python:3.7.7

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="python 3.7.7" \
  description="python 3.7.7 with pysam, sctools, requests, and a basic science stack"

COPY requirements.txt .

RUN apt-get update && apt-get install -y patch && apt-get install -y libhdf5-dev && apt-get install -y vim

RUN pip3 install -r requirements.txt

RUN mkdir /sctools/

COPY .  /sctools 

RUN  pip3 install /sctools

ARG libStatGen_version="1.0.14"


RUN cd /sctools/fastqpreprocessing &&\
    wget https://github.com/statgen/libStatGen/archive/v${libStatGen_version}.tar.gz &&\
    tar -zxvf v${libStatGen_version}.tar.gz &&\
    mv libStatGen-${libStatGen_version} libStatGen &&\
    patch libStatGen/fastq/FastQFile.cpp patches/FastQFile.cpp.patch &&\
    patch libStatGen/Makefile patches/Makefile.patch &&\
    patch libStatGen/general/Makefile patches/general.Makefile.patch &&\
    make -C libStatGen &&\
    mkdir src/obj &&\
    make -C src/ 

RUN cp /sctools/fastqpreprocessing/src/fastqprocess /usr/local/bin/

WORKDIR usr/local/bin/sctools


