FROM python:3.7.7

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="python 3.7.7" \
  description="python 3.7.7 with pysam, sctools, requests, and a basic science stack"

COPY requirements.txt .

RUN apt-get update && apt-get install -y patch && apt-get install -y libhdf5-dev && apt-get install -y vim
RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt

RUN mkdir /sctools/

COPY . /sctools 

RUN  pip3 install /sctools

ARG libStatGen_version="1.0.14"
ARG htslib_version="1.13"

RUN cd /sctools/fastqpreprocessing &&\
    wget https://github.com/statgen/libStatGen/archive/v${libStatGen_version}.tar.gz &&\
    wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 &&\
    tar -zxvf v${libStatGen_version}.tar.gz &&\
    tar -jxvf htslib-${htslib_version}.tar.bz2 &&\
    mv libStatGen-${libStatGen_version} libStatGen 

RUN cd /sctools/fastqpreprocessing &&\
    wget http://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz &&\
    tar -zxvf gzstream.tgz 

RUN cd /sctools/fastqpreprocessing &&\
    patch -f libStatGen/fastq/FastQFile.cpp patches/FastQFile.cpp.patch &&\
    patch -f libStatGen/general/BgzfFileType.cpp patches/BgzfFileType.cpp.patch &&\  
    patch libStatGen/Makefile patches/Makefile.patch &&\
    patch libStatGen/general/Makefile patches/general.Makefile.patch &&\
    make -C libStatGen 

RUN cd /sctools/fastqpreprocessing && make -C htslib-${htslib_version}/ && make -C gzstream

RUN cd /sctools/fastqpreprocessing && mkdir bin src/obj  &&  make -C src/ install

RUN cp /sctools/fastqpreprocessing/bin/* /usr/local/bin/

WORKDIR usr/local/bin/sctools


