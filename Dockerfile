FROM python:3.7.7

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="python 3.7.7" \
  description="python 3.7.7 with pysam, sctools, requests, and a basic science stack"

COPY requirements.txt .
RUN pip3 install -r requirements.txt

RUN mkdir /sctools/

COPY .  /sctools 

RUN  pip3 install /sctools

WORKDIR usr/local/bin/sctools


