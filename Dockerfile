FROM python:3.6.2

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="python 3.6.2" \
  description="python 3.6.2 with pysam, sctools, requests, and a basic science stack"

COPY requirements.txt .
RUN pip3 install -r requirements.txt

WORKDIR usr/local/bin/sctools

COPY . .