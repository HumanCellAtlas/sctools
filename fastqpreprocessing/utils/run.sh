#!/bin/bash


#  --tool=memcheck \
#   --leak-check=full  \
#  --log-file=valgrind-out.txt  \

valgrind  \
   --tool=massif \
    --time-unit=B \
    ./fastqproc a_R1.fastq.gz a_I1.fastq.gz a_R2.fastq.gz \
    b_R1.fastq.gz b_I1.fastq.gz b_R2.fastq.gz
