#!/bin/bash

./fastqprocess --verbose \
 --bam-size 0.001 \
 --barcode-length 16 \
 --umi-length 10 \
 --sample-id L8TX \
 --white-list ../../../data/L8TX/737K-august-2016.txt \
 --I1 ../../../data/L8TX/A_I1.fastq.gz \
 --R1 ../../../data/L8TX/A_R1.fastq.gz \
 --R2 ../../../data/L8TX/A_R2.fastq.gz \
 --I1 ../../../data/L8TX/B_I1.fastq.gz \
 --R1 ../../../data/L8TX/B_R1.fastq.gz \
 --R2 ../../../data/L8TX/B_R2.fastq.gz \
