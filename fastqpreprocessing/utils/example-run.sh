#!/bin/bash
rm -rf *.bam

src/fastqprocess --verbose \
 --bam-size 0.001 \
 --barcode-length 16 \
 --umi-length 10 \
 --sample-id TEST \
 --white-list ../L8TX/737K-august-2016.txt \
 --I1 ../L8TX/A_I1.fastq.gz \
 --I1 ../L8TX/B_I1.fastq.gz \
 --R1 ../L8TX/A_R1.fastq.gz \
 --R1 ../L8TX/B_R1.fastq.gz \
 --R2 ../L8TX/A_R2.fastq.gz \
 --R2 ../L8TX/B_R2.fastq.gz




# 1 sample
# --I1 ../L8TX/L8TX_171026_01_A04_I1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_A04_R1.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_A04_R2.fastq.gz 


# 7 samples
# --I1 ../L8TX/L8TX_171026_01_A04_I1.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_B04_I1.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_F03_I1.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_G03_I1.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_H03_I1.fastq.gz \
# --I1 ../L8TX/L8TX_180221_01_E12_I1.fastq.gz \
# --I1 ../L8TX/L8TX_180221_01_F12_I1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_A04_R1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_B04_R1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_F03_R1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_G03_R1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_H03_R1.fastq.gz \
# --R1 ../L8TX/L8TX_180221_01_E12_R1.fastq.gz \
# --R1 ../L8TX/L8TX_180221_01_F12_R1.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_A04_R2.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_B04_R2.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_F03_R2.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_G03_R2.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_H03_R2.fastq.gz \
# --R2 ../L8TX/L8TX_180221_01_E12_R2.fastq.gz \
# --R2 ../L8TX/L8TX_180221_01_F12_R2.fastq.gz 

# full set
# --I1 ../L8TX/L8TX_171026_01_A04_I1.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_B04_I1.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_F03_I1.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_G03_I1.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_H03_I1.fastq.gz \
# --I1 ../L8TX/L8TX_180221_01_E12_I1.fastq.gz \
# --I1 ../L8TX/L8TX_180221_01_F12_I1.fastq.gz \
# --I1 ../L8TX/L8TX_180221_01_G12_I1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_A04_R1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_B04_R1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_F03_R1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_G03_R1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_H03_R1.fastq.gz \
# --R1 ../L8TX/L8TX_180221_01_E12_R1.fastq.gz \
# --R1 ../L8TX/L8TX_180221_01_F12_R1.fastq.gz \
# --R1 ../L8TX/L8TX_180221_01_G12_R1.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_A04_R2.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_B04_R2.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_F03_R2.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_G03_R2.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_H03_R2.fastq.gz \
# --R2 ../L8TX/L8TX_180221_01_E12_R2.fastq.gz \
# --R2 ../L8TX/L8TX_180221_01_F12_R2.fastq.gz \
# --R2 ../L8TX/L8TX_180221_01_G12_R2.fastq.gz 

# 
# --I1 ../L8TX/L8TX_171026_01_A04_I1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_A04_R1.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_A04_R2.fastq.gz \
# --I1 ../L8TX/L8TX_171026_01_B04_I1.fastq.gz \
# --R1 ../L8TX/L8TX_171026_01_B04_R1.fastq.gz \
# --R2 ../L8TX/L8TX_171026_01_B04_R2.fastq.gz 
