zcat ../../L8TX/L8TX_180221_01_F12_R2.fastq.gz | head -n 4000000 > a_R1.fastq

gzip a_R1.fastq

cp a_R1.fastq.gz  b_R2.fastq.gz
cp a_R1.fastq.gz  b_I1.fastq.gz
cp a_R1.fastq.gz  b_R1.fastq.gz

cp a_R1.fastq.gz  a_R2.fastq.gz
cp a_R1.fastq.gz  a_I1.fastq.gz
