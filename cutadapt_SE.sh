#!/usr/bin

R1=(*.fastq.gz)

suffix=".fastq.gz"

for ((i=0;i<${#R1[@]};i++)); do
	sample=${R1[i]%$suffix}
	cutadapt --trim-n -m 35 -u 12 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC -o $sample.trimmed.fastq.gz ${R1[i]} > $sample.trimming.log
done
