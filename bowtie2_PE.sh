#!/usr/bin

while getopts ":1:2:" opt; do
  case $opt in
    1) read_1="$OPTARG"
    ;;
    2) read_2="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

## HG38 MAPPING
## ============
hg38=~/db/hg38/index/hg38.fa

R1=(*$read_1)
R2=(*$read_2)

suffix="$read_1"

for ((i=0;i<${#R1[@]};i++)); do
        sample=${R1[i]%$suffix}
        (bowtie2 -N 1 -p 24 -x $hg38 -1 ${R1[i]} -2 ${R2[i]} -S $sample.hg38.sam)2>$sample.hg38.log
done

for i in *.hg38.sam ; do samtools view -Sb $i > $i.bam ; done

for i in *.hg38.sam.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.hg38.sam.bam ; do samtools sort -T . -@ 24 -O bam -o $i.sorted.bam $i ; done

## R1 unmapped, R2 mapped
for i in *.hg38.sam.bam.sorted.bam ; do samtools view -u -b -f 4 -F 264 $i > $i.L.bam ; done

## R1 mapped, R2 unmapped
for i in *.hg38.sam.bam.sorted.bam ; do samtools view -u -b -f 8 -F 260 $i > $i.R.bam ; done

## R1 & R2 unmapped
for i in *.hg38.sam.bam.sorted.bam ; do samtools view -u -b -f 12 -F 256 $i > $i.B.bam ; done

## Merge unmapped reads
for i in *.hg38.sam.bam.sorted.bam ; do samtools merge $i.unmapped.bam -u $i.L.bam $i.R.bam $i.B.bam ; done

for i in *.hg38.sam.bam.sorted.bam.unmapped.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.hg38.sam.bam.sorted.bam.unmapped.bam ; do samtools index $i -@ 24 ; done

for i in *.hg38.sam.bam.sorted.bam.unmapped.bam ; do bam2fastq -q $i -o $i#.fq ; done
