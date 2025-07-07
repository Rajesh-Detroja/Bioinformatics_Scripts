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
        (bwa mem -t 24 $hg38 ${R1[i]} ${R2[i]} > $sample.sam)2>$sample.log
done

for i in *.sam ; do samtools view -Sb $i > $i.bam ; done

for i in *.sam.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.sam.bam ; do samtools sort -T . -@ 24 -O bam -o $i.sorted.bam $i ; done

for i in *.sam.bam.sorted.bam ; do samtools index $i -@ 24 ; done
