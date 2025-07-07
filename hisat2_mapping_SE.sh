#!/usr/bin

while getopts ":1:" opt; do
  case $opt in
    1) read_1="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


## Reference MAPPING
## =================
ref=~/db/Mouse/Mus_musculus_GRCm38.fa 

R1=(*$read_1)

suffix="$read_1"

for ((i=0;i<${#R1[@]};i++)); do
        sample=${R1[i]%$suffix}
	(hisat2 -p 30 -x $ref -U ${R1[i]} -S $sample.sam)2>$sample.log
done

for i in *.sam ; do samtools view -Sb $i > $i.bam ; done

for i in *.sam.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.sam.bam ; do samtools sort -T . -@ 30 -O bam -o $i.sorted.bam $i ; done

rm -rf *.sam


## featureCounts
## =============
featureCounts -O -T 30 -p -t exon -g gene_id -a ~/db/Mouse/Mus_musculus_GRCm38.gtf -o counts.txt *.sam.bam.sorted.bam
