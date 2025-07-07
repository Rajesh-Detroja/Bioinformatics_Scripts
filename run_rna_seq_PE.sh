## HG38 MAPPING
## ============

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


hg38=~/db/Rat/Rattus_norvegicus.fa
hg38_gtf=~/db/Rat/Rattus_norvegicus.gtf
R1=(*$read_1)
R2=(*$read_2)

suffix="$read_1"


for ((i=0;i<${#R1[@]};i++)); do
        sample=${R1[i]%$suffix}
        (hisat2 -p 30 --dta -x $hg38 -1 ${R1[i]} -2 ${R2[i]} -S $sample.sam)2>$sample.log
done

for i in *.sam ; do samtools view -Sb $i > $i.bam ; done

for i in *.sam.bam ; do samtools sort -T . -@ 30 -O bam -o $i.sorted.bam $i ; done

for i in *.sam.bam.sorted.bam ; do stringtie -p 30 -G $hg38_gtf -o ${i%???????????????????}.gtf -l ${i%???????????????????} $i ; done

stringtie --merge -p 30 -G $hg38_gtf -o stringtie_merged.gtf *.gtf

## Optional
#gffcompare –r $hg38_gtf –G –o merged stringtie_merged.gtf

for i in *.sam.bam.sorted.bam ; do stringtie -e -B -p 30 -G stringtie_merged.gtf -o ballgown/${i%???????????????????}/${i%???????????????????}.gtf $i ; done
