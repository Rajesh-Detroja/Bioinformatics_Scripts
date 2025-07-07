#!/usr/bin

while getopts ":1:c:" opt; do
  case $opt in
    1) read_1="$OPTARG"
    ;;
    c) config="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

## Source config file
source $config

## HG38 MAPPING
## ============
hg38=$human_index

R1=(*$read_1)

suffix="$read_1"

for ((i=0;i<${#R1[@]};i++)); do
        sample=${R1[i]%$suffix}
#       (bowtie2 -N 1 -p $threads -x $hg38 -U ${R1[i]} -S $sample.hg38.sam)2>$sample.hg38.log
#       (bwa mem -t $threads $hg38 ${R1[i]} > $sample.hg38.sam)2>$sample.hg38.log
        (STAR --genomeDir $hg38 --readFilesIn ${R1[i]} --outSAMunmapped Within --outFileNamePrefix $sample.hg38. --runThreadN $threads)2>$sample.hg38.log
done

## Options - STAR
## --readFilesCommand zcat

## Rename suffix for "STAR"
for i in *.hg38.Aligned.out.sam ; do mv $i ${i%????????????????}.sam ; done

## SAM to BAM
for i in *.hg38.sam ; do samtools view -Sb $i > $i.bam ; done

for i in *.hg38.sam.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.hg38.sam.bam ; do samtools sort -T . -@ $threads -O bam -o $i.sorted.bam $i ; done

## Extract unmapped reads
for i in *.hg38.sam.bam.sorted.bam ; do samtools view -u -b -f 4 $i > $i.unmapped.bam ; done

for i in *.hg38.sam.bam.sorted.bam.unmapped.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.hg38.sam.bam.sorted.bam.unmapped.bam ; do samtools index $i -@ $threads ; done

for i in *.hg38.sam.bam.sorted.bam.unmapped.bam ; do bam2fastq -q $i -o $i#.fq ; done



## CHIMERA MAPPING
## ===============

X=$chimera_star_index

R1=(*.hg38.sam.bam.sorted.bam.unmapped.bam_M.fq)

suffix=".hg38.sam.bam.sorted.bam.unmapped.bam_M.fq"

for ((i=0;i<${#R1[@]};i++)); do
	sample=${R1[i]%$suffix}
        (STAR --genomeDir $X --readFilesIn ${R1[i]} --outFilterMismatchNmax 2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --seedSearchStartLmax 20 --outSAMunmapped Within --outFileNamePrefix $sample.chimeras. --runThreadN $threads)2>$sample.chimera.log
done

for i in *.chimeras.Aligned.out.sam ; do samtools view -Sb $i > $i.bam ; done

for i in *.chimeras.Aligned.out.sam.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.chimeras.Aligned.out.sam.bam ; do samtools sort -T . -@ $threads -O bam -o $i.sorted.bam $i ; done

for i in *.chimeras.Aligned.out.sam.bam.sorted.bam ; do samtools index $i ; done

for i in *.chimeras.Aligned.out.sam.bam.sorted.bam ; do bedtools coverage -f 1.00 -a $chimera_bed -b $i | awk -v depth=$depth '($7 >= 1 && $4 >= depth)' | awk 'BEGIN { OFS="\t"; print "chimera", "start", "end", "depth", "covered_bp", "Junction_length", "coverage" } { print $0, "" }' > $i.coverage ; done

for i in *.chimeras.Aligned.out.sam.bam.sorted.bam.coverage ; do if [[ $(wc -l < $i ) -lt 2 ]] ; then awk 'BEGIN { OFS="\t"; print "NO_Chimera", "NA", "NA", "NA", "NA", "NA", "NA" }' >> $i ; fi ; done

for i in *.chimeras.Aligned.out.sam.bam.sorted.bam.coverage ; do awk 'NR>1 {print $1,"\t",$4}' $i > ${i%?????????}.junction.coverage ; done

for i in *.chimeras.Aligned.out.sam.bam.sorted.bam.junction.coverage ; do file=TheFileName ; sed -i '1{h;s/.*/chimera '"${i%??????????????????????????????????????????????????????????}"'/;G}' "$i" ; done

merge -k -e "0" *.chimeras.Aligned.out.sam.bam.sorted.bam.junction.coverage > chimeras.tsv

head -n 1 chimeras.tsv | awk '{$1=$1 "\t" "chimera_type" "\t" "gene1" "\t" "strand1" "\t" "gene2" "\t" "strand2" "\t" "junction_id" "\t" "Length"}1' OFS="\t" > all_chimeras.txt

join $chimera_ANN <(sort -k1,1 chimeras.tsv) | awk -v OFS="\t" '$1=$1' | sort -k2,2 >> all_chimeras.txt

grep -e "^chimera" -e "exon-exon" all_chimeras.txt | awk 'NR == 1; NR > 1 {print $0 | "sort -k9nr"}' > exon_exon_chimeras.txt

## get junction and chimera sequences
cut -f7 exon_exon_chimeras.txt | grep -v "^junction_id" | sed 's/:/\t/g' | sed 's/-/\t/g' | bedtools getfasta -fi $chimera_fa -bed - > exon_exon_chimeras_junctions.fa

cut -f7 exon_exon_chimeras.txt | grep -v "^junction_id" | cut -d":" -f1 | while read id ; do samtools faidx $chimera_fa ${id} >> exon_exon_chimeras_sequence.fa ; done


## Predict ORF / Protein Sequence

chimera_out=exon_exon_chimeras_sequence.fa

date && time getorf -sequence $chimera_out -outseq ${chimera_out%???}".fnn" -minsize 75 -find 1

perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${chimera_out%???}".fnn" > $chimera_out".1_line.fnn"

cut -f7 exon_exon_chimeras.txt | grep -v "^junction_id" | cut -d":" -f1 | while read id ; do sed -n -e '/^>'${id}'_/{$!N;p;}' -e h $chimera_out".1_line.fnn" > ${id}".fnn" ; done

cut -f7 exon_exon_chimeras.txt | grep -v "^junction_id" | cut -d":" -f1 | while read id ; do awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${id}".fnn" | awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' | sort -k1,1nr | cut -f 2- | tr "\t" "\n" | head -n 2 > ${id}".1.fnn" ; done

cat *.1.fnn | sed '/^$/d' > exon_exon_chimeras_protein.fnn

cut -f7 exon_exon_chimeras.txt | grep -v "^junction_id" | cut -d":" -f1 | while read id ; do rm -rf ${id}".fnn" ; done
rm -rf *.1_line.fnn rm -rf *.1.fnn

#Rscript ~/scripts/RPKM_TPM.R chimeras.txt

rm -rf *.hg38.log *.hg38.sam *.hg38.sam.bam *.chimeras.Aligned.out.sam *.chimeras.Aligned.out.sam.bam *.chimeras.Aligned.out.sam.bam.sorted.bam.coverage *.chimeras.Aligned.out.sam.bam.sorted.bam.junction.coverage *.chimera.log *.Log.* *.out.tab chimeras.tsv exon_exon_chimeras_sequence.fnn
