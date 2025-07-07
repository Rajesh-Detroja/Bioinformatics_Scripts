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
hg38=/home/morgensternlab/detrojar/db/hg38/index/hg38.fa

R1=(*$read_1)
R2=(*$read_2)

suffix="$read_1"

for ((i=0;i<${#R1[@]};i++)); do
        sample=${R1[i]%$suffix}
        (bwa mem -t 30 $hg38 ${R1[i]} ${R2[i]} > $sample.sam)2>$sample.log
done

for i in *.sam ; do samtools view -Sb $i > $i.bam ; done

for i in *.sam.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.sam.bam ; do samtools index $i ; done

# To activate this environment, use:
# > source activate gatk
#
# To deactivate an active environment, use:
# > source deactivate

## Add Read Group and sort
for i in *.sam.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" AddOrReplaceReadGroups -I $i \
-O ${i%????????}.rg.bam \
-RGLB lib1 \
-RGPL illumina \
-RGPU unit1 \
-RGSM ${i%????????} \
-SO coordinate \
; done

## Mark Duplicate
for i in *.rg.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" MarkDuplicates -I $i \
-O ${i%????}.marked.bam \
-M ${i%????}.marked.txt \
-CREATE_INDEX true \
-VALIDATION_STRINGENCY SILENT \
; done

## Base Recalibration
know_sites_path="/data/morgensternlab/detrojar/db/hg38/known_sites/hg38_changed_chr/filtered"
for i in *.rg.marked.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" BaseRecalibrator -I $i \
-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
-known-sites $know_sites_path/1000G_phase1.snps.high_confidence.hg38.vcf \
-known-sites $know_sites_path/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-known-sites $know_sites_path/dbsnp_138.hg38.vcf \
-O ${i%????}.recal-data.table \
; done

## Apply base quality score recalibration
for i in *.rg.marked.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" ApplyBQSR \
-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
-I $i \
-bqsr-recal-file ${i%????}.recal-data.table \
-O ${i%????}.recal.bam \
; done

## Index final BAM file
for i in *.rg.marked.recal.bam ; do samtools index $i ; done

## Running HaplotypeCaller
for i in *.rg.marked.recal.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" HaplotypeCaller -I $i \
-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
-O ${i%????????????????????}.vcf \
; done

## Activate GATKTools
source activate gatk

## Apply a Convolutional Neural Net to filter annotated variants
for i in *.rg.marked.recal.bam ; do gatk CNNScoreVariants -I $i \
-V ${i%????????????????????}.vcf \
-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
-O ${i%????????????????????}.ann.vcf \
-tensor-type read_tensor \
; done

## Apply tranche filtering
for i in *.ann.vcf ; do gatk FilterVariantTranches -V $i \
-resource $know_sites_path/hapmap_3.3.hg38.vcf \
-resource $know_sites_path/1000G_omni2.5.hg38.vcf \
-resource $know_sites_path/1000G_phase1.snps.high_confidence.hg38.vcf \
-resource $know_sites_path/dbsnp_138.hg38.vcf \
-resource $know_sites_path/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-info-key CNN_2D \
-snp-tranche 99.95 \
-indel-tranche 99.4 \
-invalidate-previous-filters \
-O ${i%????}.filtered.vcf \
; done

## Functional Annotator using VEP
for i in *.ann.filtered.vcf ; do cat $i | grep -v "^##" | cut -f1,2,3,4,5,6 > ${i%????}.tsv ; done

## Annotate VCF
mkdir vep

for i in *.ann.filtered.tsv ; do vep -i $i --tab --everything --cache --force_overwrite --filter_common --fork 30 ; mkdir vep/${i%?????????????????} ; mv variant_effect_output* vep/${i%?????????????????}/ ; done

## Filter CDS affects
for i in vep/*/*.txt ; do awk '$7 ~/stop_gained/ || $7 ~/start_lost/ || $7 ~/missense_variant/ || $7 ~/incomplete_terminal_codon_variant/ || $7 ~/synonymous_variant/ || $7 ~/coding_sequence_variant/ {print $0}' $i > ${i%????}.CDS.affect.txt ; done
