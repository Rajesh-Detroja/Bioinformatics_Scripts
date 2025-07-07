## Add Read Group and sort
for i in *.Aligned.out.sam.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" AddOrReplaceReadGroups -I=$i \
-O=${i%????????????????????}.rg.bam \
-RGLB=lib1 \
-RGPL=illumina \
-RGPU=unit1 \
-RGSM=${i%????????????????????} \
-SO=coordinate \
; done


## Mark Duplicate
for i in *.rg.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" MarkDuplicates -I=$i \
-O=${i%????}.marked.bam \
-M=${i%????}.marked.txt \
-CREATE_INDEX=true \
-VALIDATION_STRINGENCY=SILENT \
; done


## Split'N'Trim and reassign mapping qualities (Learn about it!)
for i in *.rg.marked.bam ; do gatk --java-options "-Xmx500G" SplitNCigarReads -I=$i \
-O=${i%????}.split.bam \
-R=/home/morgensternlab/detrojar/db/hg38/index/star/hg38.fa \
-tmp-dir=. \
; done


## Base Recalibration
know_sites_path="/data/morgensternlab/detrojar/db/hg38/known_sites/hg38_changed_chr/filtered"
for i in *.rg.marked.split.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" BaseRecalibrator -I=$i \
-R=/data/morgensternlab/detrojar/db/hg38/index/star/hg38.fa \
-known-sites=$know_sites_path/1000G_phase1.snps.high_confidence.hg38.vcf \
-known-sites=$know_sites_path/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-known-sites=$know_sites_path/dbsnp_138.hg38.vcf \
-O=${i%????}.recal-data.table \
; done

## Index final BAM file
for i in *.rg.marked.split.bam ; do samtools index $i ; done

## Variant calling
for i in *.rg.marked.split.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" HaplotypeCaller -I=$i \
-R=/data/morgensternlab/detrojar/db/hg38/index/star/hg38.fa \
-O=${i%????????????????????}.vcf \
-dont-use-soft-clipped-bases=true \
-stand-call-conf=20.0 \
; done


## Variant filtering
for i in *.vcf ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" VariantFiltration -V=$i \
-R=/data/morgensternlab/detrojar/db/hg38/index/star/hg38.fa \
-O=${i%????????????????????}.filtered.vcf \
-window=30 \
-cluster=3 \
-filter-name="FS" \
-filter-expression="FS > 30.0" \
-filter-name="QD" \
-filter-expression="QD < 2.0" \
-filter-name="DP" \
-filter-expression="DP > 2" \
; done

## Final Variant filtering
for i in *.filtered.vcf ; do grep -e "^#CHROM" -e "DP;FS;QD;SnpCluster" $i > ${i%????}.final.vcf ; done
for i in *.filtered.final.vcf ; do cut -f1,2,3,4,5,6 $i > ${i%????}.tsv ; done

## Annotate VCF
mkdir vep

for i in *.filtered.final.tsv ; do vep -i $i --tab --everything --cache --force_overwrite --filter_common --fork 30 ; mkdir vep/${i%???????????????????????} ; mv variant_effect_output* vep/${i%???????????????????????}/ ; done
