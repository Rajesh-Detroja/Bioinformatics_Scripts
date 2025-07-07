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

## Running HaplotypeCaller
for i in *.rg.marked.bam ; do gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" HaplotypeCaller -I $i \
-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
-O ${i%??????????????}.g.vcf \
-ERC GVCF \
; done

## Consolidating GVCFs
gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" GenomicsDBImport \
-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
-V TS_105.g.vcf \
-V TS_LID.g.vcf \
-V AC3_WB.g.vcf \
-L /home/morgensternlab/detrojar/db/hg38/hg38_chr.bed \
--genomicsdb-workspace-path gvcfs_db

# Joint calling with GenotypeGVCFs
gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" GenotypeGVCFs \
-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
-V gendb://gvcfs_db \
-O raw.SNPs.vcf

## [Optional] View comblined gVCF files 
#gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" SelectVariants \
#-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
#-V gendb://gvcfs_db \
#-O combined.g.vcf

## VariantRecalibrator builds the Gaussian mixture model
gatk --java-options "-Xmx500G -Djava.io.tmpdir=`pwd`" VariantRecalibrator \
-V raw.SNPs.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $know_sites_path/hapmap_3.3.hg38.vcf \
-resource:omni,known=false,training=true,truth=false,prior=12.0 $know_sites_path/1000G_omni2.5.hg38.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 $know_sites_path/1000G_phase1.snps.high_confidence.hg38.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $know_sites_path/dbsnp_138.hg38.vcf \
-an DP \
-an QD \
-an FS \
-an MQRankSum \
-mode SNP \
-max-gaussians 4 \
-O raw.SNPs.recal \
-tranches-file raw.SNPs.tranches \
-rscript-file recal.plots.R

## ApplyVQSR applies the filtering threshold
gatk ApplyVQSR \
-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
-V raw.SNPs.vcf \
-mode SNP \
-recal-file raw.SNPs.recal \
-tranches-file raw.SNPs.tranches \
-O recal.SNPs.vcf \
-ts-filter-level 99.0

## Genotype Refinement
#gatk CalculateGenotypePosteriors \
#-R reference.fasta \
#-V input.vcf \
#-ped family.ped \
#-supporting population.vcf \
#-O output.vcf

## Filter low confidence GQs
#gatk VariantFiltration \
#-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
#-V recal.SNPs.vcf \
#-genotype-filter-expression "GQ<20" \
#-genotype-filter-name "lowGQ" \
#-O recal.SNPs.filtered.vcf

#gatk VariantAnnotator \
#-R /home/morgensternlab/detrojar/db/hg38/index/hg38.fa \
#-V recal.SNPs.filtered.vcf \
#-A PossibleDeNovo \
#-O recal.SNPs.filtered.ann.vcf

## Variant Level Callset Evaluation
#gatk VariantEval \
#-R reference.b37.fasta \
#-eval callset.vcf \
#--comp truthset.vcf \
#-o results.eval.grp
