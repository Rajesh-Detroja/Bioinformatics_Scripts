## Discover known microbiomes
## References: https://gatkforums.broadinstitute.org/gatk/discussion/13502/no-outputs-from-pathseq

#!/usr/bin

PATHSEQ="/home/morgensternlab/detrojar/db/pathseq"

R1=(*.trimmed.sam.bam.sorted.bam.unmapped.bam)
suffix=".trimmed.sam.bam.sorted.bam.unmapped.bam"

for ((i=0;i<${#R1[@]};i++)); do

	sample=${R1[i]%$suffix}

        (gatk --java-options "-Xmx500G" PathSeqPipelineSpark  \
              --input ${R1[i]} \
              --kmer-file $PATHSEQ/pathseq_host.bfi \
              --filter-bwa-image $PATHSEQ/pathseq_host.fa.img \
              --microbe-bwa-image $PATHSEQ/pathseq_microbe.fa.img \
              --microbe-fasta $PATHSEQ/pathseq_microbe.fa \
              --taxonomy-file $PATHSEQ/pathseq_taxonomy.db \
              --min-clipped-read-length 31 \
              --min-score-identity 0.80 \
              --identity-margin 0.02 \
	      --disable-read-filter WellformedReadFilter \
              --scores-output $sample.scores.txt \
              --output $sample.output.bam \
             --filter-metrics $sample.filter_metrics.txt \
             --score-metrics $sample.score_metrics.txt
       )

#for i in *.scores.txt ; do cat $i | awk '($10 > 0 && $9 >= 20)' > $i.filtered.txt ; done

done

## Discover novel microbiomes
## samtools view -c output.pathseq.bam > output.pathseq.novel.txt
