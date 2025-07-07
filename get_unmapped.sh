## R1 unmapped, R2 mapped
for i in *.Aligned.out.sam.bam.sorted.bam ; do samtools view -u -b -f 4 -F 264 $i > $i.L.bam ; done

## R1 mapped, R2 unmapped
for i in *.Aligned.out.sam.bam.sorted.bam ; do samtools view -u -b -f 8 -F 260 $i > $i.R.bam ; done

## R1 & R2 unmapped
for i in *.Aligned.out.sam.bam.sorted.bam ; do samtools view -u -b -f 12 -F 256 $i > $i.B.bam ; done

## Merge unmapped reads
for i in *.Aligned.out.sam.bam.sorted.bam ; do samtools merge $i.unmapped.bam -u $i.L.bam $i.R.bam $i.B.bam ; done

for i in *.Aligned.out.sam.bam.sorted.bam.unmapped.bam ; do bamtools-2.4.1 stats -in $i > $i.stats.txt ; done

for i in *.Aligned.out.sam.bam.sorted.bam.unmapped.bam ; do samtools index $i -@ 30 ; done

for i in *.Aligned.out.sam.bam.sorted.bam.unmapped.bam ; do bam2fastq -q $i -o $i#.fq ; done

rm -rf *.L.bam *.R.bam *.B.bam
