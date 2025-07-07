#!/usr/bin

for i in *.sam.bam.sorted.bam ; do bamCoverage -p 30 --ignoreDuplicates -of bedgraph --blackListFileName ~/db/hg38/blacklist/hg38.blacklist.bed -b $i -o $i.bedgraph ; done

#for i in *.sam.bam.marked.bam.bedgraph ; do cat $i | awk '$4 >= 5' | sed 's/^/chr/' | sed 's/chrMT/chrM/' > $i.filtered.bedgraph ; done
