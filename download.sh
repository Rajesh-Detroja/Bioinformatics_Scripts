#while read id ; do fastq-dump --split-files $id ; done < SRR.txt
cat SRR.txt | parallel fastq-dump {}
