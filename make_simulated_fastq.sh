date && time randomreads.sh ref=AM491361.fa out=AM491361.fastq length=150 reads=30 paired=t

date && time reformat.sh in=AM491361.fastq out1=AM491361_R1.fastq out2=AM491361_R2.fastq

## java -Xmx10024m -jar igv.jar
