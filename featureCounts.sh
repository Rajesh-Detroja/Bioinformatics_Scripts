#!/usr/bin

while getopts ":1:" opt; do
  case $opt in
    1) bam="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

date && time featureCounts -p -T 30 -a /home/morgensternlab/detrojar/db/hg38/index/star/hg38.gtf -o counts.txt *$bam
