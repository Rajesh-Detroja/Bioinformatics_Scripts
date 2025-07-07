#!/usr/bin

while getopts ":A:B:" opt; do
  case $opt in
    A) sample_A="$OPTARG"
    ;;
    B) sample_B="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# paired-end data
date && time danpos dpos $sample_A/:$sample_B/ -o $sample_A.$sample_B -m 1

# single-end data
#date && time danpos dpos $sample_A/:$sample_B/ -o $sample_A.$sampleB
