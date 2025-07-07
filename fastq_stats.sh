#!/bin/sh
#fastq-stats

(( c1=0 ))
for i in $@
	do
		if [[ ( -n "$i" && "${i: -9}" == ".fastq.gz") || ( -n "$i" && "${i: -6}" == ".fq.gz")|| ( -n "$i" && "${i: -6}" == ".fastq") || ( -n "$i" && "${i: -3}" == ".fq") ]]
			then
				if [ $c1 -eq 0 ]
					then
						printf "File Name\tTotal Reads\tTotal Bases\tMin Length\tMax Length\tMean Length\n"
				fi
			
			(( c1++ ))
			if [ $c1 -gt 1 ]
				then
					echo ""
			fi
			
			echo -n $i
			fastq-stats $i | \
			grep -e "reads" -e "total bases" -e "len min" -e "len	" -e "len mean" | \
			sed 's/reads/A/g; s/total bases/B/g ; s/len min/C/g ; s/len\t/D\t/g ; s/len mean/E/g' | \
			sort | sed 's/A//g; s/B//g ; s/C//g ; s/D//g ; s/E//g'| tr "\n" "\t" | tr -s '\t' '\t'
		else
			exit 1
		fi
done
printf "\n"
