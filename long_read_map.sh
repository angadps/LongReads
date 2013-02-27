#!/bin/sh
#$ -cwd
### 6G,8::

base=$1
region=$2

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi

prefix=${base}_${region}

# One-time calling
# bwa index -a bwtsw bcm_hg18.fasta

# Defaults are: -b 3 -q 5 -r 2
# Expt range: -b [4-6] -q [2-4] -r [1-2]

bwa bwasw -t 8 -b 4 -q 4 -r 2 -w 500 reference/bcm_hg18.fasta ${prefix}.fq  > ${prefix}_temp.samd

awk -v reg="chr${region}" 'function rem(cigar) {num=0;for(it=1;it<=length(cigar);it++) {sb=substr(cigar,it,1); if(sb ~ /[0-9]/) {num=10*num+sb;} else { if(sb!="S") return 0; else return num;}}} {if((NF==3)||!match($10,"N")&&$4>0&&$3==reg&&$6!="*"&&rem($6)<50) print}' ${prefix}_temp.samd > ${prefix}_temp.sam
samtools calmd -S ${prefix}_temp.sam reference/bcm_hg18.fasta > ${prefix}.sam

rm ${prefix}_temp.samd ${prefix}_temp.sam
samtools view -bS ${prefix}.sam | samtools sort -m 2000000000 - ${prefix}.sorted
samtools index ${prefix}.sorted.bam
samtools view -h ${prefix}.sorted.bam > ${prefix}.sorted.sam

grep -v ^@ ${prefix}.sorted.sam | cut -f 1,4 | cut -d ':' -f 3 > ${prefix}.ibs

rm ${prefix}.fq
rm ${prefix}.sam

