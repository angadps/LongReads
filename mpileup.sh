#!/bin/sh
#$ -cwd
#$ -V
#$ -o mpileup.out
#$ -e mpileup.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

rbase=$1
region=$2
vbase=$3
snprate=$4
i=$SGE_TASK_ID

lent=`samtools view -H ${rbase}.sorted.bam | grep -w "SN:chr${region}" | cut -d':' -f 3`
let "overlap=$lent/200"

if [[ $i -ne 1 ]]
then
	let "start=(${i}-1)*$lent/10+1-$overlap"
else
	let "start=(${i}-1)*$lent/10+1"
fi

if [[ $i -ne 10 ]]
then
	let "end=${i}*$lent/10+$overlap"
else
	let "end=${i}*$lent/10"
fi

echo $i $lent $start $end

p=`echo $snprate | awk '{print -10*log($0)/log(10);}'`
# samtools 0.1.18
samtools mpileup -E -Q 0 -o $p -e 50 -uIf reference/hg19.fasta -d 300 -gDS -r "chr${region}:${start}-${end}" ${rbase}.sorted.bam > ${vbase}.bcf.$i # -C 1 -6

bcftools view -gvN -e -i 1.0 -p 1.1 ${vbase}.bcf.${i} > ${vbase}.${i}.temp.vcf

