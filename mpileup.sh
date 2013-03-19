#!/bin/sh
#$ -cwd

rbase=$1
region=$2
vbase=$3
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

samtools mpileup -uIf reference/bcm_hg18.fasta -d 300 -q 0 -Q 0 -gDS -r "chr${region}:${start}-${end}" ${rbase}.sorted.bam > ${vbase}.bcf.$i # -C 1 -6 -E

bcftools view -gvN -e -i 1.0 -p 1.1 ${vbase}.bcf.${i} > ${vbase}.${i}.temp.vcf

