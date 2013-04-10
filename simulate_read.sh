#!/bin/sh
#$ -cwd

snprate=$1
errrate=$2
coverage=$3
readlen=$4
region=$5
base=$6
som_base=$7

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi

prefix=${base}_${region}

echo -out_base=$prefix -coverage=$coverage -snp_rate=$snprate -err_rate=$errrate -read_len=$readlen -region=$region

./simulate -out_base=$prefix -coverage=$coverage -snp_rate=$snprate -err_rate=$errrate -read_len=$readlen -region=$region -som_file=${som_base}_${region}.list

# Temporary workaround in lieu of the bug in simulation
simlen=`wc -l ${prefix}.fq | cut -d' ' -f 1`
let "nsimreads=${simlen}/4"
let "actlen=${nsimreads}*4"
head -${actlen} ${prefix}.fq > ${prefix}_new.fq
mv ${prefix}_new.fq ${prefix}.fq


