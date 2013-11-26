#!/bin/sh
#$ -cwd
#$ -V
#$ -o simulate.out
#$ -e simulate.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

snprate=$1
errrate=$2
coverage=$3
readlen=$4
region=$5
base=$6
som_base=$7
freq=$8

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi

prefix=${base}_${region}

./simulate -out_base=$prefix -coverage=$coverage -snp_rate=$snprate -err_rate=$errrate -read_len=$readlen -region=$region -som_file=${som_base}_${region}.list -freq=$freq > simulate_${snprate}_${errrate}_${coverage}_${readlen}_${freq}_${region}.log

# Temporary workaround in lieu of the bug in simulation
simlen=`wc -l ${prefix}.fq | cut -d' ' -f 1`
let "nsimreads=${simlen}/4"
let "actlen=${nsimreads}*4"
head -${actlen} ${prefix}.fq > ${prefix}_new.fq
mv ${prefix}_new.fq ${prefix}.fq

# set args -out_base=/home/singhann/HapMut/reads/0.01_0.01_10_1000_21 -coverage=10 -snp_rate=0.01 -err_rate=0.01 -read_len=1000 -region=21 -som_file=/home/singhann/HapMut/snps/som_21.list

