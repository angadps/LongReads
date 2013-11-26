#!/bin/sh
#$ -cwd
#$ -V
#$ -o program.out
#$ -e program.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

prefix=$1
suffix=$2
region=$3
iteration=$4
som=$5
out=$6
pre=$7
coverage=$8
errrate=$9
freq=$10

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi

rbase=${prefix}_${region}.sorted
vbase=${suffix}_${region}.vcf
sbase=${som}_${region}.list
obase=${out}_${region}

# make program

./program $rbase $vbase $region $sbase $obase $coverage $errrate > program.${pre}_${region}.log

# set args /home/singhann/HapMut/reads/0.01_0.01_50_2000_21.sorted /home/singhann/HapMut/vcfs/0.01_0.01_50_2000_21.vcf 21 /home/singhann/HapMut/snps/som_0.01_0.01_50_2000_21.list /home/singhann/HapMut/output/0.01_0.01_50_2000_21 50 0.01
