#!/bin/sh
#$ -cwd

prefix=$1
suffix=$2
region=$3
iteration=$4
out=$5
pre=$6
errrate=$7

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi

rbase=${prefix}_${region}.sorted
vbase=${suffix}_${region}.vcf
obase=${out}_${region}

#make program

# ./program $iteration > program.${prefix}_${iteration}_${region}

./program $rbase $vbase $region $iteration $obase $errrate > program.${pre}_${region}_${iteration}.log

