#!/bin/sh
#$ -cwd

prefix=$1
suffix=$2
region=$3
iteration=$4
som=$5
out=$6
pre=$7
coverage=$8
errrate=$9

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi

rbase=${prefix}_${region}.sorted
vbase=${suffix}_${region}.vcf
sbase=${som}_${region}.list
obase=${out}_${region}

# cd /ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/NB
# make program
# cp program /ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/HMM
# cd ../HMM

./program $rbase $vbase $region $sbase $obase $coverage $errrate > program.${pre}_${region}_${iteration}.log

# set args /ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/reads/0.02_0.02_10_4000_21.sorted /ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/vcfs/0.02_0.02_10_4000_21.vcf 21 /ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/HMM/output/0.02_0.02_10_4000_1_21 0.02
