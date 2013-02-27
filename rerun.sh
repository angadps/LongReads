#!/bin/sh
#$ -cwd

prefix=$1
suffix=$2
region=$3
iteration=$4
out=$5
pre=$6

rbase=${prefix}_${region}.sorted
vbase=${suffix}_${region}.vcf
obase=${out}_${region}

#./simulate -out_base=/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/0.01_0.01_10_1000_${iter} -coverage=10 -snp_rate=0.01 -err_rate=0.01 -read_len=1000 -region=${iter}

cd HMM
./program $rbase $vbase $region $iteration $obase > program.${pre}_${region}_${iteration}.log
