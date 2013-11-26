#!/bin/sh
#$ -cwd

s1="0.005_0.005_5_500"
s2="0.01_0.01_5_500"
s3="0.01_0.01_10_500"
s4="0.01_0.01_5_1000"
s5="0.015_0.015_5_1000"
s6="0.015_0.015_10_1000"
s7="0.015_0.015_5_2000"
s8="0.015_0.015_10_2000"
s9="0.02_0.02_10_2000"
s10="0.015_0.015_5_4000"
s11="0.02_0.02_5_4000"
s12="0.02_0.02_10_4000"

som1="0.005_0.005_5_500"
som2="0.005_0.005_5_600"
som3="0.005_0.005_10_500"
som4="0.005_0.005_10_600"
som5="0.02_0.02_15_2000"
som6="0.02_0.02_15_900"

for i in `echo $s1 $s2 $s3 $s4 $s5 $s6 $s7 $s8 $s9 $s10 $s11 $s12`
do
#	i=0.0045_0.0005_30_600
	echo $i
	sr=`echo $i | cut -d'_' -f1`
	ir=`echo $i | cut -d'_' -f2`
	cov=`echo $i | cut -d'_' -f3`
	rl=`echo $i | cut -d'_' -f4`
	Rscript plot_novel.R $sr $ir $cov $rl
done

for i in `echo $som1 $som3 $som5`
do
	echo $i
	sr=`echo $i | cut -d'_' -f1`
	ir=`echo $i | cut -d'_' -f2`
	cov=`echo $i | cut -d'_' -f3`
	rl=`echo $i | cut -d'_' -f4`
	prev=`echo "scale=2; (100-2*$cov)/100" | bc -l`
	Rscript plot_som.R $sr $ir $cov $rl $prev
done

for i in `echo $som2 $som4 $som6`
do
	echo $i
	sr=`echo $i | cut -d'_' -f1`
	ir=`echo $i | cut -d'_' -f2`
	cov=`echo $i | cut -d'_' -f3`
	rl=`echo $i | cut -d'_' -f4`
	prev=`echo "scale=2; (80-2*$cov)/100" | bc -l`
	Rscript plot_som.R $sr $ir $cov $rl $prev
done

