#!/bin/sh
#$ -cwd

snprate=0.02
errrate=0.02
coverage=15
readlen=2000
region=$1
iteration=1

snprate=$1
errrate=$2
coverage=$3
readlen=$4
region=$5

base=$PWD
pre=${snprate}_${errrate}_${coverage}_${readlen}
read_base=${base}/reads/${pre}
vcf_base=${base}/vcfs/${pre}
som_base=${base}/snps/som_${pre}

#gcc -g simulate_read.c -o simulate_read 

if [[ $region -ne 0 ]]
then
	qsub -sync y -l mem=10G,time=4:: ./simulate_read.sh $snprate $errrate $coverage $readlen $region ${read_base} ${som_base}
	echo simulation complete

	let "qtime=3*$coverage"
	qsub -sync y -l mem=8G,time=52:: ./long_read_map.sh ${read_base} $region
	echo alignment complete

	let "qtime=3*$coverage"
	qsub -sync y -l mem=1G,time=16:: ./snp_calling.sh ${read_base} ${vcf_base} $region
	echo snp calling complete

	cd HMM
	output_base=${base}/HMM/output/${pre}
	qsub -sync y -l mem=4G,time=32:: ./program.sh ${read_base} ${vcf_base} $region $iteration ${som_base} ${output_base} $pre $coverage $errrate
	echo HMM complete
	qsub -l mem=1G,time=2:: ./analysis.sh ${output_base}_${region} ${read_base}_${region} $coverage
else
	let "qmem=400*$coverage"
	qsub -sync y -t 1-22 -l mem=10G,time=4:: ./simulate_read.sh $snprate $errrate $coverage $readlen $region ${read_base} ${som_base}
	echo simulation complete

	let "qtime=3*$coverage"
	qsub -sync y -t 1-22 -l mem=8G,time=32:: ./long_read_map.sh ${read_base} $region

	rm ${read_base}.ibs
	for i in `seq 1 22`
	do
		echo ${read_base}_${i}.ibs
		cat ${read_base}_${i}.ibs >> ${read_base}.ibs
		rm ${read_base}_${i}.ibs
	done
	echo alignment complete

	let "qtime=3*$coverage"
	qsub -sync y -t 1-22 -l mem=1G,time=16:: ./snp_calling.sh ${read_base} ${vcf_base} $region
	echo snp calling complete

	cd HMM
	output_base=${base}/HMM/output/${pre}_${iteration}
	qsub -sync y -t 1-22 -l mem=4G,time=32:: ./program.sh ${read_base} ${vcf_base} $region $iteration ${som_base} ${output_base} $pre $coverage $errrate

	rm ${output_base}.obs
	rm ${output_base}.sta
	for i in `seq 1 22`
	do
		echo $i
		cat ${output_base}_${i}.obs >> ${output_base}.obs
		cat ${output_base}_${i}.sta >> ${output_base}.sta
		rm ${output_base}_${i}.obs
		rm ${output_base}_${i}.sta
	done
	echo HMM complete

	qsub -l mem=1G,time=4:: ./analysis.sh ${output_base} ${read_base} $coverage
fi

cd ..

echo Analysis complete

