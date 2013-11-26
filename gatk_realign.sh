#!/bin/sh
#$ -cwd
#$ -V
#$ -o realign.out
#$ -e realign.err
#$ -S /bin/bash
. /usr/prog/softenv/softenv.sh

JAVA="java -Xmx8000m -Djava.io.tmpdir=./tmp"
GATK="$JAVA -jar /usr/prog/onc/seqtools/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar"

PREFIX=$1
CHR="chr$2"

INP=${PREFIX}.rg.bam
REF=reference/hg19.fasta
INDELVCF=reference/dbsnp137.vcf

$GATK \
    -L $CHR \
    -T RealignerTargetCreator \
    -I $INP \
    -R $REF \
    -known $INDELVCF \
    -o ${PREFIX}.forRealigner.intervals

$GATK \
    -L $CHR \
    -I $INP \
    -R $REF \
    -T IndelRealigner \
    -known $INDELVCF \
    -compress 5 \
    -targetIntervals ${PREFIX}.forRealigner.intervals \
    -o ${PREFIX}.realigned.bam

