#!/bin/bash

if [ $# -lt 1 ];then
    echo "Need 1 parameters! <sample name>"
    exit
fi

R=$1

# QC
fastqc -o $QCWD -f fastq $FASTQWD/${R}.fastq.gz

# Mapping
module load bowtie2
bowtie2 -p 30 -x $INDEX -U ${FASTQWD}/${R}.fastq.gz -S ${SAMWD}/${R}.sam   2>&1 >>/dev/null |tee -a ${SAMWD}/${R}.bowtie2.SE.out

# Samtools
module load samtools
samtools view -q 30 -bht ${LEN} ${SAMWD}/${R}.sam -o ${BAMWD}/${R}.q30.bam

# Bamtools
module load bamtools
bamtools sort -in ${BAMWD}/${R}.q30.bam -out ${BAMWD}/${R}.q30.sorted.bam
bamtools coverage -in ${BAMWD}/${R}.q30.sorted.bam -out $BWWD/${R}.bw

# use sicer to call broad peaks
sicer -t ${BAMWD}/${R}.q30.sorted.bam -s hg38 -o $SICERWD
wigToBigWig $SICERWD/${R}.q30.sorted-W200-normalized.wig $LEN $SICERWD/${R}.q30.sorted-W200-normalized.bw
mv $SICERWD/${R}.q30.sorted-W200-G600.scoreisland $SICERWD/${R}.bed
