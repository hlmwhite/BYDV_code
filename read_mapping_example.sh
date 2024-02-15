#!/bin/bash

##### 
# Author: Mark Whitehead, hlmwhite@liverpool.ac.uk
#
# Example code for read mapping and bam processing
#####



cat ../reference/*/data/GC*/*fna ../reference/S.avenae-primary-V1.0.fna | awk '{print $1}' | oneline_script.sh - | fold -w 80 > ref.fasta
 
bwa index ref.fasta
 
for dir in ../reads/*SA*
do
    SAMPLE=$(basename $dir)
    mkdir $SAMPLE
    R1_FILE=$(readlink -f $dir/*R1*gz)
    R2_FILE=$(readlink -f $dir/*R2*gz)
    printf "bwa mem -t 32 ref.fasta $R1_FILE $R2_FILE | samtools view -Sb - > $SAMPLE/aln.bam\n"
done > bwa.comms
 
parallel --bar -j 4 < bwa.comms
 
 
for dir in ../reads/*SA*
do
    SAMPLE=$(basename $dir)
    samtools sort -o $SAMPLE/sort.bam $SAMPLE/aln.bam
done
 
for dir in ../reads/*SA*
do
    SAMPLE=$(basename $dir)
    samtools index $SAMPLE/sort.bam &
done
 
 
for dir in ../reads/*SA*
do
    SAMPLE=$(basename $dir)
    java -Xmx10G -jar picard-2.8.2.jar MarkDuplicates I=$SAMPLE/sort.bam O=$SAMPLE/marked_duplicates.bam M=$SAMPLE/marked_dup_metrics.txt &
done
