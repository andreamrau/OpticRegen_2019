#!/usr/bin/env bash

## Convert SAM to BAM
## Remove mitochondrial, duplicate reads, and multiple mapped reads using samtools

mt=MT
aligndir=~/Data/Optic-Regen/results/ATACseq/ALIGNED_TRANS/
peaksdir=~/Data/Optic-Regen/results/ATACseq/PEAKS_TRANS/

sample=( 0ATAC-1 0ATAC-2 0ATAC-3 12ATAC-1 12ATAC-2 12ATAC-3 2ATAC-1 2ATAC-2 2ATAC-3 4ATAC-1 4ATAC-2 4ATAC-3 7ATAC-2 7ATAC-3 )
for s in "${sample[@]}"
do

saminput=${aligndir}${s}_trans.sam
baminput=${aligndir}${s}.bam
bamsortedinput=${aligndir}${s}.sorted.bam
# nomt=${aligndir}${s}.noMT.bam
nodup=${aligndir}${s}.nodup.bam
nodupunique=${aligndir}${s}.nodup.unique.bam

sbatch --mem-per-cpu=60000 --wrap="source /etc/bashrc
module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3

echo Convert SAM to BAM, sort and index...
samtools view -S -b ${saminput} > ${baminput}
samtools sort ${baminput} -o ${bamsortedinput}
samtools index ${bamsortedinput}
samtools idxstats ${bamsortedinput} > ${bamsortedinput}.idxstats
samtools flagstat ${bamsortedinput} > ${bamsortedinput}.flagstat

# I REMOVED THIS SINCE IT DOESNT LOOK LIKE THERE ARE ANY ALIGNED TO MT READS?
# echo Removing MT reads...
# samtools idxstats $bamsortedinput | awk '{print $1}' | grep -v $mt | xargs samtools view -bh $baminput 
# samtools index $nomt
# samtools idxstats $nomt
# samtools flagstat $nomt

echo Removing duplicate reads...
samtools sort -n -o ${bamsortedinput}.namesort ${bamsortedinput} 
samtools fixmate -m ${bamsortedinput}.namesort ${bamsortedinput}.fixmate
samtools sort -o ${bamsortedinput}.fixmate.sort ${bamsortedinput}.fixmate 
samtools markdup -r -s ${bamsortedinput}.fixmate.sort ${nodup}
samtools index ${nodup}
samtools idxstats $nodup > ${nodup}.idxstats
samtools flagstat $nodup > ${nodup}.flagstat

echo Removing multiple mapped reads...
samtools view -bq 1 ${nodup} > ${nodupunique}
samtools index ${nodupunique}
samtools idxstats ${nodupunique} > ${nodupunique}.idxstats
samtools flagstat ${nodupunique} > ${nodupunique}.flagstat
"

done

