#!/usr/bin/env bash

homedir='~/Optic-Regen/results/ATACseq/'
sample=( 0ATAC-1 0ATAC-2 0ATAC-3 12ATAC-1 12ATAC-2 12ATAC-3 2ATAC-1 2ATAC-2 2ATAC-3 4ATAC-1 4ATAC-2 4ATAC-3 7ATAC-2 7ATAC-3 )
for s in "${sample[@]}"
do
sbatch --wrap="source /etc/bashrc
module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3

fastqc ${homedir}TRIMMED/${s}_R1_val_1.fq.gz -o ${homedir}TRIMMED/QC/
fastqc ${homedir}TRIMMED/${s}_R2_val_2.fq.gz -o ${homedir}TRIMMED/QC/

fastqc ${homedir}TRIMMED_trim1/${s}_R1_val_1.fq.gz -o ${homedir}TRIMMED_trim1/QC/
fastqc ${homedir}TRIMMED_trim1/${s}_R2_val_2.fq.gz -o ${homedir}TRIMMED_trim1/QC/

fastqc ${homedir}MERGED_LANES/${s}_R1.fastq.gz -o ${homedir}MERGED_LANES/QC/
fastqc ${homedir}MERGED_LANES/${s}_R2.fastq.gz -o ${homedir}MERGED_LANES/QC/"
done




