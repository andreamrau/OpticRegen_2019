#!/usr/bin/env bash

homedir='~/Optic-Regen/results/RNAseq/'
sample=( 0RNA1 0RNA2 0RNA3 12RNA1 12RNA2 12RNA3 2RNA1 2RNA2 2RNA3 4RNA1 4RNA2 4RNA3 7RNA1 7RNA2 7RNA3 )
for s in "${sample[@]}"
do
sbatch --wrap="source /etc/bashrc
module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3

fastqc ${homedir}MERGED_LANES/${s}_R1.fastq.gz -o ${homedir}MERGED_LANES/QC/
fastqc ${homedir}MERGED_LANES/${s}_R2.fastq.gz -o ${homedir}MERGED_LANES/QC/

fastqc ${homedir}TRIMMED/${s}_R1_val_1.fq.gz -o ${homedir}TRIMMED/QC/
fastqc ${homedir}TRIMMED/${s}_R2_val_2.fq.gz -o ${homedir}TRIMMED/QC/

fastqc ${homedir}TRIMMED_trim1/${s}_R1_val_1.fq.gz -o ${homedir}TRIMMED_trim1/QC/
fastqc ${homedir}TRIMMED_trim1/${s}_R2_val_2.fq.gz -o ${homedir}TRIMMED_trim1/QC/
"
done
