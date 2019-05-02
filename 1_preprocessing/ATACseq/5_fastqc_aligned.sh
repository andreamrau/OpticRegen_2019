#!/usr/bin/env bash

homedir='~/Optic-Regen/results/ATACseq/'
sample=( 0ATAC-1 0ATAC-2 0ATAC-3 12ATAC-1 12ATAC-2 12ATAC-3 2ATAC-1 2ATAC-2 2ATAC-3 4ATAC-1 4ATAC-2 4ATAC-3 7ATAC-2 7ATAC-3 )
for s in "${sample[@]}"
do
sbatch --mem-per-cpu=60000 --wrap="source /etc/bashrc
module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3

fastqc ${homedir}ALIGNED/${s}.sam -o ${homedir}ALIGNED_QC/"

done




