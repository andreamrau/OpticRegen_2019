#!/usr/bin/env bash

## In TRIMMED, we have the trimmed reads (paired end mode)
## In TRIMMED_trim1, we trim an additional 1bp from the 3' end to avoid alignment problems with Bowtie1
sample=( 0ATAC-1 0ATAC-2 0ATAC-3 12ATAC-1 12ATAC-2 12ATAC-3 2ATAC-1 2ATAC-2 2ATAC-3 4ATAC-1 4ATAC-2 4ATAC-3 7ATAC-2 7ATAC-3 )
for s in "${sample[@]}"
do
sbatch --wrap="source /etc/bashrc
module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3
trim_galore --stringency 3 -q 20 --paired -o ~/Optic-Regen/results/ATACseq/TRIMMED/ \
  ~/Optic-Regen/results/ATACseq/MERGED_LANES/${s}_R1.fastq.gz \
  ~/Optic-Regen/results/ATACseq/MERGED_LANES/${s}_R2.fastq.gz > ~/Optic-Regen/results/ATACseq/TRIMMED/${s}_trash.out

trim_galore --stringency 3 -q 20 --paired --trim1 -o ~/Optic-Regen/results/ATACseq/TRIMMED_trim1/ \
  ~/Optic-Regen/results/ATACseq/MERGED_LANES/${s}_R1.fastq.gz \
  ~/Optic-Regen/results/ATACseq/MERGED_LANES/${s}_R2.fastq.gz > ~/Optic-Regen/results/ATACseq/TRIMMED_trim1/${s}_trash.out"
done


