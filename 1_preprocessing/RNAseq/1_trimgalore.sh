#!/usr/bin/env bash

## In TRIMMED, we have the trimmed reads (paired end mode)
## In TRIMMED_trim1, we trim an additional 1bp from the 3' end to avoid alignment problems with Bowtie1
homedir='~/Optic-Regen/results/RNAseq/'
sample=( 0RNA1 0RNA2 0RNA3 12RNA1 12RNA2 12RNA3 2RNA1 2RNA2 2RNA3 4RNA1 4RNA2 4RNA3 7RNA1 7RNA2 7RNA3 )
for s in "${sample[@]}"
do
sbatch --wrap="source /etc/bashrc
module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3
trim_galore --stringency 3 -q 20 --paired -o ~/Optic-Regen/results/RNAseq/TRIMMED/ ${homedir}MERGED_LANES/${s}_R1.fastq.gz \
${homedir}MERGED_LANES/${s}_R2.fastq.gz > ~/Optic-Regen/results/RNAseq/TRIMMED/${s}_trash.out

trim_galore --stringency 3 -q 20 --paired --trim1 -o ~/Optic-Regen/results/RNAseq/TRIMMED_trim1/ ${homedir}MERGED_LANES/${s}_R1.fastq.gz \
${homedir}MERGED_LANES/${s}_R2.fastq.gz > ~/Optic-Regen/results/RNAseq/TRIMMED_trim1/${s}_trash.out"
done
