#!/usr/bin/env bash

## Run kallisto with 500 bootstraps
sample=( 0RNA1 0RNA2 0RNA3 12RNA1 12RNA2 12RNA3 2RNA1 2RNA2 2RNA3 4RNA1 4RNA2 4RNA3 7RNA1 7RNA2 7RNA3 )
for s in "${sample[@]}"
do
sbatch --mem-per-cpu=64000 --job-name=${s} --wrap="source /etc/bashrc
module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3
mkdir ~/Data/Optic-Regen/results/RNAseq/KALLISTO_hdf5/KALLISTO_v2_${s}
kallisto quant -i ~/Data/Optic-Regen/tools/REF_DATA/danio_transcripts_transgene_kallisto_index -o ~/Data/Optic-Regen/results/RNAseq/KALLISTO_hdf5/KALLISTO_v2_${s} -b 500  \
 ~/Data/Optic-Regen/results/RNAseq/TRIMMED/${s}_R1_val_1.fq.gz ~/Data/Optic-Regen/results/RNAseq/TRIMMED/${s}_R2_val_2.fq.gz"
done

