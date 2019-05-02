!/usr/bin/env bash                                                                                                                                                                                         

sample=( 0RNA1 0RNA2 0RNA3 12RNA1 12RNA2 12RNA3 2RNA1 2RNA2 2RNA3 4RNA1 4RNA2 4RNA3 7RNA1 7RNA2 7RNA3 )
for s in "${sample[@]}"
do

mkdir /Users/raua/Desktop/Optic-Regen_results/RNAseq/KALLISTO_PSEUDOSAM/KALLISTO_v2_${s}

cd /Users/raua/Desktop/Optic-Regen_results/RNAseq/KALLISTO_PSEUDOSAM/TRIMMED

kallisto quant -i /Users/raua/Desktop/Optic-Regen_results/RNAseq/KALLISTO_PSEUDOSAM/danio_transcripts_transgene_simplified_kallisto_index -o /Users/raua/Desktop/Optic-Regen_results/RNAseq/KALLISTO_PSEUDOSAM/KALLISTO_v2_${s} --genomebam --gtf /Users/raua/Desktop/Optic-Regen_results/RNAseq/KALLISTO_PSEUDOSAM/Danio_rerio.GRCz10.84.gtf --chromosomes /Users/raua/Desktop/Optic-Regen_results/RNAseq/KALLISTO_PSEUDOSAM/danRer10.chrom.sizes.reformat ${s}_R1_val_1.fq.gz ${s}_R2_val_2.fq.gz

done