#!/usr/bin/env bash

## Make peaks from PE ATAC-seq aligned (BWA) reads using MACS2
## Then we can quantify peaks using featureCounts

aligndir=~/Optic-Regen/results/ATACseq/ALIGNED_TRANS/
peaksdir=~/Optic-Regen/results/ATACseq/PEAKS_TRANS_PEAKLETS_ALLTIMES/

echo Peak calling, merging replicates across all time points, calling peaklets...

sbatch --mem-per-cpu=60000 --wrap="source /etc/bashrc
module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3
echo *0 time pt
macs2 callpeak --nomodel -t ${aligndir}0ATAC-1.nodup.unique.bam ${aligndir}0ATAC-2.nodup.unique.bam ${aligndir}0ATAC-3.nodup.unique.bam ${aligndir}12ATAC-1.nodup.unique.bam ${aligndir}12ATAC-2.nodup.unique.bam ${aligndir}12ATAC-3.nodup.unique.bam ${aligndir}2ATAC-1.nodup.unique.bam ${aligndir}2ATAC-2.nodup.unique.bam ${aligndir}2ATAC-3.nodup.unique.bam ${aligndir}4ATAC-1.nodup.unique.bam ${aligndir}4ATAC-2.nodup.unique.bam ${aligndir}4ATAC-3.nodup.unique.bam ${aligndir}7ATAC-2.nodup.unique.bam ${aligndir}7ATAC-3.nodup.unique.bam -f BAMPE -g 1.37e+09 --call-summits -n ATAC.nodup.unique.macs.peaklets --keep-dup all --outdir ${peaksdir}
"


