#!/usr/bin/env bash

module load pkgsrc/2018Q1

bedtools getfasta -fi ~/Optic-Regen/tools/REF_DATA/danio_rerio.fa \
-bed ~/Optic-Regen/results/sigpeaks_2vs0_peaklets_500bp_down.bed \
-fo ~/Optic-Regen/results/sigpeaks_2vs0_peaklets_500bp_down.seq
