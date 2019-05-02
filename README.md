### Source code for Dhara et al. (2019)

This repository contains the following source code files used to analyze the RNA-seq and ATAC-seq data in Dhara et al. (2019).

It is organized as follows:

- **1_preprocessing/**
    - **ATACseq**

        - `1_merge_lanes.txt`: concatenate sequencing lanes containing technical replicates into a single file
        - `2_trimgalore.sh`: trim adaptor sequences
        - `3_fastqc.sh`: validate sequence quality
        - `4_bwa_alignment_danio_w_transgene.sh`: align reads to the zebrafish genome (GRCz10) and transgene sequence
        - `5_fastqc_aligned.sh`: validate aligned sequence quality
        - `6_macs_peaklets_alltimes.sh`: call peaks from aligned reads and identify summits of deconvoluted subpeaks were identified
        - `7_remove_MTduplicates.sh`: remove duplicates and multiple mapped reads
        - `get_fasta_from_bed.sh`: script to obtain fasta file of sequence corresponding to a provided bed file


    - **RNAseq**

        - `1_trimgalore.sh`: trim adaptor sequences 
        - `2_fastqc.sh`: validate sequence quality
        - `3_kallisto_hdf5.sh`: quantify transcriptome abundances using HDF5 data format

- **2_differential-analysis/**

    - **ATACseq**  
        - `1_differential-analysis_diffBind_all-vs-0.R`: perform differential analysis of quantified peaklets for time points 2, 4, 7, and 12 versus 0 using *diffBind*  
        - `2_sigresults-as-bed-file.R`: output bed files of significant peaks  
        - `4timepts_vs0_results/`: folder containing result files (bed and csv) for each comparison

    - **RNAseq**
        - `0_my_sleuth_functions.R`: internal function slightly modifying the *sleuth* code to enable target IDs to be updated  
        - `1_differential-analysis_sleuth.R`: perform differential analysis of expression data for time points 2, 4, 7, and 12 versus 0 using *sleuth*  
        - `4timepts_vs0_results`: folder containing result files (csv) for each comparison  
        - `format_kallisto_as_pseudosam/`: folder containing helper scripts to convert *kallisto* quantifications to a pseudosam file for visualization in IGV


- **3_peak-gene-annotation/**

    - `1_make_MEME_format.R`: script to automate making a file in MEME format based on Jaspar and Lambert IDs
    - `2_AME_reformat.R`: script to reformat AME results for downstream analysis
    - `3_TF_cooccurrence-graphs.R`: produce visualizations of TF motif co-occurrence (network plots, UpsetR plots, bar graphs, heatmaps)
    - `4_identify_distal-proximal_peaks.R`: identify annotations for peaks (e.g. distal and proximal to genes, exonic/intronic/overlapping TSS)
    - `peaklets/`: folder containing MACS output files for peaks (`ATAC.nodup.unique.macs.peaklets_peaks.narrowPeak`) and peak summits (`ATAC.nodup.unique.macs.peaklets_summits.bed`) as well as a bed file of peaklet locations (`ALLMERGED_ATAC.nodup.unique.macs.peaklets_peaks.pvalsort.narrowPeak_500bp.bed`)

- **misc/**
  - `checking_final_numbers.R`: script to check final numbers for the paper
   - `danRer10_alias.tab`: alias file providing correspondence between UCSC and Ensembl chromosome naming schemes