### Source code for Dhara et al. (2019)

This repository contains the following source code files used to analyze the RNA-seq and ATAC-seq data in Dhara et al. (2019).

It is organized as follows:

- **1_preprocessing/**
    - **ATACseq**

        - `1_merge_lanes.txt` 
        - `2_trimgalore.sh`
        - `3_fastqc.sh`
        - `4_bwa_alignment_danio_w_transgene.sh`
        - `5_fastqc_aligned.sh`
        - `6_macs_peaklets_alltimes.sh`
        - `7_remove_MTduplicates.sh`
        - `get_fasta_from_bed.sh`


    - **RNAseq**

        - `1_trimgalore.sh` 
        - `2_fastqc.sh`
        - `3_kallisto_hdf5.sh`

- **2_differential-analysis/**

    - **ATACseq**

        - `1_differential-analysis_diffBind_all-vs-0.R`  
        - `2_sigresults-as-bed-file.R`  
        - `4timepts_vs0_results`

    - **RNAseq**

        - `0_my_sleuth_functions.R`  
        - `1_differential-analysis_sleuth.R`  
        - `4timepts_vs0_results`
        - `format_kallisto_as_pseudosam`


- **3_peak-gene-annotation**  

        - `1_make_MEME_format.R`  
        - `2_AME_reformat.R` 
        - `3_TF_cooccurrence-graphs.R`
        - `4_identify_distal-proximal_peaks.R`
        - `peaklets`


- **misc**

        - `checking_final_numbers.R`  
        - `danRer10_alias.tab` 