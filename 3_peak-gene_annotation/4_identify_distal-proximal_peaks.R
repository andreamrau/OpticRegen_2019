
## Identify distal and proximal (non-exonic) peaks for each gene
##   - Proximal: peaks within +/- 1 kb window to all genes but not overlapping exons
##   - Distal: peaks within (1kb, 100kb) window around TSS but not overlapping exons
##   - Differential non-exonic peaks
##   - Non-differential non-exonic peaks

# Grcz10 (danRer10)
library(ChIPpeakAnno)
library(org.Dr.eg.db)
library(GenomicFeatures)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(rtracklayer)
library(dplyr)
library(AnnotationHub)
library(forcats)
library(tidyverse)
library(readxl)

## Load Danio Rerio UCSC annotation to anntate peaks with features
ucsc.dr10.refGene <- genes(TxDb.Drerio.UCSC.danRer10.refGene)
ucsc.dr10.exons <- exons(TxDb.Drerio.UCSC.danRer10.refGene)
ucsc.dr10.tss <- promoters(TxDb.Drerio.UCSC.danRer10.refGene, upstream = 0, downstream = 0)
ucsc.dr10.transcripts <- transcripts(TxDb.Drerio.UCSC.danRer10.refGene)


## Import the MACS output (skip the transgene), convert chromosome names using IGV alias
# macs <- read.table(paste0("PEAKS_TRANS/", 
#                           "ALLMERGED_ATAC.nodup.unique.macs_peaks.pvalsort.narrowPeak.bed"), 
#                    sep="\t", skip=1)[,1:3]
# macs <- read.table(paste0("PEAKS_TRANS/", 
#                           "ALLMERGED_ATAC.nodup.unique.macs_peaks.pvalsort.narrowPeak.bed"), 
#                    sep="\t", skip=1)[,1:3]
## AR update: use peaklets (skip first three lines for transgene)
## AR update: use 500bp peaklets (skip first three lines for trnasgene)
macs <- read.table(paste0("PEAKS_TRANS_PEAKLETS_ALLTIMES/", 
                          "ALLMERGED_ATAC.nodup.unique.macs.peaklets_peaks.pvalsort.narrowPeak_500bp.bed"),
                   sep="\t", skip = 3, stringsAsFactors = FALSE)[,1:3]
colnames(macs) <- c("chr_old", "start", "end")
alias <- read.table("danRer10_alias.tab", stringsAsFactors = FALSE)
colnames(alias) <- c("chr_old", "chr")
macs_alias <- left_join(macs, alias, by="chr_old") %>%
  dplyr::select(chr, start, end) %>%
  dplyr::mutate(seqnames=chr) %>%
  dplyr::select(seqnames, start, end)
macsOutput <- GRanges(macs_alias)


# Peaks overlapping exons, TSS, 1kb window, 100kb window------------------------------------------

# peaks_overlapping_exons <- 
#   annotatePeakInBatch(macsOutput, 
#                       PeakLocForDistance = "middle",
#                       AnnotationData=ucsc.dr10.exons, 
#                       output="overlapping", minOverlap = 0L) 
# peaks_overlapping_exons <- unique(peaks_overlapping_exons[which(!is.na(peaks_overlapping_exons$feature))])

peaks_intersect_exons <- findOverlapsOfPeaks(macsOutput, ucsc.dr10.exons, 
                                               minoverlap = 50)$peaklist$`macsOutput///ucsc.dr10.exons`
peaks_overlapping_exons <- 
  annotatePeakInBatch(macsOutput,
                      output="overlapping",
                      AnnotationData=peaks_intersect_exons)
peaks_overlapping_exons <- unique(peaks_overlapping_exons[which(!is.na(peaks_overlapping_exons$feature))])

peaks_overlapping_tss <- 
  annotatePeakInBatch(macsOutput, 
                      PeakLocForDistance = "middle",
                      FeatureLocForDistance = "start",
                      AnnotationData=ucsc.dr10.tss, 
                      output="overlapping", gap=-1L) 
peaks_overlapping_tss <- unique(peaks_overlapping_tss[which(!is.na(peaks_overlapping_tss$feature))])
peaks_in_1kb <- 
 annotatePeakInBatch(macsOutput, 
                     PeakLocForDistance = "middle",
                     FeatureLocForDistance = "TSS",
                     AnnotationData=ucsc.dr10.refGene, 
                     output="overlapping", bindingRegion=c(-1000,1000))
peaks_in_1kb <- unique(peaks_in_1kb[which(!is.na(peaks_in_1kb$feature))])
peaks_in_100kb <- 
  annotatePeakInBatch(macsOutput, 
                      PeakLocForDistance = "middle",
                      FeatureLocForDistance = "TSS",
                      AnnotationData=ucsc.dr10.refGene, 
                      output="overlapping", bindingRegion=c(-100000,100000))
peaks_in_100kb <- unique(peaks_in_100kb[which(!is.na(peaks_in_100kb$feature))])

## AR: change the DE results
## AR update: change the DE results again
sigpeaks_2vs0 <- readRDS("DIFFERENTIAL_vs0/sigpeaks_pvalsort_order_2vs0_peaklets_500bp.rds")
sigpeaks_4vs0 <- readRDS("DIFFERENTIAL_vs0/sigpeaks_pvalsort_order_4vs0_peaklets_500bp.rds")
sigpeaks_12vs0 <- readRDS("DIFFERENTIAL_vs0/sigpeaks_pvalsort_order_12vs0_peaklets_500bp.rds")
sigpeaks <- bind_rows(sigpeaks_2vs0, sigpeaks_4vs0, sigpeaks_12vs0) %>%
   dplyr::mutate(chr_old = seqnames)
# sigpeaks <- bind_rows(sigpeaks_2vs0, sigpeaks_12vs0) %>%
#   dplyr::mutate(chr_old = seqnames)
sigpeaks_alias <- left_join(sigpeaks, alias, by="chr_old") %>%
  dplyr::select(chr, start, end) %>%
  dplyr::mutate(seqnames=chr) %>%
  dplyr::select(seqnames, start, end) %>%
  arrange(seqnames, start)
sigpeaksOutput <- GRanges(sigpeaks_alias) %>% unique()  ## 182

# Proximal peaks = overlapping TSS and/or within 1kb +/- TSS but not overlapping an exon by 50bp  ----------
# Distal peaks = in (1kb, 100kb) +/- TSS but not overlapping an exon by 50bp 
proximal_peaks <- c(granges(setdiff(peaks_in_1kb, peaks_overlapping_exons)), granges(peaks_overlapping_tss)) %>%
  unique() %>% sort() ## 8961
distal_peaks <- setdiff(setdiff(peaks_in_100kb, peaks_in_1kb), peaks_overlapping_exons) %>%
  unique() %>% sort() ## 20712
nondifferential_nonexonic_peaks <- setdiff(setdiff(macsOutput, peaks_overlapping_exons), sigpeaksOutput) ## 39673
differential_nonexonic_peaks <- intersect(setdiff(macsOutput, peaks_overlapping_exons), sigpeaksOutput) ## 174

## How many nondifferential, non-exonic peaks are distal versus proximal?
length(intersect(nondifferential_nonexonic_peaks, proximal_peaks))
length(intersect(nondifferential_nonexonic_peaks, distal_peaks))


write.table(data.frame(proximal_peaks)[,1:3], 
            "ANNOTATION/proximal-1kbtss_nonexonic-50bp_peaklets_500bp.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(data.frame(distal_peaks)[,1:3], 
            "ANNOTATION/distal-1kbto100kb_nonexonic-50bp_peaklets_500bp.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(data.frame(nondifferential_nonexonic_peaks)[,1:3], 
            "ANNOTATION/nondifferential_nonexonic-50bp_peaklets_500bp.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(data.frame(differential_nonexonic_peaks)[,1:3], 
            "ANNOTATION/differential_nonexonic-50bp_peaklets_500bp.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


# DE peaks that are proximal or distal ------------------------------------------
sigpeaks_proximal <- intersect(proximal_peaks, sigpeaksOutput)
sigpeaks_distal <- intersect(distal_peaks, sigpeaksOutput)
sigpeaks_other <- setdiff(sigpeaksOutput, 
                          c(sigpeaks_proximal, 
                            GRanges(data.frame(sigpeaks_distal)[,1:3]))) ## 98 peaks

## Now retrieve genes 
sigpeaks_proximal_withGene <-
  annotatePeakInBatch(sigpeaks_proximal, 
                      PeakLocForDistance = "middle",
                      FeatureLocForDistance = "TSS",
                      AnnotationData=ucsc.dr10.refGene, 
                      output="overlapping", bindingRegion=c(-1000,1000)) %>%
  addGeneIDs(orgAnn="org.Dr.eg.db", 
             feature_id_type="entrez_id",
             IDs2Add="symbol")
sigpeaks_distal_withGene <-
  annotatePeakInBatch(sigpeaks_distal, 
                      PeakLocForDistance = "middle",
                      FeatureLocForDistance = "TSS",
                      AnnotationData=ucsc.dr10.refGene, 
                      output="overlapping", bindingRegion=c(-100000,100000)) %>%
  addGeneIDs(orgAnn="org.Dr.eg.db", 
             feature_id_type="entrez_id",
             IDs2Add="symbol")

sigpeaks_proximal_withGene <- data.frame(sigpeaks_proximal_withGene) %>%
  select(chr=seqnames, start, end, width, strand, distance, insideFeature, gene=symbol) %>%
  unite("new_seqnames", chr, start, sep = ":", remove=FALSE) %>% 
  unite("new_seqnames", new_seqnames, end, sep = "-", remove=FALSE) %>%
  select(chr, start, end, width, strand, distance, insideFeature, gene, new_seqnames)
sigpeaks_distal_withGene <- data.frame(sigpeaks_distal_withGene) %>%
  select(chr=seqnames, start, end, width, strand, distance, insideFeature, gene=symbol) %>%
  unite("new_seqnames", chr, start, sep = ":", remove=FALSE) %>% 
  unite("new_seqnames", new_seqnames, end, sep = "-", remove=FALSE) %>%
  select(chr, start, end, width, strand, distance, insideFeature, gene, new_seqnames)

write.table(data.frame(sigpeaks_proximal_withGene), 
            "ANNOTATION/differential-peaklets_proximal_withGenes_500bp.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(data.frame(sigpeaks_distal_withGene), 
            "ANNOTATION/differential-peaklets_distal_withGenes_500bp.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


# Overlap of significant peaks of with DREME ------------------------------------------

## Distal
dreme5_distpeaks <- read.table("ANNOTATION/DREME5_distpeaks.tsv", header=TRUE, fill=TRUE) %>%
  select(motif_alt_id) %>%
  separate(motif_alt_id, into = c("chr_old","position"), sep=":") %>%
  separate(position, into = c("start", "end"), sep="-") %>%
  left_join(., alias, by="chr_old") %>%
  select(chr, start, end)
dreme5_distpeaksOutput <- GRanges(dreme5_distpeaks)
## Proximal 
dreme5_proxpeaks <- read.table("ANNOTATION/DREME5_proxpeaks.tsv", header=TRUE, fill=TRUE) %>%
  select(motif_alt_id) %>%
  separate(motif_alt_id, into = c("chr_old","position"), sep=":") %>%
  separate(position, into = c("start", "end"), sep="-") %>%
  left_join(., alias, by="chr_old") %>%
  select(chr, start, end)
dreme5_proxpeaksOutput <- GRanges(dreme5_proxpeaks)

## Now match up DREME5 with distal/proximal sig peaks
sigpeaks_distal_dreme5 <-
  annotatePeakInBatch(sigpeaks_distal, 
                      PeakLocForDistance = "middle",
                      FeatureLocForDistance = "TSS",
                      AnnotationData=dreme5_distpeaksOutput, 
                      output="overlapping")
no_match <- grep("NA", names(sigpeaks_distal_dreme5))
sigpeaks_distal_dreme5 <- sigpeaks_distal_dreme5[-no_match] %>% unique %>%
  as.data.frame() %>% select(seqnames, start, end)  ## 37
sigpeaks_proximal_dreme5 <-
  annotatePeakInBatch(sigpeaks_proximal, 
                      PeakLocForDistance = "middle",
                      FeatureLocForDistance = "TSS",
                      AnnotationData=dreme5_proxpeaksOutput, 
                      output="overlapping")
no_match <- grep("NA", names(sigpeaks_proximal_dreme5))
sigpeaks_proximal_dreme5 <- sigpeaks_proximal_dreme5[-no_match] %>% unique %>%
  as.data.frame() %>% select(seqnames, start, end)  ## 17

sigpeaks_distal_dreme5$new_seqnames <- paste0(sigpeaks_distal_dreme5$seqnames, ":", sigpeaks_distal_dreme5$start, "-",
                                              sigpeaks_distal_dreme5$end)

sigpeaks_proximal_dreme5$new_seqnames <- paste0(sigpeaks_proximal_dreme5$seqnames, ":", sigpeaks_proximal_dreme5$start, "-",
                                              sigpeaks_proximal_dreme5$end)

write.table(data.frame(sigpeaks_distal_dreme5), 
            "ANNOTATION/differential-peaklets_distal_DREME5_500bp.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(data.frame(sigpeaks_proximal_dreme5), 
            "ANNOTATION/differential-peaklets_proximal_DREME5_500bp.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")



## Put the original chromosome names back on --------------------------------------
proximal_peaks <- data.frame(proximal_peaks)[,1:3]
colnames(proximal_peaks) <- c("chr", "start", "end")
proximal_peaks <- left_join(proximal_peaks, alias, by = "chr") %>%
  select(chr=chr_old, start, end)
distal_peaks <- data.frame(distal_peaks)[,1:3]
colnames(distal_peaks) <- c("chr", "start", "end")
distal_peaks <- left_join(distal_peaks, alias, by = "chr") %>%
  select(chr=chr_old, start, end)
nondifferential_nonexonic_peaks <- data.frame(nondifferential_nonexonic_peaks)[,1:3]
colnames(nondifferential_nonexonic_peaks) <- c("chr", "start", "end")
nondifferential_nonexonic_peaks <- left_join(nondifferential_nonexonic_peaks, alias, by = "chr") %>%
  select(chr=chr_old, start, end)
differential_nonexonic_peaks <- data.frame(differential_nonexonic_peaks)[,1:3]
colnames(differential_nonexonic_peaks) <- c("chr", "start", "end")
differential_nonexonic_peaks <- left_join(differential_nonexonic_peaks, alias, by = "chr") %>%
  select(chr=chr_old, start, end)

write.table(proximal_peaks, 
            "ANNOTATION/proximal-1kbtss_nonexonic-50bp_peaklets_500bp_origchr.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(distal_peaks, 
            "ANNOTATION/distal-1kbto100kb_nonexonic-50bp_peaklets_500bp_origchr.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(nondifferential_nonexonic_peaks, 
            "ANNOTATION/nondifferential_nonexonic-50bp_peaklets_500bp_origchr.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(differential_nonexonic_peaks, 
            "ANNOTATION/differential_nonexonic-50bp_peaklets_500bp_origchr.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")



## Annotate peaks by differential status, TF --------------------------------------
proximal_peaks <- c(granges(setdiff(peaks_in_1kb, peaks_overlapping_exons)), granges(peaks_overlapping_tss)) %>%
  unique() %>% sort() ## 8961
distal_peaks <- setdiff(setdiff(peaks_in_100kb, peaks_in_1kb), peaks_overlapping_exons) %>%
  unique() %>% sort() ## 20712

detf <- read_excel("../2018_1025_DETF_IDs.xlsx") %>%
  select(zebrafish_gene) %>%
  unlist()

rna_res_all <- c()
for(time in c(2,4,7,12)) {
  for(type in c("up", "down")) {
    tmp <- read.csv("../RNAseq/4_vs0_model/results_2dpi-0dpi.csv",
                    stringsAsFactors = FALSE)
    if(type == "up") {
      tmp <- tmp  %>%
        filter(qval < 0.05 & b > 0 ) %>% select(ext_gene) %>% unlist() %>% na.omit()
    } else {
      tmp <- tmp  %>%
        filter(qval < 0.05 & b < 0 ) %>% select(ext_gene) %>% unlist() %>% na.omit()
    }
    assign(paste0("rna_res_", time, "vs0_", type), tmp)
    rna_res_all <- c(rna_res_all, tmp)
  }
}
rna_res_all <- unique(rna_res_all)

peak_res_all <- vector("list", 8)
names(peak_res_all) <- c("2vs0_up_proximal", "2vs0_down_proximal", "12vs0_up_proximal", 
                         "12vs0_down_proximal", "2vs0_up_distal", "2vs0_down_distal", 
                         "12vs0_up_distal", "12vs0_down_distal")
for(time in c(2,12)) {
  for(type in c("up", "down")) {
    res <- get(paste0("sigpeaks_", time, "vs0"))
    if(type == "up") {
      peak_res <- filter(res, log2FoldChange > 0) %>% 
        mutate(chr_old=seqnames)
    } else {
      peak_res <- filter(res, log2FoldChange < 0) %>%
        mutate(chr_old=seqnames)
    }
    peak_res <- peak_res %>%
       left_join(., alias, by="chr_old") %>%
       dplyr::select(chr, start, end, width, strand) %>%
       dplyr::mutate(seqnames=chr) %>%
       dplyr::select(seqnames, start, end) %>%
       arrange(seqnames, start)
    peak_res_output <- GRanges(peak_res) %>% unique()
    
    for(distance_type in c("proximal", "distal")) {
      peak_choice <- intersect(get(paste0(distance_type, "_peaks")), 
                                       peak_res_output)
      if(!length(peak_choice)) next;
      ## Retrieve gene names
      if(distance_type == "proximal") {
        br <- c(-1000,1000)
      } else {
        br <- c(-100000,100000)
      }
      peak_choice_withGene <- 
        annotatePeakInBatch(peak_choice, 
                            PeakLocForDistance = "middle",
                            FeatureLocForDistance = "TSS",
                            AnnotationData=ucsc.dr10.refGene, 
                            output="overlapping", bindingRegion=br) %>%
        addGeneIDs(orgAnn="org.Dr.eg.db", 
                   feature_id_type="entrez_id",
                   IDs2Add="symbol")
      peak_choice_withGene <- data.frame(peak_choice_withGene) %>%
        select(chr=seqnames, start, end, width, strand, distance, insideFeature, gene=symbol) %>%
        unite("new_seqnames", chr, start, sep = ":", remove=FALSE) %>% 
        unite("new_seqnames", new_seqnames, end, sep = "-", remove=FALSE) %>%
        select(chr, start, end, width, strand, distance, insideFeature, gene, new_seqnames) %>%
        mutate(overall_DE = ifelse(gene %in% rna_res_all, TRUE, FALSE)) %>%
        mutate(!!paste0(time, "_vs0_", type, "_DE") := 
                  ifelse(gene %in% get(paste0("rna_res_", time, "vs0_", type)), TRUE, FALSE)) %>%
        mutate(is_TF = ifelse(gene %in% detf, TRUE, FALSE))

      peak_res_all[[paste0(time, "vs0_", type, "_", distance_type)]] <- peak_choice_withGene
    }
  }
}

for (i in c(1:8)){
  if(i == 4) next;
  write_excel_csv(peak_res_all[i][[1]], 
             path=paste0("ANNOTATION/differential_peaklets_500bp_withGenes_proximal-distal_",
            names(peak_res_all)[i], ".csv"))
} 

## Write out bed files so that I can output sequences
alias <- read.table("danRer10_alias.tab", stringsAsFactors = FALSE)
colnames(alias) <- c("chr_old", "chr")
for (i in c("2vs0_up_proximal", "2vs0_up_distal", "12vs0_down_distal")) {
  tmp <- read.csv(paste0("ANNOTATION/differential_peaklets_500bp_withGenes_proximal-distal_",
                            i, ".csv")) %>%
    select(chr, start, end) %>%
    unique() %>%
    left_join(., alias, by="chr") %>%
    select(-chr) %>%
    select(chr=chr_old, start, end)
  print(dim(tmp))
  write.table(tmp, file=paste0("ANNOTATION/differential_peaklets_500bp_withGenes_proximal-distal_",
                            i, ".bed"), row.names = FALSE, col.names=FALSE,
              quote=FALSE, sep="\t")
} 

## Annotate peaks by differential status with closest DE TF or DE gene --------------------------------------

mcols(ucsc.dr10.refGene)$feature <- mcols(ucsc.dr10.refGene)$gene_id
ucsc.dr10.refGene_symbol <- addGeneIDs(ucsc.dr10.refGene, orgAnn="org.Dr.eg.db", 
                                       feature_id_type="entrez_id",
                                       IDs2Add="symbol")

peak_res_all_closest <- vector("list", 4)
names(peak_res_all_closest) <- c("2vs0_up", "2vs0_down", "12vs0_up",  "12vs0_down")

for(time in c(2,12)) {
  for(type in c("up", "down")) {
    res <- get(paste0("sigpeaks_", time, "vs0"))
    if(type == "up") {
      peak_res <- filter(res, log2FoldChange > 0) %>% 
        mutate(chr_old=seqnames)
    } else {
      peak_res <- filter(res, log2FoldChange < 0) %>%
        mutate(chr_old=seqnames)
    }
    peak_res <- peak_res %>%
      left_join(., alias, by="chr_old") %>%
      dplyr::select(chr, start, end, width, strand) %>%
      dplyr::mutate(seqnames=chr) %>%
      dplyr::select(seqnames, start, end) %>%
      arrange(seqnames, start)
    peak_res_output <- GRanges(peak_res) %>% unique()
    
    peak_res_df <- data.frame(peak_res_output, stringsAsFactors=FALSE) 
    
    for(gene_type in c("overall", "overall_DE", "matched_DE", "TF", "overall_TF_DE", "matched_TF_DE")) {
     if(gene_type == "overall") {
       ## Strictly closest gene
       gene_choice <- ucsc.dr10.refGene_symbol
     } else if(gene_type == "overall_DE") {
       ## Strictly closest DE gene
       gene_choice <- ucsc.dr10.refGene_symbol[which(ucsc.dr10.refGene_symbol$symbol %in% rna_res_all)]
     } else if(gene_type == "matched_DE") {
       ## Strictly closest matched DE gene
       gene_choice <- ucsc.dr10.refGene_symbol[which(ucsc.dr10.refGene_symbol$symbol %in% 
                                                       get(paste0("rna_res_", time, "vs0_", type)))]
     } else if(gene_type == "TF") {
       ## Strictly closest TF
       gene_choice <- ucsc.dr10.refGene_symbol[which(ucsc.dr10.refGene_symbol$symbol %in% detf)]
     } else if(gene_type == "overall_TF_DE") {
       ## Strictly closest DE TF
       diff_TF <- detf[which(detf %in% rna_res_all)]
       gene_choice <- ucsc.dr10.refGene_symbol[which(ucsc.dr10.refGene_symbol$symbol %in% diff_TF)]
     } else {
       ## Strictly closest matched DE TF
       diff_TF_matched <- detf[which(detf %in% get(paste0("rna_res_", time, "vs0_", type)))]
       gene_choice <- ucsc.dr10.refGene_symbol[which(ucsc.dr10.refGene_symbol$symbol %in% diff_TF_matched)]
     }
     peak_res_output_closest <- annotatePeakInBatch(peak_res_output, 
                                 PeakLocForDistance = "middle",
                                 FeatureLocForDistance = "TSS", 
                                 AnnotationData=gene_choice, 
                                 output="nearestLocation", select="first") %>%
       addGeneIDs(orgAnn="org.Dr.eg.db", feature_id_type="entrez_id",IDs2Add="symbol") %>%
       data.frame(stringsAsFactors = FALSE) %>%
       select(seqnames, start, end, width, strand, distancetoFeature, gene=symbol) 
     peak_res_df <- left_join(peak_res_df, peak_res_output_closest, by = c("seqnames", "start", "end", "width", "strand")) %>%
       mutate(!!paste0(gene_type, "_distance") := distancetoFeature) %>%
       mutate(!!paste0(gene_type, "_gene") := gene) %>%
       select(-distancetoFeature, -gene)
    }
    peak_res_df <- peak_res_df %>%
      unite("new_seqnames", seqnames, start, sep = ":", remove=FALSE) %>% 
      unite("new_seqnames", new_seqnames, end, sep = "-", remove=FALSE) %>%
      arrange(seqnames) %>%
      select(new_seqnames, chr=seqnames, everything())
    peak_res_all_closest[[paste0(time, "vs0_", type)]] <- peak_res_df
  }
}

options(scipen = 20)
for (i in c(1:4)){
  write_excel_csv(peak_res_all_closest[i][[1]], 
                  path=paste0("ANNOTATION/differential_peaklets_500bp_withGenes_closest_",
                              names(peak_res_all_closest)[i], ".csv"))
} 




# 
# tmp_tss <- data.frame(ucsc.dr10.tss) %>% 
#   arrange(seqnames, start, end) %>%
#   mutate(merge_col = start)
# tmp_gene <- data.frame(ucsc.dr10.refGene_symbol) %>%
#   arrange(seqnames, start, end) %>%
#   mutate(merge_col = ifelse(strand == "-", end+1, start)) %>%
#   left_join(., tmp_tss, by=c("seqnames", "merge_col"))
