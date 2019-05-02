library(tidyverse)
library(GenomicRanges)
library(ChIPpeakAnno)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(org.Dr.eg.db)
library(AnnotationHub)

## Load Danio Rerio Ensembl annotation (EnsDb) to annotate peaks with features
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Danio Rerio", "GRCz10", 89)) # Query for zebrafish
Drerio <- ahDb[["AH56671"]]  ## Ensembl 89 EnsDb for Danio Rerio  
annoData <- toGRanges(Drerio, feature="transcript")
annoData[1:2]

# Read in RNA-seq differential results -------------------------------------
allgenes <- read.csv(file = "RNAseq/4_vs0_model/results_LFC.csv",
                     check.names=FALSE)
allgenes <- allgenes[,-1]
sigtranscripts <- allgenes[which(rowSums(abs(allgenes[,-c(1:2)])) > 0),] %>%
  unite(col="pattern", 3:6, sep=",") %>%
  dplyr::select(target_id, pattern) %>%
  dplyr::filter(target_id != "TRANSGENE")
sigtranscripts$target_id <- as.character(sigtranscripts$target_id)

## There are 11 missing transcripts here (NB: they are from Ensembl 85)
sigtranscripts_annot <- annoData[which(mcols(annoData)$tx_id %in% sigtranscripts$target_id)] 
sigtranscripts[which(!sigtranscripts$target_id %in% mcols(annoData)$tx_id),]
sigtranscripts_annot2 <- sigtranscripts_annot %>% as.data.frame() %>% dplyr::select(seqnames, start, end, target_id = tx_id) %>%
  left_join(., sigtranscripts, by="target_id") %>% dplyr::select(seqnames, start, end, pattern)
write.table(sigtranscripts_annot2, file="IGV/sigtranscripts_vs0.bed", col.names=FALSE, quote=FALSE, row.names=FALSE)

# Read in ATAC-seq data ----------------------------------------------------
alias <- read.table("ATACseq/danRer10_alias.tab")
colnames(alias) <- c("chr_old", "chr")
allpeaks_2vs0 <- readRDS("ATACseq/DIFFERENTIAL_vs0/allpeaks_pvalsort_2vs0.rds")
allpeaks_4vs0 <- readRDS("ATACseq/DIFFERENTIAL_vs0/allpeaks_pvalsort_4vs0.rds")
## No differential peaks at 7 vs 0
allpeaks_12vs0 <- readRDS("ATACseq/DIFFERENTIAL_vs0/allpeaks_pvalsort_12vs0.rds")
allpeaks <- data.frame(seqnames=allpeaks_2vs0$seqnames, start=allpeaks_2vs0$start, end=allpeaks_2vs0$end,
                       `2pi-0dpi`=ifelse(allpeaks_2vs0$padj < 0.05, sign(allpeaks_2vs0$log2FoldChange), 0),
                       `4pi-0dpi`=ifelse(allpeaks_4vs0$padj < 0.05, sign(allpeaks_4vs0$log2FoldChange), 0),
                       `7pi-0dpi`= 0,
                       `12pi-0dpi`=ifelse(allpeaks_12vs0$padj < 0.05, sign(allpeaks_12vs0$log2FoldChange), 0),
                       check.names=FALSE)
allpeaks[is.na(allpeaks)] <- 0
sigpeaks <- allpeaks %>% 
  unite(col="pattern", 4:7, sep=",") %>%
  dplyr::select(seqnames, start, end, pattern) %>%
  dplyr::filter(seqnames != "99999") %>%
  dplyr::filter(pattern != "0,0,0,0")
write.table(sigpeaks, file="IGV/sigpeaks_vs0.bed", col.names=FALSE, quote=FALSE, row.names=FALSE)

