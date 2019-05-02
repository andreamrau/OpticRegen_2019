library(viridis)
library(cowplot)
library(scales)
# Grcz10 (danRer10)
library(ChIPpeakAnno)
library(org.Dr.eg.db)
library(GenomicFeatures)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(BSgenome.Drerio.UCSC.danRer10)
library(rtracklayer)
library(AnnotationHub)
library(forcats)
library(tidyverse)
library(readxl)
library(AnnotationDbi)
library(tidyverse)


##-------------------------------------------------------------------
## Checking final numbers
##-------------------------------------------------------------------

## Total transcripts represented ------------------------------------

load("RNAseq/4_vs0_model.RData")

## sleuth:  The default is to filter out any features that do not have at least 5 
##    estimated counts in at least 47% of samples
results_table_vs0 <- vector("list", 4)
names(results_table_vs0) <- c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")
for(i in c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")) {
  results_table_vs0[[i]] <- read.csv(paste0("RNAseq/4_vs0_model/results_", i, ".csv"))
}

## Total number of transcripts = 57867
lapply(results_table_vs0, nrow)

## Total number of transcripts that pass the sleuth filter = 38494
lapply(results_table_vs0, function(x) sum(!is.na(x$qval)))

## 1. Number of expressed transcripts that are not DE at any time point = 34358
res <- full_join(results_table_vs0[["2dpi-0dpi"]], results_table_vs0[["4dpi-0dpi"]], by = "target_id") %>%
  full_join(., results_table_vs0[["7dpi-0dpi"]], by = "target_id") %>%
  full_join(., results_table_vs0[["12dpi-0dpi"]], by = "target_id") %>%
  select(target_id, qval.2dpi = qval.x, b.2dpi = b.x, qval.4dpi = qval.y, b.4dpi = b.y,
         qval.7dpi = qval.x.x, b.7dpi = b.x.x, qval.12dpi = qval.y.y, b.12dpi = b.y.y)

res %>% filter(qval.2dpi >= 0.05, qval.4dpi >= 0.05, qval.7dpi >= 0.05, qval.12dpi >= 0.05) %>% nrow()

## Output their ids
write.table(res %>% filter(qval.2dpi >= 0.05, qval.4dpi >= 0.05, qval.7dpi >= 0.05, qval.12dpi >= 0.05) %>% select(target_id),
            "expressed_NDE_transcripts.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

## 2. Recheck DE up, down, no change numbers at each time point wrz

##   => 2dpi-0dpi
length(which(results_table_vs0[["2dpi-0dpi"]]$qval < 0.05))
length(which(results_table_vs0[["2dpi-0dpi"]]$qval < 0.05 & 
               results_table_vs0[["2dpi-0dpi"]]$b > 0))
length(which(results_table_vs0[["2dpi-0dpi"]]$qval < 0.05 & 
               results_table_vs0[["2dpi-0dpi"]]$b < 0))
length(which(results_table_vs0[["2dpi-0dpi"]]$qval >= 0.05))

##   => 4dpi-0dpi
length(which(results_table_vs0[["4dpi-0dpi"]]$qval < 0.05))
length(which(results_table_vs0[["4dpi-0dpi"]]$qval < 0.05 & 
               results_table_vs0[["4dpi-0dpi"]]$b > 0))
length(which(results_table_vs0[["4dpi-0dpi"]]$qval < 0.05 & 
               results_table_vs0[["4dpi-0dpi"]]$b < 0))
length(which(results_table_vs0[["4dpi-0dpi"]]$qval >= 0.05))

##   => 7dpi-0dpi
length(which(results_table_vs0[["7dpi-0dpi"]]$qval < 0.05))
length(which(results_table_vs0[["7dpi-0dpi"]]$qval < 0.05 & 
               results_table_vs0[["7dpi-0dpi"]]$b > 0))
length(which(results_table_vs0[["7dpi-0dpi"]]$qval < 0.05 & 
               results_table_vs0[["7dpi-0dpi"]]$b < 0))
length(which(results_table_vs0[["7dpi-0dpi"]]$qval >= 0.05))

##   => 12dpi-0dpi
length(which(results_table_vs0[["12dpi-0dpi"]]$qval < 0.05))
length(which(results_table_vs0[["12dpi-0dpi"]]$qval < 0.05 & 
               results_table_vs0[["12dpi-0dpi"]]$b > 0))
length(which(results_table_vs0[["12dpi-0dpi"]]$qval < 0.05 & 
               results_table_vs0[["12dpi-0dpi"]]$b < 0))
length(which(results_table_vs0[["12dpi-0dpi"]]$qval >= 0.05))

##   => DE at any time
res %>% filter(qval.2dpi < 0.05 | qval.4dpi < 0.05 | qval.7dpi < 0.05 | qval.12dpi < 0.05) %>% nrow()
res %>% filter(qval.2dpi < 0.01 | qval.4dpi < 0.01 | qval.7dpi < 0.01 | qval.12dpi < 0.01) %>% nrow()


## Check the LRT factor overall value
load("RNAseq/2_splines_model.RData")
length(which(results_table_factor$qval < 0.01))

## Numbers of ATAC peaklets ------------------------------------
##    (all peaklets including both differential and nondifferential)

## Total number of peaklets: 215000
macsBed <- read.delim(paste0("ATACseq/PEAKS_TRANS_PEAKLETS_ALLTIMES/", 
                             "ATAC.nodup.unique.macs.peaklets_peaks.narrowPeak"), 
                      header=FALSE)
nrow(macsBed)

## Total number of peaklets with p-value < 10^-10: 44429
length(which(macsBed$V8 > -log(10^-10)))

## Total number of peaklets after merging overlaps: 42198 (includes 3 on transgene, so 42195 without transgene)
load("ATACseq/DIFFERENTIAL_vs0/differential-analysis_diffBind_pvalsort_trans_peaklets_500bp.RData")

## Total number of differential peaklets: 233
## Total number of non-differential peaklets: 41965

## Load Danio Rerio UCSC annotation to anntate peaks with features
ucsc.dr10.refGene <- genes(TxDb.Drerio.UCSC.danRer10.refGene)
ucsc.dr10.exons <- exons(TxDb.Drerio.UCSC.danRer10.refGene)
ucsc.dr10.tss <- promoters(TxDb.Drerio.UCSC.danRer10.refGene, upstream = 0, downstream = 0)
ucsc.dr10.transcripts <- transcripts(TxDb.Drerio.UCSC.danRer10.refGene)
ucsc.dr10.fiveUTRs <- fiveUTRsByTranscript(TxDb.Drerio.UCSC.danRer10.refGene)
ucsc.dr10.introns <- intronsByTranscript(TxDb.Drerio.UCSC.danRer10.refGene)

## Full genome size and relative size of each thing
genome <- BSgenome.Drerio.UCSC.danRer10
total_size <- sum(seqlengths(genome))
utr_size <- width(GenomicRanges::reduce(ucsc.dr10.fiveUTRs)) %>% unlist() %>% sum
exon_size <- width(setdiff(as.data.frame(ucsc.dr10.exons) %>%
                             dplyr::select(seqnames, start, end) %>%
                             unique() %>%
                             GRanges, 
                           as.data.frame(ucsc.dr10.fiveUTRs) %>%
                             dplyr::select(seqnames, start,end) %>%
                             unique() %>%
                             GRanges)) %>% unlist() %>% sum
intron_size <- setdiff(setdiff(ucsc.dr10.introns %>% as.data.frame() %>%
  dplyr::select(seqnames, start, end, width) %>%
  unique() %>% GRanges(), as.data.frame(ucsc.dr10.exons) %>%
    dplyr::select(seqnames, start, end) %>%
    unique() %>%
    GRanges), as.data.frame(ucsc.dr10.fiveUTRs) %>%
    dplyr::select(seqnames, start,end) %>%
    unique() %>%
    GRanges) %>%
  width() %>%
  unlist() %>%
  sum()
intergenic_size <- total_size - utr_size - exon_size - intron_size

  
alias <- read.table("RNAseq/danRer10_alias.tab", stringsAsFactors = FALSE)
colnames(alias) <- c("Chr", "chr")
peaks <- left_join(rowData_pvalsort, alias, by = "Chr") %>%
  dplyr::select(chr, start = Start, end = End) %>% na.omit() %>% GRanges()

peaks_overlapping_tss <- findOverlaps(peaks, ucsc.dr10.tss) %>% 
  as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()
peaks_overlapping_fiveUTRs <-  findOverlaps(peaks, ucsc.dr10.fiveUTRs) %>% 
  as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()
peaks_overlapping_exons <-  findOverlaps(peaks, ucsc.dr10.exons) %>% 
  as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()
peaks_overlapping_genes <-  findOverlaps(peaks, ucsc.dr10.refGene) %>% 
  as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()

peaks_in_1kb <- findOverlaps(peaks, ucsc.dr10.tss, maxgap = 1000) %>% 
  as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()
peaks_in_100kb <- findOverlaps(peaks, ucsc.dr10.tss, maxgap = 100000) %>% 
  as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()


peaks_in_5kb <- findOverlaps(peaks, ucsc.dr10.tss, maxgap = 5000) %>% 
  as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()
peaks_in_10kb <- findOverlaps(peaks, ucsc.dr10.tss, maxgap = 10000) %>% 
  as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()


  
## 1. Overlapping UTR: 7928
length(unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs)))
##    TSS: 7626
length(peaks_overlapping_tss)
##    nonTSS: 302
length(setdiff(peaks_overlapping_fiveUTRs, peaks_overlapping_tss))

## 2. Overlapping coding exon: 2188
length(setdiff(peaks_overlapping_exons, unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs))))

## 3. Overlapping intron: 7376
length(setdiff(peaks_overlapping_genes, 
               unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons))))


## 4. Intergenic: 24703
length(setdiff(1:length(peaks), 
               unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                        peaks_overlapping_genes))))
##    within 0-1kb from TSS: 787
length(setdiff(peaks_in_1kb, 
               unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                        peaks_overlapping_genes))))
##    1-100 kb from TSS: 16925
length(setdiff(peaks_in_100kb, 
               unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                        peaks_overlapping_genes, peaks_in_1kb))))
##    >100 kb from TSS: 6991
length(setdiff(1:length(peaks), 
               unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                        peaks_overlapping_genes, peaks_in_1kb, peaks_in_100kb))))

## 1 kb-5kb from TSS: 1753
length(setdiff(setdiff(peaks_in_5kb,
               unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                                     peaks_overlapping_genes))),
  setdiff(peaks_in_1kb, 
               unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                        peaks_overlapping_genes)))))
## 5kb-10kb from TSS: 1709
length(setdiff(setdiff(peaks_in_10kb,
                       unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                                peaks_overlapping_genes))),
               setdiff(peaks_in_5kb, 
                       unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                                peaks_overlapping_genes)))))
## 10kb-100kb from TSS: 13463
length(setdiff(setdiff(peaks_in_100kb,
                       unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                                peaks_overlapping_genes))),
               setdiff(peaks_in_10kb, 
                       unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons,
                                peaks_overlapping_genes)))))


## Can you send me the list of these peaks, whether they overlap UTR, exon, or intron, what gene, and whether the gene has a DE transcript or not?

## Cue up differential results
DE_results <- left_join(dplyr::select(results_table_vs0[["2dpi-0dpi"]], target_id, qval2=qval, ext_gene), 
                        dplyr::select(results_table_vs0[["4dpi-0dpi"]], target_id, qval4=qval, ext_gene),
                        by = c("target_id", "ext_gene")) %>%
  left_join(.,  dplyr::select(results_table_vs0[["7dpi-0dpi"]], target_id, qval7=qval, ext_gene),
            by = c("target_id", "ext_gene")) %>%
  left_join(.,  dplyr::select(results_table_vs0[["12dpi-0dpi"]], target_id, qval12=qval, ext_gene),
            by = c("target_id", "ext_gene")) %>%
  mutate(DE = ifelse(qval2 < 0.05 | qval4 < 0.05 | qval7 < 0.05 | qval12 < 0.05, 1, 0)) %>%
  dplyr::select(ext_gene, DE) %>%
  na.omit() %>%
  group_by(ext_gene) %>%
  mutate(total_DE = sum(DE)) %>%
  dplyr::select(gene= ext_gene, total_DE) %>%
  unique() %>%
  mutate(total_DE = sign(total_DE)) %>%
  ungroup()

## Match up gene names to TSS and UTR data 
ucsc.dr10.tss_df <- data.frame(promoters(TxDb.Drerio.UCSC.danRer10.refGene, upstream = 0, downstream = 0, columns = "gene_id"))
ucsc.dr10.tss_df$gene_id <- as.character(ucsc.dr10.tss_df$gene_id)
tss_df <- mapIds(org.Dr.eg.db, unlist(ucsc.dr10.tss_df$gene_id), "SYMBOL", "ENTREZID")
tss_df <- left_join(ucsc.dr10.tss_df, data.frame(gene_id = as.character(names(tss_df)), gene = tss_df, stringsAsFactors = FALSE), by = "gene_id") %>%
  dplyr::select(seqnames, start, end, gene)
ucsc.dr10.utr_df <- data.frame(fiveUTRsByTranscript(TxDb.Drerio.UCSC.danRer10.refGene, use.names=TRUE)) %>%
  mutate(REFSEQ = group_name)
utr_df <- mapIds(org.Dr.eg.db, unlist(ucsc.dr10.utr_df$REFSEQ), "SYMBOL", "REFSEQ")
utr_df <- left_join(ucsc.dr10.utr_df, data.frame(REFSEQ = as.character(names(utr_df)), gene = utr_df, stringsAsFactors = FALSE), by = "REFSEQ") %>%
  dplyr::select(seqnames, start, end, gene)
utr_tss_df <- bind_rows(tss_df, utr_df) %>% unique() %>% GRanges()
tmp <- findOverlaps(peaks, utr_tss_df) %>% as.data.frame()
utr_tss_withGenes <- data.frame(data.frame(peaks, row.names=NULL)[tmp$queryHits,], gene=utr_tss_df$gene[tmp$subjectHits], row.names=NULL) %>%
  unique() %>%
  left_join(., DE_results, by = "gene")

## Match up gene names to exon data
ucsc.dr10.exon_df <- data.frame(exons(TxDb.Drerio.UCSC.danRer10.refGene, columns = "gene_id"))
ucsc.dr10.exon_df[grep("c", as.character(ucsc.dr10.exon_df$gene_id)),"gene_id"] <- 
  unlist(lapply(as.list(ucsc.dr10.exon_df$gene_id[grep("c", as.character(ucsc.dr10.exon_df$gene_id))]), 
                            function(x) paste0(x, collapse=",")))
ucsc.dr10.exon_df$gene_id <- as.character(ucsc.dr10.exon_df$gene_id)
ucsc.dr10.exon_df <- ucsc.dr10.exon_df %>%
  separate_rows(gene_id, sep=",")
exon_df <- mapIds(org.Dr.eg.db, unlist(ucsc.dr10.exon_df$gene_id), "SYMBOL", "ENTREZID")
exon_df <- left_join(ucsc.dr10.exon_df, data.frame(gene_id = as.character(names(exon_df)), gene = exon_df, 
                                                    stringsAsFactors = FALSE), by = "gene_id") %>%
  dplyr::select(seqnames, start, end, gene) %>% unique() %>% GRanges()
exon_peaks <- peaks[setdiff(peaks_overlapping_exons, unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs)))]

tmp <- findOverlaps(exon_peaks, exon_df) %>% as.data.frame()
exon_withGenes <- data.frame(data.frame(exon_peaks, row.names=NULL)[tmp$queryHits,], 
                             gene=exon_df$gene[tmp$subjectHits], row.names=NULL) %>%
  unique() %>%
  left_join(., DE_results, by = "gene")

## Match up gene names with intron data
ucsc.dr10.gene_df <- data.frame(ucsc.dr10.refGene)
ucsc.dr10.gene_df$gene_id <- as.character(ucsc.dr10.gene_df$gene_id)
gene_df <- mapIds(org.Dr.eg.db, unlist(ucsc.dr10.gene_df$gene_id), "SYMBOL", "ENTREZID")
gene_df <- left_join(ucsc.dr10.gene_df, data.frame(gene_id = as.character(names(gene_df)), gene = gene_df, 
                                                   stringsAsFactors = FALSE), by = "gene_id") %>%
  dplyr::select(seqnames, start, end, gene) %>% unique() %>% GRanges()
intron_peaks <- peaks[setdiff(peaks_overlapping_genes, 
                              unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs, peaks_overlapping_exons)))]

tmp <- findOverlaps(intron_peaks, gene_df) %>% as.data.frame()
intron_withGenes <- data.frame(data.frame(intron_peaks, row.names=NULL)[tmp$queryHits,], 
                             gene=gene_df$gene[tmp$subjectHits], row.names=NULL) %>%
  unique() %>%
  left_join(., DE_results, by = "gene")


## Write out these tables
write.table(utr_tss_withGenes %>% dplyr::select(-width, -strand) %>% mutate(total_DE = ifelse(total_DE == 1, "DE", "NDE")), 
            "utr_peaks_withGenes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(exon_withGenes %>% dplyr::select(-width, -strand) %>% mutate(total_DE = ifelse(total_DE == 1, "DE", "NDE")), "exon_peaks_withGenes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(intron_withGenes %>% dplyr::select(-width, -strand) %>% mutate(total_DE = ifelse(total_DE == 1, "DE", "NDE")), "intron_peaks_withGenes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


## Histogram of distance to nearest TSS ------------------------------------

# distance <- distanceToNearest(peaks, ucsc.dr10.tss) %>% as.data.frame()
# nr <- nearest(peaks, ucsc.dr10.tss)
# dt <- distance(peaks[which(!is.na(nr))], ucsc.dr10.tss[na.omit(nr)])

diff_peaks <- read.table("ATACseq/DIFFERENTIAL_vs0/sigpeaks_all_peaklets_500bp.bed")
colnames(diff_peaks) <- c("Chr", "start", "end") 
diff_peaks <- diff_peaks %>%
  left_join(., alias, by = "Chr") %>%
  dplyr::select(seqnames = chr, start, end) %>%
  mutate(status = "DO")

gr.anno <- annotatePeakInBatch(peaks, AnnotationData=ucsc.dr10.tss, 
                               output="nearestLocation", ignore.strand=FALSE) 
gr.anno2 <- gr.anno %>%
  as.data.frame() %>%
  dplyr::select(seqnames, start, end, distancetoFeature, start_position) %>%
  mutate(distancetoFeature2 = distancetoFeature) %>%
#           ifelse(start_position >= start & start_position <= end, 0, distancetoFeature)) %>%
  left_join(., diff_peaks, by = c("seqnames", "start", "end")) %>%
  mutate(status = ifelse(is.na(status), "NDO", status)) %>%
  filter(!is.na(distancetoFeature2)) 



fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

main.plot <- ggplot(gr.anno2, aes(x=distancetoFeature2, color = status, fill = status)) +
  geom_density(alpha = 0.5) +
  geom_rug() +
  scale_x_continuous(name = "Distance to closest TSS", labels = comma) +
  scale_y_continuous(labels = fancy_scientific) +
  scale_fill_viridis(discrete=TRUE, begin = 0, end = 0.8) +
  scale_color_viridis(discrete=TRUE, begin = 0, end = 0.8) +
  theme_classic() +
  guides(color=FALSE, fill=FALSE)+
  facet_grid(~status,  labeller=labeller(status = c(DO = "Differentially open peaks", 
                                                    NDO = "Open peaks")))


inset.plot <- ggplot(filter(gr.anno2, status == "DO"), 
                     aes(x=distancetoFeature2, color = status, fill = status)) +
  geom_density(alpha = 0.5) +
  geom_rug() +
  scale_x_continuous(name = "", labels = comma, limits = c(-200000, 200000)) +
  scale_y_continuous(name = "", labels = fancy_scientific) +
  scale_fill_viridis(discrete=TRUE, begin = 0, end = 0.8) +
  scale_color_viridis(discrete=TRUE, begin = 0, end = 0.8) +
  theme_bw() +
  guides(color=FALSE, fill=FALSE) + 
  theme(axis.title.y=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 8))


inset.plot2 <- ggplot(filter(gr.anno2, status != "DO"), 
                     aes(x=distancetoFeature2, color = status, fill = status)) +
  geom_density(alpha = 0.5) +
  geom_rug() +
  scale_x_continuous(name = "", labels = comma, limits = c(-200000, 200000)) +
  scale_y_continuous(name = "", labels = fancy_scientific) +
  scale_fill_viridis(discrete=TRUE, begin = 0.8, end = 0) +
  scale_color_viridis(discrete=TRUE, begin = 0.8, end = 0) +
  theme_bw() +
  guides(color=FALSE, fill=FALSE) + 
  theme(axis.title.y=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 8))

plot.with.inset <-
  ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = 0.1, y = .5, width = .25, height = .4) +
  draw_plot(inset.plot2, x = 0.5, y = .5, width = .25, height = .4)
plot.with.inset


## Count overlap with PEGASUS predicted enhancers ------------------------------------

## 54872 predicted enhancers (average length = 128bp, minimum length = 10bp, median length = 90bp, maximum length = 1628bp)
pegasus <- read.table("PEGASUS_predicted_enhancers/danRer10_liftOver_CNEs_PEGASUS.data.bed", sep = "\t")
colnames(pegasus) <- c("chr", "start", "end", "enhancer_name")
quantile(pegasus$end - pegasus$start + 1)
mean(pegasus$end - pegasus$start + 1)

## 3109 peaks (7%) overlapping 4859 predicted enhancers (9%)
peaks_overlapping_pegasus <- findOverlaps(peaks, GRanges(pegasus)) %>% 
  as.data.frame
length(unlist(unique(peaks_overlapping_pegasus$queryHits)))
length(unlist(unique(peaks_overlapping_pegasus$subjectHits)))

## 4339 peaks (10%) within 1000bp of 7603 predicted enhancers (14%)
peaks_closeto_pegasus <- findOverlaps(peaks, GRanges(pegasus), maxgap = 1000) %>% 
  as.data.frame
length(unlist(unique(peaks_closeto_pegasus$queryHits)))
length(unlist(unique(peaks_closeto_pegasus$subjectHits)))

## 7 DO peaks (3%) overlapping 8 predicted enhancers
DOpeaks_overlapping_pegasus <- findOverlaps(GRanges(diff_peaks), GRanges(pegasus)) %>% 
  as.data.frame
length(unlist(unique(DOpeaks_overlapping_pegasus$queryHits)))
length(unlist(unique(DOpeaks_overlapping_pegasus$subjectHits)))

## 9 DO peaks (4%) within 1000bp of 15 predicted enhancers
DOpeaks_closeto_pegasus <- findOverlaps(GRanges(diff_peaks), GRanges(pegasus), maxgap = 1000) %>% 
  as.data.frame
length(unlist(unique(DOpeaks_closeto_pegasus$queryHits)))
length(unlist(unique(DOpeaks_closeto_pegasus$subjectHits)))



## Categorize TSS peaks wrt genes ------------------------------------

tss_peaks <- peaks[unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs))] %>%
  as.data.frame() 
tss_peaks$type <- "TSS_5primeUTR"
kb1_peaks <- peaks[setdiff(peaks_in_1kb, 
                           unique(c(peaks_overlapping_tss, peaks_overlapping_fiveUTRs)))] %>%
  as.data.frame()
kb1_peaks$type <- "0-1kb"

transcripts_to_genes_map <- readRDS("Shiny_app/data/transcripts_to_genes_map.rds")
peaks_near_tss <- bind_rows(tss_peaks, kb1_peaks) %>%
  select(-width, -strand) %>%
  left_join(., diff_peaks, by = c("seqnames", "start", "end")) %>%
  mutate(status = ifelse(is.na(status), "NDO", status)) %>%
  GRanges() %>%
  annotatePeakInBatch(., AnnotationData=GRanges(tss_df), 
                      output="nearestLocation", ignore.strand=FALSE) %>%
  as.data.frame() %>%
  mutate(feature = as.numeric(feature))
peaks_near_tss$ext_gene <- tss_df$gene[peaks_near_tss$feature]
peaks_near_tss <- peaks_near_tss %>%
  select(-fromOverlappingOrNearest, -shortestDistance, -insideFeature,
         -feature_strand, -feature, -peak, -strand, -start_position, -end_position,
         -width) %>%
  left_join(., transcripts_to_genes_map, by = "ext_gene") %>%
  select(-ens_gene)
res_DE_genes <- res %>%
  mutate(overall_DE = ifelse(is.na(qval.2dpi) , "not transcribed",
    ifelse(qval.2dpi < 0.05 | qval.4dpi < 0.05 | qval.7dpi < 0.05 | qval.12dpi < 0.05,
                             "differential", "expressed"))) %>%
  select(target_id, overall_DE)
peaks_near_tss_withDE <- left_join(peaks_near_tss, res_DE_genes, by = "target_id") %>%
  group_by(seqnames, start, end, ext_gene, type, status, distancetoFeature) %>%
  summarize(DE_status = toString(overall_DE), target_ids = toString(target_id)) %>%
  ungroup()
DE_index <- grep("differential", peaks_near_tss_withDE$DE_status)
NDE_index <- setdiff(grep("expressed", peaks_near_tss_withDE$DE_status),
                     DE_index)
notexpressed_index <- setdiff(c(grep("NA", peaks_near_tss_withDE$DE_status),
                        grep("not transcribed", peaks_near_tss_withDE$DE_status)),
                        c(DE_index, NDE_index))
peaks_near_tss_withDE$overall_DE <- rep("NA", nrow(peaks_near_tss_withDE))
peaks_near_tss_withDE$overall_DE[DE_index] <- "differential"
peaks_near_tss_withDE$overall_DE[NDE_index] <- "expressed"
peaks_near_tss_withDE$overall_DE[notexpressed_index] <- "not transcribed"
peaks_near_tss_withDE <- peaks_near_tss_withDE %>%
  mutate(peak_ID = paste0(">", seqnames, ":", start, "-", end)) %>%
  select(peak_ID, seqnames, start, end, status, type, ext_gene, overall_DE,
         distancetoFeature,
         target_ids, DE_status)


write.table(peaks_near_tss_withDE, "peaks_near_tss_withDE.txt", 
            sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)

