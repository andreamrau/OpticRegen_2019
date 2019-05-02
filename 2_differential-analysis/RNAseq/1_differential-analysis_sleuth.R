library(sleuth)
library(dplyr)
library(tidyr)
library(data.table)
library(splines)
library(biomaRt)
library(VennDiagram)
source("my_sleuth_functions.R")
#http://www.nxn.se/valent/timecourse-analysis-with-sleuth

# Read in data, construct sleuth object ----------------------------------------
base_dir <- "~/Desktop/Optic-Regen_results/RNAseq"
sample_id <- dir(file.path(base_dir,"KALLISTO_hdf5"))
## Only keep v2 samples (corresponds to when we added transgene to reference genome)
sample_id <- sample_id[grep("v2", sample_id)]
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "KALLISTO_hdf5", id))
kal_dirs
targetIds_to_transcripts <- read.table("danio_transcripts_transgene.fa.tlst",
                                       stringsAsFactors = FALSE, fill=TRUE)
colnames(targetIds_to_transcripts) <- c("target_id", "transcript", 
                                        "chr", "pos")

tmp <- strsplit(sample_id, split="_", fixed=TRUE) %>%
  lapply(function(x) x[3]) %>%
  unlist
timept <- strsplit(tmp, split="RNA") %>%
  lapply(function(x) x[1]) %>%
  unlist
replicate <- strsplit(tmp, split="RNA") %>%
  lapply(function(x) x[2]) %>%
  unlist
expdesign <- data.frame(sample=sample_id, 
                        time=factor(timept, 
                                           levels=c(0,2,4,7,12)),
                        time_numeric=as.numeric(timept),
                        replicate=replicate,
                        path=kal_dirs,
                        stringsAsFactors=FALSE) 
## Reorder samples
expdesign <- expdesign[c(1:3,7:15,4:6),]

# Insert zebrafish gene names ---------------------------------------------
run_biomart <- FALSE
if(run_biomart) {
# ensembl <- useMart("ensembl")
#  ensembl_dt <- useDataset("drerio_gene_ensembl", mart=ensembl)
# host="useast.ensembl.org", "uswest.ensembl.org"
  date <- "aug2017" ## Ensembl 90
  ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                     host=paste(date,"archive.ensembl.org",sep="."),
                     dataset="drerio_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                          "external_gene_name"), 
                        mart = ensembl)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)  
  saveRDS(t2g, "transcripts_to_genes_map.rds")
} else {
  t2g <- readRDS("transcripts_to_genes_map.rds")
}

gene_name_map <- read.table("gene_name_map.txt", header=FALSE)

# Natural spline model ---------------------------------------------
## Design matrix
full_design <- model.matrix(formula(~ ns(expdesign$time_numeric, df = 3)))
full_design
## Load parameters
so <- my_sleuth_prep(sample_to_covariates = expdesign, 
                     full_model= full_design,
                     num_cores = 1, 
                     target_mapping = t2g,
                     targetIds_to_transcripts = targetIds_to_transcripts,
                     extra_bootstrap_summary = TRUE,
                     read_bootstrap_tpm = TRUE,
                     max_bootstrap = 500)
## Fit full model
so <- sleuth_fit(so)
## Fit reduced model
so <- sleuth_fit(so, ~1, 'reduced')
## Perform the test
so <- sleuth_lrt(so, 'reduced', 'full')
results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
sleuth_significant <- dplyr::filter(results_table, qval <= 0.05)
head(sleuth_significant, 20)
dim(sleuth_significant)

sleuth_live(so)
sleuth_deploy(so, 
              base_dir = "/Users/raua/Desktop/Optic-Regen_results/RNAseq/2_splines_model/",
              overwrite=TRUE)
write.csv(results_table, file = paste0("2_splines_model/results_splines_model.csv"),
            row.names=FALSE)

# Time factor model with LRT ---------------------------------------------
full_design_factor <- model.matrix(formula(~ expdesign$time))
colnames(full_design_factor) <- c("0dpi", "2dpi-0dpi", "4dpi-0dpi", "7dpi-0dpi", "12dpi-0dpi")
full_design_factor
## Load parameters
so_factor <- my_sleuth_prep(sample_to_covariates = expdesign, 
                            full_model= full_design_factor,
                            num_cores = 1, 
                            target_mapping = t2g,
                            targetIds_to_transcripts = targetIds_to_transcripts,
                            extra_bootstrap_summary = TRUE,
                            read_bootstrap_tpm = TRUE,
                            max_bootstrap = 500)
## Fit full model
so_factor <- sleuth_fit(so_factor, fit_name="full")
## Likelihood ratio test for all comparisons 
so_factor <- sleuth_fit(so_factor, ~1, 'reduced')
so_factor <- sleuth_lrt(so_factor, 'reduced', 'full')
results_table_factor <- sleuth_results(so_factor, 'reduced:full', test_type = 'lrt')
sleuth_significant_factor <- dplyr::filter(results_table_factor, qval <= 0.05)
head(sleuth_significant_factor, 20)
dim(sleuth_significant_factor)

## Venn diagram of results between spline and factor model?
venn.diagram(list(factor = sleuth_significant_factor$target_id,
                  spline = sleuth_significant$target_id), imagetype="png", filename="factor_vs_spline_DE.png",
             lwd = 2,
             fill = c("cornflowerblue", "darkorchid1"),
             alpha = 0.75,
             label.col = "white",
             fontfamily = "serif",
             fontface = "bold",
             cat.col = c("cornflowerblue", "darkorchid1"),
             cat.fontfamily = "serif",
             cat.fontface = "bold",
             cat.dist = c(0.03, 0.03))

## Save results
save.image("2_splines_model.RData")

#################################################
## ADDED May 4, 2018
#################################################

# Time factor model, with all comparisons back to 0 ---------------------------------------------
## Load parameters
so_factor_vs0 <- my_sleuth_prep(sample_to_covariates = expdesign,
                     full_model= full_design_factor,
                     num_cores = 1,
                     target_mapping = t2g,
                     targetIds_to_transcripts = targetIds_to_transcripts,
                     extra_bootstrap_summary = TRUE,
                     read_bootstrap_tpm = TRUE,
                     max_bootstrap = 500)
## Fit full model
so_factor_vs0 <- sleuth_fit(so_factor_vs0, fit_name="full")
## Wald test for each comparison
so_factor_vs0 <- sleuth_wt(so_factor_vs0, "4dpi-0dpi", which_model="full")
so_factor_vs0 <- sleuth_wt(so_factor_vs0, "2dpi-0dpi", which_model="full")
so_factor_vs0 <- sleuth_wt(so_factor_vs0, "7dpi-0dpi", which_model="full")
so_factor_vs0 <- sleuth_wt(so_factor_vs0, "12dpi-0dpi", which_model="full")

for(i in c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")) {
  results_table_vs0 <- sleuth_results(so_factor_vs0, i, test_type = 'wt')
  write.csv(results_table_vs0, file = paste0("4_vs0_model/results_", i, ".csv"),
            row.names=FALSE)
  results_table_vs0_select <- results_table_vs0[which(results_table_vs0$qval < 0.05),]
  results_table_vs0_select <- results_table_vs0_select[,c("target_id", "ext_gene", 
                                                          "b", "qval")]
  colnames(results_table_vs0_select) <- c("target_id", "ext_gene", "LFC", "padj")
  write.csv(results_table_vs0_select, file = paste0("4_vs0_model/results_", i, "_IPA.csv"),
            row.names=FALSE)
}

save.image("4_vs0_model.RData")


#################################################
## ADDED May 31, 2018
#################################################

load("4_vs0_model.RData")

LFC_vs0 <- vector("list", length = 4)
names(LFC_vs0) <- c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")
for(i in names(LFC_vs0)) {
  results_table_vs0 <- sleuth_results(so_factor_vs0, i, test_type = 'wt')
  tmp <- data.frame(target_id=results_table_vs0$target_id, 
                    ext_gene = results_table_vs0$ext_gene,
                    logFC=results_table_vs0$b, padj=results_table_vs0$qval)
  tmp$logFC <- ifelse(tmp$padj > 0.05, 0, ifelse(tmp$logFC > 0, 1, -1))
  tmp$padj <- NULL
  tmp$logFC[is.na(tmp$logFC)] <- 0
  colnames(tmp) <- c("target_id", "ext_gene", i)
  LFC_vs0[[i]] <- tmp
}
LFC_vs0 <- full_join(LFC_vs0[["2dpi-0dpi"]], LFC_vs0[["4dpi-0dpi"]], by = c("target_id", "ext_gene")) %>%
  full_join(., LFC_vs0[["7dpi-0dpi"]], by=c("target_id", "ext_gene")) %>%
  full_join(., LFC_vs0[["12dpi-0dpi"]], by=c("target_id", "ext_gene"))
write.csv(LFC_vs0, file = paste0("4_vs0_model/results_LFC.csv"))



#################################################
## ADDED June 5, 2018
#################################################

load("4_vs0_model.RData")
load("2_splines_model.RData")

res <- vector("list", 4)
names(res) <- c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")
for(i in c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")) {
  res[[i]] <- sleuth_results(so_factor_vs0, i, test_type = 'wt')
}

library(venn)
x <- list(factor = results_table_factor[which(results_table_factor$qval < 0.05),1],
          vs0 = unique(c(res[[1]][which(res[[1]]$qval < 0.05),1],
                         res[[2]][which(res[[2]]$qval < 0.05),1],
                         res[[3]][which(res[[3]]$qval < 0.05),1],
                         res[[4]][which(res[[4]]$qval < 0.05),1])))
png("factor_vs_vs0_Venn.png")
venn(x, ilabels=TRUE, zcolor="style")
dev.off()


# !!! OBSOLETE !!!
# ## Gene level analysis: load parameters
# so_gene <- my_sleuth_prep(sample_to_covariates = expdesign,
#                      full_model= full_design,
#                      num_cores = 1,
#                      target_mapping = t2g,
#                      targetIds_to_transcripts = targetIds_to_transcripts,
#                      extra_bootstrap_summary = TRUE,
#                      read_bootstrap_tpm = TRUE,
#                      max_bootstrap = 500,
#                      aggregation_column = 'ext_gene')
# ## Fit full model
# so_gene <- sleuth_fit(so_gene, fit_name="full")
# ## Wald test for each comparison
# so_gene <- sleuth_wt(so_gene, "4dpi-0dpi", which_model="full")
# so_gene <- sleuth_wt(so_gene, "2dpi-0dpi", which_model="full")
# so_gene <- sleuth_wt(so_gene, "7dpi-0dpi", which_model="full")
# so_gene <- sleuth_wt(so_gene, "12dpi-0dpi", which_model="full")
# ## Likelihood ratio test for any differences from the 0 time pt
# so_gene <- sleuth_fit(so_gene, ~1, 'reduced')
# so_gene <- sleuth_lrt(so_gene, 'reduced', 'full')
# ## Live results
# sleuth_live(so_gene)
# sleuth_deploy(so_gene, 
#               base_dir = "/Users/raua/Desktop/Optic-Regen_results/RNAseq/LIVE_initial_model_aggregation")
