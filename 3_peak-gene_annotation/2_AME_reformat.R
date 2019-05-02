library(readxl)
library(tidyverse)
library(writexl)
library(stringr)
library(data.table) 

#------------------------------------------------------------------------------------
# motif ID = unique identifier
# name = the name I've assigned
# cluster order = assigned 1-6 based on which cluster they appear in the DETF heatmap
# unique TFID = I've assigned each of the DETFs that we have motifs for a unique TFID
#------------------------------------------------------------------------------------
motifs_cluster_order <- read_excel("AME_reformat/2019_0110_TF_Motif_cluster_family_assignments_update.xlsx",
                                   sheet = "motif-cluster-symbol (2)") %>%
  select(contains("__1")) %>%
  select(-"X__1") %>%
  mutate_all(funs(as.character)) %>%
  mutate_all(funs(toupper))
colnames(motifs_cluster_order) <- unlist(strsplit(colnames(motifs_cluster_order), split = "__1", fixed=TRUE))
colnames(motifs_cluster_order) <- sub("cluster", "clust", colnames(motifs_cluster_order))

res <- vector("list", 14)
names(res) <- paste0(rep(c("1clust5", "3clust2", "7clust1", "2clust4", "5clust7", "6clust3", "4clust6"), 
                         each = 2), rep(c("_prox", "_dist"), times = 7))
                     
for(clus in c("1clust5", "3clust2", "7clust1", "2clust4", "5clust7", "6clust3", "4clust6")) {
  for(type in c("prox", "dist")) {
    ## Read in files
    if(!clus %in% c("7clust1", "5clust7", "6clust3", "4clust6")) {
      dat_orig <- read_excel(paste0("AME_reformat/2019_0108_AME_", clus, "_JASPAR-LambertDETF-revised.xlsx"), 
                             sheet=paste0("DETF_", clus, type), skip=13)
      seq_orig <- read_excel(paste0("AME_reformat/2019_0108_AME_", clus, "_JASPAR-LambertDETF-revised.xlsx"), 
                             sheet=paste0("sequences_", clus, type))
    } else{
      dat_orig <- read_excel(paste0("AME_reformat/2019_0108_AME_", clus, "_JASPAR-LambertDETF-revised.xlsx"), 
                             sheet=paste0("DETF_", type), skip=13)
      seq_orig <- read_excel(paste0("AME_reformat/2019_0108_AME_", clus, "_JASPAR-LambertDETF-revised.xlsx"), 
                             sheet=paste0("sequences_", type))
    }

    ## Follow Ava's procedure ------------------------------------------------------------
    dat <- dat_orig %>% 
      filter(!str_detect(rank, "#")) %>%
      ## Remove unnecessary columns
      select(-rank, -motif_DB, -consensus, -`p-value`, -`adj_p-value`, -`E-value`, -tests,
             -FASTA_max, -pos, -neg, -PWM_min) %>%
      ## Calculate TP/FP ratio
      mutate(ratio = `%TP` / `%FP`) %>%
      separate_rows(motif_alt_ID, sep = ",") %>%
      separate_rows(motif_alt_ID, sep = "::") %>%
      mutate(motif_alt_ID = toupper(motif_alt_ID))
    
    dat$motif_alt_ID <- unlist( lapply(strsplit(dat$motif_alt_ID, split = "(", fixed=TRUE), function(x) x[1]))
    dat <- dat %>% filter(motif_alt_ID %in% unlist(na.omit(motifs_cluster_order[,clus])))
      
      # ## Merge JASPAR IDs with DETF cluster assignments
      # left_join(., motifs_cluster_order, by = "motif_ID") %>%
      # rename(DETF_cluster = `cluster order`) %>%
      # ## Split multiple DETF clusters into separate rows
      # separate_rows(DETF_cluster, sep=",") %>%
      # ## Sort by DETF cluster
      # arrange(DETF_cluster)
      ## Choose only the relevant TFs based on cluster
 #   if(clus == "1clust5") {
 #      dat <- filter(dat, DETF_cluster == 1)
 #    } else if(clus == "2clust4") {
 #      dat <- filter(dat, DETF_cluster %in% c(1,2)) 
 #    } else if(clus == "3clust2") {
 # #     dat <- filter(dat, DETF_cluster %in% c(2,3))
 #      dat <- filter(dat, DETF_cluster %in% c(3))
 #    } else if(clus == "4clust6") {
 #      dat <- filter(dat, DETF_cluster == 4)
 #    } else if(clus == "5clust7") {
 #      dat <- filter(dat, DETF_cluster == 5)
 #    } else {
 #      dat <- filter(dat, DETF_cluster == 6)
 #    }
 #            
    TF_choice <- unlist(dat$motif_alt_ID)                                 
    
    ## AR: the length of unique seq_IDs is the number to divide 
    seq <- seq_orig %>% 
      filter(!str_detect(motif_DB, "#")) %>%
      ## Remove fp sequences
      filter(class != "fp") %>%
      ## Remove unnecessary columns
      select(-motif_DB, -FASTA_score, -PWM_score, -class) %>%
      left_join(., dat, by = "motif_ID") %>%
      select(motif_ID, seq_ID, Name=motif_alt_ID)  %>%
      separate_rows(Name, sep=",") %>%
      ## Filter to include only chosen TFs based on cluster 
      filter(Name %in% TF_choice) %>%
      arrange(Name) %>%
      ## Remove duplicates within each TF
      select(-motif_ID) %>%
      group_by(Name) %>%
      unique() %>%
      ungroup()
    
    seq_list <- vector("list", length = length(unique(seq$Name)))
    names(seq_list) <- unique(seq$Name)
    for(j in unique(seq$Name)) {
      seq_list[[j]] <- seq %>% filter(Name == j) %>% select(seq_ID) %>% unlist()
    }

    tmp <- t(rbindlist(lapply(seq_list, function(x) data.table(t(x))), fill=TRUE)) %>%
      as.data.frame(row.names, stringsAsFactors = FALSE)
    rownames(tmp) <- NULL
    colnames(tmp) <- unique(seq$Name)
    
    res[[paste0(clus, "_", type)]] <- tmp
  }
}


write_xlsx(res,
            paste0("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx"))

