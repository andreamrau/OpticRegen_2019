library(readxl)
library(tidyverse)
library(ggraph)
library(widyr)
library(igraph)
library(tidygraph)
library(viridis)
library(UpSetR)
library(RColorBrewer)
library(forcats)

TF_families <-  read_excel("AME_reformat/2019_0110_TF_Motif_cluster_family_assignments_update.xlsx",
                           sheet = "symbol-family-cluster") %>%
  select(TF=`HGNC symbol`, DBD)
plot_data <- read_excel("AME_reformat/2019_0110_TF_Motif_cluster_family_assignments_update.xlsx",
                        sheet = "motif-cluster-symbol (2)") %>% 
  select(13:19) 
colnames(plot_data) <- sub("cluster", "clust", colnames(plot_data))
colnames(plot_data) <- strsplit(colnames(plot_data), split = "..", fixed=TRUE) %>%
  lapply(., function(x) x[1]) %>% unlist()

## Subset TFs
select_TFs <- read_excel("AME_reformat/2019_0117_DEtranscript_surroundingpeak_summary.xlsx",
                         sheet = "summary") %>%
  select(14:15)
  # select("X__9", "X__10")
select_TFs <- select_TFs[-which(rowSums(is.na(select_TFs)) == 2),]
tmp <- grep("10% cutoff", select_TFs$..14)
sel_TFs <- vector("list", 7)
names(sel_TFs) <- c("1clust5",  "2clust4", 
                    "3clust2", "4clust6",
                    "5clust7", "6clust3",
                    "7clust1")
tmp <- c(tmp, nrow(select_TFs)+1)
for(i in 1:(length(tmp)-1)) {
  sel_TFs[[i]] <- select_TFs[(tmp[i]+1):(tmp[i+1]-1),]
}

for(i in c("1clust5_prox", "1clust5_dist", "2clust4_prox", "2clust4_dist", 
           "3clust2_prox", "3clust2_dist", "4clust6_prox", "4clust6_dist",
           "5clust7_prox", "5clust7_dist", "6clust3_prox", "6clust3_dist",
           "7clust1_prox", "7clust1_dist")) {
  
  dat0 <- read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
                     sheet=i)
  dat <- matrix(0, nrow=length(unique(unlist(dat0))), ncol = ncol(dat0))
  rownames(dat) <- unique(unlist(dat0))
  colnames(dat) <- colnames(dat0)
  for(j in 1:nrow(dat)) {
    for(k in 1:ncol(dat)) {
      if(rownames(dat)[j] %in% unlist(dat0[,k])) dat[j,k] <- 1
    }
  }
  dat <- as.data.frame(dat)
  
  ## Remove unneeded TFs
  to_include <- sel_TFs[[unlist(strsplit(i, split = "_", fixed=TRUE))[1]]]$..14
  remove_col <- which(!colnames(dat) %in% to_include)
  if(length(remove_col)) {
    cat(i, "remove following TFs:", colnames(dat)[remove_col], "\n")
    dat <- dat[,-remove_col]    
  }
  
  dat$ID <- 1:nrow(dat)
  dat_counts <- dat %>% gather(key = TF, value = count, -ID) %>%
    filter(count == 1)
  new_dat <- dat_counts %>% pairwise_count(TF, ID, sort=TRUE, upper=FALSE)
  vertex_df <- data.frame(TF=unique(dat_counts$TF), TF_size = colSums(select(dat, -ID)),
                          row.names=NULL, check.names = FALSE)
  vertex_df <- left_join(vertex_df, TF_families, by = "TF")
  
  if(i == "1clust5_prox") {
    df_dat <- new_dat
    df_dat$choice <- i
    df_vertex <- vertex_df
    df_vertex$choice <- i
  } else{
    new_dat$choice <- i
    vertex_df$choice <- i
    df_dat <- bind_rows(df_dat, new_dat)
    df_vertex <- bind_rows(df_vertex, vertex_df)
  }
}

col <- c(brewer.pal(11, "Set3")[-c(2, 9)],  ## Remove yellow and grey
         colorRampPalette(c("grey10", "grey90"))(7))
names(col) <- unique(df_vertex$DBD)[c(3,2,1,9,7,10,8,11,4,5,6,16,13,14,15,12)]

####################################################################################
## Graph plots, separate for proximal and distal
####################################################################################

node_sizes <- matrix(NA, nrow=0, ncol = 4)
colnames(node_sizes) <- c("TF", "TF_size", "DBD", "choice")
node_sizes <- as.data.frame(node_sizes)
for(i in c("1clust5_prox", "1clust5_dist", "2clust4_prox", "2clust4_dist", 
           "3clust2_prox", "3clust2_dist", "4clust6_prox", "4clust6_dist",
           "5clust7_prox", "5clust7_dist", "6clust3_prox", "6clust3_dist",
           "7clust1_prox", "7clust1_dist")) {
  set.seed(12345)
  tmp_df_dat <- filter(df_dat, choice == i)
  tmp_df_vertex <- filter(df_vertex, choice == i)
  ## Make it into percentages using the total number of edges
  tmp_df_dat$n <- tmp_df_dat$n / sum(tmp_df_dat$n) * 100
  ## Make it into percentages using the total number of peaks (see AME_reformat.R file!)
  clus <- strsplit(i, split = "_", fixed=TRUE) %>% lapply(., function(x) x[1]) %>% unlist()
  type <- strsplit(i, split = "_", fixed=TRUE) %>% lapply(., function(x) x[2]) %>% unlist() 
  if(!clus %in% c("7clust1", "5clust7", "6clust3", "4clust6")) {
    seq_orig <- read_excel(paste0("AME_reformat/2019_0108_AME_", clus, "_JASPAR-LambertDETF-revised.xlsx"), 
                           sheet=paste0("sequences_", clus, type))
  } else{
    seq_orig <- read_excel(paste0("AME_reformat/2019_0108_AME_", clus, "_JASPAR-LambertDETF-revised.xlsx"), 
                           sheet=paste0("sequences_", type))
  }
  total_peaks <- seq_orig %>% 
    dplyr::select(seq_ID) %>% 
    unlist() %>% 
    strsplit(., split = "_shuf") %>% 
    lapply(., function(x) x[1]) %>% 
    unlist() %>% 
    unique() %>% 
    na.omit() %>% 
    length()
#  print(c(i, total_peaks))
  tmp_df_vertex$TF_size <- tmp_df_vertex$TF_size / total_peaks * 100
  print(c(i, max(tmp_df_vertex$TF_size), max(tmp_df_dat$n)))
  g <- tmp_df_dat %>% 
    graph_from_data_frame(directed=FALSE, 
                          vertices = tmp_df_vertex) %>%
    ggraph() +
    scale_color_manual(drop = FALSE, values = col) +
    geom_edge_link(aes(edge_width = n, edge_color = n), edge_alpha = 0.65) +
    scale_edge_color_continuous(low = "grey95", high = "black", guide = "none") +
    geom_node_point(aes(size = TF_size, color=DBD)) +
    theme_void(base_size = 15) +
    geom_node_text(aes(label = name), size = 10, repel = TRUE) + 
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.margin = margin(18,18,18,18)) +
    scale_edge_width(limits = c(0, 50), guide = "none",
                     breaks = c(1, 10, 25, 50)) +
    scale_size(range = c(0, 25),
               guide = 'none', breaks = c(20,40,60,80,100)) +
    guides(color=guide_legend(title="DBD", override.aes = list(size=5)), 
           size = guide_legend(title="% peaks\nenriched\nfor motif"),
           edge_width= guide_legend(title="% shared\nenriched\npeaks")) 
  ggsave(g, filename = paste0("AME_reformat/TF_cooccurrence_plots/graphs/", i, "_graph.pdf"), 
         width=unit(8.9, "cm"), height=unit(11.5, "cm"))
  node_sizes <- rbind(node_sizes, tmp_df_vertex)
}

library(writexl)
write_xlsx(node_sizes, "AME_reformat/TF_cooccurrence_plots/graphs/node_sizes.xlsx")

## Remake to change colors
col2 <- brewer.pal(8, "Set3")[-c(2, 7)]
names(col2) <- unique(filter(df_vertex, 
                             choice %in% c("1clust5_prox", "1clust5_dist")) %>%
                               select(DBD) %>% unlist)

for(i in c("1clust5_prox", "1clust5_dist")) {
  set.seed(12345)
  tmp_df_dat <- filter(df_dat, choice == i)
  tmp_df_vertex <- filter(df_vertex, choice == i)
  ## Make it into percentages using the total number of edges
  tmp_df_dat$n <- tmp_df_dat$n / sum(tmp_df_dat$n) * 100
  ## Make it into percentages using the total number of peaks (see AME_reformat.R file!)
  clus <- strsplit(i, split = "_", fixed=TRUE) %>% lapply(., function(x) x[1]) %>% unlist()
  type <- strsplit(i, split = "_", fixed=TRUE) %>% lapply(., function(x) x[2]) %>% unlist() 
  if(!clus %in% c("7clust1", "5clust7", "6clust3", "4clust6")) {
    seq_orig <- read_excel(paste0("AME_reformat/2019_0108_AME_", clus, "_JASPAR-LambertDETF-revised.xlsx"), 
                           sheet=paste0("sequences_", clus, type))
  } else{
    seq_orig <- read_excel(paste0("AME_reformat/2019_0108_AME_", clus, "_JASPAR-LambertDETF-revised.xlsx"), 
                           sheet=paste0("sequences_", type))
  }
  total_peaks <- seq_orig %>% 
    dplyr::select(seq_ID) %>% 
    unlist() %>% 
    strsplit(., split = "_shuf") %>% 
    lapply(., function(x) x[1]) %>% 
    unlist() %>% 
    unique() %>% 
    na.omit() %>% 
    length()
  #  print(c(i, total_peaks))
  tmp_df_vertex$TF_size <- tmp_df_vertex$TF_size / total_peaks * 100
  print(c(i, max(tmp_df_vertex$TF_size), max(tmp_df_dat$n)))
  g <- tmp_df_dat %>% 
    graph_from_data_frame(directed=FALSE, 
                          vertices = tmp_df_vertex) %>%
    ggraph() +
    scale_color_manual(drop = FALSE, values = col2) +
    geom_edge_link(aes(edge_width = n, edge_color = n), edge_alpha = 0.65) +
    scale_edge_color_continuous(low = "grey95", high = "black", guide = "none") +
    geom_node_point(aes(size = TF_size, color=DBD)) +
    theme_void(base_size = 15) +
    geom_node_text(aes(label = name), size = 10, repel = TRUE) + 
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.margin = margin(18,18,18,18)) +
    scale_edge_width(limits = c(0, 50), guide = "none",
                     breaks = c(1, 10, 25, 50)) +
    scale_size(range = c(0, 25),
               guide = 'none', breaks = c(20,40,60,80,100)) +
    guides(color=guide_legend(title="DBD", override.aes = list(size=5)), 
           size = guide_legend(title="% peaks\nenriched\nfor motif"),
           edge_width= guide_legend(title="% shared\nenriched\npeaks")) 
  ggsave(g, filename = paste0("AME_reformat/TF_cooccurrence_plots/graphs/", i, "_graph_newcolor.pdf"), 
         width=unit(8.9, "cm"), height=unit(11.5, "cm"))
}


####################################################################################
## Bar graphs of number of peaks enriched for each TF
####################################################################################

for(i in c("1clust5",  "2clust4", 
           "3clust2",  "4clust6", 
           "5clust7",  "6clust3",
           "7clust1")) {
  
  df <- filter(df_vertex, choice %in% c(paste0(i, "_prox"), paste0(i, "_dist"))) %>%
    separate(choice, into = c("cluster", "type"), sep = "_")
  df$TF <- factor(df$TF) %>%
    fct_expand(LETTERS[1:(23 - length(unique(df$TF)))])
#  TFnames <- unique(df$TF)
#  df$TF <- factor(df$TF, levels = c(unique(df$TF), LETTERS[1:(23 - length(unique(df$TF)))]))
  g <- ggplot(df) +
    geom_bar(aes(x=reorder(TF, -(df %>% group_by(TF) %>% mutate(n = sum(TF_size)) %>% ungroup() %>% select(n) %>% unlist())), 
                 y=TF_size, fill=type), stat="identity", width = 0.9) +
    scale_fill_manual(values = c("gray23", "gray88"), name = "", 
                      labels = c("distal", "proximal"), drop=FALSE) +
    scale_x_discrete(drop=FALSE,
                     labels = c("A" = "", "B" = "", "C" = "", "D" = "", "E" = "", "F" = "",
                                "G" = "", "H" = "", "I" = "", "J" = "", "K" = "", "L" = "",
                                "M" = "", "N" = "", "O" = "", "P" = "", "Q" = "", "R" = "",
                                "S" = "", "T" = "", "U" = "", "V" = "", "W" = "", "X" = "")) +
#                     labels = c(TFnames, rep("", nlevels(df$TF) - length(TFnames)))) +
    xlab("") +
    ylab("Number of peaks with enriched motifs") +
    ylim(c(0,4250)) +
    theme_classic() +
    theme(legend.position = c(0.8,0.8), text = element_text(size = 12),
          axis.text.x = element_text(angle =90, size = 14, vjust=0.5, hjust = 1),
          axis.text.y = element_text(size = 12))
  ggsave(g, filename = paste0("AME_reformat/TF_cooccurrence_plots/barplots/", i, "_barplot.pdf"),
         width=unit(8.9, "cm"), height=unit(5, "cm"))
}


####################################################################################
## UpSetR plots
####################################################################################

# df_col <- df_vertex %>% group_by(DBD) %>% count() %>% arrange(desc(n))
# df_col$color <- c(brewer.pal(12, "Set3"),
#                   colorRampPalette(c("grey0", "grey90"))(11))
# col <- df_col$color
# names(col) <- df_col$DBD
col <- c(col, "grey95")
names(col)[length(col)] <- "CSD"

## UpsetR plots combining proximal and distal
cutoff <- 10
Myfunc <- function(row, type) {
  data <- (row["type"] %in% type)
}

for(i in c("1clust5",  "2clust4", 
           "3clust2", "4clust6",
           "5clust7", "6clust3",
           "7clust1")) {
  
  ## Proximal
  dat0a <- read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
                     sheet=paste0(i, "_prox")) %>% as.list() %>% lapply(na.omit)
  ## Remove unneeded TFs
  to_include <- sel_TFs[[unlist(strsplit(i, split = "_", fixed=TRUE))[1]]]$X__9
  dat0a <- dat0a[which(names(dat0a) %in% to_include)]
  mda <- data.frame(TF = colnames(read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
                                            sheet=paste0(i, "_prox") ))) %>%
    left_join(., TF_families, by = "TF") %>%
    select(sets = TF, DBD)
  ## Remove intersections of size < cutoff
  d0 <- fromList(dat0a)
  d0 <- d0 %>%
    unite_("ID", colnames(d0), remove=FALSE) %>%
    group_by(ID) %>% count %>%
    filter(n >= cutoff) %>%
    select(ID) %>% unlist()
  d0b <- fromList(dat0a)
  d0b <- d0b %>% 
    unite_("ID", colnames(d0b), remove=FALSE) %>%
    filter(ID %in% d0) %>%
    select(-ID)
  d0prox <- d0b
  if(nrow(d0prox) > 0) {
    d0prox$type <- "proximal"
  }

  ## Distal
  dat0b <- read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
                      sheet=paste0(i, "_dist")) %>% as.list() %>% lapply(na.omit)
  ## Remove unneeded TFs
  to_include <- sel_TFs[[unlist(strsplit(i, split = "_", fixed=TRUE))[1]]]$X__9
  dat0b <- dat0b[which(names(dat0b) %in% to_include)]
  mdb <- data.frame(TF = colnames(read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
                                             sheet=paste0(i, "_dist") ))) %>%
    left_join(., TF_families, by = "TF") %>%
    select(sets = TF, DBD)
  ## Remove intersections of size < cutoff
  d0 <- fromList(dat0b)
  d0 <- d0 %>%
    unite_("ID", colnames(d0), remove=FALSE) %>%
    group_by(ID) %>% count %>%
    filter(n >= cutoff) %>%
    select(ID) %>% unlist()
  d0b <- fromList(dat0b)
  d0b <- d0b %>% 
    unite_("ID", colnames(d0b), remove=FALSE) %>%
    filter(ID %in% d0) %>%
    select(-ID)
  d0dist <- d0b
  d0dist$type <- "distal"
  
  md <- bind_rows(mda, mdb) %>% unique()
  
  if(nrow(d0prox) > 0) {
    d0 <- bind_rows(d0prox, d0dist)
    d0[is.na(d0)] <- 0
    d0 <- d0 %>% select(type, everything())
    queries <- list( list(query = Myfunc, params = list(c("distal")), color = "gray23", active = T,
                         query.name = "Distal"),
                    list(query = Myfunc, params = list(c("proximal")), color = "grey75", active = T,
                         query.name = "Proximal"))
    if(min(colSums(d0[,-1])) == 0) {
      d0 <- d0[,-c(which(colSums(d0[,-1]) == 0)+1)]
    }
  } else {
    d0 <- d0dist
    d0[is.na(d0)] <- 0
    d0 <- d0 %>% select(type, everything())
    queries <- list( list(query = Myfunc, params = list(c("distal")), color = "gray23", active = T,
                          query.name = "Distal"))
    if(min(colSums(d0[,-1])) == 0) {
      d0 <- d0[,-c(which(colSums(d0[,-1]) == 0)+1)]
    }
  }
  sets <- (arrange(md, DBD) %>% select(sets) %>% unlist)
  sets <- sets[sets %in% colnames(d0)]

  pdf( paste0("AME_reformat/TF_cooccurrence_plots/upsetr_new/", i, "_upsetr_combined.pdf"), 
       onefile = FALSE, width = 10, height = 7)
  upset(d0, nsets = ncol(d0), 
        order.by=c("freq"), 
        nintersects = NA,
        set.metadata = list(data = md, plots = list(
          list(type = "matrix_rows", 
               column = "DBD", alpha = 0.5,
               colors = col[which(names(col) %in% md$DBD)],
               text.scale = 0.5))),
        query.legend = "bottom",
        queries = queries,
        keep.order = TRUE,
        decreasing = c(TRUE, TRUE),
        sets = sets)
#,
#        text.scale=0.4,
#        point.size=0.35,
#        line.size = 0.2)
  dev.off()
  
}


  


####################################################################################
## Heatmaps of proportion of solo, double, etc enrichment in each cluster
####################################################################################

enrich_counts <- vector("list", 14)
names(enrich_counts) <- c("1clust5_prox", "1clust5_dist", "2clust4_prox", "2clust4_dist", 
                          "3clust2_prox", "3clust2_dist", "4clust6_prox", "4clust6_dist",
                          "5clust7_prox", "5clust7_dist", "6clust3_prox", "6clust3_dist",
                          "7clust1_prox", "7clust1_dist")

for(i in c("1clust5_prox", "1clust5_dist", "2clust4_prox", "2clust4_dist", 
           "3clust2_prox", "3clust2_dist", "4clust6_prox", "4clust6_dist",
           "5clust7_prox", "5clust7_dist", "6clust3_prox", "6clust3_dist",
           "7clust1_prox", "7clust1_dist")) {
  cat("***", i, "\n")
  ## Read in data
  dat0a <- read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
                      sheet=paste0(i)) 
  ## Remove unneeded TFs
  to_include <- sel_TFs[[unlist(strsplit(i, split = "_", fixed=TRUE))[1]]]$..14
  dat0a <- dat0a[which(colnames(dat0a) %in% to_include)]
  mda <- data.frame(TF = colnames(read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
                                             sheet=paste0(i) ))) %>%
    left_join(., TF_families, by = "TF") %>%
    select(sets = TF, DBD)
  
  tmp <- gather(dat0a, key="TF", value="peak") %>%
    group_by(peak) %>%
    mutate(TF = paste0(TF, collapse=",")) %>%
    ungroup() %>%
    unique() %>%
    mutate(n = unlist(lapply(strsplit(TF, split = ","), length))) %>%
    filter(!is.na(peak))
  
  mat <- matrix(0, nrow=ncol(dat0a), ncol=ncol(dat0a))
  rownames(mat) <- colnames(dat0a)
  colnames(mat) <- paste0("n", 1:ncol(mat))
  for(j in 1:nrow(mat)) {
    cat(j, "\n")
    tmp_filter <- tmp %>% slice(grep(rownames(mat)[j], tmp$TF))
    tab <- table(tmp_filter$n)
    mat[j, as.numeric(names(tab))] <- tab
  }
  
  # mat <- matrix(0, nrow=ncol(dat0a), ncol=ncol(dat0a))
  # rownames(mat) <- colnames(dat0a)
  # colnames(mat) <- paste0("n", 1:ncol(mat))
  # for(j in 1:nrow(mat)) {
  #   cat(j, "\n")
  #   tmp <- unlist(na.omit(dat0a[,j]))
  #   for(k in 1:length(tmp)) {
  #     mat[j,length(grep(tmp[k], dat0a))] <- mat[j,length(grep(tmp[k], dat0a))] + 1
  #   }
  # }
  enrich_counts[[i]] <- mat
}

for(i in 1:14) {
  tmp1 <- enrich_counts[[i]]
  colnames(tmp1) <- 1:ncol(tmp1)
  # if(i %in% c(1,3,5,7,9,11,13)) {
  #   name = "# of\nproximal\npeaks"
  # } else {
  #   name = "# of\ndistal\npeaks"
  #   
  # }
  # pdf(paste0("AME_reformat/TF_cooccurrence_plots/heatmaps/", names(enrich_counts)[i], "_count_heatmap.pdf"),
  #     width = 7, height = 5)
  # h1 <- Heatmap(tmp1, cluster_columns = FALSE, name = name,
  #         column_title = "Number of shared TF motifs enriched in peaks",
  #         column_title_side = "bottom", col=viridis(100),
  #         column_names_gp = gpar(fontsize = 18),
  #         row_names_gp = gpar(fontsize = 18))
  # print(h1)
  # dev.off()
  
  if(i %in% c(1,3,5,7,9,11,13)) {
    name = "% of\nproximal\npeaks"
  } else {
    name = "% of\ndistal\npeaks"
    
  }
  pdf(paste0("AME_reformat/TF_cooccurrence_plots/heatmaps/", names(enrich_counts)[i], "_prop_heatmap.pdf"),
      width = 7, height = 5)
  h1 <- Heatmap(tmp1 / rowSums(tmp1), cluster_columns = FALSE, name = name,
                column_title = "Number of shared TF motifs enriched in peaks",
                column_title_side = "bottom", col=viridis(100),
                column_names_gp = gpar(fontsize = 18),
                row_names_gp = gpar(fontsize = 18))
  print(h1)
  dev.off()
}




## Separate UpsetR plots for proximal and distal
# for(i in c("1clust5_prox", "1clust5_dist", "2clust4_prox", "2clust4_dist", 
#            "3clust2_prox", "3clust2_dist", "4clust6_prox", "4clust6_dist",
#            "5clust7_prox", "5clust7_dist", "6clust3_prox", "6clust3_dist",
#            "7clust1_prox", "7clust1_dist")) {
#   
#   dat0 <- read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
#                      sheet=i) %>% as.list() %>% lapply(na.omit)
#   md <- data.frame(TF = colnames(read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
#                                             sheet=i) )) %>%
#     left_join(., TF_families, by = "TF") %>%
#     select(sets = TF, DBD)
#   ## Remove intersections of size 1
#   d0 <- fromList(dat0)
#   d0 <- d0 %>%
#     unite_("ID", colnames(d0), remove=FALSE) %>%
#     group_by(ID) %>% count %>%
#     filter(n > 1) %>%
#     select(ID) %>% unlist()
#   d0b <- fromList(dat0)
#   d0b <- d0b %>% 
#     unite_("ID", colnames(d0b), remove=FALSE) %>%
#     filter(ID %in% d0) %>%
#     select(-ID)
#   pdf( paste0("AME_reformat/TF_cooccurrence_plots/upsetr_new/", i, "_upsetr.pdf"), 
#        width=15, height = 7, 
#        onefile = FALSE)
#   upset(d0b, nsets = ncol(d0b), 
#         order.by=c("freq", "degree"), 
#         nintersects = NA,
#         set.metadata = list(data = md, plots = list(
#                                                     list(type = "matrix_rows", 
#                                                          column = "DBD", alpha = 0.5,
#                                                          colors = col[which(names(col) %in% md$DBD)],
#                                                          text.scale = 0.5))),
#         keep.order = TRUE,
#         decreasing = c(TRUE, TRUE),
#         sets = arrange(md, DBD) %>% select(sets) %>% unlist %>% rev)#,
#         # text.scale=0.4,
#         # point.size=0.35,
#         # line.size = 0.2)
#   dev.off()
# }

####################################################################################
## Graph plots, separate for proximal and distal

# TF_families <-  read_excel("AME_reformat/2019_0110_TF_Motif_cluster_family_assignments_update.xlsx",
#                            sheet = "symbol-family-cluster") %>%
#   select(TF=`HGNC symbol`, DBD)
# plot_data <- read_excel("AME_reformat/2019_0110_TF_Motif_cluster_family_assignments_update.xlsx",
#                         sheet = "motif-cluster-symbol (2)") %>%
#   select(contains("cluster")) %>%
#   select(contains("__"))
# colnames(plot_data) <- sub("cluster", "clust", colnames(plot_data))
# colnames(plot_data) <- sub("__1", "", colnames(plot_data))
# 
# 
# for(i in c("1clust5_prox", "1clust5_dist", "2clust4_prox", "2clust4_dist",
#            "3clust2_prox", "3clust2_dist", "4clust6_prox", "4clust6_dist",
#            "5clust7_prox", "5clust7_dist", "6clust3_prox", "6clust3_dist",
#            "7clust1_prox", "7clust1_dist")) {
# 
#   dat0 <- read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
#                            sheet=i)
#   dat <- matrix(0, nrow=length(unique(unlist(dat0))), ncol = ncol(dat0))
#   rownames(dat) <- unique(unlist(dat0))
#   colnames(dat) <- colnames(dat0)
#   for(j in 1:nrow(dat)) {
#     for(k in 1:ncol(dat)) {
#       if(rownames(dat)[j] %in% unlist(dat0[,k])) dat[j,k] <- 1
#     }
#   }
#   dat <- as.data.frame(dat)
# 
#   ## Remove unneeded TFs
#   to_include <- unique(unlist(strsplit(unlist(na.omit(plot_data[,unlist(strsplit(i, split = "_", fixed=TRUE))[1]])),
#                                        split = ",", fixed=TRUE)))
#   remove_col <- which(!colnames(dat) %in% to_include)
#   if(length(remove_col)) {
#     cat(i, "remove following TFs:", colnames(dat)[remove_col], "\n")
#     dat <- dat[,-remove_col]
#   }
# 
# 
#   dat$ID <- 1:nrow(dat)
#   dat_counts <- dat %>% gather(key = TF, value = count, -ID) %>%
#     filter(count == 1)
#   new_dat <- dat_counts %>% pairwise_count(TF, ID, sort=TRUE, upper=FALSE)
#   vertex_df <- data.frame(TF=unique(dat_counts$TF), TF_size = colSums(select(dat, -ID)),
#                           row.names=NULL, check.names = FALSE)
#   vertex_df <- left_join(vertex_df, TF_families, by = "TF")
# 
#   if(i == "1clust5_prox") {
#     df_dat <- new_dat
#     df_dat$choice <- i
#     df_vertex <- vertex_df
#     df_vertex$choice <- i
#   } else{
#     new_dat$choice <- i
#     vertex_df$choice <- i
#     df_dat <- bind_rows(df_dat, new_dat)
#     df_vertex <- bind_rows(df_vertex, vertex_df)
#   }
# }
# 
# col <- viridis(length(unique(df_vertex$DBD)))
# names(col) <- unique(df_vertex$DBD)
# for(i in c("1clust5_prox", "1clust5_dist", "2clust4_prox", "2clust4_dist",
#            "3clust2_prox", "3clust2_dist", "4clust6_prox", "4clust6_dist",
#            "5clust7_prox", "5clust7_dist", "6clust3_prox", "6clust3_dist",
#            "7clust1_prox", "7clust1_dist")) {
# set.seed(12345)
# g <- filter(df_dat, choice == i) %>%
#   graph_from_data_frame(directed=FALSE,
#                         vertices = filter(df_vertex, choice == i)) %>%
#   ggraph() +
#   scale_color_manual(drop = FALSE, values = col) +
#   geom_edge_link(aes(edge_width = n, edge_color = n), edge_alpha = 0.25) +
#   scale_edge_color_continuous(low = "white", high = "black", guide = "none") +
#   geom_node_point(aes(size = TF_size, color=DBD)) +
#   theme_void() +
#   geom_node_text(aes(label = name), repel = TRUE) +
#   theme(legend.position = "bottom", legend.direction = "vertical",
#         legend.margin = margin(18,18,18,18)) +
#   scale_edge_width(limits = c(0, max(df_dat$n)), guide = "none") +
#   scale_size(limits = c(0, max(df_vertex$TF_size)), guide = 'none') +
# #  guides(color="none") +
#   ggtitle(i)+
#   guides(color=guide_legend(title="DBD"),
#          size = guide_legend(title="# peaks\nenriched\nfor motif"),
#          edge_width= guide_legend(title="# shared\nenriched\npeaks"))
# ggsave(g, filename = paste0("AME_reformat/", i, "_graph.pdf"), width=unit(8.9, "cm"), height=unit(7, "cm"))
# }
# 
#
#
#
#
# for(i in c("1clust5_prox", "1clust5_dist", "2clust4_prox", "2clust4_dist",
#            "3clust2_prox", "3clust2_dist", "4clust6_prox", "4clust6_dist",
#            "5clust7_prox", "5clust7_dist", "6clust3_prox", "6clust3_dist",
#            "7clust1_prox", "7clust1_dist")) {
#   
#   dat0 <- read_excel("AME_reformat/2019_0115_refined_DETF_list_motif_clusters_summary.xlsx",
#                      sheet=i)
#   dat <- matrix(0, nrow=length(unique(unlist(dat0))), ncol = ncol(dat0))
#   rownames(dat) <- unique(unlist(dat0))
#   colnames(dat) <- colnames(dat0)
#   for(j in 1:nrow(dat)) {
#     for(k in 1:ncol(dat)) {
#       if(rownames(dat)[j] %in% unlist(dat0[,k])) dat[j,k] <- 1
#     }
#   }
#   dat <- as.data.frame(dat)
#   
#   ## Remove unneeded TFs
#   to_include <- unique(unlist(strsplit(unlist(na.omit(plot_data[,unlist(strsplit(i, split = "_", fixed=TRUE))[1]])),
#                                        split = ",", fixed=TRUE)))
#   remove_col <- which(!colnames(dat) %in% to_include)
#   if(length(remove_col)) {
#     cat(i, "remove following TFs:", colnames(dat)[remove_col], "\n")
#     dat <- dat[,-remove_col]
#   }
#   
#   
#   dat$ID <- 1:nrow(dat)
#   dat_counts <- dat %>% gather(key = TF, value = count, -ID) %>%
#     filter(count == 1)
#   new_dat <- dat_counts %>% pairwise_count(TF, ID, sort=TRUE, upper=FALSE)
#   vertex_df <- data.frame(TF=unique(dat_counts$TF), TF_size = colSums(select(dat, -ID)),
#                           row.names=NULL, check.names = FALSE)
#   vertex_df <- left_join(vertex_df, TF_families, by = "TF")
#   
#   if(i == "1clust5_prox") {
#     df_dat <- new_dat
#     df_dat$choice <- i
#     df_vertex <- vertex_df
#     df_vertex$choice <- i
#   } else{
#     new_dat$choice <- i
#     vertex_df$choice <- i
#     df_dat <- bind_rows(df_dat, new_dat)
#     df_vertex <- bind_rows(df_vertex, vertex_df)
#   }
# }


####################################################################################
## Graph plots, proximal and distal together
####################################################################################
# 
# for(i in c("1clust5",  "2clust4", 
#            "3clust2",  "4clust6", 
#            "5clust7",  "6clust3",
#            "7clust1")) {
#   set.seed(12345)
#   tmp <- filter(df_dat, choice %in% c(paste0(i, "_prox"), paste0(i, "_dist"))) %>% 
#     separate(choice, sep = "_", into=c("misc", "choice")) %>%
#     group_by(item1, item2) %>%
#     summarise(n = sum(n), choice = paste0(choice, collapse=".")) %>%
#     mutate(choice = ifelse(choice == "dist", "Distal", ifelse(choice == "prox", "Proximal", "Distal")))
#   g <- tmp %>%
#     graph_from_data_frame(directed=FALSE, 
#                           vertices = df_vertex %>% 
#                             filter(choice  %in% c(paste0(i, "_prox"), paste0(i, "_dist"))) %>%
#                             select(-choice) %>% group_by(TF) %>% 
#                             summarise(TF_size = sum(TF_size)) %>% 
#                             ungroup() %>%
#                             left_join(., unique(select(df_vertex, -TF_size, -choice)), by = "TF")) %>%
#     ggraph(layout = "nicely") +
#     scale_color_manual(drop = FALSE, values = col) +
#     geom_edge_fan(aes(edge_width = n, edge_linetype = choice), edge_alpha = 0.25) +
# #    scale_edge_color_manual(values = c("grey70", "black")) +
#     geom_node_point(aes(size = TF_size, color=DBD)) +
#     theme_void(base_size = 15) +
#     geom_node_text(aes(label = name), size = 6, repel = TRUE) + 
#     theme(legend.position = "bottom", legend.direction = "vertical",
#           legend.margin = margin(18,18,18,18)) +
#     scale_edge_linetype_manual(values=c("solid", "11"), name = "", limits = c("Distal and/or proximal", "Uniquely proximal")) +
#     scale_edge_width(limits = c(0, max(df_dat$n)), guide = "none") +
#     scale_size(range = c(0, 20), limits = c(0, max(df_vertex$TF_size)),
#                guide = 'none') +
#     #  guides(color="none") + 
#     #   ggtitle(i)+
#     guides(color=guide_legend(title="DBD", override.aes = list(size=5)), 
#            size = guide_legend(title="# peaks\nenriched\nfor motif"),
#            edge_width= guide_legend(title="# shared\nenriched\npeaks")) 
#   ggsave(g, filename = paste0("AME_reformat/TF_cooccurrence_plots/graphs_combined/", i, "_graph.pdf"), 
#          width=unit(8.9, "cm"), height=unit(10, "cm"))
# }

