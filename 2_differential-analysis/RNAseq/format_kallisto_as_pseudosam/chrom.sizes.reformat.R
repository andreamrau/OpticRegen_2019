library(tidyverse)
chrom.sizes <- read.table("danRer10.chrom.sizes", stringsAsFactors = FALSE)
colnames(chrom.sizes) <- c("chr", "size")
# alias <- read.table("../../ATACseq/danRer10_alias.tab")
# colnames(alias) <- c("chr", "chr_old")
# chrom.sizes.reformat <- inner_join(chrom.sizes, alias, by="chr_old") %>%
#   select(chr, size)

chrom.sizes.reformat <- chrom.sizes
chrom.sizes.reformat$chr_new <-  substr(chrom.sizes$chr, 4, 20)
chrom.sizes.reformat$chr_new <- ifelse(substr(chrom.sizes.reformat$chr_new, 1, 3) == "Un_",
                                       substr(chrom.sizes.reformat$chr_new, 4, 20),
                                       chrom.sizes.reformat$chr_new)

chrom.sizes.reformat$chr_new <-
  unlist(lapply(strsplit(chrom.sizes.reformat$chr_new, split="v"), 
                paste, collapse = "."))
chrom.sizes.reformat <- chrom.sizes.reformat %>%
  select(chr_new, size)
  
write.table(chrom.sizes.reformat, "danRer10.chrom.sizes.reformat", col.names=FALSE,
            row.names=FALSE, quote=FALSE, sep="\t")
