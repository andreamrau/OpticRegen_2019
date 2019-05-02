library(readxl)
file_name <- "DETFs_JASPAR_Lambert_revised.meme"

## Read in Excel file 
JASPAR <- read_excel("DETFs_JASPAR_Lambert_IDs.xlsx", sheet="DETF_JASPAR_motif_IDs")
Lambert <- read_excel("2018_0108_revised_DETFs_JASPAR_Lambert_IDs.xlsx", sheet="DETF_Lambert_uniquemotif_IDs")

memeLines <- c()
## Overall header
memeLines <- c(memeLines, "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n")
  
## Add in JASPAR motifs
for(i in 1:nrow(JASPAR)) {
  id <- unlist(JASPAR[i,1])
  motif <- readLines(paste0("JASPAR2018_CORE_vertebrates_pfm_MEME/", id, ".meme"))[-c(1:9)]
  motif <- paste0(motif[-length(motif)], collapse="\n")
  memeLines <- c(memeLines, paste0(motif, "\n"))
}

## Add in Lambert motifs
for(i in 1:nrow(Lambert)) {
  id <- unlist(Lambert[i,1])
  name <- unlist(Lambert[i,2])
  motif <- readLines(paste0("LambertPWMs/", id, ".txt"))
  w <- length(motif) - 1
  motif <- paste0(unlist(lapply(strsplit(motif[-1], split = "\t", fixed=TRUE), function(x) paste0("  ", paste0(x[-1], 
                                                                                           collapse="  "), 
                                                                                    collapse=""))),
                  collapse = "\n")
  
  
  memeLines <- c(memeLines, paste0("MOTIF ", id, " ", name, "\nletter-probability matrix: alength= ", 4,
                                   " w= ", w, " nsites= ", 20, " E=0\n", motif, 
                                   "\nURL http://humantfs.ccbr.utoronto.ca/download.php\n"))
  
}

fileConn <- file(file_name)
writeLines(memeLines, fileConn)
close(fileConn)


