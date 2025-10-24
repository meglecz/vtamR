## load
library(vtamR)
library(dplyr)
library(ggplot2)
library(rRDP)
library(rRDPData)

remove.packages("rRDP")
remove.packages("rRDPData")


setwd("/home/meglecz/vtamR/")
library("devtools")
library("roxygen2")
load_all(".")
roxygenise()
usethis::use_roxygen_md()


asv <- "/home/meglecz/vtamR/local/TAS_dryad/out_16S_pooled/2_Clustered_ASV_16S.csv"





### set up
rdp_path <- "~/rdp_classifier_2.14/dist/classifier.jar"
rdp_db <- "~/rdp_classifier_2.14/src/data/classifier/16srrna/rRNAClassifier.properties"
sep <- ","
outdir <- "~/vtamR_rdp_out"

asv <- "/home/meglecz/vtamR/local/TAS_dryad/out_16S_pooled/2_Clustered_ASV_16S.csv"
    

TaxAsssignRDP_local <- function(
      read_count,
      rdp_path = "classifier.jar",
      rdp_db,
      max_memory=8,
      sep=","
    ){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  ### make fasta file with all ASVs
  fasta <- file.path(tempdir(), "asv.fasta")
  outfile <- file.path(tempdir(), "rdp_output.tsv")
  write_fasta_df(read_count_df, fasta)
  
  
  ### Use local installation of rdp_classifier_2.14, downlodable from https://sourceforge.net/projects/rdp-classifier/
  memory = paste("-Xmx", max_memory, "g", sep="")
  args <- c(
    memory,
    "-jar", rdp_path,
    "classify",
    "-t", rdp_db,
    "-f", "fixrank",
    "-o", outfile,
    fasta
  )
  run_system2("java", args, quiet=FALSE)
  
}


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DECIPHER", "DECIPHERData"))

BiocManager::install("DECIPHER", force=TRUE)
BiocManager::install("DECIPHERData", force=TRUE)

library(DECIPHER)
browseVignettes("DECIPHER")
library(DECIPHERData)


BiocManager::version()


# select columns
rdp_assignments <- read.table(outfile, sep="\t", header=FALSE, fill = TRUE) %>%
  select(1,3,5,6,8,9,11,12,14,15,17,18,20)

colnames(rdp_assignments) <- c("asv_id", "domain", "domain_bootstrap", "phylum", "phylum_bootstrap", "class", "class_bootstrap", "order", "order_bootstrap", "family", "family_bootstrap", "genus", "genus_bootstrap")

# delete Chloroplast
rdp_assignments <- rdp_assignments %>%
  filter(class != "Chloroplast")

for(i in seq(2, ncol(rdp_assignments), by=2)){
  rdp_assignments[[i]][rdp_assignments[[i+1]] < 0.8]  <- NA # delete taxa with < 0.8 bootstrap
  rdp_assignments[[i+1]][rdp_assignments[[i+1]] < 0.8]  <- NA
  rdp_assignments[[i+1]][grepl("^Gp[0-9]+$", rdp_assignments[[i]])]  <- NA # delete taxa Gpxx
  rdp_assignments[[i]][grepl("^Gp[0-9]+$", rdp_assignments[[i]])]  <- NA
}

asv_unique <- clustered_df %>%
  select(asv_id, asv) %>%
  distinct()

# add ASVs 
rdp_assignments <- left_join(rdp_assignments, asv_unique, by=c("asv_id"))

outfile <- paste(outdir, "3_TaxAssign_TAS_",my_marker,"_rdp_clean.csv", sep="")
write.csv(rdp_assignments, outfile, row.names = FALSE)
