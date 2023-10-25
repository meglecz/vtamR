install.packages("devtools")
install.packages("roxygen2")
install.packages("seqinr")
install.packages("tidyr")
install.packages("data.table")

library("devtools")
library("roxygen2")
library("seqinr")
library("dplyr")
library("tidyr")
library("utils") # to handle zipped files
library("data.table") # reading files quickly
#library("Biostrings")
setwd("/home/meglecz")
#filename <- "vtamR_out_large_dataset/sorted/Sea18-IIICBR-Sea18_CA20_1-1.fasta"
#filename <- "/home/meglecz/vtamR_large_files/out/merged/MFZR1_S4_L001_R1_001.fasta" # 345Mb
#filename <- "/home/meglecz/vtamR_large_files/out/sorted/run1-MFZR-14Mon07-3.fasta" #15Mb, modified to single lines
filename <- "/home/meglecz/vtamR_large_files/out/sorted/run1-MFZR-TnegExt1_prerun-1.fasta" # small file some seq single line others are multi

# #7Mb, modified to single lines
system.time(df <- fread_fasta(file=filename, ids=F, dereplicate=F)) # only single line fasta
#       0.084       0.004       0.055 
system.time(df <- read.fasta(file=filename))
#       0.856       0.000       0.854 
system.time(df <- read_fasta_seq(filename=filename, dereplicate=F))
#       0.476       0.000       0.479 
system.time(df <- fread_fasta_multiline(filename=filename, dereplicate=F))
#      0.464       0.000       0.430 


# 345Mb includes multiple lines
system.time(df_rf <- read.fasta(file=filename))
#     64.260       0.468      64.720
system.time(df_rmf <- read_fasta_seq(filename=filename, dereplicate=F))
#       32.092       0.452      32.541
system.time(df_frmf <- fread_fasta_multiline(filename=filename, dereplicate=F))
#      29.756       0.144      29.181 

dim(df_rf)
dim(df_rmf)
dim(df_frmf)

fread_fasta_multiline <- function(filename=filename, dereplicate=F){
  # can deal with sequences in multiple lines
  # data.table package
  data <- fread(file=filename, header=F)

  colnames(data) <- c("asv")
  data$asv <- gsub(">.*", ">", data$asv)
  
  data <- do.call(paste, c(as.list(data$asv), sep = ""))
  data <- as.data.frame(strsplit(data, ">"))
  colnames(data) <- c("asv")
  data <- data %>%
    filter(!(asv==""))
  

  if(dereplicate){
    data <- data %>%
      group_by(asv) %>%
      summarize(read_count = length(asv))
  }
  
  return(df)
}

read_fasta_seq <- function(filename=filename, dereplicate=F){
  # can deal with sequences in multiple lines
  # only R-base
  file_connection <- file(filename, "r")
  data <- readLines(file_connection, n = -1)
  close(file_connection)
  
  data <- gsub(" ", "_", data)
  data <- gsub(">[^ ]+", ">", data, fixed=F, perl=T)
  data <- do.call(paste, c(as.list(data), sep = ""))
  data <- as.data.frame(strsplit(data, ">"))
  colnames(data) <- c("asv")
  data <- data %>%
    filter(!(asv==""))
  
  if(dereplicate){
    data <- data %>%
      group_by(asv) %>%
      summarize(read_count = length(asv))
  }
  return(data)
}

fread_fasta <- function(file=filename, ids=F, dereplicate=F){
  # cannot deal with sequences in multiple lines
  # data.table package
# quicker than read.fasta from seqinR, but can deal with fasta, where each sequence is on a single line
  df <- fread(file=filename, header=F)
  colnames(df) <- c("asv")
  df <- df %>%
    filter(!(asv==""))
  n <- nrow(df)
  
  if(ids){
    
    colnames(df) <- c("seq_ID")
    # make df with seq_IDs
    df_id <- df %>%
      filter(startsWith(seq_ID, ">"))
    n2 <- nrow(df_id)
    if(n2 * 2 != n){
      stop("Probably sequences are printed on several lines. This function can deal with fasta file where each sequences is written in a single line")
    }     
    # delete ids from df
    colnames(df) <- c("asv")
    df<- df %>%
      filter(!(startsWith(asv, ">")))
    df <- cbind(df_id, df)
    df$seq_ID <- gsub(">", "", df$seq_ID, perl=T)
    rm(df_id)
  }else{
    df <- df %>%
      filter(!(startsWith(asv, ">")))
    n2 <- nrow(df)
    if(n2 * 2 != n){
      stop("Probably sequences are prineted on several lines. This function can deal with fasta file where each sequences is written in a single line")
    }
  }
  
  
  if(dereplicate){
    df <- df %>%
      group_by(asv) %>%
      summarize(read_count = length(asv))
  }
  
  return(df)
}




