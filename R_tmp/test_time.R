

## load
library(vtamR)
library(dplyr)
#library(ggplot2)
#library(rRDP)
#library(rRDPData)

setwd("/home/meglecz/vtamR/")
setwd("C:/Users/emese/vtamR")
library("devtools")
library("roxygen2")
load_all(".")
roxygenise()
usethis::use_roxygen_md()

### set up
cutadapt_path <- "~/miniconda3/envs/vtam/bin/cutadapt"
vsearch_path <- "~/miniconda3/envs/vtam/bin/vsearch"
blast_path <- "~/miniconda3/envs/vtam/bin/blastn"
swarm_path <- "swarm"
pigz_path <- "pigz"
sep <- ","
stat_dir <- "/home/meglecz/vtamR/R_tmp/time"
syst <- "Bombyx"


### Win
cutadapt_path <- "C:/Users/Public/cutadapt"
vsearch_path <- "C:/Users/Public/vsearch-2.23.0-win-x86_64/bin/vsearch"
blast_path <- "C:/Users/Public/blast-2.16.0+/bin/blastn"
swarm_path <- "C:/Users/Public/swarm-3.1.5-win-x86_64/bin/swarm"
pigz_path <- "C:/Users/Public/pigz-win32/pigz"
num_threads <- 0
sep <- ","
stat_dir <- "C:/data/EPI09/out/time_stat"
syst <- "windows"


################# Demo data
fastq_dir <- system.file("extdata/demo/fastq", package = "vtamR")
fastqinfo <-  system.file("extdata/demo/fastqinfo.csv", package = "vtamR")
mock_composition <-  system.file("extdata/demo/mock_composition.csv", package = "vtamR")
asv_list <-  system.file("extdata/demo/asv_list.csv", package = "vtamR")
taxonomy <- system.file("extdata/db_test/taxonomy_reduced.tsv", package = "vtamR")
blast_db <- system.file("extdata/db_test", package = "vtamR")
blast_db <- file.path(blast_db, "COInr_reduced")
# Bombyx
outdir <- "/home/meglecz/vtamR_demo_out"
# win
outdir <- ""


################# Big Data
# Bombyx
outdir <- "/home/meglecz/vtamR_test_EPI09_COI/test_Merge_SortReads"
fastq_dir <- "/home/meglecz/vtamR_large_files/EPI09"
fastqinfo <- "/home/meglecz/vtamR_large_files/EPI09/metainfo/mock_composition_EPI09_COI.csv"
ds <- "EPI9"
# Win
outdir <- "C:/data/EPI09/out"
fastq_dir <- "C:/data/EPI09"
fastqinfo <- "C:/data/EPI09/metainfo/fastqinfo_Epi09_COI.csv"
ds <- "EPI9"



time_df <- data.frame(
  Step = character(),
  user = numeric(),
  system = numeric(),
  elapsed = numeric(),
  stringsAsFactors = FALSE)

###### Test merge 
# input compressed
# output: compressed/uncompressed
# compress_methode: R/pigz 

### Merge  output: uncompressed; compress_methode: R
t1 <- proc.time()
merged_dir_uncompress <- file.path(outdir, "Merge_uncompress")
fastainfo_df_uncompress <- Merge(fastqinfo, 
                                 fastq_dir=fastq_dir, 
                                 vsearch_path=vsearch_path, 
                                 compress_method="R",
                                 outdir=merged_dir_uncompress,
                                 compress=FALSE)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "Merge_uncompress_R",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))

# Merge  output: compressed; compress_methode: R
t1 <- proc.time()
merged_dir_compress <- file.path(outdir, "Merge_compress")
fastainfo_df_compress <- Merge(fastqinfo, 
                         fastq_dir=fastq_dir, 
                         vsearch_path=vsearch_path, 
                         compress_method="R",
                         outdir=merged_dir_compress,
                         compress=TRUE)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "Merge_compress_R",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


# Merge  output: uncompressed; compress_methode: pigz
t1 <- proc.time()
merged_dir_uncompress <- file.path(outdir, "Merge_uncompress")
fastainfo_df_uncompress <- Merge(fastqinfo, 
                                 fastq_dir=fastq_dir, 
                                 vsearch_path=vsearch_path, 
                                 compress_method="pigz",
                                 outdir=merged_dir_uncompress,
                                 compress=FALSE)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "Merge_uncompress_pigz",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))

# Merge  output: compressed; compress_methode: pigz
t1 <- proc.time()
merged_dir_compress <- file.path(outdir, "Merge_compress")
fastainfo_df_compress <- Merge(fastqinfo, 
                               fastq_dir=fastq_dir, 
                               vsearch_path=vsearch_path, 
                               compress_method="pigz",
                               outdir=merged_dir_compress,
                               compress=TRUE)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "Merge_compress_pigz",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))

file_stat <- paste("Merge_", ds, "_", syst, ".csv", sep="")
file_stat <- file.path(stat_dir, file_stat)
write.csv(time_df, file=file_stat)


###### Test RandomSample 
# input: uncompressed
# output: uncompressed
# compress_methode: R => irrelevant, since I/O are uncompressed 
# use_vsearch: TRUE/FALSE

time_df <- data.frame(
  Step = character(),
  user = numeric(),
  system = numeric(),
  elapsed = numeric(),
  stringsAsFactors = FALSE)

## Vsearch
t1 <- proc.time()
RandomSeqdir <- file.path(outdir, "RandomSeq_uncompress_vsearch_R")
fastainfo_random_uncompress <- RandomSeq(fastainfo=fastainfo_df_uncompress, 
                                       n = 40000,
                                        fasta_dir=merged_dir_uncompress,
                                        outdir=RandomSeqdir, 
                                        use_vsearch=TRUE,
                                        vsearch_path=vsearch_path,
                                        randseed=123, 
                                        compress_method="R",
                                        pigz_path=pigz_path,
                                        num_threads=0,
                                        compress=F, 
                                        quiet=TRUE)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "RandomSeq_uncompress_vsearch_R",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))



time_df <- data.frame(
  Step = character(),
  user = numeric(),
  system = numeric(),
  elapsed = numeric(),
  stringsAsFactors = FALSE)

## No Vsearch 
fastainfo_df_uncompress <- read.csv("C:/data/EPI09/out/Merge_uncompress/fastainfo.csv")
t1 <- proc.time()
RandomSeqdir <- file.path(outdir, "RandomSeq_uncompress_NOvsearch_R")
fastainfo_random_uncompress <- RandomSeq(fastainfo=fastainfo_df_uncompress, 
                                         n = 10000000,
                                         fasta_dir=merged_dir_uncompress,
                                         outdir=RandomSeqdir, 
                                         use_vsearch=FALSE,
                                         vsearch_path=vsearch_path,
                                         randseed=123, 
                                         compress_method="R",
                                         pigz_path=pigz_path,
                                         num_threads=0,
                                         compress=F, 
                                         quiet=TRUE)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "RandomSeq_uncompress_NOvsearch_R",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))

file_stat <- paste("RandomSeq_", ds, "_", syst, ".csv", sep="")
file_stat <- file.path(stat_dir, file_stat)
write.csv(time_df, file=file_stat)


###### Test SortReads 
# input: uncompressed Merge/RandomSeq
# output: compressed/uncompressed => check if difference, and will be useful for Dereplicate
# compress_method: R => irrelevant, if check_reverse == FALSE (most cases)


time_df <- data.frame(
  Step = character(),
  user = numeric(),
  system = numeric(),
  elapsed = numeric(),
  stringsAsFactors = FALSE)
# input Merge output uncompress
fastainfo_df_uncompress <- read.csv(file.path(outdir, "Merge_uncompress/fastainfo.csv"))
fasta_dir <- file.path(outdir, "Merge_uncompress")
sorted_dir <- file.path(outdir, "SortReads_Merge_out_uncompress")

t1 <- proc.time()
sampleinfo_uncompress <- SortReads(fastainfo=fastainfo_df_uncompress,
                                    fasta_dir=fasta_dir,
                                    outdir=sorted_dir, 
                                    cutadapt_path=cutadapt_path,
                                    vsearch_path=vsearch_path, 
                                    compress_method="R",
                                    check_reverse=F, 
                                    tag_to_end=T, 
                                    primer_to_end=T, 
                                    compress=F,
                                    quiet=T)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "SortReads_Merge_out_uncompress",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


# input Merge output compress
sorted_dir <- file.path(outdir, "SortReads_Merge_out_compress")
t1 <- proc.time()
sampleinfo_uncompress <- SortReads(fastainfo=fastainfo_df_uncompress,
                                   fasta_dir=fasta_dir,
                                   outdir=sorted_dir, 
                                   cutadapt_path=cutadapt_path,
                                   vsearch_path=vsearch_path, 
                                   compress_method="R",
                                   check_reverse=F, 
                                   tag_to_end=T, 
                                   primer_to_end=T, 
                                   compress=T,
                                   quiet=T)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "SortReads_Merge_out_compress",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


# input Randomseq - output uncompress
fastainfo_df_uncompress <- read.csv(file.path(outdir, "RandomSeq_uncompress_NOvsearch_R/fastainfo.csv"))
fasta_dir <- file.path(outdir, "RandomSeq_uncompress_NOvsearch_R")
sorted_dir <- file.path(outdir, "SortReads_RandomSeq_out_uncompress")
t1 <- proc.time()
sampleinfo_uncompress <- SortReads(fastainfo=fastainfo_df_uncompress,
                                   fasta_dir=fasta_dir,
                                   outdir=sorted_dir, 
                                   cutadapt_path=cutadapt_path,
                                   vsearch_path=vsearch_path, 
                                   compress_method="R",
                                   check_reverse=F, 
                                   tag_to_end=T, 
                                   primer_to_end=T, 
                                   compress=F,
                                   quiet=T)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "SortReads_RandomSeq_out_uncompress",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


# input Randomseq - output uncompress
sorted_dir <- file.path(outdir, "SortReads_RandomSeq_out_compress")
t1 <- proc.time()
sampleinfo_uncompress <- SortReads(fastainfo=fastainfo_df_uncompress,
                                   fasta_dir=fasta_dir,
                                   outdir=sorted_dir, 
                                   cutadapt_path=cutadapt_path,
                                   vsearch_path=vsearch_path, 
                                   compress_method="R",
                                   check_reverse=F, 
                                   tag_to_end=T, 
                                   primer_to_end=T, 
                                   compress=T,
                                   quiet=T)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "SortReads_RandomSeq_out_compress",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))

file_stat <- paste("Sortreads_", ds, "_", syst, ".csv", sep="")
file_stat <- file.path(stat_dir, file_stat)
write.csv(time_df, file=file_stat)


###### Test Dereplicate
# input: uncompressed/compressed RandomSeq

time_df <- data.frame(
  Step = character(),
  user = numeric(),
  system = numeric(),
  elapsed = numeric(),
  stringsAsFactors = FALSE)

sampleinfo <- read.csv(file.path(outdir, "SortReads_RandomSeq_out_uncompress/sampleinfo.csv"))
sorted_dir <- file.path(outdir, "SortReads_RandomSeq_out_uncompress")
out <- file.path(outdir, "Dereplicate/Dereplicted_RandomSeq_input_uncompress.csv")
t1 <- proc.time()
read_count_df <- Dereplicate(sampleinfo=sampleinfo, 
                             dir=sorted_dir,
                             outfile=out,
                             quiet=T)
  
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "Dereplicted_RandomSeq_input_uncompress",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


  
sampleinfo <- read.csv(file.path(outdir, "SortReads_RandomSeq_out_compress/sampleinfo.csv"))
sorted_dir <- file.path(outdir, "SortReads_RandomSeq_out_compress")
out <- file.path(outdir, "Dereplicate/Dereplicted_RandomSeq_input_compress.csv")
t1 <- proc.time()
read_count_df <- Dereplicate(sampleinfo=sampleinfo, 
                             dir=sorted_dir,
                             outfile=out,
                             quiet=T)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "Dereplicted_RandomSeq_input_compress",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


file_stat <- paste("Dereplicate_", ds, "_", syst, ".csv", sep="")
file_stat <- file.path(stat_dir, file_stat)
write.csv(time_df, file=file_stat)











### R
fastainfo_df <- RandomSeq2(fastainfo="/home/meglecz/vtamR_demo_out/merged_uncompress/fastainfo.csv", 
                       n=40000,
                       fasta_dir="/home/meglecz/vtamR_demo_out/merged_uncompress",
                       outdir="/home/meglecz/vtamR_demo_out/RandomSeq_uncompress", 
                       use_vsearch=FALSE,
                       vsearch_path=vsearch_path,
                       compress=F, 
                       quiet=TRUE)


fastainfo_df <- RandomSeq2(fastainfo="/home/meglecz/vtamR_demo_out/merged_uncompress/fastainfo.csv", 
                           n=40000,
                           fasta_dir="/home/meglecz/vtamR_demo_out/merged_uncompress",
                           outdir="/home/meglecz/vtamR_demo_out/RandomSeq_gz", 
                           use_vsearch=FALSE,
                           vsearch_path=vsearch_path,
                           compress=T, 
                           quiet=TRUE)


fastainfo_df <- RandomSeq2(fastainfo="/home/meglecz/vtamR_demo_out/merged_gz/fastainfo.csv", 
                           n=40000,
                           fasta_dir="/home/meglecz/vtamR_demo_out/merged_gz",
                           outdir="/home/meglecz/vtamR_demo_out/RandomSeq_gz_uncompress", 
                           use_vsearch=FALSE,
                           vsearch_path=vsearch_path,
                           compress=F, 
                           quiet=TRUE)


fastainfo_df <- RandomSeq2(fastainfo="/home/meglecz/vtamR_demo_out/merged_gz/fastainfo.csv", 
                           n=40000,
                           fasta_dir="/home/meglecz/vtamR_demo_out/merged_gz",
                           outdir="/home/meglecz/vtamR_demo_out/RandomSeq_gz_gz", 
                           use_vsearch=FALSE,
                           vsearch_path=vsearch_path,
                           compress=T, 
                           quiet=TRUE)

### vsearch
fastainfo_df <- RandomSeq2(fastainfo="/home/meglecz/vtamR_demo_out/merged_uncompress/fastainfo.csv", 
                           n=40000,
                           fasta_dir="/home/meglecz/vtamR_demo_out/merged_uncompress",
                           outdir="/home/meglecz/vtamR_demo_out/RandomSeq_uncompress_uncompress", 
                           use_vsearch=TRUE,
                           vsearch_path=vsearch_path,
                           compress=F, 
                           quiet=F)


fastainfo_df <- RandomSeq2(fastainfo="/home/meglecz/vtamR_demo_out/merged_uncompress/fastainfo.csv", 
                           n=40000,
                           fasta_dir="/home/meglecz/vtamR_demo_out/merged_uncompress",
                           outdir="/home/meglecz/vtamR_demo_out/RandomSeq_uncompress_gz", 
                           use_vsearch=TRUE,
                           vsearch_path=vsearch_path,
                           compress=T, 
                           quiet=TRUE)


fastainfo_df <- RandomSeq2(fastainfo="/home/meglecz/vtamR_demo_out/merged_gz/fastainfo.csv", 
                           n=40000,
                           fasta_dir="/home/meglecz/vtamR_demo_out/merged_gz",
                           outdir="/home/meglecz/vtamR_demo_out/RandomSeq_gz_uncompress", 
                           use_vsearch=TRUE,
                           vsearch_path=vsearch_path,
                           compress=F, 
                           quiet=TRUE)


fastainfo_df <- RandomSeq2(fastainfo="/home/meglecz/vtamR_demo_out/merged_gz/fastainfo.csv", 
                           n=40000,
                           fasta_dir="/home/meglecz/vtamR_demo_out/merged_gz",
                           outdir="/home/meglecz/vtamR_demo_out/RandomSeq_gz_gz", 
                           use_vsearch=TRUE,
                           vsearch_path=vsearch_path,
                           compress=T, 
                           quiet=TRUE)


#' Select random sequences
#' 
#' Random select n sequences from each input fasta file. 
#'   
#' Do not work on Windows! If using Windows, please, use RandomSeqWindows
#'  
#' @param fastainfo Data frame or csv file with a fasta column containing input 
#' fasta file names. Files can be gzip compressed.
#' @param n Positive integer: the number of randomly selected sequences from each input file.
#' @param fasta_dir Character string: directory that contains the 
#' input fasta files
#' @param outdir Character string: output directory.
#' @param vsearch_path Character string: path to vsearch executables. 
#' @param pigz Logical: Use pigz for compression/decompression. If FALSE uses `R.utils`.
#' @param pigz_path Character string: Path to `pigz` executable.
#' @param num_threads Positive integer: Number of CPUs. If 0, use all available CPUs.
#' @param randseed Positive integer: seed for random sampling. 
#' 0 (default value) means to use a pseudo-random seed. 
#' A given non-zero seed produces always the same result.
#' @param compress logical: Compress output files to gzip format.
#' @param sep Field separator character in input and output csv files.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns The input data frame with updated file names and read counts.
#' @examples
#' \dontrun{
#' fastainfo_df <- RandomSeq(fastainfo, 
#'    fasta_dir="data/fasta", 
#'    outdir="data/randomseq"
#'    )
#' }
#' @export
#' 
RandomSeq <- function(fastainfo,
                      n, 
                      fasta_dir="",
                      outdir="", 
                      vsearch_path="vsearch",
                      pigz= FALSE,
                      pigz_path= "pigz",
                      num_threads = 0,
                      randseed=0,
                      compress=F,
                      sep=",", 
                      quiet=T){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  
  # if run on non linux-like system, use RandomSeqWindows
  if(!is_linux()){
    fastainfo_df <- RandomSeqWindows(fastainfo, 
                                     fasta_dir=fasta_dir, 
                                     outdir=outdir, 
                                     n, 
                                     randseed=randseed, 
                                     compress=compress, 
                                     sep=sep)
    return(fastainfo_df)
  }
  
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  #  CheckFileinfo(file=fastainfo_df, dir=fasta_dir, file_type="fastainfo", sep=sep, quiet=TRUE)
  
  # quite fast for uncompressed and gz files
  check_dir(fasta_dir)
  check_dir(outdir)
  
  fastainfo_df$new_file <- NA
  fastainfo_df$read_count <- NA
  
  unique_fasta <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(unique_fasta)){ # go through all fasta files
    input_fasta <- unique_fasta[i]
    print(input_fasta)
    # stop if zip file
    if(endsWith(input_fasta, ".zip")){
      stop("Zip files are not supported")
    }
    # adjust output filename in function of the compression
    output_fasta <- input_fasta
    if(compress && !endsWith(output_fasta, ".gz")){
      output_fasta <- paste(output_fasta, ".gz", sep="")
    }
    if(!compress && endsWith(output_fasta, ".gz")){
      output_fasta <- sub(".gz$", "", output_fasta)
    }
    input_fasta_p <- file.path(fasta_dir, input_fasta)
    output_fasta_p <- file.path(outdir, output_fasta)
    
    seq_n <- count_seq(file=input_fasta_p)
    if(!quiet){
      msg <- paste("Number of sequences in:", seq_n)
      print(msg)
    }    
    if(n > seq_n ){ # not enough seq
      # print msg
      msg <- paste("WARNING:", input_fasta_p,"has",
                   seq_n,"sequences, which is lower than", n,
                   ". The file is copied to the",outdir,
                   "directory without subsampling", 
                   sep=" "
      )
      print(msg)
      # copy and compress/decompress input file according to the need
      if(compress && !endsWith(input_fasta_p, ".gz")){ # input uncompressed, and compress = T
        #        output <- compress_file(filename=input_fasta_p, remove_input=F) # compress input, keep original
        output <- smart_gzip(file=input_fasta_p,
                             remove = FALSE,
                             pigz=pigz,
                             pigz_path = pigz_path,
                             num_threads = num_threads,
                             quiet = quiet,
                             compress = TRUE)        
        file.rename(output, output_fasta_p) # move compressed file to the target location
      }else{
        if(!compress && endsWith(input_fasta_p, ".gz")){ # input compressed, and compress = F
          #         output <- decompress_file(filename=input_fasta_p, remove_input=F) # compress input, keep original
          output <- smart_gzip(file=input_fasta_p,
                               remove = FALSE,
                               pigz=pigz,
                               pigz_path = pigz_path,
                               num_threads = num_threads,
                               quiet = quiet,
                               compress = FALSE)
          file.rename(output, output_fasta_p) # move compressed file to the target location
        }else{ # output, input same compression
          file.copy(input_fasta_p, output_fasta_p)
        }
      }
      fastainfo_df$new_file[which(fastainfo_df$fasta==input_fasta)] <- output_fasta
      fastainfo_df$read_count[which(fastainfo_df$fasta==input_fasta)] <- seq_n
      next()
    } # not enough seq
    
    # enough seq => resample
    # do not transform large numbers to scientific forms, since it would lead to an error in vsearch
    options(scipen=100)
    output_fasta_p <- gsub(".gz", "", output_fasta_p) # vsearch makes decompressed files
    
    ##### run vsearch
    # Build argument vector
    args <- c("--fastx_subsample", input_fasta_p,
              "--fastaout", output_fasta_p,
              "--sample_size",  n,
              "--randseed", randseed
    )
    if(num_threads > 0){
      args <- append(args, c("--threads", num_threads), after = 2)
    }
    run_system2(vsearch_path, args, quiet=quiet)
    
    options(scipen=0)
    
    if(compress){ # compress the output file 
      #      output_fasta_p <- compress_file(filename=output_fasta_p, remove_input=T)
      output_fasta_p <- smart_gzip(file=output_fasta_p,
                                   remove = TRUE,
                                   pigz=pigz,
                                   pigz_path = pigz_path,
                                   num_threads = num_threads,
                                   quiet = quiet,
                                   compress = TRUE)
    }
    
    fastainfo_df$new_file[which(fastainfo_df$fasta==input_fasta)] <- output_fasta
    fastainfo_df$read_count[which(fastainfo_df$fasta==input_fasta)] <- n
  }# all files
  
  fastainfo_df$fasta <- fastainfo_df$new_file
  fastainfo_df <- fastainfo_df %>%
    select(-new_file)
  new_fastainfo <- file.path(outdir, "fastainfo.csv")
  write.table(fastainfo_df, file = new_fastainfo,  row.names = F, sep=sep)
  return(fastainfo_df)
}














### demultiplex
demultiplexed_dir <- file.path(outdir, "demultiplexed")
sampleinfo_df <- SortReads(fastainfo_df, 
                           fasta_dir=merged_dir, 
                           outdir=demultiplexed_dir, 
                           check_reverse=TRUE, 
                           cutadapt_path=cutadapt_path, 
                           vsearch_path=vsearch_path,
                           tag_to_end = T,
                           primer_to_end=T)

###############
### dereplicate
###############
updated_asv_list <- file.path(outdir, "updated_asv_list.tsv")
demultiplexed_dir <- file.path(outdir, "demultiplexed")
sampleinfo <- file.path(demultiplexed_dir, "sampleinfo.csv")
read_count_df <- Dereplicate(sampleinfo, 
                             dir=demultiplexed_dir, 
                             input_asv_list=asv_list,
                             output_asv_list = updated_asv_list)

### stat
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())

stat_df <- GetStat(read_count_df, stat_df, stage="input_sample_replicate", params=NA)

#### swarm
by_sample <- TRUE
d=1
fastidious= TRUE
quiet=TRUE
outfile <- file.path(outdir, "filter", "2_Swarm_by_sample.csv")

read_count_df <- ClusterASV(read_count_df,
                            method = "swarm",
                            swarm_d=d,
                            fastidious=fastidious,
                            by_sample=by_sample,
                            group = TRUE,
                            path=swarm_path
                            )


stat_df <- GetStat(read_count_df, stat_df, stage="Swarm", params=by_sample)

#### LFNglobalReadCount
global_read_count_cutoff = 2

read_count_df <- LFNglobalReadCount(read_count_df, 
                                    cutoff=global_read_count_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNglobalReadCount", params=NA)

#### FilterIndel
read_count_df <- FilterIndel(read_count_df)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterIndel", params=NA)

### FilterCodonStop
genetic_code = 5
read_count_df <- FilterCodonStop(read_count_df, 
                                 genetic_code=genetic_code)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterCodonStop", params=NA)


read_count_df_backup <- read_count_df

### FilterExternalContaminant
conta_file <- file.path(outdir, "tmp", "external_contamination.csv")
read_count_df <- FilterExternalContaminant(read_count_df, 
                          sample_types=sampleinfo, 
                          conta_file=conta_file)

stat_df <- GetStat(read_count_df, stat_df, stage="FilterExternalContaminant", params=NA)

### FilterChimera
abskew=2
by_sample = T
sample_prop = 0.8
read_count_df <- FilterChimera(read_count_df, 
                               vsearch_path=vsearch_path, 
                               by_sample=by_sample, 
                               sample_prop=sample_prop, 
                               abskew=abskew)

stat_df <- GetStat(read_count_df, stat_df, stage="FilterChimera", params=NA)

#### FilterRenkonen
cutoff <- 0.4
read_count_df <- FilterRenkonen(read_count_df, 
                                cutoff=cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterRenkonen", params=cutoff)


fas <- "/home/meglecz/vtamR/inst/extdata/demo/mock_ncbi.fasta"
outdir <- "/home/meglecz/vtamR/inst/extdata/demo/mock"

taxonomy <- taxonomy_COInr

MakeMockCompositionLTG(read_count_df, 
                      fas=fas,
                      taxonomy=taxonomy_COInr,
                      blast_path = blast_path,
                      sampleinfo = sampleinfo_df,
                      outdir= outdir)







make_taxid_file(file=fas, outfile=taxids)









### FilterPCRerror

mock_composition_df <- read.csv(mock_composition)

mock_composition_df[5,"asv"] <- "TTTTTTTTTTTT"
mock_composition_df[6,"asv"] <- "AAAAAAAAAa"

dir_opt = file.path(outdir, "OptimizeLFNreadCountLFNvariant") 
opt_rc_var <- OptimizeLFNreadCountLFNvariant(read_count_df, 
                                           outdir = dir_opt,
                                           known_occurrences = NULL, 
                                           mock_composition = mock_composition_df,
                                           sampleinfo = sampleinfo_df,
                                           habitat_proportion = 0.5,
                                           sep=",",
                                           min_lfn_read_count_cutoff=10, 
                                           max_lfn_read_count_cutoff=100, 
                                           increment_lfn_read_count_cutoff=5, 
                                           min_lnf_variant_cutoff=0.001, 
                                           max_lnf_variant_cutoff=0.01, 
                                           increment_lnf_variant_cutoff=0.001, 
                                           by_replicate=FALSE, 
                                           min_replicate_number=2, 
                                           quiet=T)

opt_rc_var <- OptimizeLFNreadCountLFNvariant(read_count_df, 
                                             outdir = dir_opt,
                                             known_occurrences = "/home/meglecz/vtamR_demo_out/OptimizeLFNreadCountLFNvariant/known_occurrences.csv", 
                                             sep=",",
                                             min_lfn_read_count_cutoff=10, 
                                             max_lfn_read_count_cutoff=100, 
                                             increment_lfn_read_count_cutoff=5, 
                                             min_lnf_variant_cutoff=0.001, 
                                             max_lnf_variant_cutoff=0.01, 
                                             increment_lnf_variant_cutoff=0.001, 
                                             by_replicate=FALSE, 
                                             min_replicate_number=2, 
                                             quiet=T)




optPCR <- OptimizePCRerror(read_count=read_count_df, 
                 mock_composition=mock_composition_df, 
                 vsearch_path=vsearch_path, 
                 max_mismatch=2, 
                 min_read_count=5)

optLFN_sample <- OptimizeLFNsampleReplicate(read_count_df, mock_composition=mock_composition_df)

optLFN_var <- OptimizeLFNreadCountLFNvariant(read_count_df, 
                                           known_occurrences, 
                                           sep=",",
                                           outfile="", 
                                           min_lfn_read_count_cutoff=10, 
                                           max_lfn_read_count_cutoff=100, 
                                           increment_lfn_read_count_cutoff=5, 
                                           min_lnf_variant_cutoff=0.001, 
                                           max_lnf_variant_cutoff=0.01, 
                                           increment_lnf_variant_cutoff=0.001, 
                                           by_replicate=FALSE, 
                                           min_replicate_number=2, 
                                           quiet=T
)


pcr_error_var_prop <- 0.05
max_mismatch <- 2
read_count_df <- FilterPCRerror(read_count_df, 
                                vsearch_path=vsearch_path, 
                                pcr_error_var_prop=pcr_error_var_prop, 
                                max_mismatch=max_mismatch)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterPCRerror", params=NA)

### LFNsampleReplicate
lfn_sample_replicate_cutoff <- 0.004
read_count_df <- LFNsampleReplicate(read_count_df, 
                                    cutoff=lfn_sample_replicate_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNsampleReplicate", params=NA)

### FilterMinReplicate
min_replicate_number <- 2
read_count_df <- FilterMinReplicate(read_count_df, 
                                    cutoff=min_replicate_number)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterMinReplicate", params=NA)


### LFNvariant
lnf_variant_cutoff = 0.001
read_count_df_lnf_variant <- LFNvariant(read_count_df, 
                                        cutoff=lnf_variant_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNvariant", params=NA)

### LFNreadCount
lfn_read_count_cutoff <- 10
read_count_df_lfn_read_count <- LFNreadCount(read_count_df, 
                                             cutoff=lfn_read_count_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNreadCount", params=NA)




### Combine results
read_count_df <- PoolFilters(read_count_df_lfn_read_count, 
                             read_count_df_lnf_variant)
stat_df <- GetStat(read_count_df, stat_df, stage="Combine results", params=NA)

# delete temporary data frames
rm(read_count_df_lfn_read_count)
rm(read_count_df_lnf_variant)

### FilterMinReplicate
min_replicate_number <- 2
read_count_df <- FilterMinReplicate(read_count_df, 
                                    cutoff=min_replicate_number)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterMinReplicate", params=NA)

read_count_df_backup <- read_count_df

plot_swarm <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      outfile="13_FilterMinReplicate.png")






ASVspecificCutoff_df_round <- ASVspecificCutoff(read_count_df,  mock_composition=mock_composition,
                              by_replicate=FALSE, 
                              outfile="tmp/ASVspecificCutoff_by_replicate_false.csv")






tmp2 <- LFNvariant2(read_count_df, 
                        asv_specific_cutoff = NULL,
                        cutoff=0.01,
                        by_replicate=FALSE, 
                        outfile="", 
                        sep=",", 
                        min_read_count_prop=0.7)

tmp <- LFNvariant(read_count_df, 
                    cutoff=0.01,
                    by_replicate=FALSE, 
                    outfile="", 
                    sep=",", 
                    min_read_count_prop=0.7)









### MakeKnownOccurrences performance_metrics
results <- MakeKnownOccurrences(read_count_samples_df, 
                                sortedinfo=sortedinfo_df, 
                                mock_composition=mock_composition)


updated_asv_list <- file.path(outdir, "updated_asv_list_end.tsv")
UpdateASVlist(asv_list1 = read_count_samples_df,
              asv_list2 =asv_list, 
              outfile=updated_asv_list
)

### TaxAssign
asv_tax <- TaxAssign(asv=read_count_samples_df, 
                     taxonomy=taxonomy, 
                     blast_db=blast_db, 
                     blast_path=blast_path, 
                     num_threads=num_threads)


plot <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      num_threads=num_threads,
                                      outfile="2_swarm.png")


### WriteASVtable
outfile=file.path(outdir, "Final_asvtable_with_TaxAssign.csv")
asv_table_df <- WriteASVtable(read_count_samples_df, 
                              outfile=outfile, 
                              asv_tax=asv_tax, 
                              sortedinfo=sortedinfo_df, 
                              add_empty_samples=T, 
                              add_sums_by_sample=T, 
                              add_sums_by_asv=T, 
                              add_expected_asv=T, 
                              mock_composition=mock_composition)


#####################
#####################
#####################
# make mOTU

#####################
### mOTU with swarm

stat_df <- GetStat(read_count_samples_df, stat_df, stage="PoolReplicates", params=NA)

by_sample <- FALSE
d = 7
read_count_df_swarm_motu <- Swarm(read_count_samples_df, 
                       swarm_path=swarm_path, 
                       swarm_d=d,
                       fastidious=FALSE,
                       num_threads=num_threads, 
                       by_sample=by_sample)

stat_df <- GetStat(read_count_df_swarm_motu, stat_df, stage="swarm_motu_7", params=d)

### mOTU with ClusterSize
identity <- 0.97 
read_count_df_clustersize_motu <- ClusterSize(read_count_samples_df,
                                     id=identity, 
                                     vsearch_path=vsearch_path,
                                     num_threads=num_threads,
                                     by_sample=FALSE)

stat_df <- GetStat(read_count_df_clustersize_motu, stat_df, stage="clustersize_motu_7", params=identity)



### ClusterSize
identity <- 0.97 
read_count_samples_df <- ClusterSize(read_count_samples_df,
                                     id=identity, 
                                     vsearch_path=vsearch_path,
                                     by_sample=FALSE)

### TaxAssign
asv_tax <- TaxAssign(asv=read_count_samples_df, 
                     taxonomy=taxonomy, 
                     blast_db=blast_db, 
                     blast_path=blast_path, 
                     num_threads=num_threads)

### WriteASVtable
outfile=file.path(outdir, "Final_asvtable_with_TaxAssign.csv")
asv_table_df <- WriteASVtable(read_count_samples_df, 
                              outfile=outfile, 
                              asv_tax=asv_tax, 
                              sortedinfo=sortedinfo_df, 
                              add_empty_samples=T, 
                              add_sums_by_sample=T, 
                              add_sums_by_asv=T, 
                              add_expected_asv=T, 
                              mock_composition=mock_composition)




#################
#################
################
# test ClusterSize after dereplicate
################

#### ClusterSize by sample
## sample_replicate
by_sample <- TRUE
read_count_df_ClusterSize_sample_replicate <- ClusterSize(read_count_df, 
                                                          id= 0.97,
                                                          vsearch_path=vsearch_path, 
                                                          num_threads=num_threads, 
                                                          by_sample=by_sample)


stat_df <- GetStat(read_count_df_ClusterSize_sample_replicate, stat_df, stage="ClusterSize_sample_replicate", params=by_sample)

## sample
read_count_df_ClusterSize_sample <- ClusterSize(read_count_sample_df, 
                                          id= 0.97,
                                          vsearch_path=vsearch_path, 
                                          num_threads=num_threads, 
                                          by_sample=by_sample)


stat_df <- GetStat(read_count_df_ClusterSize_sample, stat_df, stage="ClusterSize_sample", params=by_sample)

#### ClusterSize by sample=FALSE

## sample_replicate
by_sample <- FALSE
read_count_df_ClusterSize_all_sample_replicate <- ClusterSize(read_count_df, 
                                                              id= 0.97,
                                                              vsearch_path=vsearch_path, 
                                                              num_threads=num_threads, 
                                                              by_sample=by_sample)

stat_df <- GetStat(read_count_df_ClusterSize_all_sample_replicate, stat_df, stage="ClusterSize_all_sample_replicate", params=by_sample)

## sample
read_count_df_ClusterSize_all_sample <- ClusterSize(read_count_sample_df, 
                                             id= 0.97,
                                             vsearch_path=vsearch_path, 
                                       num_threads=num_threads, 
                                       by_sample=by_sample)


stat_df <- GetStat(read_count_df_ClusterSize_all_sample, stat_df, stage="ClusterSize_all_sample", params=by_sample)





################
# test swarm after dereplicate
################
#### Swarm by sample
## sample_replicate
by_sample <- TRUE
read_count_df_swam_sample_replicate <- Swarm(read_count_df, 
                       swarm_path=swarm_path, 
                       num_threads=num_threads, 
                       by_sample=by_sample
)
stat_df <- GetStat(read_count_df_swam_sample_replicate, stat_df, stage="Swarm_sample_replicate", params=by_sample)

## sample
read_count_df_swam_sample <- Swarm(read_count_sample_df, 
                                   swarm_path=swarm_path, 
                                   num_threads=num_threads, 
                                   by_sample=by_sample)


stat_df <- GetStat(read_count_df_swam_sample, stat_df, stage="Swarm_sample", params=by_sample)

#### Swarm by sample=FALSE

## sample_replicate
by_sample <- FALSE
read_count_df_swam_all_sample_replicate <- Swarm(read_count_df, 
                                             swarm_path=swarm_path, 
                                             num_threads=num_threads, 
                                             by_sample=by_sample
)
stat_df <- GetStat(read_count_df_swam_all_sample_replicate, stat_df, stage="Swarm_all_sample_replicate", params=by_sample)

## sample
read_count_df_swam_all_sample <- Swarm(read_count_sample_df, 
                                   swarm_path=swarm_path, 
                                   num_threads=num_threads, 
                                   by_sample=by_sample)


stat_df <- GetStat(read_count_df_swam_all_sample, stat_df, stage="Swarm_all_sample", params=by_sample)



plot <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      num_threads=num_threads,
                                      outfile="density_plot_1_15_3.png")

tmp <- PairwiseIdentity(read_count_df, 
                             min_id = 0.8, 
                             num_threads=num_threads,
                             vsearch_path=vsearch_path)


###################################
###################################
# test TaxAssignRDP
file <- "/home/meglecz/vtamR/tmp/test_rdp/9_PoolReplicates_TAS2_16S.csv"
read_count_16S <- read.csv(file)

taxa <- TaxAsssignRDP(asv=read_count_16S, confidence=0.7, max_memory=8, rm_chloroplast=FALSE)

plot_vsearch <- PairwiseIdentityPlotPerClusterIdentityThreshold(read_count_16S,
                                                                identity_min=0.9,
                                                                identity_max=0.99,
                                                                identity_increment=0.01,
                                                                min_id = 0.8, 
                                                                vsearch_path=vsearch_path)

scatterplot_vsearch <- PlotClusterClasstification(
  read_count=read_count_16S,
  taxa=taxa,
  clustering_method="vsearch",
  cluster_params=c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99),
  vsearch_path=vsearch_path,
  taxlevels= c("species", "genus"))



read_count_samples_df_ClusterSize <- ClusterASV(read_count_16S,
                                                method = "vsearch",
                                                identity = 0.97,
                                                group = TRUE,
                                                by_sample=FALSE,
                                                path=vsearch_path)

dim(read_count_16S)
dim(read_count_samples_df_ClusterSize)

meta_info <- "/home/meglecz/vtamR/tmp/test_rdp/fastqinfo.csv"
mock <- "/home/meglecz/vtamR/tmp/test_rdp/mock_composition.csv"
out <- "/home/meglecz/vtamR/tmp/test_rdp/asv_table.csv"
asv_table_df <- WriteASVtable(read_count_16S, 
                              asv_tax=taxa, 
                              sortedinfo=meta_info, 
                              pool_replicates=TRUE,
                              add_empty_samples=T, 
                              add_sums_by_sample=T, 
                              add_sums_by_asv=T, 
                              add_expected_asv=T, 
                              outfile=out,
                              mock_composition=mock
)




##########################################################################
##########################################################################
##########################################################################
# DONE

#################"smart_gzip
#### compress with smart_gzip Input file 19 Gb
t1 <- proc.time()
smart_gzip(file = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/Epi09-COI-R2_S17_R1_001.fasta",
           outfile = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/R.fasta.gz",
           remove = F,
           method = "R",
           num_threads = 0,
           quiet = F,
           compress = T)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "smart_gzip_R",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


t1 <- proc.time()
smart_gzip(file = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/Epi09-COI-R2_S17_R1_001.fasta",
           outfile = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/gzip.fasta.gz",
           remove = F,
           method = "gzip",
           num_threads = 0,
           quiet = F,
           compress = T)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "smart_gzip_gzip",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


t1 <- proc.time()
smart_gzip(file = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/Epi09-COI-R2_S17_R1_001.fasta",
           outfile = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/pigz.fasta.gz",
           remove = F,
           method = "pigz",
           num_threads = 0,
           quiet = F,
           compress = T)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "smart_gzip_pigz",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))



#### uncompress with smart_gzip
t1 <- proc.time()
smart_gzip(file = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/R.fasta.gz",
           outfile = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/R.fasta",
           remove = F,
           method = "R",
           num_threads = 0,
           quiet = F,
           compress = F)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "smart_gzip_R_unzip",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


t1 <- proc.time()
smart_gzip(file = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/R.fasta.gz",
           outfile = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/gzip.fasta",
           remove = F,
           method = "gzip",
           num_threads = 0,
           quiet = F,
           compress = F)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "smart_gzip_gzip_unzip",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))


t1 <- proc.time()
smart_gzip(file = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/R.fasta.gz",
           outfile = "/home/meglecz/vtamR_test_EPI09_COI/test_zip/pigz.fasta",
           remove = F,
           method = "pigz",
           num_threads = 0,
           quiet = F,
           compress = F)
t <- proc.time() - t1
time_df <- rbind(time_df, data.frame(Step = "smart_gzip_pigz_unzip",user = t["user.self"],system = t["sys.self"], elapsed = t["elapsed"], stringsAsFactors = FALSE))

write.csv(time_df, file="/home/meglecz/vtamR_test_EPI09_COI/time_tests/smart_gzip_linux.csv", row.names=FALSE)

#################"
# RandomSampleFastaLinux on demo dataset
#### RandomSampleFastaLinux total < n
RandomSampleFastaLinux(fasta="/home/meglecz/vtamR_demo_out/Merge_uncompress/mfzr_1_fw.fasta",
                       outfile="/home/meglecz/vtamR_demo_out/tmp/mfzr_1_fw_random_seq.fasta",
                       n=40000,
                       randseed = NULL, 
                       quiet=TRUE) 

RandomSampleFastaLinux(fasta="/home/meglecz/vtamR_demo_out/Merge_uncompress/mfzr_1_fw.fasta",
                       outfile="/home/meglecz/vtamR_demo_out/tmp/mfzr_1_fw_random_seq.fasta.gz",
                       n=40000,
                       randseed = NULL, 
                       quiet=TRUE) 

RandomSampleFastaLinux(fasta="/home/meglecz/vtamR_demo_out/Merge_uncompress/mfzr_1_fw.fasta",
                       outfile="/home/meglecz/vtamR_demo_out/tmp/mfzr_1_fw_random_seq.fasta",
                       n=40000,
                       randseed = NULL, 
                       quiet=TRUE)


RandomSampleFastaLinux(fasta="/home/meglecz/vtamR_demo_out/Merge_compress/mfzr_1_fw.fasta.gz",
                       outfile="/home/meglecz/vtamR_demo_out/tmp/mfzr_1_fw_random_seq.fasta.gz",
                       n=40000,
                       randseed = NULL, 
                       quiet=TRUE)

#### RandomSampleFastaLinux total > n
RandomSampleFastaLinux(fasta="/home/meglecz/vtamR_demo_out/Merge_uncompress/mfzr_2_fw.fasta",
                       outfile="/home/meglecz/vtamR_demo_out/tmp/mfzr_2_fw_random_seq.fasta",
                       n=40000,
                       randseed = NULL, 
                       quiet=TRUE) 

RandomSampleFastaLinux(fasta="/home/meglecz/vtamR_demo_out/Merge_uncompress/mfzr_2_fw.fasta",
                       outfile="/home/meglecz/vtamR_demo_out/tmp/mfzr_2_fw_random_seq.fasta.gz",
                       n=40000,
                       randseed = NULL, 
                       quiet=TRUE) 

RandomSampleFastaLinux(fasta="/home/meglecz/vtamR_demo_out/Merge_uncompress/mfzr_2_fw.fasta",
                       outfile="/home/meglecz/vtamR_demo_out/tmp/mfzr_2_fw_random_seq.fasta",
                       n=40000,
                       randseed = NULL, 
                       quiet=TRUE)


RandomSampleFastaLinux(fasta="/home/meglecz/vtamR_demo_out/Merge_compress/mfzr_2_fw.fasta.gz",
                       outfile="/home/meglecz/vtamR_demo_out/tmp/mfzr_2_fw_random_seq.fasta.gz",
                       n=40000,
                       randseed = NULL, 
                       quiet=TRUE)

