install.packages("devtools")
install.packages("roxygen2")

library("seqinr")
library("dplyr")
library("tidyr")

setwd("~/vtamR")
####
# define input filenames
fastqdir <- "small_test"
fileinfo <- "user_input/fileinfo_small.csv"
fastqdir <- "mfzr/sorted/"
fileinfo <- "user_input/fileinfo_mfzr.csv"

fastqdir <- "/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/"
fileinfo <-"/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/fileinfo_vtamr.csv"

fastqdir <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/"
fileinfo <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/fileinfo_vtamr.csv"

outdir <- "~/vtamR/out/bat"
outdir <- chek_dir_name(outdir)
dir.create(outdir, recursive =TRUE)
#####################################################
# Functions
#####################################################
#######################################
#######################################

chek_dir_name <- function(dir){
  
  bool <- endsWith(dir, "/")
  if(!bool){
    dir <- paste(dir, '/', sep="")
  }
  return(dir)
}
#######################################
#######################################
get_stat <- function(read_count_df, stat_df, stage, params){
  #define a temporary data frame
  df <- data.frame(parameters=character(),
                   asv_count=integer(),
                   read_count=integer(),
                   sample_count=integer(),
                   sample_replicate_count=integer())
  # get 4 different counts and place then into a data fram
  sample_repl <- paste(read_count_df$sample, read_count_df$replicate, sep="-")
  df[1,"parameters"] <-params
  df[1,"asv_count"] <-length(unique(read_count_df$asv))
  df[1,"read_count"] <-  sum(read_count_df$read_count)
  df[1,"sample_count"] <-length(unique(read_count_df$sample))
  df[1,"sample_replicate_count"] <-length(unique(sample_repl))
  # add rowname
  rownames(df) <- c(stage)
  # add new data to stat_df
  rbind(stat_df, df)
}
#######################################
#######################################
# read input only sequences from input fasta file
# transform all nt to upper case
# make a unique list of reads and count them
# return a tibble with read and nb_reads as vars
read_fasta_count_reads <- function (filename) {
  fas <- read.fasta(filename, seqonly = T)
  # transform lust to matrix (1 column)
  fas <-as.matrix(fas)
  # all nt to upper case letters
  fas <-toupper(fas)
  # transform matrix to df. If I want to do it directly from list, it will put each sequence to one separate column
  fas <- as.data.frame(fas)
  colnames(fas) <- c("asv")
  
  # make unique list of sequences and count them
  read_count_df <- fas %>%
    group_by(asv) %>%
    summarize(read_count=length(asv))
  return(read_count_df)
}
#######################################
#######################################
read_fastas_from_fileinfo <- function (fileinfo, fastqdir) {
  # read all fasta files in fileino to a read_count_df
  # read fileinfo to df
  fileinfo_df <- read.csv(fileinfo, header=T, sep="\t")
  fastqdir <- chek_dir_name(fastqdir)
  # define empty read_count_df to pool the results of variables
  read_count_df <- data.frame(asv=character(),
                              read_count=integer(),
                              plate=character(),
                              marker=character(),
                              sample=character(),
                              replicate=character())
  # read all fasta files in fileinfo and count the reads
  for(i in 1:length(fileinfo_df$file)){
    filename <- paste(fastqdir, fileinfo_df$file[i], sep = "")
    print(filename)
    if(file.size(filename) < 50){ # an empty file gzipped has 42 size => skip these files, An unzipped fasta with 80 nt is bigger than 50
      next
    }

    read_count_df_tmp <- read_fasta_count_reads(filename)
    read_count_df_tmp$plate <- fileinfo_df[i,"plate"]
    read_count_df_tmp$marker <- fileinfo_df[i,"marker"]
    read_count_df_tmp$sample <- fileinfo_df[i,"sample"]
    read_count_df_tmp$replicate <- fileinfo_df[i,"replicate"]
    read_count_df <- rbind(read_count_df, read_count_df_tmp)
  }
  rm(read_count_df_tmp)
  # reorder columns
  read_count_df <- read_count_df[, c("asv", "plate", "marker", "sample", "replicate", "read_count")]
  return(read_count_df)
}
#######################################
#######################################
# delete asv with less than min_read_count_plate in the whole dataset
delete_low_read_count_run <- function (read_count_df, min_read_count_plate) {
  df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count=sum(read_count)) %>%
    filter(read_count > min_read_count_plate)
  read_count_df <- filter(read_count_df, (asv %in% df$asv))
  return(read_count_df)
}
#######################################
#######################################
# filter_lnf_read_count
filter_lnf_read_count <- function (read_count_df, min_read_count_occurrence) {
  read_count_df <- filter(read_count_df,  (read_count >= min_read_count_occurrence))
  return(read_count_df)
}
#######################################
#######################################
# lfn_sample_replicate
filter_lnf_sample_replicate <- function (read_count_df, lnf_sample_replicate_prop) {
  
  #Make wide format from de df
  wide_read_count_df <- as.data.frame(pivot_wider(read_count_df, names_from = c(plate, marker, sample, replicate), values_from = read_count, values_fill=0, names_sep = ";"))
  # when using sweep non-numeric values make a problem => keep asv for later
  asv_list <- wide_read_count_df$asv
  wide_read_count_df$asv <- NULL
  # divide all values by the sum of the columns
  wide_read_count_df_prop <- sweep(wide_read_count_df, MARGIN = 2, STATS = colSums(wide_read_count_df), FUN = "/")
  # 0 to values under lnf_sample_replicate_prop
  wide_read_count_df_prop[wide_read_count_df_prop < lnf_sample_replicate_prop] <- 0
  # NA to values of the df with counts, where the prop has been put to zero
  wide_read_count_df[wide_read_count_df_prop == 0] <- NA
  # add asv column
  wide_read_count_df$asv <- asv_list
  # retransform to long format; avoid lines with Z values
  read_count_df_lnf_sample_replicate <- pivot_longer(wide_read_count_df, cols = -asv, names_to = c("plate", "marker", "sample", "replicate"), values_to = "read_count", names_sep = ";", values_drop_na = TRUE)
  return(read_count_df_lnf_sample_replicate)
}

#######################################
#######################################
# filter_lnf_variant
filter_lnf_variant_all <- function (read_count_df, lnf_sample_variant_prop) {
  
  #Make wide format from de df
  wide_read_count_df <- as.data.frame(pivot_wider(read_count_df, names_from = c(plate, marker, sample, replicate), values_from = read_count, values_fill=0, names_sep = ";"))
  # when using sweep non-numeric values make a problem => keep asv for later
  asv_list <- wide_read_count_df$asv
  wide_read_count_df$asv <- NULL
  # divide all values by the sum of the columns
  wide_read_count_df_prop <- sweep(wide_read_count_df, MARGIN = 1, STATS = rowSums(wide_read_count_df), FUN = "/")
  # 0 to values under lnf_sample_replicate_prop
  wide_read_count_df_prop[wide_read_count_df_prop < lnf_sample_variant_prop] <- 0
  # NA to values of the df with counts, where the prop has been put to zero
  wide_read_count_df[wide_read_count_df_prop == 0] <- NA
  # add asv column
  wide_read_count_df$asv <- asv_list
  # retransform to long format; avoid lines with Z values
  read_count_df_lnf_sample_replicate <- pivot_longer(wide_read_count_df, cols = -asv, names_to = c("plate", "marker", "sample", "replicate"), values_to = "read_count", names_sep = ";", values_drop_na = TRUE)
  return(read_count_df_lnf_sample_replicate)
}
#######################################
#######################################
# filter_lnf_variant_replicate
filter_lnf_variant_replicate <- function (read_count_df, lnf_sample_variant_replicate_prop) {
  
  #get the list of replicates
  replicate_list <- unique(read_count_df$replicate)
  # defile an empty dataframe
  read_count_df_lnf_variant_replicate  <- data.frame(asv=character(),
                                                     read_count=integer(),
                                                     plate=character(),
                                                     marker=character(),
                                                     sample=character(),
                                                     replicate=character())
  # make one dataframe for each replicate and filter them as in filter_lnf_variant
  for(i in replicate_list){
    replicate_df <- filter(read_count_df,  (replicate == i))
    replicate_df_filtered <- filter_lnf_variant(replicate_df, lnf_sample_variant_replicate_prop)
    read_count_df_lnf_variant_replicate <- rbind(read_count_df_lnf_variant_replicate, replicate_df_filtered)
  }
  return(read_count_df_lnf_variant_replicate)
}
#######################################
#######################################
# use filter_lnf_variant_replicate or filter_lnf_variant_all in function of the by_replicate parameter (default 1)
filter_lnf_variant <- function (read_count_df, lnf_variant_prop, by_replicate) {
  if(by_replicate == T){
    read_count_df_lnf_variant <- filter_lnf_variant_all(read_count_df, lnf_variant_prop)
  }
  else{
    read_count_df_lnf_variant <- filter_lnf_variant_replicate(read_count_df, lnf_variant_prop)
  }
  return(read_count_df_lnf_variant)
}
#######################################
#######################################
# pool_lfn # takes 2 df (same columns) and keeps only lines if non-zero in both
pool_2_lfn_df <- function (df1, df2) {
  merged_df <- merge(df1, df2, by = c("asv", "plate", "marker" ,"sample", "replicate"),  all = FALSE)
  test <- merged_df$read_count.x == merged_df$read_count.y
  if(FALSE %in% test){
    print("ERROR: Not all read counts are identical in data.frame comparison")
  }else{
    merged_df$read_count <- merged_df$read_count.x
    merged_df$read_count.y <- NULL    
    merged_df$read_count.x <- NULL    
  }
  return(merged_df)
}
#######################################
#######################################
# pool all count_read_df_xxx, and keep only occurrences present in all filters
pool_lfn <- function (...) {
  df_list <- list(...)
  merged <-  df_list[[1]]
  for(i in 2:length(df_list)){
    merged <- pool_2_lfn_df(merged, df_list[[i]])
  }
  return(merged)
}
#######################################
#######################################
#####################################################
# End functions
#####################################################

# Measure runtime using system.time()
start_time <- Sys.time()  # Record the start time
# define stat dataframe that will be completed with counts after each step
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())


# read input fasta files in fileinfo, demultiplex and count the number of reads in each plate-sample-replicate
read_count_df <- read_fastas_from_fileinfo(fileinfo, fastqdir)
# make stat counts
stat_df <- get_stat(read_count_df, stat_df, "input", NA)


# Eliminate variants with less then min_read_count_plate reads in the dataser (plate)
min_read_count_plate = 2
read_count_df <- delete_low_read_count_run(read_count_df, min_read_count_plate)
stat_df <- get_stat(read_count_df, stat_df, "min_read_count_plate", min_read_count_plate)
#dim(read_count_df)

###
### LFN filters
###
# lfn absolute read count
min_read_count_occurrence = 2
read_count_df_lfn_read_count <- filter_lnf_read_count(read_count_df, min_read_count_occurrence)
stat_df <- get_stat(read_count_df_lfn_read_count, stat_df, "min_read_count_occurrence", min_read_count_occurrence)
#dim(read_count_df_lfn_read_count)

# lfn sample_replicate (by column)
lnf_sample_replicate_prop = 0.1
read_count_df_lnf_sample_replicate <- filter_lnf_sample_replicate(read_count_df, lnf_sample_replicate_prop)
stat_df <- get_stat(read_count_df_lnf_sample_replicate, stat_df, "lnf_sample_replicate_prop", lnf_sample_replicate_prop)
#dim(read_count_df_lnf_sample_replicate)

# lfn sample_variant (by line)
lnf_variant_prop = 0.1
by_replicate = TRUE
param_values <- paste(lnf_variant_prop, by_replicate, sep=";")
read_count_df_lnf_variant <- filter_lnf_variant(read_count_df, lnf_variant_prop, by_replicate)
stat_df <- get_stat(read_count_df_lnf_variant, stat_df, "lnf_variant_prop", param_values)
#dim(read_count_df_lnf_variant)

# write the readcut file before lfn filtering for documentation
write.csv(read_count_df, file = paste(outdir, "read_count_after_min_read_count_plate.csv", sep=""))
# pool the results of the different filter to one dataframe; keep only occurrences that passed all filters
read_count_df <- pool_lfn(read_count_df_lfn_read_count, read_count_df_lnf_variant, read_count_df_lnf_sample_replicate)
stat_df <- get_stat(read_count_df, stat_df, "lfn_all", NA)
# delete temporary data frames
read_count_df_lfn_read_count <- NULL
read_count_df_lnf_variant <- NULL
read_count_df_lnf_sample_replicate <- NULL

write.csv(stat_df, file = paste(outdir, "count_stat.csv", sep=""))

end_time <- Sys.time()  # Record the end time
runtime <- end_time - start_time  # Calculate the runtime
print(runtime)
