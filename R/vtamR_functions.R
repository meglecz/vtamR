#' Check directory
#' 
#' Create dir if does not exists.
#' Add slash to the end of the directory name.
#' 
#' @param dir directory name
#' @export
check_dir <- function(dir){
  if(!(endsWith(dir, "/"))){
    dir <- paste(dir, '/', sep="")
  }
  if(!dir.exists(dir)){
    dir.create(dir, recursive =TRUE)
  }
  return(dir)
}

#' Get read, variant, sample and replicate counts
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param stat_df data frame with the following variables: parameters, asv_count, read_count, read_count, sample_count, sample_replicate_count
#' @param stage name of the filtering step; it is used for labeling the line in the stat_df
#' @param params parameter value (or their concatenation if more than one) used for the filtering step
#' @export
#' 
get_stat <- function(read_count_df, stat_df, stage, params=NA){
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

#' Read a fasta file, demultiplex reads and count them
#' 
#' Read only sequences from input fasta file.
#' Transform all nt to upper case.
#' Demultiplex reads to ASVs.
#' Count the number of reads of each ASVs.
#' Return a tibble with asv and nb_reads.
#' 
#' @param file tsv file with columns: plate	marker	sample	replicate	file
#' @export
read_fasta_count_reads <- function (file) {
  fas <- read.fasta(file, seqonly = T)
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

#' Read all fasta files in fileinfo file to a data frame
#' 
#' Read all fasta file in fileinfo file (plate	marker	sample	replicate	filename).
#' Demultiplex reads to ASVs.
#' Count the number of reads of each ASVs in each sample-replicate.
#' Returns a df ("asv", "plate", "marker", "sample", "replicate", "read_count").
#' 
#' @param file tsv file with columns: plate	marker	sample	replicate	file
#' @param dir name of the directory with fasta file 
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
read_fastas_from_fileinfo <- function (file, dir, write_csv=F, outdir=NA) {
  # read all fasta files in fileinfo to a read_count_df
  # read fileinfo to df
  fileinfo_df <- read.csv(file, header=T, sep="\t")
  dir <- check_dir(dir)
  # define empty read_count_df to pool the results of variables
  read_count_df <- data.frame(asv=character(),
                              read_count=integer(),
                              plate=character(),
                              marker=character(),
                              sample=character(),
                              replicate=character())
  # read all fasta files in fileinfo and count the reads
  for(i in 1:length(fileinfo_df$filename)){
    fas <- paste(dir, fileinfo_df$filename[i], sep = "")
    print(fas)
    if(file.size(fas) < 50){ #' an empty file gzipped has 42 size => skip these files, An unzipped fasta with 80 nt is bigger than 50
      next
    }

    read_count_df_tmp <- read_fasta_count_reads(fas)
    read_count_df_tmp$plate <- fileinfo_df[i,"plate"]
    read_count_df_tmp$marker <- fileinfo_df[i,"marker"]
    read_count_df_tmp$sample <- fileinfo_df[i,"sample"]
    read_count_df_tmp$replicate <- fileinfo_df[i,"replicate"]
    read_count_df <- rbind(read_count_df, read_count_df_tmp)
  }
  rm(read_count_df_tmp)
  # reorder columns
  read_count_df <- read_count_df[, c("asv", "plate", "marker", "sample", "replicate", "read_count")]
  
  if(write_csv){
     write.csv(read_count_df, file = paste(outdir, "Input.csv", sep=""))
  }
  return(read_count_df)
}

#' LFN_global_read_count
#' 
#' Eliminate ASVs with less than cutoff reads in the dataset
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param cutoff minimum number of reads for an ASV; default=10
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
LFN_global_read_count <- function (read_count_df, cutoff=10, write_csv=F, outdir=NA) {
  df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count=sum(read_count)) %>%
    filter(read_count > cutoff)
  read_count_df <- filter(read_count_df, (asv %in% df$asv))
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "LFN_global_read_count.csv", sep=""))
  }
  return(read_count_df)
}

#' LFN_read_count
#' 
#' Eliminate occurrences with less than cutoff reads
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param cutoff minimum number of reads for an occurrence; default=10; default=10
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
LFN_read_count <- function (read_count_df, cutoff=10, write_csv=F, outdir=NA) {
  read_count_df <- filter(read_count_df,  (read_count >= cutoff))
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "LFN_read_count.csv", sep=""))
  }
  return(read_count_df)
}

#' LFN_sample_replicate
#' 
#' Eliminate occurrences where the readcount/sum(reacount of the sample-replicate) is less than cutoff
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param cutoff minimum proportion for an occurrence within a sample-replicate; default=0.001
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
LFN_sample_replicate <- function (read_count_df, cutoff=0.001, write_csv=F, outdir=NA) {
  
  #Make wide format from de df
  wide_read_count_df <- as.data.frame(pivot_wider(read_count_df, names_from = c(plate, marker, sample, replicate), values_from = read_count, values_fill=0, names_sep = ";"))
  # when using sweep non-numeric values make a problem => keep asv for later
  asv_list <- wide_read_count_df$asv
  wide_read_count_df$asv <- NULL
  # divide all values by the sum of the columns
  wide_read_count_df_prop <- sweep(wide_read_count_df, MARGIN = 2, STATS = colSums(wide_read_count_df), FUN = "/")
  # 0 to values under cutoff
  wide_read_count_df_prop[wide_read_count_df_prop < cutoff] <- 0
  # NA to values of the df with counts, where the prop has been put to zero
  wide_read_count_df[wide_read_count_df_prop == 0] <- NA
  # add asv column
  wide_read_count_df$asv <- asv_list
  # re-transform to long format; avoid lines with O values
  read_count_df_lnf_sample_replicate <- pivot_longer(wide_read_count_df, cols = -asv, names_to = c("plate", "marker", "sample", "replicate"), values_to = "read_count", names_sep = ";", values_drop_na = TRUE)
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "LFN_sample_replicate.csv", sep=""))
  }
  return(read_count_df_lnf_sample_replicate)
}
#' LFN_variant_all
#' 
#' Eliminate occurrences where the read_count/sum(read_count of the asv) is less than cutoff 
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param cutoff minimum proportion for an occurrence within an asv; default=0.001
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
LFN_variant_all <- function (read_count_df, cutoff=0.001, write_csv=F, outdir=NA) {
  
  #Make wide format from de df
  wide_read_count_df <- as.data.frame(pivot_wider(read_count_df, names_from = c(plate, marker, sample, replicate), values_from = read_count, values_fill=0, names_sep = ";"))
  # when using sweep non-numeric values make a problem => keep asv for later
  asv_list <- wide_read_count_df$asv
  wide_read_count_df$asv <- NULL
  # divide all values by the sum of the columns
  wide_read_count_df_prop <- sweep(wide_read_count_df, MARGIN = 1, STATS = rowSums(wide_read_count_df), FUN = "/")
  # 0 to values under cutoff
  wide_read_count_df_prop[wide_read_count_df_prop < cutoff] <- 0
  # NA to values of the df with counts, where the prop has been put to zero
  wide_read_count_df[wide_read_count_df_prop == 0] <- NA
  # add asv column
  wide_read_count_df$asv <- asv_list
  # retransform to long format; avoid lines with 0 values
  read_count_df_lnf_sample_replicate <- pivot_longer(wide_read_count_df, cols = -asv, names_to = c("plate", "marker", "sample", "replicate"), values_to = "read_count", names_sep = ";", values_drop_na = TRUE)
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "LFN_variant_all.csv", sep=""))
  }
  return(read_count_df_lnf_sample_replicate)
}
#' LFN_variant_replicate
#' 
#' Eliminate occurrences where the read_count/sum(read_count of all asvs in replicate) is less than cutoff 
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param cutoff minimum proportion for an occurrence within a asv-replicate; default=0.001
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory 
#' @export
#' 
LFN_variant_replicate <- function (read_count_df, cutoff=0.001, write_csv=F, outdir=NA) {
  
  #get the list of replicates
  replicate_list <- unique(read_count_df$replicate)
  # defile an empty dataframe
  read_count_df_lnf_variant_replicate  <- data.frame(asv=character(),
                                                     read_count=integer(),
                                                     plate=character(),
                                                     marker=character(),
                                                     sample=character(),
                                                     replicate=character())
  # make one data frame for each replicate and filter them as in filter_lnf_variant
  for(i in replicate_list){
    replicate_df <- filter(read_count_df,  (replicate == i))
    replicate_df_filtered <- LFN_variant_all(replicate_df, cutoff, write_csv, outdir)
    read_count_df_lnf_variant_replicate <- rbind(read_count_df_lnf_variant_replicate, replicate_df_filtered)
  }
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "LFN_variant_replicate.csv", sep=""))
  }
  return(read_count_df_lnf_variant_replicate)
}

#' LFN_variant
#' 
#' If by_replicate=F: Eliminate occurrences where the read_count/sum(read_count of the asv) is less than cutoff 
#' If by_replicate=T: Eliminate occurrences where the read_count/sum(read_count of the asv with its replicate) is less than cutoff 
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param cutoff minimum proportion for an occurrence within an asv or an asv-replicate; default=0.001
#' @param by_replicate T/F; default=FALSE
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
LFN_variant <- function (read_count_df, cutoff=0.001, by_replicate=FALSE, write_csv=F, outdir=NA) {
  if(by_replicate == T){
    read_count_df_lnf_variant <- LFN_variant_replicate(read_count_df, cutoff, write_csv, outdir)
  }
  else{
    read_count_df_lnf_variant <- LFN_variant_all(read_count_df, cutoff, write_csv, outdir)
  }
  return(read_count_df_lnf_variant)
}

#' pool_2df
#' 
#' pool 2 count_read_df data frames, and keep only occurrences present in all filters
#' 
#' @param df1 data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param df2 data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @export
pool_2df <- function (df1, df2) {
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

#' pool_LFN
#' 
#' pool all count_read_df data frames, and keep only occurrences present in all filters
#' 
#' @param ... a list of data frames
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
#' 
pool_LFN <- function (... , write_csv=F, outdir=NA) {
  df_list <- list(...)
  merged <-  df_list[[1]]
  for(i in 2:length(df_list)){
    merged <- pool_2df(merged, df_list[[i]])
  }
  
  if(write_csv){
    write.csv(merged, file = paste(outdir, "LFN_pooled.csv", sep=""))
  }
  return(merged)
}

#' FilterMinReplicateNumber
#' 
#' Filter out all occurrences where the asv in not present in cutoff replicates
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param cutoff Minimum number of replicates; default=2
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
#'
FilterMinReplicateNumber <- function(read_count_df, cutoff=2, write_csv=F, outdir=NA){
  # add a temporary column with asv and sample concatenated
  read_count_df$tmp <- paste(read_count_df$asv,  read_count_df$sample, sep="-")
  # make a df_tmp containing the number of replicates for each asv-sample combination
  df_tmp <- read_count_df  %>%
    group_by(tmp) %>%
    summarize(repl_number=length(tmp))  %>%
    filter(repl_number >= cutoff)
  # keep only asv-sample if present at least in min_replicate_number replicates
  read_count_df <- filter(read_count_df, (read_count_df$tmp %in% df_tmp$tmp))
  read_count_df$tmp <- NULL
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "FilterMinReplicateNumber.csv", sep=""))
  }
  return(read_count_df)
}

#' FilterIndel
#' 
#' Filter out all ASVs, if the modulo 3 of their length is not the same as that of the majority of the ASVs
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
#' 
FilterIndel <- function(read_count_df, write_csv=F, outdir){
  # add a column with the modulo 3 of the length of the asvs
  read_count_df$mod3 <- nchar(read_count_df$asv) %% 3
  # make a tibble with modulo3 of the length of the ASVs and their count 
  # ASVs are counted as many times as they occur, so the most frequent ASV have higher weight
  # read_counts are not taken into account
  tmp <- read_count_df %>%
    group_by(mod3) %>%
    summarize(length_modulo=length(mod3)) %>%
    arrange(desc(length_modulo))
  # get the modulo 3 the most frequent
  my_modulo3 <- as.integer(tmp[1,"mod3"])
  # select only the lines with asv length compatible with the most frequent modulo3
  read_count_df <- read_count_df %>%
    filter(mod3 == my_modulo3)
  
  # delete the temporary column
  read_count_df$mod3 <- NULL
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "FilerIndel.csv", sep=""))
  }
  return(read_count_df)
}

#' codon_stops_from_genetic_code
#' 
#' Returns a vector of codon stops corresponding the the genetic code number in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
#'  
#' @param genetic_code genetic code number from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
#' @export
#' 
codon_stops_from_genetic_code <- function(genetic_code=5){
  if(genetic_code == 1){
    return(c("TAA","TAG","TGA"))
  }
  else if(genetic_code == 2){
    return(c("TAA","TAG","AGA", "AGG"))
  }
  else if(genetic_code == 3 || genetic_code == 4 || genetic_code == 5 || genetic_code == 9 || genetic_code == 10 || genetic_code == 13 || genetic_code == 21 || genetic_code == 24 || genetic_code == 25 || genetic_code == 31){
    return(c("TAA","TAG"))
  }
  else if(genetic_code == 6 || genetic_code == 14 || genetic_code == 33){
    return(c("TAG"))
  }
  else if(genetic_code == 16){
    return(c("TAA","TGA"))
  }
  else if(genetic_code == 11 || genetic_code == 12 || genetic_code == 26 || genetic_code == 28 ){
    return(c("TAA","TAG","TGA"))
  }
  else if(genetic_code == 22){
    return(c("TCA","TAA", "TGA"))
  }
  else if(genetic_code == 23){
    return(c("TTA","TAA", "TAG", "TGA"))
  }
  else if(genetic_code == 27 || genetic_code == 29 || genetic_code == 30){
    return( c("TGA"))
  }
  else{
    return(c())
  }
}


#' FilterCodonStop
#' 
#' Filter out all ASVs, if there is a codon stop in all three frames of the direct strand
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @param genetic_code genetic code number from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
#' @export
#' 
FilterCodonStop <- function (read_count_df, write_csv=F, outdir=outdir, genetic_code=5){
  
  codon_stops <- codon_stops_from_genetic_code(genetic_code=genetic_code)
  if(length(codon_stops) == 0){
    print("WARNING: The Genetic Code Number provided does not correspond to any code in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes\nThe FilterCodonStop step is skipped")
    return(read_count_df)
  }
  
  # get unique list of variants to a df
  unique_asv_df <- data.frame(asv = unique(read_count_df$asv))
  
  # define a function to transform each sequence to a list of nucleotides
  seq_to_list_of_nt <- function(seq){
    return(strsplit(seq, "")[[1]])
  }
  # define a function to group a list of nucleotides to pieces of three (codons) starting at start_pos
  pool_nt_to_codons <- function(seq, start_pos){
    return(splitseq(seq, frame = start_pos, word = 3))
  }
  # check if a codon is present among a list of codons
  check_codon_stops <- function(codon_list, codon_stops=codon_stops){
    if(any(codon_stops %in% codon_list)){
      return(as.numeric(1))
    }
    else{
      return(as.numeric(0))
    }
  }
  # apply the function seq_to_list_of_nt to each unique sequence
  seqs1 <- lapply(unique_asv_df$asv, seq_to_list_of_nt)
  
  # Go through the 3 reading frames and check if there is a codon sotp in each of them
  for(j in 0:2){
    # get codons in reading frame j
    seqs2 <- lapply(seqs1, pool_nt_to_codons, j)
    # check if there is at least one codon stop among the codons; retiurn 1 if CodonStop, 0 if not
    unique_asv_df[[j+2]] <- as.vector(lapply(seqs2, check_codon_stops, codon_stops))
  }
  
  # Convert columns 2 to 4 to numeric and make a column with the sum of the 3 frames
  unique_asv_df[, 2:4] <- sapply(unique_asv_df[, 2:4], as.numeric)
  unique_asv_df$CodonStop <- rowSums(unique_asv_df[,2:4])
  
  # Keep only sequences where there is at least one frame without codon stop
  unique_asv_df <-unique_asv_df %>%
    filter(CodonStop < 3)
  # filter out ASV from read_count_df, where here is a codon stop in all reading frame
  read_count_df <- filter(read_count_df, (asv %in% unique_asv_df$asv))
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "Input.csv", sep=""))
  }
  return(read_count_df)
}

#' write_fasta_seq_as_id
#' 
#' Write a fasta file using the sequences as ID
#'  
#' @param sequences list of sequences
#' @param filename output file name, including path
#' @export
#' 
write_fasta_seq_as_id <- function(sequences, filename) {
  # Open the file for writing
  file <- file(filename, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(sequences)) {
    header <- paste0(">", sequences[[i]])
    writeLines(c(header, sequences[[i]], ""), file)
  }
  # Close the file
  close(file)
}

#' flagPCRerror_vsearch
#' 
#' Filter out all ASVs if they are very similar (max_mismatch) to another more frequent ASV (pcr_error_var_prop) in a dataframe
#' The whole plate can be analysed et once (by_sampe=F)
#'  
#' @param unique_asv_df data frame with the following variables: asv, read_count; ASVs should be unique
#' @param pcr_error_var_prop if the proportion of read_counts of two similar ASVs is bellow pcr_error_var_prop, the less abundant is flagged as a PCR error
#' @param max_mismatch maximum number of mismatches (gaps included) to consider two ASVs as similar
#' @param outdir name of the output directory for temporary files
#' @param vsearch_path path to vsearch executable; can be empty if vsearch in the the PATH
#' @export
#' 
flagPCRerror_vsearch <- function(unique_asv_df, outdir="", vsearch_path="", pcr_error_var_prop=0.1, max_mismatch=1){
  
  # no ASV in the unique_asv_df => return a dataframe with 0 for all ASVs in PCRerror column
  if(length(unique_asv_df$asv) == 0){ 
    unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
    return(unique_asv_df)
  }
  
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste(outdir, 'tmp_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  outdir_tmp <- check_dir(outdir_tmp)
  
  # make fasta file with unique reads; use sequences as ids
  fas <- paste(outdir_tmp, 'unique.fas', sep="")
  write_fasta_seq_as_id(unique_asv_df$asv, fas)
  # vsearch --usearch_global to find highly similar sequence pairs
  vsearch_out <- paste(outdir_tmp, 'unique_vsearch_out.out', sep="")
  vsearch <- paste(vsearch_path, "vsearch --usearch_global ", fas, " --db ", fas, " --iddef 1 --self --id 0.90 --maxaccepts 0 --maxrejects 0 --userfields 'query+target+ids+aln' --userout ", vsearch_out, sep="")
  #https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/system
  system(vsearch)

  # no vsearch hit => return unique_asv_df completed with a PCRerror, with 0 for all ASVs
  if(file.size(vsearch_out) == 0){
    unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
    # Delete the temp directory
    unlink(outdir_tmp, recursive = TRUE)
    return(unique_asv_df)
  }
  
  # read vsearch results
  results_vsearch<- read.csv(vsearch_out, header = FALSE, sep="\t")
  colnames(results_vsearch) <- c("query","target","nb_ids","aln")
  # none of the values easily outputted by vsearch take into the extenal gaps as a diff => correct this, based on the alnlen and the number of identities
  results_vsearch$nb_diff <- nchar(results_vsearch$aln) - results_vsearch$nb_ids
  # delete unnecessary columns
  results_vsearch <- select(results_vsearch, -c(nb_ids, aln))
  # keep only pairs with max_mismatch differences 
  results_vsearch <- results_vsearch %>%
    filter(nb_diff <= max_mismatch)
  
  # rename columns and add read counts to query and target ASVs in results_vsearch from unique_asv_df
  results_vsearch <- rename(results_vsearch, asv = query)
  results_vsearch <- left_join(results_vsearch, unique_asv_df, by="asv")
  results_vsearch <- rename(results_vsearch, qasv = asv)
  results_vsearch <- rename(results_vsearch, asv = target)
  results_vsearch <- rename(results_vsearch, qread_count = read_count)
  results_vsearch <- left_join(results_vsearch, unique_asv_df, by="asv")  
  results_vsearch <- rename(results_vsearch, tread_count = read_count)
  results_vsearch <- rename(results_vsearch, tasv = asv)
  
  # flag target ASV as a PCR error, if low read_count compared to query ASV
  results_vsearch$PCRerror_target <- ((results_vsearch$qread_count * pcr_error_var_prop) > results_vsearch$tread_count)
  # keep only one column (tasv) with unique ASVs, that were flagged as PCRerror
  results_vsearch <- results_vsearch %>%
    filter(PCRerror_target==TRUE) %>%
    group_by(tasv) %>%
    select(tasv)
  # complete unique_asv_df with a PCRerror column
  unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
  unique_asv_df$PCRerror[unique_asv_df$asv %in% results_vsearch$tasv] <- 1
  
  # Delete the temp directory
  unlink(outdir_tmp, recursive = TRUE)
  
  return(unique_asv_df)
}

#' FilerPCRerror
#' 
#' Filter out all ASVs if they very similar (max_mismatch) to another more frequent ASV (pcr_error_var_prop)
#' The whole plate can be analyzed et once (by_sampe=F) or sample by sample
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param pcr_error_var_prop if the proportion of read_counts of two similar ASVs is bellow pcr_error_var_prop, the less abundant is flagged as a PCR error
#' @param max_mismatch maximum number of mismatches (gaps included) to consider two ASVs as similar
#' @param by_sample T/F, if T ASVs are flagged as an PCR error separately for each sample
#' @param sample_prop if by_sample=T, the ASV must be flagged as an PCRerror in sample_prop of the cases to be eliminated
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @param vsearch_path path to vsearch executable; can be empty if vsearch in the the PATH
#' @export
#' 
FilterPCRerror <- function(read_count_df, write_csv=F, outdir="", vsearch_path="", pcr_error_var_prop=0.1, max_mismatch=1, by_sample=T, sample_prop=0.8){
  
  # get unique list of ASVs with their total read_count in the run
  unique_asv_df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count = sum(read_count)) %>%
    arrange(desc(read_count))
  
  if(by_sample){ # sample by sample
    sample_list <- unique(read_count_df$sample)
    # loop over samples
    for(sample_loc in sample_list){
      # get unique list of ASVs with their total read_count in the sample
      unique_asv_df_sample <- read_count_df %>%
        filter(sample == sample_loc) %>%
        group_by(asv)%>%
        summarize(read_count = sum(read_count))%>%
        arrange(desc(read_count))
      
      # flag PCR errors; add one column to unique_asv_df for each sample with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df_sample <- flagPCRerror_vsearch(unique_asv_df_sample, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch)
      
      # remove read_count column
      unique_asv_df_sample <- unique_asv_df_sample[, -which(names(unique_asv_df_sample) == "read_count")]
      # add a column for for each sample to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df <- left_join(unique_asv_df, unique_asv_df_sample, by = "asv")
    }
  }
  else{ # whole dataset
    # add a column to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
    unique_asv_df <- flagPCRerror_vsearch(unique_asv_df, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch)
  }
  
  # count the number of times each ASV has been flagged and when it has not. Ignore NA, when the ASV is not present in the sample
  unique_asv_df$yes <- rowSums(unique_asv_df[3:ncol(unique_asv_df)] == 1, na.rm = TRUE)
  unique_asv_df$no <- rowSums(unique_asv_df[3:(ncol(unique_asv_df)-1)] == 0, na.rm = TRUE)
  # keep only ASVs, that are not flagged in sample_prop proportion of the samples where they are present  
  unique_asv_df <- unique_asv_df %>%
    filter(no/(yes+no) >= (1-sample_prop))
  
  # eliminate potential PCRerrors from read_count_df
  read_count_df <- read_count_df %>%
    filter(asv %in% unique_asv_df$asv)
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "Input.csv", sep=""))
  }
  return(read_count_df)
}

