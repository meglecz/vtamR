make_known_occurrences <- function(read_count_samples_df, fileinfo="", mock_composition="", sep=",", out=""){
  
  # read info on samples types and keep only relevant info
  fileinfo_df <- read.csv(fileinfo, header=T, sep=sep) %>%
    select(plate, marker, sample, sample_type, habitat)
  # make a pms (plate.marker.sample) column
  fileinfo_df$pms <- paste(fileinfo_df$plate, fileinfo_df$marker, fileinfo_df$sample, sep=".")
  # get unique lines to avoid replicates
  fileinfo_df <- unique(fileinfo_df)
  
  # define data frame for known occurrences
  occurrence_df <- read_count_samples_df
  occurrence_df$pms <- paste(occurrence_df$plate, read_count_samples_df$marker, read_count_samples_df$sample, sep=".")
  occurrence_df$pmsasv <- paste(occurrence_df$pms, read_count_samples_df$asv, sep=".")
  occurrence_df$action <- rep(NA, nrow(occurrence_df))
  
  # flag occurrences in negative control samples as delete
  occurrence_df <- flag_from_negative_controls(occurrence_df, fileinfo_df)
  
  # flag all expected occurrences in mock samples as "keep", NA for tolerate, and delete for all others
  occurrence_df <- flag_from_mock(occurrence_df, mock_composition, fileinfo_df, sep=sep)
  
  # flag occurrences as keep with low read count in habitat, compared to the others habitats
  occurrence_df <- flag_from_habitat(occurrence_df, fileinfo_df) 
  
  # keep nly relevant colums and line, sort data
  occurrence_df <- occurrence_df %>%
    select(plate,marker,sample,action,asv) %>%
    filter(!is.na(action)) %>%
    arrange(plate, marker, sample, action)
  # write to outfile
  write.table(occurrence_df, file=out, row.names = F, sep=sep)
}

flag_from_negative_controls <- function(occurrence_df, fileinfo_df){
  # keep only negative controls in fileinfo_df
  fileinfo_df <- fileinfo_df %>%
    filter(sample_type=="negative")
  # get unique list of pms with negative controls
  neg <- unique(fileinfo_df$pms)
  # set all occurrences to keep if in negative control
  occurrence_df$action <- ifelse(occurrence_df$pms %in% neg, "delete", occurrence_df$action)
  return(occurrence_df)
}

flag_from_mock <- function(occurrence_df, mock_composition, fileinfo_df, sep=","){
  
  mock_composition_df <- read.csv(mock_composition, header=T, sep=sep)
  mock_composition_df$pms <- paste(mock_composition_df$plate, mock_composition_df$marker, mock_composition_df$sample, sep=".")
  mock_composition_df$pmsasv <- paste(mock_composition_df$pms, mock_composition_df$asv, sep=".")
  
  # set temporarily all occurrences in mock as delete
  mock_samples <- unique(mock_composition_df$pms)
  occurrence_df$action <- ifelse(occurrence_df$pms %in% mock_samples, "delete", occurrence_df$action)
  
  # reset expected occurrences in mocks as keep
  keep_occurrences <- mock_composition_df %>%
    filter(action=="keep")
  occurrence_df$action <- ifelse(occurrence_df$pmsasv %in% keep_occurrences$pmsasv, "keep", occurrence_df$action)
  
  # reset tolerate occurrences in mocks as NA
  tolerate_occurrences <- mock_composition_df %>%
    filter(action=="tolerate")
  occurrence_df$action <- ifelse(occurrence_df$pmsasv %in% tolerate_occurrences$pmsasv, NA, occurrence_df$action)
  
  return(occurrence_df)
}


flag_from_habitat <- function(occurrence_df, fileinfo_df){
  # add habitat and sample_type to occurrence_df
  occurrence_df <- left_join(occurrence_df, fileinfo_df)
  # group by asv and habitat and count the total number of reads for each habitat-asv combination
  tmp <- occurrence_df %>%
    group_by(habitat, asv) %>%
    summarize(habitat_read_count=sum(mean_read_count)) %>%
    filter(!is.na(habitat))
  # count the number of habitats for each asv and keep only the ones present in at least two different habitats
  tmp2 <- tmp %>%
    group_by(asv) %>%
    summarize(nb_habitat=length(asv)) %>%
    filter(nb_habitat>1)
  # keep only selected asvs in tmp
  tmp <- tmp[tmp$asv %in% tmp2$asv, ]
  # get the total readcount for each asv in tmp
  tmp3 <- tmp %>%
    group_by(asv) %>%
    summarize(sum_read_count = sum(habitat_read_count))
  # add total readcount of asv to tmp and keep only lines where asv ih babitats where it is less frequent than in the others
  tmp <- left_join(tmp, tmp3, by="asv")
  tmp <- tmp[tmp$habitat_read_count/tmp$sum_read_count < 0.5, ]
  # keep only pertinent columns in tmp and add hab_action column with "delete"
  tmp <- tmp %>%
    select(habitat, asv)
  tmp$hab_action <- rep("delete", nrow(tmp))
  
  occurrence_df <- left_join(occurrence_df, tmp, by=c("habitat", "asv"))
  occurrence_df$action[which(occurrence_df$hab_action=="delete")] <- "delete"
  
  occurrence_df <- occurrence_df %>%
    select(-hab_action) %>%
    select(-sample_type) %>%
    select(-habitat)
  
  return(occurrence_df)
}





####################################################
####################################################
####################################################
####################################################

# info divers


decompress<- function(file=input_fasta, dir=fasta_dir){
  backup_wd <- getwd()
  setwd(dir)
  
  decompressed_file <- ""
  
  if(endsWith(file, ".zip")){
    unzip(file)
    decompressed_file <- sub(".zip", "", file)
  }else if(endsWith(file, ".gz")){
    # Open a connection to the gzipped file
    file_connection <- gzfile(file, "rb")  # "rb" stands for read binary mode
    
    # Read the contents of the gzipped file
    file_contents <- readLines(file_connection)
    
    # Close the file connection
    close(file_connection)
    
  }else{
    setwd(backup_wd)
    stop("Only gz and zip files are supported")
  }
  
  
  # reset wd
  setwd(backup_wd)
}


compressed_file <- compress_file("/home/meglecz/vtamR/local/out/sorted/plate1-MFZR-14ben01-1.fasta", compress="gz", remove_uncompressed=T)
compress_file <- function(file, compress="gz", remove_uncompressed=T){
  #### !!!!! zip is not working correctly in minux, since it produces and archive with the embedded directory structure in the path => correct it when using windows
  outfile <- file
  if(compress == "gz"){
    # Specify the path for the gzipped output file
    file_gz <- paste(file, ".gz", sep="")
    # Open the existing uncompressed file for reading
    file_content <- readBin(file, "raw", file.info(file)$size)
    # Create a gzipped copy of the file
    gz <- gzfile(file_gz, "wb")
    writeBin(file_content, gz)
    close(gz)
    if(remove_uncompressed){
      file.remove(file)
    }
    outfile <- file_gz
  }
  if(compress == "zip"){
    file_zip <- paste(file, ".zip", sep="")
    zip(file_zip, file, recursive = FALSE)
    if(remove_uncompressed){
      file.remove(file)
    }
    outfile <- file_zip
  }
  
  return(outfile)
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



### run blast with seqinr
library("seqinr")
# Load the query sequence from a file or define it directly
query_sequence <- read.fasta(file = "/home/meglecz/vtamR/local/small_test/small_test1_1.fas")
# Perform local sequence search using blast()
results <- blast(query = query_sequence, database = "/home/meglecz/vtamR/local/small_test/small_test1_1.fas", type = "DNA")
# Print the results
print(results)


if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")
library("DECIPHER")

# Load the query sequence from a file or define it directly
query_sequence <- readDNAStringSet("/home/meglecz/vtamR/local/small_test/small_test1_1.fas")

# Load the reference sequences from a file or define them directly
reference_sequences <- readDNAStringSet("database.fasta")

# Perform local sequence search using blastn()
results <- blastn(query_sequence, reference_sequences)

# Print the results
print(results)

### get the list of function in a package
functions_seqinr <- ls("package:seqinr", all.names = TRUE)
print(functions_seqinr)


### Define the command to run the third-party program
command <- "path/to/program"
args <- c("arg1", "arg2")
# Run the command using system2()
output <- system2(command, args, stdout = TRUE, stderr = TRUE)
# Print the output
print(output)


command <- "path/to/program arg1 arg2"
# Run the command using system()
system(command)



# start blast from R
myPipe <- pipe( "blastall -p blastp -i text.fasta -d data.fasta" )
results <- read.table( myPipe )
colnames( blastResults ) <- c( "QueryID",  "SubjectID", "Perc.Ident",
                               "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                               "S.start", "S.end", "E", "Bits" )


# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html



# Create a connection to an external command (e.g., "ls" command on Unix/Linux)
cmd <- pipe("ls", "r")  # Open the "ls" command for reading

# Read the output from the command
output <- readLines(cmd)

# Close the connection
close(cmd)

# Print the output
print(output)



makeblastdb -in local/out/small/tmp_1687945245/test.fas -dbtype nucl
blastn –db local/out/small/tmp_1687945245/test.fas –query local/out/small/tmp_1687945245/test.fas –outfmt 6 –out local/out/small/tmp_1687945245/test_blastout.out

pipe_vsearch <- pipe(vsearch)
results <- read.table( pipe_vsearch )

blast <- "blastn -db local/out/small/tmp_1687945245/test.fas -query local/out/small/tmp_1687945245/test.fas -outfmt 6"
myPipe <- pipe(blast)
results <- read.table( myPipe )
colnames( results ) <- c( "QueryID",  "SubjectID", "Perc.Ident",
                          "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                          "S.start", "S.end", "E", "Bits" )


#' swarm_from_sortreads
#' 
#' Run swarm (https://github.com/torognes/swarm) on input read_count_df data frame, pool variants of the same cluster sum reads of the underlying ASVs
#' Return a data frame with the same structure as the input
#' Swarm can be run sample by sample (by_sample=T) or for the whole data set in ine go
#' 
#' @param read_count data frame or csv file with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param outfile name of the output file with asv_id, sample, replicate, read_count and asv columns; Optional; If empty the results are not written to a file
#' @param swarm_path Path to th swarm executable (Default: TRUE)
#' @param by_sample [T/F], if TRUE, swarm run separately for each sample
#' @param num_threads Number of CPUs
#' @param swarm_d positive integer, d parameter for swarm (1 by default); maximum number of differences allowed between two amplicons, meaning that two amplicons will be grouped if they have d (or less) differences.
#' @param fastidious [T/F] when working with d = 1, perform a second clustering pass to reduce the number of small clusters (Default: TRUE)
#' @param sep separator for the output file
#' @export
#' 

swarm_from_sortreads <- function(sortedinfo, dir="", outfile="", swarm_path="", num_threads=1, swarm_d=1, fastidious=T, sep=",", by_sample=T){
  # !!!!!!!!!!!!!!!!!!!!!!! not yet finished
  # can accept df or file as an input
  if(is.character(sortedinfo)){
    # read known occurrences
    sortedinfo_df <- read.csv(sortedinfo, header=T, sep=sep)
  }else{
    sortedinfo_df <- sortedinfo
  }
  
  if(by_sample){
    # get list of samples 
    sample_list <- unique(sortedinfo_df$sample)
    
    # run swarm for each sample
    for(s in sample_list){
      print(s)
      # select occurrences for sample
      files <- sortedinfo_df %>%
        filter(sample==s) %>%
        select(filename) %>%
        distinct()
      files$filename <- paste(dir, files$filename, sep="")
      df_sample <- dereplicate_reads(files$filename)
      # run swarm
      df_sample <- run_swarm(df_sample, swarm_path=swarm_path, num_threads=num_threads, swarm_d=swarm_d, fastidious=fastidious)
      # add output of the sample to the total data frame
      out_df <- rbind(out_df, df_sample)
    }
  }else{ # run swarm for all samples together
    files <- sortedinfo_df %>%
      select(filename) %>%
      distinct()
    files$filename <- paste(dir, files$filename, sep="")
    df <- dereplicate_reads(files$filename)
    out_df <- run_swarm(df, swarm_path=swarm_path, num_threads=num_threads, swarm_d=swarm_d, fastidious=fastidious)
  }
  
  if(outfile != ""){
    write.table(out_df, file = outfile,  row.names = F, sep=sep)
  }
  return(out_df)
  
}

#' dereplicate_reads
#' 
#' reads all input file in the input vector containing file names
#' demultiplex reads and count them
#' returns a data frame with asv, asv_id (arbitrary, valid only within the data frame) and read_count columns
#' 
#' @param files vector with file names
#' @export
#' 
dereplicate_reads <- function(files){
  
  df <- data.frame("asv"=as.character(),
                   "read_count"=as.integer())
  
  for(file in files){
    tmp <- read_fasta_to_df(file) %>%
      group_by(sequence) %>%
      summarize(read_count=n()) %>%
      select(asv=sequence, read_count)
    
    df <- rbind(df, tmp) %>%
      group_by(asv) %>%
      summarize(read_count=sum(read_count))
  }
  
  df$asv_id <- rownames(df)
  return(df)
}



#' OptimizeLFNReadCountAndLFNvariant_with_prerunning_other_filters
#' 
#' Suggest optimal parametres for lfn_read_count_cutoff and lnf_variant_cutoff 
#' The the LFN_read_count and LFN_variant is run for a series of parameter value combinations followed by FilterMinReplicateNumber. 
#' For each parameter combination the number of FN, TP, and FP is reported. 
#' 
#' This script can run LFN_sample_replicate, FilterPCRerror and FilterMinReplicateNumber before optimizing LFN_variant and LFN_variant. 
#' This is helpful if known_occurrences has been established before running these filters, 
#' and thus avoid of pushing too high the cutoff for LFN_read_count and LFN_variant, due to FP occurrences that will be eliminated by other filters.
#' However the most forward solution is to run LFN_sample_replicate, FilterPCRerror, FilterMinReplicateNumber first, than make known_occurrences.
#' 
#' on the input using parameters set by the user (ideally optimized ones), to eliminate part of the noise. 

#' The results are written to data frame and to an outfile if the filename is given.
#' Users should chose the parameter setting that minimizes, FN and FP.
#'  
#' @param read_count data frame or csv file with the following variables: asv_id, sample, replicate, read_count, asv
#' @param known_occurrences  data frame or file produced by make_known_occurrences function, with known FP and TP
#' @param sep separator used in csv files
#' @param outfile name of the output file; Optional; If empty the results are not written to a file
#' @param min_lfn_read_count_cutoff the lowest cutoff value for LFN_read_count function (10 by default). 
#' @param max_lfn_read_count_cutoff the highest cutoff value for LFN_read_count function (100 by default). 
#' @param increment_lfn_read_count_cutoff values from min_lfn_read_count_cutoff to max_lfn_read_count_cutoff are tested by increment_lfn_read_count_cutoff increment (5 by default). 
#' @param min_lnf_variant_cutoff the lowest cutoff value for LFN_variant function (0.001 by default). 
#' @param max_lnf_variant_cutoff the highest value for LFN_variant function (0.01 by default).
#' @param increment_lnf_variant_cutoff values from min_lnf_variant_cutoff to max_lnf_variant_cutoff are tested by increment_lnf_variant_cutoff increment (0.001 by default). 
#' @param by_replicate T/F (False by default); see LFN_variant function
#' @param lfn_sample_replicate_cutoff cutoff value for LFN_sample_replicate (see LFN_sample_replicate function; default NA)
#' @param pcr_error_var_prop cutoff value for FilterPCRerror (see FilterPCRerror function; default NA)
#' @param vsearch_path path to vsearch executables
#' @param max_mismatch parameter for FilterPCRerror (see FilterPCRerror function; default 1)
#' @param by_sample parameter for FilterPCRerror (see FilterPCRerror function; default T)
#' @param sample_prop for FilterPCRerror (see FilterPCRerror function; default 0.8)
#' @param min_replicate_number for FilterMinReplicateNumber (see FilterMinReplicateNumber function; default 1)
#' @param verbose [T/F]  if TRUE print out progress; T by default
#' @export
#'

OptimizeLFNReadCountAndLFNvariant_with_prerunning_other_filters <- function(read_count, known_occurrences="", sep=",", outfile="", min_lfn_read_count_cutoff=10, max_lfn_read_count_cutoff=100, increment_lfn_read_count_cutoff=5, min_lnf_variant_cutoff=0.001, max_lnf_variant_cutoff=0.01, increment_lnf_variant_cutoff=0.001, by_replicate=FALSE, lfn_sample_replicate_cutoff=NA, pcr_error_var_prop=NA, vsearch_path="", max_mismatch=1, by_sample=T, sample_prop=0.8, min_replicate_number=2, verbose=T){
  #  read_count_df = optimize_read_count_df
  #  min_lfn_read_count_cutoff = 10
  #  min_lnf_variant_cutoff = 0.001
  #  rc_cutoff = 55
  #  var_cutoff = 0.05
  #  by_replicate = T
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # can accept df or file as an input
  if(is.character(known_occurrences)){
    # read known occurrences
    known_occurrences_df <- read.csv(known_occurrences, header=T, sep=sep)
  }else{
    known_occurrences_df <- known_occurrences
  }
  check_fileinfo(file=known_occurrences_df, file_type="known_occurrences", sep=sep)
  
  df <- read_count_df
  # filter by sample-replicate
  if(!is.na(lfn_sample_replicate_cutoff)){
    df <- LFN_sample_replicate(read_count_df, cutoff=lfn_sample_replicate_cutoff, sep=sep)
    df <- FilterMinReplicateNumber(df, cutoff=min_replicate_number, sep=sep)
  }
  # FilterPCRerror
  if(!is.na(pcr_error_var_prop)){
    df <- FilterPCRerror(df, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, sep=sep)
    df <- FilterMinReplicateNumber(df, cutoff=min_replicate_number, sep=sep)
  }
  
  
  # make a series of cutoff values for LFN_read_count
  rc_cutoff_list <- seq(from=min_lfn_read_count_cutoff, to=max_lfn_read_count_cutoff, by=increment_lfn_read_count_cutoff)
  # make a series of cutoff values for LFN_read_count
  var_cutoff_list <- seq(from=min_lnf_variant_cutoff, to=max_lnf_variant_cutoff, by=increment_lnf_variant_cutoff)
  
  out_df <- data.frame( 
    lfn_sample_replicate_cutoff =numeric(),
    pcr_error_var_prop =numeric(),
    lfn_read_count_cutoff=numeric(),
    lnf_variant_cutoff=numeric(),
    FN=numeric(),
    TP=numeric(),
    FP=numeric()
  )
  # go through all parameter combination and count the number of TP and FN
  
  for(rc_cutoff in rc_cutoff_list){
    df_tmp <- df
    #LFN_read_count
    df_tmp <- LFN_read_count(df_tmp, rc_cutoff)
    for(var_cutoff in var_cutoff_list){
      # LFN_variant
      df_tmp <- LFN_variant(df_tmp, var_cutoff, by_replicate=by_replicate)
      # FilterMinReplicateNumber
      df_tmp <- FilterMinReplicateNumber(df_tmp, min_replicate_number)
      # PoolReplicates
      df_tmp_sample <- PoolReplicates(df_tmp, digits=0)
      # pool readcount info and known occurrences info
      ko <- full_join(df_tmp_sample, known_occurrences_df, by=c("sample", "asv")) %>%
        filter(!is.na(action)) %>% # keep only lines mentioned in the known occurrences
        filter(!(is.na(read_count) & action=="delete")) # delete lines if asv is not present (read_count==NA) and the action is delete
      # get the number of FN (misssing expected occurrences) 
      missing <- ko %>%
        filter(is.na(read_count) & action=="keep")
      FN_count <- nrow(missing)
      # get the number of TP and FP
      ko <- ko %>%
        filter(!(is.na(read_count) & action=="keep")) %>%
        group_by(action) %>%
        summarise(TPFP=length(action)) %>%
        ungroup()
      
      TP_count <- 0
      if ("keep" %in% ko$action) {
        TP_count <- subset(ko, action == "keep")$TPFP
      }
      FP_count <- 0
      if ("delete" %in% ko$action) {
        FP_count <- subset(ko, action == "delete")$TPFP
      }
      new_line <- data.frame(lfn_sample_replicate_cutoff=lfn_sample_replicate_cutoff, pcr_error_var_prop=pcr_error_var_prop, lfn_read_count_cutoff=rc_cutoff, lnf_variant_cutoff=var_cutoff ,FN=FN_count, TP=TP_count, FP=FP_count)
      if(verbose){
        print(new_line)
      }
      out_df <- bind_rows(out_df, new_line )
    }
  }
  
  out_df <- out_df %>%
    arrange(FN, FP, lnf_variant_cutoff, lfn_read_count_cutoff)
  
  if(outfile != ""){
    write.table(out_df, file=outfile, sep=sep, row.names = F)
  }
  return(out_df)
}
