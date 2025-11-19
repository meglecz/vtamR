### Functions non used any more

#' RandomSeqR
#' 
#' Randomly select `n` sequences from each input FASTA file in the input data frame.
#' This is a wrapper to run `RandomSampleFastaR` on a series of FASTA files.
#'  
#' This function works on any operating system and can handle large files, 
#' but may be slow. For Linux-like systems, consider the `RandomSeq` function, 
#' which is faster.
#'  
#' @param fastainfo Data frame or CSV file containing a `fasta` column with input file names. 
#'   Files can be gzipped.
#' @param fasta_dir Character string: directory containing the input FASTA files.
#' @param n Positive integer: number of sequences to randomly select from each file.
#' @param outdir Character string: directory to write the output FASTA files.
#' @param randseed Positive integer: seed for random sampling. Default 0 uses a pseudo-random seed. 
#'   A non-zero value ensures reproducible results.
#' @param compress Logical: If TRUE, output files are compressed in gzip format.
#' @param quiet Logical: If TRUE, suppress informational messages; only warnings and errors are shown.
#' @param sep Character: field separator for input and output CSV files.
#' @return Updated input data frame with file name extensions adjusted (if needed) 
#'   and `read_counts` updated.
#' @examples
#' \dontrun{
#' RandomSeqR(
#'   fastainfo = fastainfo_df, 
#'   n = 100, 
#'   fasta_dir = "out/fasta", 
#'   outdir = "out/randomseq", 
#'   randseed = 2261, 
#'   compress = TRUE
#' )
#' }
#' 
#' @export
#' 
RandomSeqR <- function(fastainfo, 
                       n,
                       fasta_dir="",
                       outdir="", 
                       randseed=0, 
                       compress=T, 
                       sep=",",
                       quiet=TRUE){
  
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  
  check_dir(fasta_dir)
  check_dir(outdir)
  
  unique_fasta <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(unique_fasta)){ # go through all fasta files
    input_fasta <- unique_fasta[i]
    input_fasta_p <- file.path(fasta_dir, input_fasta)
    
    # outline name is the same as the input, but different folder
    # adjusted compression
    outfile <- input_fasta
    if(!compress && grepl("\\.gz$", input_fasta)){
      outfile <- sub("\\.gz$", "", outfile)
    }
    if(compress && !grepl("\\.gz$", input_fasta)){
      outfile <- paste(outfile, "gz", sep=".")
    }
    outfile_p <- file.path(outdir, outfile)
    seqn <- RandomSampleFastaR(fasta=input_fasta_p, 
                               outfile = outfile_p,
                               n=n, 
                               randseed = randseed,
                               quiet=quiet)
    
    fastainfo_df$fasta[which(fastainfo_df$fasta == input_fasta)] <- outfile
    fastainfo_df$read_count[which(fastainfo_df$fasta == outfile)] <- seqn
  } # end for
  write.table(fastainfo_df, file = file.path(outdir, "fastainfo.csv"),  row.names = F, sep=sep)
  return(fastainfo_df)
}


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
#' @param compress_method Character or logical. Compression method: `"pigz"`, `"gzip"`, or `"R"`.  
#'   `"pigz"` requires `pigz` to be installed and in the system path (or `pigz_path` specified).  
#'   `"gzip"` is Linux-only.  
#'   `"R"` uses `R.utils`, which is cross-platform but slower.
#' @param pigz_path Character string: Path to `pigz` executable. Only needed is pigz
#' is used for file compression, and it is not in the PATH.
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
RandomSeq_old <- function(fastainfo,
                      n, 
                      fasta_dir="",
                      outdir="", 
                      vsearch_path="vsearch",
                      compress_method= "R",
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
                             method = compress_method,
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
                               method = compress_method,
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
                                   method = compress_method,
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


#' Random Sequences Windows
#' 
#' Random select n sequences from each input fasta file. 
#' The output is the same compression type (if any) as the input.
#'  
#' This function can work on any operating systems, but it is relatively slow. 
#' Check out the `RandomSeq` function on linux-like systems.
#'  
#' @param fastainfo Data frame or csv file with a `fasta` column 
#' containing input file names. Files can be compressed in gz format.
#' @param fasta_dir Character string: directory that contains the input fasta files.
#' @param n Positive integer: the number of randomly selected sequences.
#' @param outdir Character string: output directory.
#' @param randseed Positive integer: seed for random sampling. 
#' 0 by default means to use a pseudo-random seed. 
#' A given non-zero seed produces always the same result.
#' @param compress logical: Compress output files to gzip format.
#' @param sep Field separator character in input and output csv files.
#' @returns Updated input data frame with file name extensions adjusted, 
#' if necessary and read_counts updated.
#' @examples
#' \dontrun{
#' RandomSeqWindows(fastainfo=fastainfo_df, 
#'     n=100, 
#'     fasta_dir="out/fasta", 
#'     outdir="out/radomseq", 
#'     randseed=2261, 
#'     compress=T)
#' }
#' @export
#' 
RandomSeqWindows <- function(fastainfo, 
                             n,
                             fasta_dir="",
                             outdir="", 
                             randseed=0, 
                             compress=T, 
                             sep=","
){
  
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  
  check_dir(fasta_dir)
  check_dir(outdir)
  
  unique_fasta <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(unique_fasta)){ # go through all fasta files
    input_fasta <- unique_fasta[i]
    input_fasta_p <- file.path(fasta_dir, input_fasta)
    
    selected_seq_df <- select_sequences(file=input_fasta_p, n, randseed=randseed)
    fastainfo_df$read_count[
      which(fastainfo_df$fasta == input_fasta)
    ] <- nrow(selected_seq_df)
    
    outfile <- file.path(outdir, input_fasta)
    # the extension of the outfile will be adjusted according to compress
    outfile <- write_df_to_fasta(selected_seq_df, out=outfile, compress=compress) 
    outfile <- sub(outdir, "", outfile)
    fastainfo_df$fasta[which(fastainfo_df$fasta == input_fasta)] <- outfile
    fastainfo_df$read_count[which(fastainfo_df$fasta == input_fasta)] <- nrow(selected_seq_df)
  } # end for
  write.table(fastainfo_df, file = file.path(outdir, "fastainfo.csv"),  row.names = F, sep=sep)
  return(fastainfo_df)
}

#' Decompress a gzipped file
#' 
#' Decompress the input gzipped file.
#' If remove_input is TRUE, deleted the compressed input file.
#' Not adapted for large files.
#'  
#' @param filename Character string: gzip compressed input file.
#' @param remove_input logical: Remove the input compressed file.
#' @returns Character string: output decompressed file.
#' @examples
#' \dontrun{
#' outfile <- decompress_file(filename="data/test.fasta.gz", remove_input=T)
#' }
#' @export 
#
decompress_file <- function(filename="", remove_input=F){
  # this version of compression can work in all systems, but might not 
  # work with very large files, since it reads the file
  
  # make output filename
  outfile <- gsub(".gz", "", filename)
  
  if(outfile == filename){
    stop("The input file must have .gz extention")
  }
  # read compressed file
  compressed_con <- gzfile(filename, "rb")
  text_content <- readLines(compressed_con)
  close(compressed_con)
  # write uncompressed file
  writeLines(text_content, outfile)
  if(remove_input){
    file.remove(filename)
  }
  return(outfile)
}



#' Compress file
#'  
#' Compress input file to gzip format.
#' 
#' This function work in all operating systems, but not adapetd to very large files. 
#'
#' @param filename Character string: uncompressed input file.
#' @param remove_input logical: Remove the input uncompressed file.
#' @returns Character string: output gz compressed file.
#' @examples
#' \dontrun{
#' outfile <- compress_file(filename="data/test.fasta", remove_input=T)
#' }
#' @export 
#' 
compress_file <- function(filename="", remove_input=F){
  
  # Specify the path for the gzipped output file
  outfile_gz <- paste(filename, ".gz", sep="")
  # Open the existing uncompressed file for reading
  file_content <- readBin(filename, "raw", file.info(filename)$size)
  # Create a gzipped copy of the file
  gz <- gzfile(outfile_gz, "wb")
  writeBin(file_content, gz)
  close(gz)
  
  if(remove_input){
    file.remove(filename)
  }
  return(outfile_gz)
}

#' Cluster all input ASV by swarm
#' 
#' Cluster all input ASV by swarm and return a data frame with asv_id, cluster_id
#' columns
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, asv, read_count.
#' @param swarm_d Positive integer: d for Swarm.
#' @param fastidious logical: when working with d = 1, perform a second 
#' clustering pass to reduce the number of small clusters.
#' @param swarm_path Character string: path to swarm executable. 
#' @param num_threads Positive integer: Number of CPUs. If 0, use all available CPUs.
#' @param outfile Character string: name of the output csv file 
#' with asv_id, cluster_id columns
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the following columns: asv_id, cluster_id
#' @examples 
#' \dontrun{
#' cluster_df <- GetClusterIdSwarm(read_count_df, 
#'                      swarm_d=7, 
#'                      swarm_path="swarm",
#'                      num_threads=8)
#' }
#' @export
GetClusterIdSwarm <- function(read_count, 
                              swarm_d=1, 
                              fastidious=FALSE,
                              swarm_path="swarm", 
                              num_threads=0, 
                              outfile="", 
                              sep=",", 
                              quiet=TRUE){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  ##### make df if read_count is file
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  # if called from a ClusterASV, without specifying the path, it has a "" value
  if(swarm_path==""){
    swarm_path<- "swarm"
  }
  
  # check coherence between swarm_d and fastidious
  if(fastidious && swarm_d != 1 ){
    stop("ERROR: The fastidious argument must be use with swarm_d 1")
  }
  
  #####
  # make a df with unique asv, asv_id and readcount (sum)
  asv_df <- read_count_df %>%
    group_by(asv, asv_id) %>%
    summarize(read_count = sum(read_count), .groups="drop")
  
  ##### make tmp dir and input files
  tmp_dir <-paste('tmp_swarm_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  tmp_dir <- file.path(tempdir(), tmp_dir)
  check_dir(tmp_dir)
  
  # make a fasta file with unique asv format adapted to swarm
  input_swarm <- file.path(tmp_dir, "swarm_input.fasta")
  writeLines(paste(">", asv_df$asv_id, "_", 
                   asv_df$read_count, "\n", 
                   asv_df$asv, 
                   sep=""), 
             input_swarm)
  
  # Outfile name. Each line is a cluster, with asv_ids separated  by space
  out_swarm <- file.path(tmp_dir, "out_swarm.txt")
  
  ##### run swarm
  # Build argument vector
  args <- c(
    "-d", swarm_d,
    "-o", out_swarm,
    input_swarm
  )
  if(num_threads > 0){
    args <- append(args, c("-t", num_threads), after = 2)
  }
  if(fastidious){
    args <- append(args, c("-f"), after = 2)
  }
  
  run_system2(swarm_path, args, quiet=quiet)
  
  print("SWARM finished")
  
  
  #####
  # make a data frame with asv_id and cluster_id columns, 
  # where asv_id has all swarm input asv_id, 
  # and  cluster_id is the name of the cluster they belong to
  
  # read.table is unpredictable, when the number of element is variable among lines. Use a more complicated, but more sure solution.
  # Read the file line by line
  lines <- readLines(out_swarm)
  # Split each line by whitespace
  split_lines <- strsplit(lines, "[[:space:]]+")
  # Determine the maximum number of fields
  max_cols <- max(sapply(split_lines, length))
  # Convert to a data frame and fill empty cells by NA
  cluster_df <- as.data.frame(
    do.call(
      rbind, lapply(
        split_lines,
        function(x) c(x, rep(NA, max_cols - length(x)))
      )
    ),
    stringsAsFactors = FALSE
  )
  
  ### Make output data frame with asv_id, cluster_id columns
  #  repeat each cluster_id as many times as columns
  # transpose and flatten to make asv_id
  # This will produce some lines with NA for asv_id. Filter them out afterwards.
  cluster_df <- data.frame(asv_id = as.vector(t(cluster_df[,])),
                           cluster_id = rep(cluster_df$V1, 
                                            each = ncol(cluster_df))
  )
  # delete lines with NA values in asv_id
  cluster_df <- cluster_df %>%
    filter(!is.na(asv_id))
  # delete read_counts from id
  cluster_df$cluster_id <- as.numeric(
    sub("_[0-9]+", "", cluster_df$cluster_id ))
  cluster_df$asv_id <- as.numeric(
    sub("_[0-9]+", "", cluster_df$asv_id ))
  
  unlink(tmp_dir, recursive = TRUE)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(cluster_df, file = outfile,  row.names = F, sep=sep)
  }
  return(cluster_df)
}