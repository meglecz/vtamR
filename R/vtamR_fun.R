#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first if_else
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym :=
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL

#' Check directory
#' 
#' Check if directory exists, and create it if not.
#' If it is a file path, get its directory path and create dir in necessary.
#' 
#' @param path Character string naming a directory or a file including path.
#' @param is_file logical; If TRUE, it is a file. Directory otherwise.
#' @examples 
#' \dontrun{
#' check_dir(path="data")
#' }
#' @export

check_dir <- function(path, is_file=FALSE){
  
  if(is_file){
    dir_to_create <- dirname(path)
  }else{
    dir_to_create <- path
  }
  
  if(!dir.exists(dir_to_create)){
      dir.create(dir_to_create, recursive =TRUE)
    }
}

#' GetStat
#' 
#' Get read, variant, sample and replicate counts.
#' Complete the stat_df with the above statistics.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param stat_df Data frame with the following variables: parameters, 
#' asv_count, read_count, read_count, sample_count, sample_replicate_count.
#' @param stage Character string: name of the filtering step. It is used as 
#' row name in the stat_df.
#' @param params Character string of concatenated parameter values used for 
#' the filtering step.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @returns Data frame (stat_df) completed with a new line.
#' @examples
#' \dontrun{
#' GetStat(read_count_df, 
#'    stat_df, 
#'    stage="LFNvariant", 
#'    params="0.002;by_replicate=TRUE"
#'    )
#' GetStat(read_count_df, 
#'    stat_df, 
#'    stage="FilterIndel", 
#'    params="0.002;by_replicate=TRUE", 
#'    outfile="out/ReadCount_stat.csv"
#'    )
#' }
#' @export
#' 
GetStat <- function(read_count, stat_df, stage="", params=NA, outfile=""){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  #define a temporary data frame
  df <- data.frame(parameters=character(),
                   asv_count=integer(),
                   read_count=integer(),
                   sample_count=integer(),
                   sample_replicate_count=integer())
  # get 4 different counts and place then into a data frame
  df[1,"parameters"] <-params
  df[1,"asv_count"] <-length(unique(read_count_df$asv))
  df[1,"read_count"] <-  sum(read_count_df$read_count)
  df[1,"sample_count"] <-length(unique(read_count_df$sample))
  if("replicate" %in% colnames(read_count_df)){
    sample_repl <- paste(read_count_df$sample, read_count_df$replicate, sep="-")
    df[1,"sample_replicate_count"] <-length(unique(sample_repl))
  }else{
    df[1,"sample_replicate_count"] <- NA
  }
  # add rowname
  rownames(df) <- c(stage)
  # add new data to stat_df
  stat_df <- rbind(stat_df, df)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(stat_df, file = outfile,  row.names = F, sep=sep)
  }
  return(stat_df)
}

#' Merge forward and reverse reads
#' 
#' Merge forward are reverse fastq read pairs and transform the files to fasta.
#' Output fasta files can be compressed if compress option is used.
#' The output fastainfo.csv file is similar to the fastqinfo file, but the 
#' fastq columns are replaced by a fasta column with the name of the output files.
#'  
#' @param fastqinfo Data frame or csv file with columns: tag_fw,primer_fw,tag_rv,
#' primer_rv,sample,sample_type(mock/negative/real),habitat(optional),replicate,
#' fastq_fw,fastq_rv  
#' @param fastq_dir Character string: directory with input fastq files 
#' (listed in fastqinfo_df$fastq_fw and fastqinfo_df$fastq_rv).
#' @param vsearch_path  Character string: path to vsearch executables. 
#' @param outdir Character string: output directory.
#' @param fastq_ascii ASCII character number (33 or 64) used as the basis for 
#' the FASTQ quality score.
#' @param fastq_maxdiffs Maximum number of non-matching nucleotides allowed in 
#' the overlapping region (positive integer).
#' @param fastq_maxee Discard sequences with more than the specified number of 
#' expected error (positive integer).
#' @param fastq_minlen Discard sequences with less than fastq_minlen bases 
#' (positive integer).
#' @param fastq_maxlen Discard sequences with more than fastq_maxlen bases 
#' (positive integer).
#' @param fastq_minmergelen Minimum length of the merged sequence 
#' (positive integer).
#' @param fastq_maxmergelen Maximum length of the merged sequence 
#' (positive integer).
#' @param fastq_maxns Discard sequences with more than fastq_maxns of N’s 
#' (positive integer).
#' @param fastq_truncqual Truncate sequences starting from the first base with 
#' fastq_truncqual base quality score or lower (positive integer).
#' @param fastq_minovlen Minimum overlap between the merged reads 
#' (positive integer).
#' @param fastq_allowmergestagger Boolean to allow to merge staggered read pairs. 
#' Staggered pairs are pairs where the 3’ end of the reverse read has an 
#' overhang to the left of the 5’ end of the forward read.
#' @param sep Field separator character in input and output csv files.
#' @param compress logical: Compress output files to gzip format.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame corresponding to the output fastainfo.csv file
#' @examples
#' \dontrun{
#' Merge(fastqinfo_df, 
#'    fastq_dir="data/fastqdir", 
#'    vsearch_path="C:/Users/Public", 
#'    outdir="data/fastadir", 
#'    compress=T, 
#'    quiet=F
#'    )
#' Merge(fastqinfo_df, 
#'    fastq_dir="data/fastqdir", 
#'    outdir="data/fastadir", 
#'    fastq_maxdiffs=5, 
#'    fastq_maxee=2, 
#'    fastq_minlen=60, 
#'    fastq_maxlen=100, 
#'    fastq_minmergelen=80, 
#'    fastq_maxmergelen=100, 
#'    fastq_maxns=1, 
#'    fastq_truncqual=20, 
#'    fastq_minovlen=20,
#'    fastq_allowmergestagger=T
#'    )
#' }
#' @export
#'

Merge <- function(fastqinfo, 
                  fastq_dir, 
                  vsearch_path="vsearch", 
                  outdir="", 
                  fastq_ascii=33, 
                  fastq_maxdiffs=10, 
                  fastq_maxee=1, 
                  fastq_minlen=50, 
                  fastq_maxlen=500, 
                  fastq_minmergelen=50, 
                  fastq_maxmergelen=1000, 
                  fastq_maxns=0, 
                  fastq_truncqual=10, 
                  fastq_minovlen=50, 
                  fastq_allowmergestagger=F, 
                  sep=",", 
                  compress=F, 
                  quiet=T){
  
  check_dir(fastq_dir)
  check_dir(outdir)
  
  # can accept df or file as an input
  if(is.character(fastqinfo)){
    # read known occurrences
    fastqinfo_df <- read.csv(fastqinfo, header=T, sep=sep)
  }else{
    fastqinfo_df <- fastqinfo
  }
  CheckFileinfo(file=fastqinfo_df, dir=fastq_dir, file_type="fastqinfo", sep=sep, quiet=TRUE)
  
  # get unique list of fastq file pairs
  tmp <- fastqinfo_df %>%
    select(fastq_fw, fastq_rv)
  tmp <- unique(tmp)
  tmp$fasta <- NA
  tmp$read_count <- NA
  
  for(i in 1:nrow(tmp)){# for each file pairs
    
    # use the name of the fw fastq file and replace extension by fasta (uncompressed)
    outfile <- sub("\\..*", ".fasta", tmp[i,1])
    tmp$fasta[i] <- outfile
    outfile <- file.path(outdir, outfile)
    # add path to input filenames
    fw_fastq <- file.path(fastq_dir, tmp[i,1])
    rv_fastq <- file.path(fastq_dir, tmp[i,2])
    
    #Decompress input files, since they are cannot be treated directly by vsearch on the OS
    if(!is_linux() && endsWith(fw_fastq, ".gz")){ 
      fw_fastq <- decompress_file(fw_fastq, remove_input=F)
      rv_fastq <- decompress_file(rv_fastq, remove_input=F)
    }
    
    if(fw_fastq == outfile){ # stop the run if input and output files have the same name
      stop("ERROR: Input and output directories for fastq and fasta files are indentical. 
           Please, give a different output directory")
    }
    
    # vsearch can accept gz files in linux, but not on windows, the outfile is always uncompressed
    vsearch <- paste(vsearch_path, 
                     " --fastq_mergepairs ", fw_fastq,
                     " --reverse ", rv_fastq ,
                     " --fastaout ", outfile,
                     " --quiet --fastq_ascii ",fastq_ascii,
                     " --fastq_maxdiffs ", fastq_maxdiffs, 
                     " --fastq_maxee ", fastq_maxee, 
                     " --fastq_minlen ", fastq_minlen, 
                     " --fastq_maxlen ",fastq_maxlen, 
                     " --fastq_minmergelen ",fastq_minmergelen,
                     " --fastq_maxmergelen ",fastq_maxmergelen,
                     " --fastq_maxns ", fastq_maxns, 
                     " --fastq_truncqual ", fastq_truncqual, 
                     " --fastq_minovlen ", fastq_minovlen, 
                     sep=""
                     )
    if(fastq_allowmergestagger){ # if reads are longer than the amplicon
      vsearch <- paste(vsearch, " --fastq_allowmergestagger", sep="")
    }
    if(!quiet){
      print(vsearch)
    }
    system(vsearch, show.output.on.console = FALSE)
    seq_n <- count_seq(outfile)
    tmp$read_count[i] <- seq_n

    # vsearch produces uncompressed files even if input is compressed => compress output file
    if(compress){
      if(is_linux()){
        cmd <- paste("gzip", outfile, sep=" ")
        system(cmd)
      }else{ # this version of compression can work in all systems, but might not work with very large files
        out <- compress_file(filename=outfile, remove_input=T)
      }
      
      # correct output filename in fastainfo if necessary
      if(!endsWith(tmp$fasta[i], ".gz")){ 
        tmp$fasta[i] <- paste(tmp$fasta[i], ".gz", sep="")
      }
    }
    
    original_fw_fastq <- file.path(fastq_dir, tmp[i,1])
    if( original_fw_fastq != fw_fastq){# the input fastq has been unzipped for vsearch => rm unzipped file to free space
      file.remove(fw_fastq)
      file.remove(rv_fastq)
    }
  }
  # make fastainfo file
  fastainfo_df <- left_join(fastqinfo_df, tmp, by=c("fastq_fw", "fastq_rv")) %>%
    select(-fastq_fw, -fastq_rv)
  write.table(fastainfo_df, file = file.path(outdir, "fastainfo.csv"),  row.names = F, sep=sep)
  
  return(fastainfo_df)
  
}

#' Test if OS is Linux-like
#' 
#' Returns TRUE if operating system is linux-like. 
#'  
#' @returns TRUE if "sysname" returned by Sys.info() starts with linux, 
#' sunos, darwin, gnu or unix
#' @examples
#' \dontrun{
#' os_linux <- is_linux()
#' }
#' @export 
#
is_linux <- function(){
  
  system_info <- Sys.info()
  os <- tolower(system_info["sysname"])
  
  # Check the operating system
  if (startsWith(os, "windows")) {
    return(FALSE)
  } else if (startsWith(os, "linux")) {
    return(TRUE)
  } else if (startsWith(os, "sunos")) {
    return(TRUE)
  } else if (startsWith(os, "darwin")) {
    return(TRUE)
  } else if (startsWith(os, "gnu")) {
    return(TRUE)
  } else if (startsWith(os, "unix")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Decompress a gzipped file
#' 
#' Decompress the input gzipped file.
#' If remove_input is TRUE, deleted the compressed input file.
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
#' This function work in all operating systems, but might not work with very large files, 
#' since it reads the file to memory.
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
                      randseed=0,
                      compress=F,
                      sep=",", 
                      quiet=T){
  
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
    output_fasta_p <- file.path(outdir, input_fasta)
    
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
        output <- compress_file(filename=input_fasta_p, remove_input=F) # compress input, keep original
        file.rename(output, output_fasta_p) # move compressed file to the target location
      }else{
        if(!compress && endsWith(input_fasta_p, ".gz")){ # input compressed, and compress = F
          output <- decompress_file(filename=input_fasta_p, remove_input=F) # compress input, keep original
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
    # do not transform large numbers to scentific forms, since it would lead to an error in vsearch
    options(scipen=100)
    output_fasta_p <- gsub(".gz", "", output_fasta_p) # vsearch makes decompressed files
    vsearch_cmd <- paste(vsearch_path, " --fastx_subsample ",
                         input_fasta_p, " --fastaout ", 
                         output_fasta_p, " --sample_size ", 
                         n, " --randseed ", 
                         randseed, 
                         sep="")
    if(!quiet){
      print(vsearch_cmd)
    }
    system(vsearch_cmd)
    options(scipen=0)
    
    if(compress){ # compress the output file 
      output_fasta_p <- compress_file(filename=output_fasta_p, remove_input=T)
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

#' Count sequences in fasta
#' 
#' Count the number of sequences in the input fasta file. 
#' 
#' Input can be uncompressed or gzip compressed file, 
#' but other compression types are not supported.
#'  
#' @param file Character string: input fasta file (including path).
#' @returns The number of sequences in the input file.
#' @examples
#' \dontrun{
#' n <- count_seq2(file="data/test.fasta")
#' }
#' @export
count_seq2 <- function(file){
  
  if(endsWith(file, '.zip') || endsWith(file, '.bz2') || endsWith(file, '.xz')){
    stop("File compression type is not supported")
  }
  
  # work on uncompressed and gz files
  if(endsWith(file, "gz")){
    file_connection <- gzfile(file, "rb")
  }else{
    file_connection <- file(file, "r")
  }
  data <- as.data.frame(readLines(file_connection, n = -1))
  close(file_connection)
  
  colnames(data) <- c("read")
  data <- data %>%
    filter(startsWith(read, ">"))
  
  count = nrow(data)
  rm(data)
  
  return(count)
}

#' Count sequences in fasta
#' 
#' Count the number of sequences in the input fasta file. 
#' 
#' Input can be uncompressed or gzip compressed file, 
#' but other compression types are not supported.
#' For linux-like systems, uses the bash commands grep and wc and it is quick. 
#' For other operating systems it uses count.fields, and it is slower.
#'  
#' @param file Character string: input fasta file (including path).
#' @returns The number of sequences in the input file.
#' @examples
#' \dontrun{
#' n <- count_seq(file="data/test.fasta")
#' }
#' @export
count_seq <- function(file) {
  
  if(endsWith(file, '.zip') || endsWith(file, '.bz2') || endsWith(file, '.xz')){
    stop("File compression type is not supported")
  }
  
  if(is_linux()){
    if(endsWith(file, '.gz')){
      cmd <- paste("zcat", file, "| grep '>' | wc -l", sep=" ")
    }else{
      cmd <- paste("grep '>' ",file, "| wc -l", sep=" ")
    }
    seq_count <- as.integer(system(cmd, intern=TRUE))
    return(seq_count)
  }else{ # non-linux
  
    if(endsWith(file, '.gz')){
      con <- gzfile(file, "rb")
    }else{
      con <- file(file, "r")
    }
    # get the number of fields per line, using '>' as separator
    field_count <- as.data.frame(count.fields(con, sep = ">")) 
    close(con)
    colnames(field_count) <- c("field_n")
    
    field_count <- field_count %>%
      filter(field_n > 1)
    
    seq_count <- nrow(field_count)
    
    return(seq_count)
  }
}

#' Trim primers
#' 
#' Trim primers from the input fasta file
#' 
#' Input fasta can be uncompressed or gzip, bz2 compressed file, 
#' but other compression types are not supported.
#'   
#' @param fasta Character string: input fasta file (icluding path). 
#' @param outfile Character string: output fasta file (icluding path).
#' @param primer_fw Character string: forward primer (IUPAC ambiguity codes are accepted).
#' @param primer_rv Character string: reverse primer (IUPAC ambiguity codes are accepted).
#' @param vsearch_path Character string: path to vsearch executables. 
#' @param cutadapt_path Character string: path to cutadapt executables. 
#' @param check_reverse logical: if TRUE, check the reverse complementary 
#' sequences of the input fasta as well.
#' @param primer_to_end logical: primers follow directly the tags 
#' (no heterogeneity spacer).
#' @param cutadapt_error_rate Real (0-1): maximum proportion of errors 
#' between primers and reads (for tags, exact match is required).
#' @param cutadapt_minimum_length Positive integer: minimum length of the 
#' trimmed sequence.
#' @param cutadapt_maximum_length Positive integer: maximum length of the 
#' trimmed sequence.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns
#' NA; Produces an output fasta file.
#' @examples
#' \dontrun{
#' TrimPrimer_OneFile(fasta="data/test.fasta", 
#'     outfile="out/test_trimmed.fasta", 
#'     primer_fw="TCCACTAATCACAARGATATTGGTAC", 
#'     primer_rv="WACTAATCAATTWCCAAATCCTCC", 
#'     check_reverse=T, 
#'     primer_to_end=T, 
#'     cutadapt_minimum_length=300,
#'     cutadapt_maximum_length=400
#'     )
#' }
#' @export
#' 
TrimPrimer_OneFile <- function(fasta, 
                               outfile, 
                               primer_fw, 
                               primer_rv, 
                               cutadapt_path="cutadapt", 
                               vsearch_path="vsearch", 
                               check_reverse=F, 
                               primer_to_end=T, 
                               cutadapt_error_rate=0.1,
                               cutadapt_minimum_length=50,
                               cutadapt_maximum_length=500, 
                               quiet=T
                               ){
  
  if(fasta == outfile){
    msg <- paste("Input and output filenames are identical:", fasta, "Please, change one of them!", sep=" ")
    stop(msg)
  }
  original_output <- outfile
  # if check_reverse, the output of vsearch --fastx_revcomp is uncompressed => 
  # make uncompressed outfiles for fw and rv, pool, then compress
  if(check_reverse && (endsWith(outfile, ".gz") || endsWith(outfile, ".bz2")) ){
    outfile <- sub( "\\.gz$", "", outfile)
    outfile <- sub( "\\.bz2$", "", outfile)
    if(fasta == outfile){
      msg <- paste("Input and output filenames are identical:", fasta, "Please, change one of them!", sep=" ")
      stop(msg)
    }
  }
  
  primer_rv_rc <- reverse_complement(primer_rv)
  if(primer_to_end){
    primer_trim_cmd <- paste(cutadapt_path, 
                             " --cores=0 --quiet -e ", cutadapt_error_rate ,
                             " --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length,
                             " --maximum-length ", cutadapt_maximum_length, 
                             " -g ^", primer_fw, "...", primer_rv_rc, "$ --output ", outfile,
                             " ", fasta, 
                             sep=""
                             )
  } else{
    primer_trim_cmd <- paste(cutadapt_path,
                             " --cores=0 --quiet -e ", cutadapt_error_rate ,
                             " --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length,
                             " --maximum-length ", cutadapt_maximum_length, 
                             ' -g "', primer_fw, 
                             ';min_overlap=',nchar(primer_fw),'...', primer_rv_rc,  
                             ';min_overlap=',nchar(primer_rv_rc),
                             '" --output ', outfile, " ", 
                             fasta,
                             sep=""
                             )
  }
  if(!quiet){
    print(primer_trim_cmd)
  }
  system(primer_trim_cmd)
  
  if(check_reverse){
    # exchange fw and rv primers
    primer_rv_rc  <- reverse_complement(primer_fw)
    primer_fw <- primer_rv
    # change output filename
    out_rv <- sub("\\.", "_rv.", outfile)
    
    if(primer_to_end){
      primer_trim_cmd <- paste(cutadapt_path, " --cores=0 --quiet -e ", cutadapt_error_rate,
                               " --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length,
                               " --maximum-length ", cutadapt_maximum_length,
                               " -g ^", primer_fw, "...", primer_rv_rc, "$ --output ", out_rv,
                               " ", fasta,
                               sep=""
                               )
    } else{
      primer_trim_cmd <- paste(cutadapt_path, " --cores=0 --quiet -e ", cutadapt_error_rate,
                               " --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length,
                               " --maximum-length ", cutadapt_maximum_length, 
                               ' -g "', primer_fw, ';min_overlap=',nchar(primer_fw),'...', primer_rv_rc,
                               ';min_overlap=',nchar(primer_rv_rc),'" --output ', out_rv, 
                               " ", fasta, 
                               sep=""
                               )
    }
    if(!quiet){
      print(primer_trim_cmd)
    }
    system(primer_trim_cmd)
    
    # reverse complement rv file and append it to the outfile
    if(file.size(out_rv) > 0){ # there are sequences in the reverse trimmed file
      # reverse complement sequences in out_rv file
      out_rv_rc <- sub("\\.", "_rc.", out_rv)
      rev_comp_cmd <- paste(vsearch_path, " --fastx_revcomp ", out_rv, 
                            " --fastaout ", out_rv_rc, 
                            " --quiet", 
                            sep=""
                            )
      if(!quiet){
        print(rev_comp_cmd)
      }
      system(rev_comp_cmd)
      # append content of minus_rc to plus file
      file.append(outfile, out_rv_rc)
      unlink(out_rv_rc)
    }
    unlink(out_rv)
    
    if(outfile != original_output ){ # Output should be compressed
      outfile <- compress_file(filename=outfile, remove_input=T)
    }
  }
}

#' Trim primers from fasta files
#' 
#' Trim primers from each fasta file in input fastainfo data frame or csv file. 
#' Keep only trimmed reads.
#' Can check both orientations.
#' Count the number of reads in the output files.
#'   
#' @param fastainfo Data frame or csv file with the following columns: 
#' tag_fw,primer_fw,tag_rv,primer_rv,sample,sample_type,habitat,replicate,fasta,read_count
#' @param fasta_dir Character string: directory of the fasta files to be trimmed
#' @param compress logical: Compress output files to gzip format.
#' @param outdir Character string: output directory for the trimmed fasta files.
#' @param cutadapt_path Character string: path to cutadapt executables. 
#' @param vsearch_path Character string: path to vsearch executables.
#' @param check_reverse logical: if TRUE, check the reverse complementary sequences 
#' of the input fasta as well.
#' @param primer_to_end logical: primers follow directly the tags (no heterogeneity spacer).
#' @param cutadapt_error_rate Real (0-1): maximum proportion of errors 
#' between primers and reads (for tags, exact match is required).
#' @param cutadapt_minimum_length Positive integer: minimum length of the 
#' trimmed sequence.
#' @param cutadapt_maximum_length Positive integer: maximum length of the 
#' trimmed sequence.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns fastainfo_df the input data frame, with updated files names and sequence counts.
#' @examples
#' \dontrun{
#' fastainfo_df <- TrimPrimer(fastainfo, 
#'     fasta_dir="data/fasta", 
#'     outdir="out", 
#'     compress=T, 
#'     check_reverse=T, 
#'     primer_to_end=F, 
#'     cutadapt_error_rate=0.1, 
#'     cutadapt_minimum_length=300, 
#'     cutadapt_maximum_length=350, 
#'     quiet=T
#'     )
#' }
#' @export
#' 
TrimPrimer <- function(fastainfo, 
                       fasta_dir="", 
                       outdir="", 
                       compress=F, 
                       cutadapt_path="cutadapt", 
                       vsearch_path="vsearch", 
                       check_reverse=F, 
                       primer_to_end=T, 
                       cutadapt_error_rate=0.1, 
                       cutadapt_minimum_length=50, 
                       cutadapt_maximum_length=500, 
                       quiet=T
                       ){
  
  check_dir(fasta_dir)
  check_dir(outdir)
  
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  CheckFileinfo(file=fastainfo_df, dir=fasta_dir, file_type="fastainfo", sep=sep, quiet=TRUE)
  
  # upper case for all primers and tags
  fastainfo_df$primer_fw <- toupper(fastainfo_df$primer_fw)
  fastainfo_df$primer_rv <- toupper(fastainfo_df$primer_rv)
  # make a column for output filenames
  fastainfo_df$filename <- NA
  
  # check dirs
  check_dir(outdir)
  check_dir(fasta_dir)
  
  for(i in 1:nrow(fastainfo_df)){ # for each input fasta
    
    # define output file name
    input <- fastainfo_df$fasta[i]
    output <- fastainfo_df$fasta[i]
    if(compress && !endsWith(output, ".gz")){ # add .gz if necessary
      output <- sub("\\..+", ".fasta.gz", output)
    }
    fastainfo_df$filename[i] <- output
    input <- file.path(fasta_dir, input)
    output <- file.path(outdir, output)
    TrimPrimer_OneFile(input, 
                       outfile=output, 
                       primer_fw=fastainfo_df$primer_fw[i], 
                       primer_rv=fastainfo_df$primer_rv[i], 
                       cutadapt_path=cutadapt_path, 
                       vsearch_path=vsearch_path, 
                       check_reverse=check_reverse, 
                       primer_to_end=primer_to_end, 
                       cutadapt_error_rate=cutadapt_error_rate, 
                       cutadapt_minimum_length=cutadapt_minimum_length, 
                       cutadapt_maximum_length=cutadapt_maximum_length, 
                       quiet=quiet
                       )
    # count reads
    seq_n <- count_seq(output)
    fastainfo_df$read_count[i] <- seq_n
  }
  fastainfo_df <- fastainfo_df %>%
    select(sample, sample_type, habitat, replicate, fasta=filename, read_count)
  
  write.table(fastainfo_df, 
              file = file.path(outdir, "sortedinfo.csv"),  
              row.names = F, 
              sep=sep
              )
  return(fastainfo_df)
}

#' Demultiplex and trim off tags and primers
#' 
#' Demultiplex each input fasta file using the tag combinations at the 
#' extremities of the merged reads.
#' Trim primers from demultiplexed reads.
#' 
#' The output sortedinfo.csv file is similar to the input fastainfo, 
#' but the but do not have tag and primer columns.
#'  
#' @param fastainfo Data frame or csv file with columns: 
#' tag_fw,primer_fw,tag_rv,primer_rv,sample,sample_type(mock/negative/real),
#' habitat(optional),replicate,fasta  
#' @param fasta_dir Character string: directory with input fasta files 
#' (listed in the fasta columns of fastainfo).
#' @param vsearch_path Character string: path to vsearch executables. 
#' @param cutadapt_path Character string: path to cutadapt executables. 
#' @param outdir Character string: output directory.
#' @param check_reverse logical: if TRUE, check the reverse complementary 
#' sequences of the input fasta as well.
#' @param tag_to_end logical: tags are at the extremity of the reads 
#' (starting at the first base).
#' @param primer_to_end logical: primers follow directly the tags 
#' (no heterogeneity spacer).
#' @param cutadapt_error_rate Real (0-1): maximum proportion of errors 
#' between primers and reads (for tags, exact match is required).
#' @param cutadapt_minimum_length Positive integer: minimum length of the 
#' trimmed sequence.
#' @param cutadapt_maximum_length Positive integer: maximum length of the 
#' trimmed sequence.
#' @param sep Field separator character in input and output csv files.
#' @param compress logical: Compress output files to gzip format.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame corresponding to the output sortedinfo.csv file 
#' and one fasta file for each tag combination for each input fasta file, 
#' containing trimmed reads.
#' @examples
#' \dontrun{
#' fastainfo_df <- SortReads(fastainfo=fastainfo_df, 
#'      fasta_dir="data/fasta", 
#'      outdir="data/sorted", 
#'      check_reverse=F, 
#'      tag_to_end=T, 
#'      primer_to_end=T, 
#'      cutadapt_minimum_length=300,
#'      cutadapt_maximum_length=350, 
#'      sep=","
#'      )
#' }
#' @export

SortReads <- function(fastainfo, 
                      fasta_dir, 
                      outdir="", 
                      cutadapt_path="cutadapt",
                      vsearch_path="vsearch", 
                      check_reverse=F, 
                      tag_to_end=T, 
                      primer_to_end=T, 
                      cutadapt_error_rate=0.1,
                      cutadapt_minimum_length=50,
                      cutadapt_maximum_length=500,
                      sep=",",
                      compress=F,
                      quiet=T
                      ){
  
  check_dir(fasta_dir)
  check_dir(outdir)
  
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  
  CheckFileinfo(file=fastainfo_df, dir=fasta_dir, file_type="fastainfo", sep=sep, quiet=TRUE)
  
  #########
  # SortReads_no_reverse does the whole demultilexing, trimming and compress on the + strand
  # If sequences are not oriented, the -strand should be checked => 
  # run SortReads_no_reverse of plus strand and on - strand after 
  # swapping fw and rev tags and primers,
  # take the reverse complement of the -strand results (vsearch)
  # pool the results of the 2 strands
  # compress if necessary
  
  # run on strand +
  if(check_reverse){
    #### use +strand, output to sorted_dir, uncompressed
    sortedinfo_df <- SortReads_no_reverse(fastainfo_df, 
                                          fasta_dir=fasta_dir, 
                                          outdir=outdir, 
                                          cutadapt_path=cutadapt_path, 
                                          tag_to_end=tag_to_end, 
                                          primer_to_end=primer_to_end, 
                                          cutadapt_error_rate=cutadapt_error_rate, 
                                          cutadapt_minimum_length=cutadapt_minimum_length,
                                          cutadapt_maximum_length=cutadapt_maximum_length,
                                          sep=sep, 
                                          compress=F,
                                          quiet=quiet
                                          )
    
    #### use - strand
    # swap fw and rv tags and primers
    fastainfo_df_tmp <- fastainfo_df %>%
      select(tag_fw_tmp = tag_rv, 
             tag_rv_tmp = tag_fw, 
             primer_fw_tmp = primer_rv, 
             primer_rv_tmp = primer_fw, 
             sample, 
             sample_type,habitat, 
             replicate, 
             fasta) %>%
      select(tag_fw = tag_fw_tmp, 
             tag_rv = tag_rv_tmp, 
             primer_fw = primer_fw_tmp, 
             primer_rv = primer_rv_tmp, 
             sample, 
             sample_type,
             habitat, 
             replicate, 
             fasta
             )
    # make temp dir 
    check_dir(outdir)
    rc_dir <- paste('rc_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    rc_dir <- file.path(tempdir(), rc_dir)
    check_dir(rc_dir)
    # run sortreads on for reverse strand
    sortedinfo_df <- SortReads_no_reverse(fastainfo_df_tmp, 
                                          fasta_dir=fasta_dir, 
                                          outdir=rc_dir, 
                                          cutadapt_path=cutadapt_path, 
                                          tag_to_end=tag_to_end, 
                                          primer_to_end=primer_to_end, 
                                          cutadapt_error_rate=cutadapt_error_rate, 
                                          cutadapt_minimum_length=cutadapt_minimum_length,
                                          cutadapt_maximum_length=cutadapt_maximum_length,
                                          sep=sep,
                                          compress=F, 
                                          quiet=quiet
                                          )
    
    ### reverse complement and pool
    # get list of files demultiplexed on - strand
    files <- list.files(path = rc_dir, pattern=".fasta")
    # Filter the files based on the motif using regular expressions
    # reverse complement sequences on the minus stand, and append info to the plus strand output
    files <- grep(pattern = "\\.fasta", x = files, value = TRUE)
    for(i in 1:length(files)){
      plus <- file.path(outdir, files[i])
      minus <- file.path(rc_dir, files[i])
      minus_rc <- paste("rc_", files[i], sep="")
      minus_rc <- file.path(rc_dir, minus_rc)
      if(file.exists(minus) && file.size(minus) > 0){
        # reverse complement sequences in minus file
        rev_comp_cmd <- paste(vsearch_path, " --fastx_revcomp ", minus, 
                              " --fastaout ", minus_rc,
                              " --quiet", 
                              sep=""
                              )
        if(!quiet){
          print(rev_comp_cmd)
        }
        system(rev_comp_cmd)
        # append content of minus_rc to plus file
        file.append(plus, minus_rc)
      }
    }
    
    # delete temporary reverse_comp dir
    unlink(rc_dir, recursive = TRUE)
    
    ### compress
    if(compress){
      
      for(i in 1:nrow(sortedinfo_df)){
        
        file <- sortedinfo_df$fasta[i]
        sortedinfo_df$fasta[i] <- paste(file, ".gz", sep="") # correct output filename
        file <- file.path(outdir, file) # add path
        file_gz <- compress_file(file, remove_input=T) # compress file
      }
    }
  }
  else{
    # check only + strand
    sortedinfo_df <- SortReads_no_reverse(fastainfo_df, 
                                          fasta_dir=fasta_dir,
                                          outdir=outdir, 
                                          cutadapt_path=cutadapt_path, 
                                          tag_to_end=tag_to_end, 
                                          primer_to_end=primer_to_end, 
                                          cutadapt_error_rate=cutadapt_error_rate, 
                                          cutadapt_minimum_length=cutadapt_minimum_length, 
                                          cutadapt_maximum_length=cutadapt_maximum_length, 
                                          sep=sep, 
                                          compress=compress, 
                                          quiet=quiet
                                          )
  }
  
  sortedinfo_df <- get_read_counts(sortedinfo_df, dir=outdir)
  write.table(sortedinfo_df, file = file.path(outdir, "sortedinfo.csv"),  row.names = F, sep=sep)
  
  return(sortedinfo_df)
}

#' Add read_count to df
#' 
#' Count the number of reads in all fatsa files in the fasta column of a df.
#' Add read_count column to df.
#'  
#' @param df Data frame with fasta column containing name of fasta files 
#' @param dir Character string: directory name with the fasta files.
#' @returns  Data frame with the read_count column added to the input.
#' @examples
#' \dontrun{
#' df <- get_read_counts(df, sorted_dir)
#' }
#' @export
#' 
get_read_counts <- function(df, dir){
  
  df$read_count <- NA
  
  fastas <- unique(df$fasta)
  
  for(file in fastas){
    file_path <- file.path(dir, file)
    read_n <- count_seq(file_path)
    df$read_count[which(df$fasta==file)] <- read_n
  }
return(df)
}

#' Demultiplex and trim off tags and primers without checking reverse strand
#' 
#' Same as SortReads, but do not check the reverse complement of the sequences.
#' Demultiplex each input fasta file using the tag combinations 
#' at the extremities of the merged reads.
#' Trim primers from demultiplexed reads.
#' 
#' The output sortedinfo.csv file is similar to the fastainfo file, 
#' but the but do not have tag and primer columns.
#'  
#' @param fastainfo Data frame or csv file with columns: 
#' tag_fw,primer_fw,tag_rv,primer_rv,sample,sample_type(mock/negative/real),habitat(optional),replicate,fasta  
#' @param fasta_dir Character string: directory with input fasta files 
#' (listed in the fasta columns of fastainfo).
#' @param cutadapt_path Character string: path to cutadapt executables.
#' @param outdir Character string: output directory.
#' @param tag_to_end logical: tags are at the extremity of the reads 
#' (starting at the first base).
#' @param primer_to_end logical: primers follow directly the tags 
#' (no heterogeneity spacer).
#' @param cutadapt_error_rate Real (0-1): maximum proportion of errors 
#' between primers and reads (for tags, exact match is required).
#' @param cutadapt_minimum_length Positive integer: minimum length of the 
#' trimmed sequence.
#' @param cutadapt_maximum_length Positive integer: maximum length of the 
#' trimmed sequence.
#' @param sep Field separator character in input and output csv files.
#' @param compress logical: Compress output files to gzip format.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns  Data frame corresponding to the output sortedinfo.csv file 
#' and one fasta file for each tag combination for each input fasta file, 
#' containing trimmed reads.
#' @examples
#' \dontrun{
#' fastainfo_df <- SortReads_no_reverse(fastainfo=fastainfo_df, 
#'     fasta_dir="data/fasta", 
#'     outdir="data/sorted", 
#'     tag_to_end=T, 
#'     primer_to_end=T, 
#'     cutadapt_minimum_length=300,
#'     cutadapt_maximum_length=350, 
#'     sep=","
#'     )
#' }
#' @export
#' 
SortReads_no_reverse <- function(fastainfo, 
                                 fasta_dir, 
                                 outdir="", 
                                 cutadapt_path="cutadapt", 
                                 tag_to_end=T, 
                                 primer_to_end=T, 
                                 cutadapt_error_rate=0.1,
                                 cutadapt_minimum_length=50,
                                 cutadapt_maximum_length=500, 
                                 sep=",",  
                                 compress=F, 
                                 quiet=T
                                 ){
  # do the complete job of demultiplexing and trimming of input file without checking the reverse sequences
  
  check_dir(fasta_dir)
  check_dir(outdir)
  
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  
  # upper case for all primers and tags
  fastainfo_df$tag_fw <- toupper(fastainfo_df$tag_fw)
  fastainfo_df$tag_rv <- toupper(fastainfo_df$tag_rv)
  fastainfo_df$primer_fw <- toupper(fastainfo_df$primer_fw)
  fastainfo_df$primer_rv <- toupper(fastainfo_df$primer_rv)
  # make a column for output filenames
  fastainfo_df$filename <- NA
  
  # get unique list of input fasta files
  fastas <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(fastas)){ # for each input fasta
    # select lines in fastainfo_df that corresponds to a given input fasta file
    fasta_file <- fastas[i]
    df <- fastainfo_df %>%
      filter(fasta==fasta_file)
    
    # Make a tmp_dir_fasta in tempdir specific to a fasta file. It will contain the tagtrimmed files.
    # This can be deleted at the end and avoid reusing tagtrimmed files created for a preious fasta file
    tmp_fasta_file <- paste(fasta_file, "_", trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    tmp_dir_fasta <- file.path(tempdir(), tmp_fasta_file)
    
    # Delete tmp_dir_fasta if exists (previous run crushed before deleting it) 
    if (dir.exists(tmp_dir_fasta)) {
      unlink(tmp_dir_fasta, recursive = TRUE)
    }
    # Create it 
    dir.create(tmp_dir_fasta)

    # make a tags.fasta file with all tag combinations of the fasta to be demultiplexed
    tag_file <- make_adapter_fasta(fastainfo_df, 
                                   fasta_file=fasta_file, 
                                   tag_to_end=tag_to_end, 
                                   outdir=tmp_dir_fasta
                                   )

     # add path
    fasta_file <- file.path(fasta_dir, fasta_file)
    # demultiplex fasta, write output to tmp file
    demultiplex_cmd = paste(cutadapt_path, 
                            " --cores=0 --quiet -e 0 --no-indels --trimmed-only -g file:",
                            tag_file," -o ", tmp_dir_fasta, "/tagtrimmed-{name}.fasta ",
                            fasta_file, 
                            sep=""
                            )
    if(!quiet){
      print(demultiplex_cmd)
    }
    system(demultiplex_cmd)
    
      # for a given marker, there is only one primer combination
      primer_fwl <- df[1,"primer_fw"]
      primer_rvl <- df[1,"primer_rv"]
      primer_rvl_rc <- reverse_complement(primer_rvl)
      
      for(f in 1:nrow(df)){# go through each de-multiplexed, tag-trimmed file and trim primers
        outfilename <- paste(df[f,"sample"], df[f,"replicate"], sep="-")
        outfilename <- paste(outfilename, ".fasta", sep="")
        if(compress){
          outfilename <- paste(outfilename, ".gz", sep="")
        }
        # complete fastainfo_df with output fasta name
        fastainfo_df$filename[
          which(fastainfo_df$sample==df[f,"sample"] & 
                  fastainfo_df$replicate==df[f,"replicate"])
          ]<- outfilename
        # add path to output file
        primer_trimmed_file <- file.path(outdir, outfilename)
        tag_trimmed_file <- paste("tagtrimmed-", 
                                  df[f,"tag_fw"], "-", 
                                  df[f,"tag_rv"], 
                                  ".fasta", 
                                  sep=""
                                  )
        tag_trimmed_file <- file.path(tmp_dir_fasta, tag_trimmed_file)
        if(primer_to_end){
          primer_trim_cmd <- paste(cutadapt_path, 
                                   " --cores=0 --quiet -e ",cutadapt_error_rate ,
                                   " --no-indels --trimmed-only --minimum-length ", 
                                   cutadapt_minimum_length,
                                   " --maximum-length ", cutadapt_maximum_length, 
                                   " -g ^", primer_fwl, "...", primer_rvl_rc, "$ --output ", 
                                   primer_trimmed_file, " ", 
                                   tag_trimmed_file, 
                                   sep=""
                                   )
        }
        else{
          primer_trim_cmd <- paste(cutadapt_path, 
                                   " --cores=0 --quiet -e ",cutadapt_error_rate ,
                                   " --no-indels --trimmed-only --minimum-length ",
                                   cutadapt_minimum_length ,
                                   " --maximum-length ", cutadapt_maximum_length, 
                                   ' -g "', primer_fwl, 
                                   ';min_overlap=',nchar(primer_fwl),'...', primer_rvl_rc,
                                   ';min_overlap=',nchar(primer_rvl_rc),
                                   '" --output ', primer_trimmed_file, " ", 
                                   tag_trimmed_file, 
                                   sep=""
                                   )
        }
        if(!quiet){
          print(primer_trim_cmd)
        }
        system(primer_trim_cmd)
      } # end tag-trimmed 
    # delete the tmp dir with the tag-trimmed files
    unlink(tmp_dir_fasta, recursive = TRUE)
  }# end fasta
  
  # make sortedinfo file
  fastainfo_df <- fastainfo_df %>%
    select(-fasta) %>%
    select(sample, sample_type, habitat, replicate, "fasta" = filename)
    
  return(fastainfo_df)
}

#' Make a fasta file with adapters
#' 
#' Make a fasta file with tag combinations using cutadapt syntax. 
#' This will be used in SortReads to demultiplex the input fasta file.
#' 
#' @param fastainfo_df Data frame with columns: tag_fw,tag_rv,fasta  
#' @param fasta_file Character string: fasta file that needs 
#' to be demultiplexed (present in the fasta column of fastainfo_df).
#' @param outdir Character string: output directory.
#' @param tag_to_end logical: tags are at the extremity of the reads.
#' @returns tags.fasta file containing the tag combinations present in fasta_file. 
#' NA if all tags are NA in fastainfo_df for fasta_file.
#' @examples 
#' \dontrun{
#' make_adapter_fasta(fastainfo_df=fastainfo_df, 
#'     fasta_file="test.fasta", 
#'     tag_to_end=F, 
#'     outdir="data/out"
#'     )
#' }
#' @export
#'
make_adapter_fasta <- function(fastainfo_df, fasta_file, tag_to_end=T, outdir=""){
  
  # select tag combinations for the fasta file
  tags <- fastainfo_df %>%
    filter(fasta==fasta_file) %>%
    select(tag_fw, tag_rv)
    
  # make unique tag combinations and add necessary columns
  tags <- unique(tags)
  
  # return NA if all tags are NA
  if(nrow(tags) == 1){ # only one tag combination
    if(is.na(tags$tag_fw[1]) && is.na(tags$tag_rv[1])){ # no tags
      return(NA)
    }
  }
  
  tags$tag_fw <- toupper(tags$tag_fw)
  tags$tag_rv <- toupper(tags$tag_rv)
  tags$tag_rv_rc <- lapply(tags$tag_rv, reverse_complement)
  tags$tag_fwl <- lapply(tags$tag_fw, nchar)
  tags$tag_rvl <- lapply(tags$tag_rv, nchar)
  
  # Specify the file path
  check_dir(outdir)
  tag_file <- file.path(outdir, "tags.fasta")
  # initialize the content of the tag_file
  text <- c()
  if(tag_to_end){
    #>tag_fw-tag_rv
    #^tcgatcacgatgt...gctgtagatcgaca$
    for(j in 1:nrow(tags)){
      title <- paste(">", tags[j,"tag_fw"], "-",   tags[j,"tag_rv"], sep="")
      seq <- paste("^", tags[j,"tag_fw"], "...",   tags[j,"tag_rv_rc"], "$",sep="")
      text <- c(text, c(title,seq))
    }
  }else{
    #>tag_fw-tag_rv
    #tcgatcacgatgt;min_overlap=13...gctgtagatcgaca;min_overlap=14
    for(j in 1:nrow(tags)){
      title <- paste(">", tags[j,"tag_fw"], "-",   tags[j,"tag_rv"], sep="")
      seq <- paste(tags[j,"tag_fw"], 
                   ";min_overlap=", 
                   tags[j,"tag_fwl"], "...",  
                   tags[j,"tag_rv_rc"], 
                   ";min_overlap=", 
                   tags[j,"tag_rvl"], 
                   sep=""
                   )
      text <- c(text, c(title,seq))
    }
  }
  # write file
  writeLines(text, tag_file)
  return(tag_file)
}


#' Reverse Complement a sequence
#' 
#' Reverse complement a DNA sequence. IUPAC codes are accepted
#' 
#' @param sequence Character string: DNA sequence.
#' @returns
#' The reverse complement of the input sequence.
#' @examples
#' \dontrun{
#' reverse_complement(sequence="AAATGCRC")
#' }
#' @export
#'

reverse_complement <- function(sequence){
  # define complementary bases
  comp <- data.frame(orig=c("A","T","C","G","R","Y","W","S","M","K","B","H","D","V","N",
                            "a","t","c","g","r","y","w","s","m","k","b","h","d","v","n"),
                     complement=c("T","A","G","C","Y","R","W","S","K","M","V","D","H","B","N",
                                  "t","a","g","c","y","r","w","s","k","m","v","d","h","b","n")
  )
  # revers, split sequences and make data frame
  sequence_df <- data.frame(reverse=rev(strsplit(sequence, NULL)[[1]]))
  # add complimentary nt
  sequence_df <- left_join(sequence_df, comp, by=c("reverse"="orig"))
  if(any(is.na(sequence_df$complement))){
    print(sequence)
    stop("ERROR: Sequence contains non-IUPAC character")
  }
  # collapse vector to string
  reverse_comp <- paste(sequence_df$complement, collapse = "")
  
  return(reverse_comp)
}

#' Read all fasta files to a data frame and dereplicate
#' 
#' Read all fasta files in the fasta columns of sortedinfo data frame (or csv file).
#' Dereplicate reads to ASVs. 
#' Count the number of reads of each ASV in each input file.
#' Add a unique asv_id to each asv. 
#' 
#' If asv_list is given (containing earlier asv and asv_id pairs), 
#' uses already existing asv_id when possible for the present data
#' and adds new unique asv_id, to new ASVs.
#' If updated_asv_list is given, writes an updated file containing 
#' all asv - asv_id pairs. 
#' 
#' @param sortedinfo Data frame or csv file with columns: 
#' sample,replicate,fasta, optional: sample_type,habitat. 
#' The fasta column contains the names of the fasta file to be dereplicated.
#' @param dir Character string naming of the directory with input fasta files. 
#' @param outfile Character string: csv for the output data frame 
#' (asv_id, sample, replicate, read_count). If empty, no file is written. 
#' @param sep Field separator character in input and output csv files.
#' @param asv_list Character string: file with asvs and asv_ids 
#' from earlier analyses. Optional. It is used to homogenize asv_ids 
#' between different data sets.
#' @param updated_asv_list Character string naming of the output file 
#' with updated asv_id - asv pairs. Optional.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the following columns:
#' asv_id,sample,replicate,read_count,asv
#' @examples
#' \dontrun{
#' Dereplicate(sortedinfo=sortedinfo, dir="data/sorted", asv_list="data/asv_list.csv")
#' }
#' @export
#' 
Dereplicate <- function(sortedinfo, 
                        dir="", 
                        outfile="", 
                        sep=",", 
                        asv_list="", 
                        updated_asv_list="", 
                        quiet=T
                        ){
  # can accept df or file as an input
  if(is.character(sortedinfo)){
    # read known occurrences
    sortedinfo_df <- read.csv(sortedinfo, header=T, sep=sep)
  }else{
    sortedinfo_df <- sortedinfo
  }
  CheckFileinfo(file=sortedinfo_df, 
                dir=dir, 
                file_type="sortedinfo", 
                sep=sep, 
                quiet=TRUE
                )
  
  # read all fasta files in sortedinfo to a read_count_df
  if(nchar(dir)>0){
    check_dir(dir)
  }
  # define empty read_count_df to pool the results of variables
  read_count_df <- data.frame(asv=character(),
                              read_count=integer(),
                              sample=character(),
                              replicate=character())
  # read all fasta files in sortedinfo and count the reads
  for(i in 1:length(sortedinfo_df$fasta)){
    fas <- file.path(dir, sortedinfo_df$fasta[i])
    if(!quiet){
      print(fas)
    }
    # an empty file gzipped has 42 size => 
    # skip these files, An unzipped fasta with 80 nt is bigger than 50
    if(file.size(fas) < 50){ 
      next
    }
    # returns data frame with asv and read_count columns
    read_count_df_tmp <- read_fasta_to_df(fas, dereplicate=T) 
    read_count_df_tmp$sample <- sortedinfo_df[i,"sample"]
    read_count_df_tmp$replicate <- sortedinfo_df[i,"replicate"]
    read_count_df <- rbind(read_count_df, read_count_df_tmp)
  }
  rm(read_count_df_tmp)
  # reorder columns
  read_count_df <- read_count_df[, c("asv", "sample", "replicate", "read_count")]
  
  # add asv_id and write updated_asv_list if filename is given
  read_count_df <- add_ids(read_count_df, 
                           asv_list=asv_list, 
                           sep=sep,  
                           updated_asv_list=updated_asv_list
                           )
  
  # write read_count table
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Add ids to ASVs
#' 
#' Add asv_ids to a data frame or csv that have an asv column. 
#' Can take into account already existing asv - asv_id pairs 
#' (from earlier data sets) present in asv_list.
#' 
#' If updated_asv_list file name is given, the input asv_list is 
#' completed by new asvs and asv_ids.
#' 
#' @param read_count Data frame or csv with columns: asv,sample,replicate,read_count
#' @param asv_list Data frame or csv with columns: file with asv - asv_id
#' Optional. It is used to homogenize asv_ids between different data sets.
#' @param updated_asv_list Character string: output file with 
#' the updated asv - asv_id pairs. Optional.
#' @param sep Field separator character in input and output csv files.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with asv_id column added to the input data frame.
#' @examples
#' \dontrun{
#' add_ids <- function(read_count_df, asv_list=asv_list)
#' }
#' @export
#' 
add_ids <- function(read_count, 
                    asv_list=NULL, 
                    updated_asv_list="", 
                    sep=",", 
                    quiet=T
                    ){
  
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  ### read earlier asv_id - asv pairs
  if(is.character(asv_list)){
    if(asv_list == ""){ # no filename
      asv_df <- data.frame("asv_id"=integer(),
                           "asv"=as.character())
    }else{  # read known asv
      asv_df <- read.csv(asv_list, header=T, sep=sep)
    }
  }else if (is.null(asv_list)){
    asv_df <- data.frame("asv_id"=integer(),
                         "asv"=as.character())
  }else{
    asv_df <- asv_list
  }
  t <- check_one_to_one_relationship(asv_df) # stop execution, if FALSE
  

  # list of unique asvs
  asv_uniq <- unique(read_count_df$asv)
  # list of unique asvs, not in the asv_df
  new_asvs <- asv_uniq[!asv_uniq %in% asv_df$asv]
  
  if(length(new_asvs) > 0){ #There are new ASVs
      
    if(nrow(asv_df)>0){
      max_id <- max(asv_df$asv_id)
    }else{
      max_id <- 0
    }
    new_ids <- seq(from =max_id+1, to = (max_id + length(new_asvs)), by=1)
    new_asvs_df <- data.frame(
      "asv_id" = new_ids, "asv"=new_asvs)
    # add new asvs to asv_df
    asv_df <- rbind(asv_df, new_asvs_df)
  }
  
  
  # add asv_id to read_count_df
  read_count_df <- left_join(read_count_df, asv_df, by="asv") %>%
    select(asv_id, sample, replicate, read_count, asv)
  
  # if asv_list should be updated, write it to a new file
  if(updated_asv_list != ""){
    write.table(asv_df, file=updated_asv_list, row.names = FALSE, sep=sep)
  }
  
  return(read_count_df)
}

#' Read sequences from fasta
#' 
#' Read sequences from a fasta file. 
#' 
#' Input fasta can be gzip compressed or uncompressed.
#' 
#' @param filename Character string: input fasta file name including full path.
#' @param dereplicate logical: If TRUE, return ASVs with read counts.
#' @returns Data frame with one read in each line (dereplicate==F), or a data frame with 
#' asv and read_count columns (dereplicate==T).
#' @examples
#' \dontrun{
#' read_df <- read_fasta_seq(filename="data/test.fasta", dereplicate=F)
#' asv_df <- read_fasta_seq(filename="data/test.fasta", dereplicate=T)
#' }
#' @export 
#' 

read_fasta_seq <- function(filename=filename, dereplicate=F){
  # can deal with sequences in multiple lines
  # only R-base
  # quicker than read.fasta from seqinr
  if(endsWith(filename, ".gz")){
    file_connection <- gzfile(filename, "rb")
  }else{
    file_connection <- file(filename, "r")
  }
  data <- readLines(file_connection, n = -1)
  close(file_connection)
  
  data <- gsub(" ", "_", data)
  data <- gsub(">[^ ]+", ">", data, fixed=F, perl=T)
  data <- do.call(paste, c(as.list(data), sep = ""))
  data <- as.data.frame(strsplit(data, ">"))
  colnames(data) <- c("read")
  data <- data %>%
    filter(!(read==""))
  
  data$read <-toupper(data$read)
  
  if(dereplicate){
    data <- data %>%
      group_by(read) %>%
      summarize(read_count = length(read)) %>%
      select(asv=read, read_count) %>%
      ungroup()
  }
  
  return(data)
}


#' Run swarm
#' 
#' Runs swarm [https://github.com/torognes/swarm](https://github.com/torognes/swarm) 
#' on input read_count data frame or csv file. 
#' Pools variants of the same cluster and sums read counts of the 
#' underlying ASVs.
#' 
#' Swarm can be run sample by sample (by_sample=T) or for the whole data set 
#' in one go.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv, plate, marker, sample, replicate (optional), read_count.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param swarm_path Character string: path to swarm executables. 
#' @param by_sample logical: run swarm separately for each sample.
#' @param num_threads Positive integer: Number of CPUs.
#' @param swarm_d Positive integer: d parameter for swarm.
#' Maximum number of differences allowed between two ASVs, 
#' meaning that two ASVs will be grouped if they have d (or less) differences.
#' @param fastidious logical: when working with d = 1, perform a second 
#' clustering pass to reduce the number of small clusters.
#' @param sep Field separator character in input and output csv files.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the same structure as the input, but ASVs of 
#' the same cluster pooled to one row.
#' @examples
#' \dontrun{
#' Swarm(read_count=read_count_df, swarm_path=swarm_path, num_threads=4, by_sample=T)
#' }
#' @export
#' 
Swarm <- function(read_count, 
                  outfile="", 
                  swarm_path="swarm", 
                  num_threads=1, 
                  swarm_d=1, 
                  fastidious=T, 
                  sep=",", 
                  by_sample=T, 
                  quiet=T
                  ){
  
   # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  if(by_sample){
    # make an empty output df, with the same columns and variable types as read_count_df
    out_df <- read_count_df %>%
      filter(asv=="")
    
    # get list of samples 
    sample_list <- unique(read_count_df$sample)
    
    # run swarm for each sample
    for(s in sample_list){
      if(!quiet){
        print(s)
      }
      
      # select occurrences for sample
      df_sample <- read_count_df %>%
        filter(sample==s)
      # run swarm
      df_sample <- run_swarm(df_sample, 
                             swarm_path=swarm_path, 
                             num_threads=num_threads, 
                             swarm_d=swarm_d, 
                             fastidious=fastidious,
                             quiet=quiet
                             )
      # add output of the sample to the total data frame
      out_df <- rbind(out_df, df_sample)
    }
  }else{ # run swarm for all samples together
    out_df <- run_swarm(read_count_df, 
                        swarm_path=swarm_path,
                        num_threads=num_threads,
                        swarm_d=swarm_d,
                        fastidious=fastidious, 
                        quiet=quiet
                        )
  }
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(out_df, file = outfile,  row.names = F, sep=sep)
  }
  return(out_df)
  
}

#' Run Swarm on one file
#' 
#' Runs swarm [https://github.com/torognes/swarm](https://github.com/torognes/swarm) 
#' on input read_count data frame or csv file. 
#' Pools variants of the same cluster and sums read counts of the 
#' underlying ASVs.
#' 
#' `run_swarm` runs swarm on the whole data set in one go, while 
#' `Swarm` can run swarm sample by sample or in one go.
#' 
#' @param read_count_df Data frame with the following variables: 
#' asv, plate, marker, sample, replicate (optional), read_count.
#' @param swarm_path Character string: path to swarm executables. 
#' @param num_threads Positive integer: Number of CPUs.
#' @param swarm_d Positive integer: d parameter for swarm.
#' Maximum number of differences allowed between two ASVs, 
#' meaning that two ASVs will be grouped if they have d (or less) differences.
#' @param fastidious logical: when working with d = 1, perform a second 
#' clustering pass to reduce the number of small clusters.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the same structure as the input, but ASVs of 
#' the same cluster pooled to one row.
#' @examples
#' \dontrun{
#' run_swarm(read_count=read_count_df, swarm_path=swarm_path, num_threads=4)
#' }
#' @export
#' 
run_swarm <- function(read_count_df, 
                      swarm_path="swarm", 
                      num_threads=1, 
                      swarm_d=1, 
                      fastidious=T,
                      quiet=T){
  
  tmp_dir <-paste('tmp_swarm_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  tmp_dir <- file.path(tempdir(), tmp_dir)
  check_dir(tmp_dir)
  
  
  ### make df with unique asv and read_count
  df_unique <- read_count_df %>%
    group_by(asv, asv_id) %>%
    summarize(sum_read_count = sum(read_count), .groups="drop_last") %>%
    ungroup() 

  ### make a fasta with dereplicated sequences  
  input_swarm <- file.path(tmp_dir, "swarm_input.fasta")
  writeLines(paste(">", df_unique$asv_id, "_", 
                   df_unique$sum_read_count, "\n", 
                   df_unique$asv, 
                   sep="" 
                   ), 
             input_swarm)
  
  df_unique <- df_unique %>%
    select(-sum_read_count)
  ### run swarm
  #  representatives <- file.path(tmp_dir, "representatives.fasta")
  clusters <- file.path(tmp_dir, "clusters.txt")
  swarm <- paste(swarm_path, 
                 " -d ",swarm_d,
                 " -t ", num_threads, 
                 " -o ", clusters, 
                 sep=""
                 )
  if(fastidious){
    swarm <- paste(swarm, "-f", sep=" ")
  }
  swarm <- paste(swarm, input_swarm, sep=" ")
  if(!quiet){
    print(swarm)
  }
  system(swarm, show.output.on.console = FALSE)

  
  ###
  # pool clusters in read_count_df
  ###
  # make a data frame with representative and clustered columns, 
  # where clustered has all swarm input sequences id, 
  # and  representative is the name of the cluster they belong to
#  print("Pooling ASVs to clusters")
  cluster_df <- read.table(clusters, fill =TRUE, strip.white=TRUE, header = FALSE)
  cluster_df <- data.frame(representative = rep(cluster_df$V1, 
                                                each = ncol(cluster_df)),
                           clustered = as.vector(t(cluster_df[,])))
  # delete line with no values in clustered
  cluster_df <- cluster_df %>%
    filter(clustered != "")
  # delete read counts from id
  cluster_df$representative <- sub("_[0-9]+", "", cluster_df$representative )
  cluster_df$representative <- as.numeric(cluster_df$representative)
  cluster_df$clustered <- sub("_[0-9]+", "", cluster_df$clustered )
  cluster_df$clustered <- as.numeric(cluster_df$clustered)
  # add representative asv
  cluster_df <- left_join(cluster_df, df_unique, by= c("representative" = "asv_id")) %>%
    select("representative_id"=representative, representative_asv=asv, "clustered_id"=clustered)

  # free space
  remove(df_unique)
  unlink(input_swarm)
  unlink(clusters)
  unlink(tmp_dir, recursive=TRUE)
  
#  print("Replace asv by representative sequences")
  # replace asv by representative sequences in read_count_df
  read_count_df <- left_join(read_count_df, cluster_df,  by= c("asv_id" = "clustered_id"))
  
  
  if("replicate" %in% colnames(read_count_df)){ # if replicates in input, keep replicates
    
  read_count_df <- read_count_df %>%
    select(-asv, -asv_id) %>%
    group_by(representative_asv, representative_id, sample, replicate) %>%
    summarize(read_count_cluster=sum(read_count), .groups="drop_last") %>%
    rename("asv" = representative_asv, 
           "read_count"=read_count_cluster, 
           "asv_id"=representative_id) %>%
    select(asv_id, sample, replicate, read_count, asv) %>%
    ungroup()
  }else{ # if no replicate column
    read_count_df <- read_count_df %>%
      select(-asv, -asv_id) %>%
      group_by(representative_asv, representative_id, sample) %>%
      summarize(read_count_cluster=sum(read_count), .groups="drop_last") %>%
      rename("asv" = representative_asv, 
             "read_count"=read_count_cluster, 
             "asv_id"=representative_id) %>%
      select(asv_id, sample, read_count, asv) %>%
      ungroup()
  }
  
  return(read_count_df)
}

#' Check one to one relationship
#' 
#' Check if there is a one to one relationship between 
#' unique ASVs and unique asv_ids in the input data frame.
#' 
#' The same asv - asv_id combination can appear more than once in the data frame.
#' 
#' @param df Data frame with the following variables: asv_id, asv 
#' (can have other columns as well).
#' @returns  TRUE if one to one relationship between asv - asv_id pairs. 
#' Otherwise stops the run with an error message.
#' @examples
#' \dontrun{
#' check_one_to_one_relationship(df=read_count_df)
#' }
#' @export
#' 
check_one_to_one_relationship <- function(df){
  
  # make unique asv-asv_id combinations
  df <- df %>%
    select(asv_id, asv) %>%
    distinct()
  
  # check if more then one asv per asv_id
  unique_asv_id <- df %>%
    group_by(asv_id) %>%
    summarize(count= length(asv)) %>%
    filter(count>1) %>%
    ungroup()
  if(nrow(unique_asv_id) > 0 ){
    print(unique_asv_id)
    stop("Some of the the asv_ids belong to multile asv")
  }
  
  # check if more then one asv_id per asv
  unique_asv <- df %>%
    group_by(asv) %>%
    summarize(count= length(asv_id)) %>%
    filter(count>1) %>%
    ungroup()
  if(nrow(unique_asv) > 0 ){
    print(unique_asv)
    stop("Some of the asv has multiple asv_ids")
  }
  
  return(TRUE)
}

#' Update ASV list
#' 
#' Pools unique asv - asv_id pairs in the input data frames or csv files. 
#' Check if there is no conflict within an between the input data.
#' Returns a complete and unique list of asv - asv_id pairs.
#' 
#' 
#' @param asv_list1 Data frame or csv file with columns:
#' asv_id,asv
#' @param asv_list2 Data frame or csv file with columns:
#' asv_id,asv
#' @param outfile Character string: output csv file name. 
#' @param sep Field separator character in input and output csv files.
#' @param return_df logical; if TRUE returns a data frame
#' @returns Updated data frame or csv file with all unique asv_id - asv pairs 
#' in the input data frame and csv file.
#' If there is a conflict within or between the input data stops with 
#' an error message.
#' @examples
#' \dontrun{
#' UpdateASVlist(asv_list1=read_count_df, 
#' asv_list2="data/asv_list.csv", 
#' outfile="out/updated_asv_list.csv"
#' )
#' }
#' @export
#' 
UpdateASVlist <- function(asv_list1, asv_list2, outfile, sep=",", return_df=FALSE){
  
  if(is.character(asv_list1)){
    df1 <- read.csv(asv_list1, header=T, sep=sep)
  }else{
    df1 <- asv_list1
  }
  # select columns and make unique
  df1 <- df1 %>%
    select(asv_id, asv) %>%
    distinct()
  t <- check_one_to_one_relationship(df1)
  
  if(is.character(asv_list2)){
    df2 <- read.csv(asv_list2, header=T, sep=sep)
  }else{
    df2 <- asv_list2
  }  
  # select columns and make unique
  df2 <- df2 %>%
    select(asv_id, asv) %>%
    distinct()
  t <- check_one_to_one_relationship(df2)
  
  # pool 
  df1 <- rbind(df1, df2) %>%
    distinct() %>%
    arrange(asv_id)
  t <- check_one_to_one_relationship(df1)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(df1, file=outfile, row.names = FALSE, sep=sep)
  }
  if(return_df){
    return(df1)
  }
}

#' LFNglobalReadCount
#' 
#' Eliminate ASVs with less than cutoff reads in the data set.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param cutoff Positive integer: minimum number of reads for an ASV in the 
#' whole data set. Bellow this cutoff, the ASV is deleted from the data set.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written. 
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- LFNglobalReadCount(read_count_df, cutoff=2)
#' }
#' @export
#' 
LFNglobalReadCount <- function (read_count, cutoff=10, outfile="", sep=",") {
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count_all=sum(read_count)) %>%
    filter(read_count_all > cutoff) %>%
    ungroup()
  read_count_df <- filter(read_count_df, (asv %in% df$asv))
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' LFNreadCount
#' 
#' Eliminate occurrences (presence of the ASV in a sample-replicate)
#' with less than cutoff reads.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param cutoff Positive integer: minimum number of reads for an occurrence.
#' Bellow this cutoff, the occurrence (presence of the ASV in a sample-replicate)
#' is deleted.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written. 
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count <- LFNreadCount(read_count_df, cutoff=20)
#' }
#' @export
#' 
LFNreadCount <- function (read_count, cutoff=10, outfile="", sep=",") {
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  read_count_df <- filter(read_count_df,  (read_count >= cutoff))
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' LFNsampleReplicate
#' 
#' Eliminate occurrences (presence of the ASV in a sample-replicate)
#' where the `read_count / sum(read_count of the sample-replicate)` is less than 
#' cutoff.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param cutoff Real (0-1): minimum proportion of the read count of an
#' occurrence within the reads of its sample-replicate. Bellow this cutoff
#'  the occurrence is deleted.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written. 
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- LFNsampleReplicate(read_count_df, cutoff=0.005)
#' }
#' @export
LFNsampleReplicate <- function (read_count, cutoff=0.001, outfile="", sep=",") {
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  sum_by_column_df <- read_count_df %>%
    group_by(sample,replicate) %>%
    summarize(sr_sum = sum(read_count), .groups="drop_last") %>%
    ungroup()
  
  read_count_df <- left_join(read_count_df, sum_by_column_df, 
                             by=c("sample","replicate")) %>%
    filter(read_count/sr_sum >= cutoff) %>%
    select(-sr_sum) %>%
    ungroup()
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' LFNvariant
#' 
#' If by_replicate is FALSE: Eliminate occurrences where the 
#' (read_count/read_count of the asv in the data set) is less than cutoff.
#' If by_replicate is TRUE: Eliminate occurrences where the 
#' (read_count/read_count of the asv in its replicate) is less than cutoff.
#' 
#' Issues a warning if the total read count of an ASV has been reduced 
#' bellow min_read_count_prop, since it can indicate a to high cutoff value.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param cutoff Real (0-1): minimum proportion of the read count of
#'  an occurrence within all reads of the asv or asv-replicate. Bellow this cutoff
#'  the occurrence is deleted.
#' @param by_replicate logical: Compare read count of the occurrence to the 
#' read counts of the ASV-replicate.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param min_read_count_prop Real (0-1): If the proportion of the read count 
#' of a variant in the output compared to the input is less then 
#' min_read_count_prop, prints out a warning, since it suggest a 
#' to high cutoff value
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- LFNvariant(read_count_df, cutoff=0.005, min_read_count_prop=0.8)
#' }
#' @export
#' 
LFNvariant <- function(read_count, 
                       cutoff=0.001, 
                       by_replicate=FALSE, 
                       outfile="", 
                       sep=",", 
                       min_read_count_prop=0.7
                       ){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # Get the list of asv and the number of samples-replicates present and the total number of reads;
  # This will be compared to values after filtering
  asvs <- read_count_df %>%
    group_by(asv_id) %>%
    summarize("sample_count_input" = length(sample), "read_count_input"=sum(read_count)) %>%
    filter(read_count_input > 10) %>%
    ungroup()
  
  
  if(by_replicate){
    sum_by_asv <- read_count_df %>%
      group_by(asv,replicate) %>%
      summarize(asv_sum = sum(read_count), .groups="drop_last") %>%
      ungroup()
    read_count_df <- left_join(read_count_df, sum_by_asv, by=c("asv", "replicate")) %>%
      filter(read_count/asv_sum >= cutoff) %>%
      select(-asv_sum)
    
  }else{
    sum_by_asv <- read_count_df %>%
      group_by(asv) %>%
      summarize(asv_sum = sum(read_count)) %>%
      ungroup()
    read_count_df <- left_join(read_count_df, sum_by_asv, by="asv") %>%
      filter(read_count/asv_sum >= cutoff) %>%
      select(-asv_sum)
  }
  
  ###
  # Check if filter do not eliminate occurrences with relatively high readcount 
  ###
  asvs_output <- read_count_df %>%
    group_by(asv_id) %>%
    summarize("sample_count_output" = length(sample), 
              "read_count_output"=sum(read_count)) %>%
    ungroup()
  # join sample and read counts before and after filtering
  asvs <- left_join(asvs, asvs_output, by="asv_id")
  asvs$sample_prop <- asvs$sample_count_output / asvs$sample_count_input
  asvs$read_count_prop <- asvs$read_count_output / asvs$read_count_input
  asvs <- asvs %>%
    filter(read_count_prop<min_read_count_prop) %>%
    arrange(read_count_prop, sample_prop) %>%
    select("asv_id", "read_count_input", "read_count_output", "read_count_prop",
           "sample_count_input", "sample_count_output", "sample_prop")
  
  if(nrow(asvs > 0)){
    print("WARNING: The following ASVs have lost a high proportion of their 
          reads during this filtering step. 
          The cutoff value of LFNvariant function might need to be reduced.")
    print(asvs)
  }
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' PoolFilters
#' 
#' Pool all count_read_df data frames, and keep only occurrences 
#' present in all filters.
#' 
#' @param ... Data frames with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- PoolFilters(read_count_df1, read_count_df2)
#' }
#' @export
#' 
PoolFilters <- function(... , outfile="", sep=","){
  df_list <- list(...)
  merged <-  df_list[[1]]
  for(i in 2:length(df_list)){
    merged <- inner_join(merged, df_list[[i]])
  }
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(merged, file = outfile,  row.names = F, sep=sep)
  }
  return(merged)
}

#' FilterMinReplicate
#' 
#' Filter out all occurrences where the asv in not present in at least 
#' `cutoff` number of replicates of the sample.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param cutoff Positive integer: minimum number of replicates.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- FilterMinReplicate(read_count_df, cutoff=3)
#' }
#' @export
#'
FilterMinReplicate <- function(read_count, cutoff=2, outfile="", sep=","){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  # read_count_df <- df
  # add a temporary column with asv and sample concatenated
  read_count_df$tmp <- paste(read_count_df$asv,  read_count_df$sample, sep="-")
  # make a df_tmp containing the number of replicates for each asv-sample combination
  df_tmp <- read_count_df  %>%
    group_by(tmp) %>%
    summarize(repl_number=length(tmp))  %>%
    filter(repl_number >= cutoff) %>%
    ungroup()
  # keep only asv-sample if present at least in min_replicate_number replicates
  read_count_df <- filter(read_count_df, (read_count_df$tmp %in% df_tmp$tmp))
  read_count_df$tmp <- NULL
  
  if(outfile !=""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' FilterIndel
#' 
#' Filter out all ASVs, if the modulo 3 of their length is not the same as 
#' that of the majority of the ASVs.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <-FilterIndel(read_count_df)
#' }
#' @export
#' 
FilterIndel <- function(read_count, outfile="", sep=","){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # add a column with the modulo 3 of the length of the asvs
  read_count_df$mod3 <- nchar(read_count_df$asv) %% 3
  # make a tibble with modulo3 of the length of the ASVs and their count 
  # ASVs are counted as many times as they occur, so the most frequent ASV have higher weight
  # read_counts are not taken into account
  tmp <- read_count_df %>%
    group_by(mod3) %>%
    summarize(length_modulo=length(mod3)) %>%
    arrange(desc(length_modulo)) %>%
    ungroup()
  # get the modulo 3 the most frequent
  my_modulo3 <- as.integer(tmp[1,"mod3"])
  # select only the lines with asv length compatible with the most frequent modulo3
  read_count_df <- read_count_df %>%
    filter(mod3 == my_modulo3)
  
  # delete the temporary column
  read_count_df$mod3 <- NULL
  
  if(outfile !=""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Codon stops from genetic code
#' 
#' Returns a vector of codon stops corresponding the genetic 
#' code number in [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes).
#'  
#' @param genetic_code Integer: genetic code number from [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes).
#' @returns Vector of codon stops.
#' @examples
#' \dontrun{
#' codon_stops_from_genetic_code(genetic_code=1)
#' }
#' @export
#' 
codon_stops_from_genetic_code <- function(genetic_code=5){
  if(genetic_code == 1){
    return(c("TAA","TAG","TGA"))
  }
  else if(genetic_code == 2){
    return(c("TAA","TAG","AGA", "AGG"))
  }
  else if(genetic_code == 3 || genetic_code == 4 || genetic_code == 5 || 
          genetic_code == 9 || genetic_code == 10 || genetic_code == 13 || 
          genetic_code == 21 || genetic_code == 24 || genetic_code == 25 || 
          genetic_code == 31){
    return(c("TAA","TAG"))
  }
  else if(genetic_code == 6 || genetic_code == 14 || genetic_code == 33){
    return(c("TAG"))
  }
  else if(genetic_code == 16){
    return(c("TAA","TGA"))
  }
  else if(genetic_code == 11 || genetic_code == 12 || 
          genetic_code == 26 || genetic_code == 28 ){
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
#' Filter out all ASVs, if there is a codon stop in all three reading 
#' frames of the direct strand.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param genetic_code Positive integer: genetic code number from [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes).
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- FilterCodonStop(read_count_df, genetic_code=5)
#' }
#' @export
#' 
FilterCodonStop <- function(read_count, outfile="", genetic_code=5, sep=","){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  codon_stops <- codon_stops_from_genetic_code(genetic_code=genetic_code)
  if(length(codon_stops) == 0){
    print("WARNING: The Genetic Code Number provided does not correspond to any code in 
          https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes\n
          The FilterCodonStop step is skipped")
    return(read_count_df)
  }
  
  # get unique list of variants to a df
  unique_asv_df <- data.frame(asv = unique(read_count_df$asv))
  
  # define a function to transform each sequence to a list of nucleotides
  seq_to_list_of_nt <- function(seq){
    return(strsplit(seq, "")[[1]])
  }
  # define a function to group a list of nucleotides to pieces of three (codons)
  # starting at start_pos
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
    # check if there is at least one codon stop among the codons; 
    # return 1 if CodonStop, 0 if not
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
  
  if(outfile !=""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Write fasta file
#' 
#' Write a vector of sequeces to a fasta file using 
#' either the sequences as ID (seq_as_id==TRUE), or arbitrary sequence IDs 
#' (seq_as_id==FALSE)
#'  
#' @param sequences Vector of sequences.
#' @param filename Character string: name of the output fasta file.
#' @param seq_as_id logical: Use sequences as seqID.
#' @returns fasta file
#' @examples
#' \dontrun{
#' write_fasta(sequences=read_count_df$asv, filename="out/seq.fasta", seq_as_id=T)
#' }
#' @export
#' 
write_fasta <- function(sequences, filename, seq_as_id=F) {
  # Open the file for writing
  file <- file(filename, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(sequences)) {
    if(seq_as_id){
      header <- paste0(">", sequences[[i]])
    }else{
      header <- paste0(">", i)
    }
    writeLines(c(header, sequences[[i]], ""), file)
  }
  # Close the file
  close(file)
}

#' Flag PCR error using vsearch
#' 
#' Identify potential PCRerrors: ASVs very similar (`max_mismatch`) to another 
#' more frequent ASV (`pcr_error_var_prop`) in the input data frame.
#' Adds a column to the input data frame, with 1 if sequence is a 
#' probable PCR error, and 0 otherwise.
#'  
#' @param unique_asv_df Data frame with the following variables: 
#' asv, read_count. ASVs must be unique.
#' @param pcr_error_var_prop Real (0-1): if the proportion of read counts 
#' of two similar ASVs is bellow `pcr_error_var_prop`, the less abundant 
#' is flagged as a PCR error.
#' @param max_mismatch Positive integer: maximum number of mismatches 
#' (gaps included) to consider two ASVs as similar.
#' @param vsearch_path Character string: path to vsearch executables.
#' @returns Input data frame completed by a PCRerror column 
#' (1 if potential PCR error, 0 otherwise).
#' @examples
#' \dontrun{
#' unique_asv_df <- read_count_df %>%
#'   group_by(asv) %>%
#'   summarize(read_count = sum(read_count))
#' unique_asv_df_flagged <- flagPCRerror_vsearch(unique_asv_df, 
#'      vsearch_path=vsearch_path, 
#'      pcr_error_var_prop=0,2, 
#'      max_mismatch=2
#'      )
#' }
#' @export
#' 
flagPCRerror_vsearch <- function(unique_asv_df,
                                 vsearch_path="vsearch", 
                                 pcr_error_var_prop=0.1, 
                                 max_mismatch=1
                                 ){
  
  # no ASV in the unique_asv_df => return a dataframe with 0 for all ASVs in PCRerror column
  if(length(unique_asv_df$asv) == 0){ 
    unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
    return(unique_asv_df)
  }
  
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste('tmp_PCRerror_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  outdir_tmp <- file.path(tempdir(), outdir_tmp)
  check_dir(outdir_tmp)
  
  # make fasta file with unique reads; use sequences as ids
  fas <- file.path(outdir_tmp, 'unique.fas')
  write_fasta(unique_asv_df$asv, fas, seq_as_id=T)
  # vsearch --usearch_global to find highly similar sequence pairs
  vsearch_out <- file.path(outdir_tmp, 'unique_vsearch_out.out')
  vsearch <- paste(vsearch_path, " --usearch_global ", fas, 
                   " --db ", fas, 
                   " --userout ",  vsearch_out,
                   " --quiet --iddef 1 --self --id 0.90 ",
                   " --maxaccepts 0 --maxrejects 0 ",
                   " --userfields ",'"query+target+ids+aln"', 
                  sep="")
  #https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/system
  system(vsearch)
  
  # no vsearch hit => return unique_asv_df completed with a PCRerror, with 0 for all ASVs
  if(!file.exists(vsearch_out) || file.size(vsearch_out) == 0){
    unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
    # Delete the temp directory
    unlink(outdir_tmp, recursive = TRUE)
    return(unique_asv_df)
  }
  
  # read vsearch results
  results_vsearch<- read.csv(vsearch_out, header = FALSE, sep="\t")
  colnames(results_vsearch) <- c("query","target","nb_ids","aln")
  # none of the values easily outputted by vsearch take into the external gaps as a diff 
  # => correct this, based on the alnlen and the number of identities
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
  results_vsearch$PCRerror_target <- 
    ((results_vsearch$qread_count * pcr_error_var_prop) >= results_vsearch$tread_count)
  # keep only one column (tasv) with unique ASVs, that were flagged as PCRerror
  results_vsearch <- results_vsearch %>%
    filter(PCRerror_target==TRUE) %>%
    group_by(tasv) %>%
    select(tasv) %>%
    ungroup()
  # complete unique_asv_df with a PCRerror column
  unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
  unique_asv_df$PCRerror[unique_asv_df$asv %in% results_vsearch$tasv] <- 1
  
  # Delete the temp directory
  unlink(outdir_tmp, recursive = TRUE)
  
  return(unique_asv_df)
}

#' Filter PCR error
#' 
#' Filter out an ASVs if it is very similar (`max_mismatch`) to another 
#' more frequent ASV (`pcr_error_var_prop`).
#' 
#' The whole data set can be analyzed at once (`by_sample=F`) 
#' or sample by sample.
#'  
#' @param read_count Data frame or csv file with the following variables:  
#' asv_id, sample, replicate, read_count, asv.
#' @param pcr_error_var_prop Real (0-1): if the proportion of read_counts of 
#' two similar ASVs is less or equal to `pcr_error_var_prop`, 
#' the less abundant is flagged as a PCR error.
#' @param max_mismatch Positive integer: maximum number of mismatches 
#' (gaps included) to consider two ASVs as similar
#' @param by_sample logical: if TRUE ASVs are flagged as an PCR error 
#' separately for each sample.
#' @param sample_prop Real (0-1): if by_sample=TRUE, the ASV must be 
#' flagged as a PCRerror in `sample_prop` of the samples to be eliminated.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param vsearch_path Character string: path to vsearch executables. 
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- FilterPCRerror(read_count_df, 
#'      vsearch_path=vsearch_path, 
#'      pcr_error_var_prop=0.2, 
#'      max_mismatch=2, 
#'      by_sample=T, 
#'      sample_prop=0.8
#'      )
#' filtered_read_count_df <- FilterPCRerror(read_count_df, 
#'     vsearch_path=vsearch_path, 
#'     pcr_error_var_prop=0.2,
#'     max_mismatch=2, 
#'     by_sample=F
#'     )
#' }
#' @export
#' 
FilterPCRerror <- function(read_count,
                           outfile="", 
                           vsearch_path="vsearch", 
                           pcr_error_var_prop=0.1,
                           max_mismatch=1, 
                           by_sample=T, 
                           sample_prop=0.8, 
                           sep=","
                           ){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  # get unique list of ASVs with their total read_count in the run
  unique_asv_df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count = sum(read_count)) %>%
    arrange(desc(read_count)) %>%
    ungroup()
  
  if(by_sample){ # sample by sample
    sample_list <- unique(read_count_df$sample)
    # loop over samples
    for(sample_loc in sample_list){
      # get unique list of ASVs with their total read_count in the sample
      unique_asv_df_sample <- read_count_df %>%
        filter(sample == sample_loc) %>%
        group_by(asv)%>%
        summarize(read_count = sum(read_count))%>%
        arrange(desc(read_count)) %>%
        ungroup()
      
      # flag PCR errors; 
      # add one column to unique_asv_df for each sample with 1 if ASV is flagged in the sample, 
      # 0 otherwise
      unique_asv_df_sample <- flagPCRerror_vsearch(unique_asv_df_sample, 
                                                   vsearch_path=vsearch_path, 
                                                   pcr_error_var_prop=pcr_error_var_prop, 
                                                   max_mismatch=max_mismatch
                                                   )
      
      # remove read_count column
      unique_asv_df_sample$read_count <- NULL
      # add a column for for each sample to unique_asv_df, 
      # with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df <- left_join(unique_asv_df, unique_asv_df_sample, by = "asv")
    }
  }
  else{ # whole dataset
    # add a column to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
    unique_asv_df <- flagPCRerror_vsearch(unique_asv_df, 
                                          vsearch_path=vsearch_path, 
                                          pcr_error_var_prop=pcr_error_var_prop, 
                                          max_mismatch=max_mismatch
                                          )
  }
  
  # count the number of times each ASV has been flagged and when it has not. 
  # Ignore NA, when the ASV is not present in the sample
  unique_asv_df$yes <- rowSums(unique_asv_df[3:ncol(unique_asv_df)] == 1, na.rm = TRUE)
  unique_asv_df$no <- rowSums(unique_asv_df[3:(ncol(unique_asv_df)-1)] == 0, na.rm = TRUE)
  # keep only ASVs, 
  # that are not flagged in sample_prop proportion of the samples where they are present  
  unique_asv_df <- unique_asv_df %>%
    filter(yes/(yes+no) >= sample_prop)
  
  # eliminate potential PCRerrors from read_count_df
  read_count_df <- read_count_df %>%
    filter(!asv %in% unique_asv_df$asv)
  
  if(outfile !=""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}


#' Flag chimeric sequences
#' 
#' Select chimeras in a data frame of unique ASVs.
#' Add a column to the input data frame with 1 if the ASV is a probable 
#' chimera and 0 otherwise.
#'  
#' @param unique_asv_df Data frame with the following variables: 
#' asv, read_count; ASVs must be unique.
#' @param abskew Positive integer: a chimera must be at least `abskew` 
#' times less frequent that the parental ASVs.
#' @param vsearch_path Character string: path to vsearch executables. 
#' @returns
#' The input data frame completed by chimera `column`. 
#' 1 if potential chimera, 0 otherwise.
#' @examples
#' \dontrun{
#' unique_asv_df <- read_count_df %>%
#'   group_by(asv) %>%
#'   summarize(read_count = sum(read_count))
#' flagChimera(unique_asv_df, vsearch_path=vsearch_path, abskew=2)
#' }
#' @export
#'
flagChimera <- function(unique_asv_df, vsearch_path="vsearch", abskew=2){

  # no ASV in the unique_asv_df => return a data frame with 0 for all ASVs in Chimera column
  if(length(unique_asv_df$asv) == 0){ 
    unique_asv_df$chimera <- rep(0, length(unique_asv_df$asv))
    return(unique_asv_df)
  }
  
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste('tmp_FilterChimera_', 
                      trunc(as.numeric(Sys.time())), 
                      sample(1:100, 1), 
                      sep=''
                      )
  outdir_tmp <- file.path(tempdir(), outdir_tmp)
  check_dir(outdir_tmp)
  
  # make fasta file with unique reads; use sequences as ids
  fas <- file.path(outdir_tmp, 'unique.fas')
  # Open the file for writing
  file <- file(fas, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(unique_asv_df$asv)) {
    header <- paste0(">",unique_asv_df$asv[i], ";size=", unique_asv_df$read_count[i], sep="")
    writeLines(c(header, unique_asv_df$asv[i], ""), file)
  }
  close(file)
  
  # vsearch --usearch_global to find highly similar sequence pairs
  vsearch_out <- file.path(outdir_tmp, 'unique_vsearch_out.out')
  vsearch <- paste(vsearch_path, " --uchime3_denovo ", 
                   fas, " --quiet --abskew ", abskew ,
                   " --uchimeout  ", vsearch_out, 
                   sep=""
                   )
  system(vsearch)
  
  # no vsearch hit => return unique_asv_df completed with a PCRerror, with 0 for all ASVs
  if(!file.exists(vsearch_out) || file.size(vsearch_out) == 0){
    unique_asv_df$chimera <- rep(0, length(unique_asv_df$asv))
    # Delete the temp directory
    unlink(outdir_tmp, recursive = TRUE)
    return(unique_asv_df)
  }
  
  # read vsearch results
  results_vsearch<- read.csv(vsearch_out, header = FALSE, sep="\t")
  # keep only pertinent columns
  results_vsearch <- select(results_vsearch, c(2, ncol(results_vsearch)))
  colnames(results_vsearch) <- c("asv", "chimera")
  results_vsearch$asv <- gsub(";size=[0-9]+", "", results_vsearch$asv)
  # keep only chimeras
  results_vsearch <- results_vsearch %>%
    filter(chimera == "Y")
  
  # complete unique_asv_df with chimara info
  unique_asv_df$chimera <- rep(0, length(unique_asv_df$asv))
  unique_asv_df$chimera[unique_asv_df$asv %in% results_vsearch$asv] <- 1
  
  # Delete the temp directory
  unlink(outdir_tmp, recursive = TRUE)
  return(unique_asv_df)
}

#' FilterChimera
#' 
#' Filter out Chimeras.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param abskew Positive integer: a chimera must be at least `abskew` 
#' times less frequent that the parental ASVs.
#' @param by_sample logical: ASVs are flagged as chimera separately 
#' for each sample.
#' @param sample_prop logical: if by_sample=TRUE, the ASV deleted if 
#' they are flagged as chimera in at least `sample_prop` of the 
#' samples among the sample they are present.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param vsearch_path Character string: path to vsearch executables. 
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- FilterChimera(read_count_df, 
#'      vsearch_path=vsearch_path, 
#'      by_sample=T, 
#'      sample_prop=0.7, 
#'      abskew=4
#'      )
#' filtered_read_count_df <- FilterChimera(read_count_df,
#'      vsearch_path=vsearch_path, 
#'      by_sample=F, 
#'      abskew=4
#'      )
#' }
#' @export
#' 
FilterChimera <- function(read_count, 
                          outfile="", 
                          vsearch_path="vsearch", 
                          by_sample=T, 
                          sample_prop=0.8, 
                          abskew=2, 
                          sep=","
                          ){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  # get unique list of ASVs with their total read_count in the run
  unique_asv_df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count = sum(read_count)) %>%
    arrange(desc(read_count)) %>%
    ungroup()
  
  if(by_sample){ # sample by sample
    sample_list <- unique(read_count_df$sample)
    # loop over samples
    for(sample_loc in sample_list){
      # get unique list of ASVs with their total read_count in the sample
      unique_asv_df_sample <- read_count_df %>%
        filter(sample == sample_loc) %>%
        group_by(asv)%>%
        summarize(read_count = sum(read_count))%>%
        arrange(desc(read_count)) %>%
        ungroup()
      
      # flag chimeras; 
      # add one column to unique_asv_df for each sample with 1 if ASV is flagged in the sample, 
      # 0 otherwise
      unique_asv_df_sample <- flagChimera(unique_asv_df_sample, 
                                          vsearch_path=vsearch_path, 
                                          abskew=abskew
                                          )
      
      # remove read_count column
      unique_asv_df_sample <- select(unique_asv_df_sample, -c("read_count"))
      # add a column for each sample to unique_asv_df, with 1 if ASV is flagged in the sample,
      # 0 otherwise
      unique_asv_df <- left_join(unique_asv_df, unique_asv_df_sample, by = "asv")
    }
  }else{ # whole dataset
    # add a column to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
    unique_asv_df <- flagChimera(unique_asv_df, vsearch_path=vsearch_path, abskew=abskew)
  }
  
  # count the number of times each ASV has been flagged and when it has not. 
  # Ignore NA, when the ASV is not present in the sample
  unique_asv_df$yes <- rowSums(unique_asv_df[3:ncol(unique_asv_df)] == 1, na.rm = TRUE)
  unique_asv_df$no <- rowSums(unique_asv_df[3:(ncol(unique_asv_df)-1)] == 0, na.rm = TRUE)
  # keep only ASVs, that are flagged in sample_prop proportion 
  # of the samples where they are present  
  unique_asv_df <- unique_asv_df %>%
    filter(yes/(yes+no) >= sample_prop)
  
  # eliminate potential Chimeras from read_count_df
  read_count_df <- read_count_df %>%
    filter(!asv %in% unique_asv_df$asv)
  
  if(outfile !=""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Calculate Renkonen Distance
#' 
#' Calculate renkonen distance between two data frames with ASVs and read counts.
#'  
#' @param df1,df2 Data frames with asv and read_count columns.
#' @returns
#' Renkonen distance (Real; 0-1) between the two date frames.
#' @examples
#' \dontrun{
#' df1 <- read_count_df %>%
#'   filter(sample=="tpos1") %>%
#'   group_by(asv) %>%
#'   summarize(read_count = sum(read_count))
#' df2 <- read_count_df %>%
#'   filter(sample=="tnegtag") %>%
#'   group_by(asv) %>%
#'   summarize(read_count = sum(read_count))
#' calculate_renkonen_dist(df1, df2)
#' }
#' @export
#'
calculate_renkonen_dist <- function(df1, df2){
  df1 <- df1 %>%
    select(asv, "read_count1"=read_count)
  df2 <- df2 %>%
    select(asv, "read_count2"=read_count)
  
  df <- full_join(df1, df2, by="asv")
  # replace NA by 0
  df <- df %>%
    mutate(read_count1 = ifelse(is.na(read_count1), 0, read_count1)) %>%
    mutate(read_count2 = ifelse(is.na(read_count2), 0, read_count2))
  # calculate  (number of reads for variant x in replicate i) / (number of reads in replicate i)
  df$read_count1 <- df$read_count1/ sum(df$read_count1)
  df$read_count2 <- df$read_count2/ sum(df$read_count2)
  # minimum of the above proportion between the 2 replicates
  df$min <- pmin(df$read_count1, df$read_count2)
  rdist <- 1- sum(df$min)
}

#' Make Renkonen Distances
#' 
#' Calculate the Renkonen distance between pairs of sample-replicates.
#' 
#' @param read_count_df Data frame with asv, sample, replicate, and read_count columns.
#' @param compare_all logical: if TRUE calculate the Renkonen distance among all pairs 
#' of sample-replicates. Only between replicates of the same sample otherwise.
#' @returns Data frame with the following columns: 
#' sample1,sample2,replicate1,replicate2,renkonen_d,sample_comp
#' (within, if sample1 and sample2 are identical, between otherwise).
#' @examples
#' \dontrun{
#' MakeRenkonenDistances(read_count_df, compare_all=FALSE)
#' }
#' @export
#' 
MakeRenkonenDistances <- function(read_count_df, compare_all=FALSE){
  
  df <- read_count_df %>%
    select(asv, sample, replicate, read_count)
  df$sr <- paste(df$sample, df$replicate, sep="-")
  # list of samples
  sample_replicate_df <- df %>%
    select(sample, replicate, sr) %>%
    unique()
  
  # make empty data frame
  renkonen_df <- data.frame("sample1" = character(),
                            "sample2" = character(),
                            "sr1"  = character(),
                            "sr2"  = character(),
                            "renkonen_d" = numeric())
  
  # loop over all pairs of sample-replicates within sample
  for(i in 1:(nrow(sample_replicate_df)-1)){
    sri <- sample_replicate_df$sr[i]
    sampi <- sample_replicate_df$sample[i]
    repli <- sample_replicate_df$replicate[i]
    
    for(j in ((i+1):nrow(sample_replicate_df))){
      srj <- sample_replicate_df$sr[j]
      sampj <- sample_replicate_df$sample[j]
      replj <- sample_replicate_df$replicate[j]
      
      
      if(sampi == sampj){ # same sample
        dfi <- filter(df, sr == sri)
        dfj <- filter(df, sr == srj)
        rdist <- calculate_renkonen_dist(dfi, dfj)
        # add line to renkonen_df
        new_line <- data.frame(sample1 = sampi, 
                               sample2 = sampj, 
                               replicate1 = repli, 
                               replicate2 = replj, 
                               renkonen_d = rdist
                               )
        renkonen_df <- rbind(renkonen_df, new_line)
      }else{ # different sample
        if(compare_all){ # calculate only if within and between sample comparison is necessary
          dfi <- filter(df, sr == sri)
          dfj <- filter(df, sr == srj)
          rdist <- calculate_renkonen_dist(dfi, dfj)
          # add line to renkonen_df
          new_line <- data.frame(sample1 = sampi, 
                                 sample2 = sampj, 
                                 replicate1 = repli, 
                                 replicate2 = replj, 
                                 renkonen_d = rdist
                                 )
          renkonen_df <- rbind(renkonen_df, new_line)
        }
      }
    }
  }
  return(renkonen_df)
}

#' FilterRenkonen
#' 
#' Filter out all replicates that have renkonen distances above `cutoff` 
#' to most other replicates of the same sample.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param cutoff Real (0-1):  Filter out all replicates that have 
#' renkonen distances above cutoff to most other replicates of the same sample.
#' @param renkonen_distance_quantile Real (0-1): if cutoff value is not given, 
#' use the `renkonen_distance_quantile` to determine cutoff value 
#' (e.g. with 0.9 as `renkonen_distance_quantile`, 90% of the distances are bellow cutoff)
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- FilterRenkonen(read_count_df, cutoff = 0.6)
#' filtered_read_count_df <- FilterRenkonen(read_count_df, renkonen_distance_quantile=0.9)
#' }
#' @export
#' 
FilterRenkonen <- function(read_count, 
                           outfile="",
                           cutoff = NA, 
                           renkonen_distance_quantile=0.9,
                           sep=","
                           ){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # calculate Renkonen distances between all pairs of replicates of within sample
  renkonen_df <- MakeRenkonenDistances(read_count_df, compare_all=FALSE) %>%
    select("sample" = sample1, replicate1, replicate2, renkonen_d) %>%
    arrange(renkonen_d)
  
  # determine the cut off renkonen distance; values > cutoff are considered as high
  if(is.na(cutoff)){  
    last_row <- floor(length(renkonen_df$renkonen_d) * renkonen_distance_quantile)
    cutoff <- renkonen_df$renkonen_d[last_row]
  }
#  msg <- paste("The cutoff value for Renkonen distances is ", cutoff)
#  print(msg)
  # get list of samples
  sample_list <- unique(renkonen_df$sample)
  # filter out replicates sample by sample
  for(samp in sample_list){
    sample_df <- renkonen_df %>%
      filter(sample == samp)
    # make a df with the replicate columns exchanged
    sample_tmp <- data.frame("sample" = sample_df$sample,
                             "replicate1"  = sample_df$replicate2,
                             "replicate2"  = sample_df$replicate1,
                             "renkonen_d" = sample_df$renkonen_d)
    # complete sample_df to include distances between repl X and Y and also Y and X
    sample_df <- rbind(sample_df, sample_tmp)
    # get unique list of replicates
    replicate_list <-  unique(sample_df$replicate1)
    # the minimum number of distances to be bellow cutoff, to keep the replicate
    min_number_of_distances_bellow_cutoff <- (length(replicate_list) -1) / 2
    
    # keep only distances above cutoff in sample_df
    sample_df <- sample_df %>%
      filter(renkonen_d > cutoff)
    
    # count for each replicate the number of distances above ctuoff 
    sample_df <- sample_df %>%
      group_by(replicate1) %>%
      summarize(n_dist=length(renkonen_d)) %>%
      filter(n_dist > min_number_of_distances_bellow_cutoff) %>%
      ungroup()
    
    # eliminate replicates with too many distances above cutoff
    read_count_df <- read_count_df %>%
      filter(!(sample == samp & replicate %in% sample_df$replicate1))
  }
  
  if(outfile !=""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)  
}

#' Pool replicates sample by sample
#' 
#' Take the mean non-zero read counts over replicates for each sample and asv.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param digits Positive integer: round the mean read counts to digits.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @returns Data frame with the following columns: 
#' asv, plate, marker, sample, read_count (over replicates)
#' @examples
#' \dontrun{
#' PoolReplicates(read_count_df)
#' }
#' @export 
#'
PoolReplicates <- function(read_count, digits=0, outfile="", sep=","){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  read_count_samples_df <- read_count_df %>%
    group_by(asv_id,sample,asv) %>%
    summarize(read_count = mean(read_count), .groups="drop_last") %>%
    select(asv_id, sample, read_count, asv) %>%
    ungroup()
  
  read_count_samples_df$read_count <- round(read_count_samples_df$read_count, digits =digits)
  
  if(outfile !=""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_samples_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_samples_df)
}


#' TaxAssign
#' 
#' Find Lowest Taxonomic Group (LTG) for each ASV in the input data frame.
#' 
#' TaxAssign uses the mkLTG algorithm described 
#' in [Meglécz, 2024](https://rdcu.be/dxABF) and  
#' [https://github.com/meglecz/mkLTG](https://github.com/meglecz/mkLTG).
#' This is a BLAST based method using a series of identity percentage of cutoff
#' values to validate BLAST hits and each cutoff is associated with
#' further parameters:
#' * pcov: Percentage of coverage
#' * phit: Perentage of the validated hits to be included in LTG
#' * taxn: Minimum number of taxa among the validated hits
#' * seqn: Minimum number of sequences among the validated hits
#' * refres: Minimum resolution of the hit to be validated
#' * ltgres: Maximum resolution of the LTG.
#'  
#' @param asv Data frame or csv file containing an asv and asv_id columns.
#' @param ltg_params Data frame or csv file with a list of 
#' percentage of identity values (pid) and associated parameters 
#' (pcov,phit,taxn,seqn,refres,ltgres).
#' @param taxonomy TSV file containing the following columns: 
#' tax_id,parent_tax_id,rank,name_txt,old_tax_id(has been merged to tax_id),
#' taxlevel (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 
#' 2: kingdom, 1: superkingdom, 0: root).
#' @param blast_db Character string naming the BLAST database.
#' @param blast_path Character string: path to BLAST executable. 
#' @param fill_lineage Boolean. If TRUE, fill in missing higher-level taxa 
#' in the lineage using the name of the next known lower-level taxon, prefixed 
#' by the current taxonomic level 
#' (e.g., No_kingdom_Chrysophyceae if kingdom is missing but class is known).
#' @param num_threads Positive integer: number of CPUs.
#' @param tax_sep Field separator character used in taxonomy file.
#' @param sep Field separator character in input and output csv files.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the following columns:
#' asv_id,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,
#' superkingdom_taxid,superkingdom,kingdom_taxid,kingdom,
#' phylum_taxid,phylum,class_taxid,class,order_taxid,order,
#' family_taxid,family,genus_taxid,genus,species_taxid,species,
#' pid,pcov,phit,taxn,seqn,refres,ltgres,asv
#' @examples
#' \dontrun{
#' TaxAssign(asv=read_count_df, taxonomy="xxxxxx", blast_db="xxxxxxxxx", num_threads=4)
#' }
#' @export
#'
#'
TaxAssign <- function(asv, 
                      taxonomy, 
                      blast_db, 
                      blast_path="blastn", 
                      ltg_params="", 
                      outfile="", 
                      num_threads=1, 
                      tax_sep="\t", 
                      sep=",",
                      quiet=T, 
                      fill_lineage=TRUE
                      ){


taxonomy <- path.expand(taxonomy)
blast_db <- path.expand(blast_db)

# can accept df or file as an input
if(is.character(asv)){
  asv_df <- read.csv(asv, header=T, sep=sep)
}else{
  asv_df <- asv
}
# get unique list
asv_df <- asv_df %>%
  select(asv_id, asv) %>%
  distinct() %>%
  arrange(asv_id)
t <- check_one_to_one_relationship(asv_df)

if (is.character(ltg_params)){ 
  if(ltg_params == ""){ # default value for ltg_params_df
    ltg_params_df = data.frame( pid=c(100,97,95,90,85,80),
                                pcov=c(70,70,70,70,70,70),
                                phit=c(70,70,70,70,70,70),
                                taxn=c(1,1,2,3,4,4),
                                seqn=c(1,1,2,3,4,4),
                                refres=c(8,8,8,7,6,6),
                                ltgres=c(8,8,8,8,7,7)
    )
  }else{ # read params from file
    ltg_params_df <- read.csv(ltg_params, header=T, sep=sep)
  }
} else{ # ltg_params is df
  ltg_params_df <- ltg_params
}


#### Read taxonomy info 
# read taxonomy file; 
# quote="" is important, since some of the taxon names have quotes and this should be ignored
tax_df <- read.delim(taxonomy, header=T, sep=tax_sep, fill=T, quote="") %>%
  select(tax_id, parent_tax_id, rank, name_txt, old_tax_id, taxlevel)

# make data frame with old taxids as line numbers and taxids in a columns
old_taxid <- tax_df %>%
  filter(!is.na(old_tax_id)) %>%
  select(tax_id, old_tax_id)
# delete old_tax_ids from tax_df and make taxids unique
tax_df <- tax_df %>%
  select(-old_tax_id)
tax_df <- unique(tax_df)

####
# create a tmp directory for temporary files using time and a random number
outdir_tmp <- paste('tmp_TaxAssign_', 
                    trunc(as.numeric(Sys.time())), 
                    sample(1:100, 1), 
                    sep=''
                    )
outdir_tmp <- file.path(tempdir(), outdir_tmp)
check_dir(outdir_tmp)

### run blast and clean/complete results
# run blast and read read results to data frame 
# (blast_res columns: "qseqid","pident","qcovhsp","staxids")
# Query seqid are the same as the asv_id
blast_res <- run_blast(asv_df,
                       blast_db=blast_db, 
                       blast_path=blast_path, 
                       outdir=outdir_tmp, 
                       qcov_hsp_perc=min(ltg_params_df$pcov),
                       perc_identity=min(ltg_params_df$pid), 
                       num_threads=num_threads, 
                       quiet=quiet
                       )
# add update old taxids to valid ones
blast_res <- update_taxids(blast_res, old_taxid)
# add taxlevel
blast_res <- left_join(blast_res, tax_df, by=c("staxids" = "tax_id")) %>%
  select(-parent_tax_id, -rank, -name_txt) # "qseqid"   "pident"   "qcovhsp"  "staxids"  "taxlevel"

### make a lineage for each taxid in blast_res
lineages <- get_lineage_ids(unique(blast_res$staxids), tax_df)

# new data frame with all asv, asv_id and NA for the other columns
taxres_df <- asv_df %>%
  mutate(ltg_taxid = NA, 
         pid=NA, 
         pcov=NA, 
         phit=NA, 
         taxn=NA, 
         seqn=NA, 
         refres=NA, 
         ltgres=NA)

for(i in 1:nrow(taxres_df)){ # go through all sequences 
  for(p in 1:nrow(ltg_params_df)){ # for each pid
    pidl <- ltg_params_df[p,"pid"]
    pcovl <- ltg_params_df[p,"pcov"]
    phitl <- ltg_params_df[p,"phit"]
    taxnl <- ltg_params_df[p,"taxn"]
    seqnl <- ltg_params_df[p,"seqn"]
    refresl <- ltg_params_df[p,"refres"]
    ltgresl <- ltg_params_df[p,"ltgres"]
    
    # filter the blastres according to qseqid,  pid, pcov, refres
    df_intern <- blast_res %>%
      filter(qseqid==taxres_df$asv_id[i] & pident>=pidl & qcovhsp>=pcovl & taxlevel>=refresl)
    
    # check if enough taxa and seq among validated hits
    tn <- length(unique(df_intern$staxids))
    if(tn >= taxnl & nrow(df_intern) >= seqnl ){
      # make ltg if all conditions are met
      ltg <- make_ltg(df_intern$staxids, lineages, phit = phitl)
      # fill out line with the ltg and the parmeters that were used to get it
      taxres_df[i, 3:(ncol(taxres_df))] <- 
        list(ltg, pidl, pcovl, phitl, taxnl, seqnl, refresl, ltgresl)
      break
    } # end if
  } # end p (pids)
} # end i (asvs)

# get the ranked lineage for each taxid in taxres_df
ranked_lineages <- get_ranked_lineages(
  unique(taxres_df$ltg_taxid), 
  tax_df, 
  fill_lineage= TRUE
  )
# add lineage to taxres_df
taxres_df <- left_join(taxres_df, ranked_lineages, by="ltg_taxid") %>%
  select(asv_id,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,domain_taxid,
         domain,kingdom_taxid,kingdom,phylum_taxid,phylum,class_taxid,class,
         order_taxid,order,family_taxid,family,genus_taxid,genus,
         species_taxid,species,pid,pcov,phit,taxn,seqn,refres,ltgres,asv)
# adjust resolution if it is higher than ltgres
taxres_df <- adjust_ltgres(taxres_df, tax_df)
# taxres_df data frame with the following columns: 
# asv_id,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,superkingdom_taxid,
# superkingdom,kingdom_taxid,kingdom,phylum_taxid,phylum,class_taxid,class,
# order_taxid,order,family_taxid,family,genus_taxid,genus,species_taxid,species,pid,
# pcov,phit,taxn,seqn,refres,ltgres,asv

# delete temporary  dir
unlink(outdir_tmp, recursive = TRUE)

if(outfile != ""){
  check_dir(outfile, is_file=TRUE)
  write.table(taxres_df, file = outfile,  row.names = F, sep=sep)
}

return(taxres_df)
}


#' Run BLAST
#' 
#' Run BLAST using ASVs of the input data fram as queries.
#'  
#' @param df Data frame containing and asv and asv_id columns.
#' @param blast_db BLAST DB including path.
#' @param blast_path Character string: path to BLAST executables. 
#' @param outdir Character string: output directory.
#' @param qcov_hsp_perc Real between 0 and 100: minimum query coverage.
#' @param perc_identity Real between 0 and 100: minimum percentage of identity.
#' @param num_threads Positive integer: number of threads.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns data frame with the BLAST results: qseqid,pident,qcovhsp,staxids
#' @examples
#' \dontrun{
#' run_blast(df=read_count_df, 
#'     blast_db="xxxxxxxx", 
#'     blast_path="blastn", 
#'     qcov_hsp_perc=80, 
#'     perc_identity=00, 
#'     num_threads=4
#'     )
#' }
#' @export
#'
run_blast <- function(df, 
                      blast_db, 
                      blast_path="blastn", 
                      outdir="", 
                      qcov_hsp_perc=70, 
                      perc_identity=70, 
                      num_threads=1, 
                      quiet=T
                      ){

  check_dir(outdir)
  
  # make fasta file with unique reads; use numbers as ids  
#  seqs <- unique(df$asv)
  fas <- file.path(outdir, 'unique.fas')
  write_fasta_df(df, outfile=fas, read_count=FALSE)
#  write_fasta(seqs, fas, seq_as_id=T)
  # define the name of the output file
  blast_out <- file.path(outdir, 'blast.out')
  
  task = "megablast"
  e = 1e-20
  dust = "yes"
  max_target_seqs=500
  
  blast <- paste(blast_path, " -task ", task, 
                 " -db ",blast_db ,
                 " -query ",fas,
                 " -evalue ",e,
                 " -out ",blast_out,
                 ' -outfmt "6 qseqid pident qcovhsp staxids" -dust ',dust,
                 " -qcov_hsp_perc ",qcov_hsp_perc,
                 " -perc_identity ",perc_identity,
                 " -num_threads ",num_threads,
                 " -max_target_seqs ",max_target_seqs, 
                 sep=""
                 )
  if(!quiet){
    print(blast)
  }
  system(blast)
  
  # read BLAST results; 
  # take care of lines where there is several taxids in the staxids column 
  # This can happen in ncbi nt
  blast_res <- read_blast_res(file=blast_out)
  return(blast_res)
}

#' Read BLAST results to data frame
#' 
#' Read BLAST result to a data frame.
#' If more than one taxid for a hit, make a separate line for each taxid.
#' 
#'  
#' @param file Character string naming the output of BLAST. 
#' Tab separated colums: qseqid,pident,qcovhsp,staxids
#' @returns Data frame with the following columns: qseqid,pident,qcovhsp,staxids
#' @examples
#' \dontrun{
#' read_blast_res(file="blastout.txt")
#' }
#' @export
#'
read_blast_res <- function(file){
  
  blast_res <- read.delim(file, header=F, sep="\t", fill=T, quote="")
  colnames(blast_res) <- c("qseqid","pident","qcovhsp","staxids") 
  
  # if BLAST against NCBI nt, 
  # the staxids can contain more than one taxids, separated by ";" 
  # => make a separate line for each
  # This is a relatively rare case, so to avoid a long loop over each line, 
  # first select lines with multiple taxids, expand them and then pool the 
  # results with the other lines
  # select lines with multiple taxids
  blast_res$staxids <- as.character(blast_res$staxids)
  df_multiple_taxids <- blast_res[grepl(";", blast_res$staxids),]
  # select lines with single taxids
  blast_res <- blast_res[!grepl(";", blast_res$staxids),]

  # make as many lines as different taxids for each input line
  df_multiple_taxids <- df_multiple_taxids %>%
    rowwise() %>%
    do(expand_rows(.))
  # change tibble to data frame
  df_multiple_taxids <- as.data.frame(df_multiple_taxids)
  
  # pool expanded and single taxid results
  blast_res <- rbind(df_multiple_taxids, blast_res) %>%
    arrange(qseqid, desc(pident))
  
  blast_res$staxids <- as.integer(blast_res$staxids)
  return(blast_res)
}

#' Expand Rows
#' 
#' Split the taxid column to different taxids in the input data frame.
#' Make one line for each taxid.
#'  
#' @param row Row of a data frame with the following columns: 
#' qseqid,pident,qcovhsp,staxids
#' @returns Data frame with the following columns: qseqid,pident,qcovhsp,staxids
#' @examples
#' \dontrun{
#' expand_rows(row)
#' }
#' @export
#'
expand_rows <- function(row){
  staxids <- as.character(row$staxids)
  taxids_list <- unlist(strsplit(staxids, ";"))
  new_rows <- data.frame(
    qseqid = rep(row$qseqid, length(taxids_list)),
    pident = rep(row$pident, length(taxids_list)),
    qcovhsp = rep(row$qcovhsp, length(taxids_list)),
    staxids = as.integer(taxids_list)
  )
  return(new_rows)
}


#' Update taxids
#' 
#' Replaces old taxids (merged to other taxIDs in NCBI taxonomy) by valid ones.
#'  
#' @param df Data frame with the following columns: qseqid,pident,qcovhsp,staxids
#' @param old_taxid Data frame with the following columns:  tax_id,old_tax_id
#' @returns Data frame with the following columns: qseqid,pident,qcovhsp,staxids,
#' where the merged taxids has been replaced by the up to date ones.
#' @examples
#' \dontrun{
#' update_taxids(df=blastout_df, old_taxid=old_taxid_df)
#' }
#' @export
#'
update_taxids <- function(df, old_taxid){
  # df is a data frame with the following columns: qseqid,pident,qcovhsp,staxids
  # old_taxid is a data frame with the following columns:  tax_id,old_tax_id
  
  # replace old taxids (if any) in df by up to date ones 
  df <- left_join(df, old_taxid, by=c("staxids" = "old_tax_id"))
  df$staxids[which(!is.na(df$tax_id))] <- df$tax_id[which(!is.na(df$tax_id))]
  # delete tax_id column since the values (if non NA were used to replace staxids)
  df <- df %>%
    select(-tax_id)
  return(df)
}

#' Get the complete taxonomic lineage using taxIDS
#' 
#' Get the complete taxonomic lineage of each taxID in the input vector.
#'  
#' @param taxids Vector of taxIDs (taxonomic IDs).
#' @param tax_df Data frame with the following columns: 
#' tax_id, parent_tax_id, rank, name_txt, taxlevel 
#' (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 2: kingdom, 1: superkingdom, 0: root).
#' @returns Data frame with each line starting by a taxID followed by a vector taxids
#' in its full lineage (starting by the highest taxonomic level).
#' @examples
#' \dontrun{
#' taxids <- c(9606, 9593)
#' tax_df <- read.delim(taxonomy="xxxxxxxx", header=T, sep="\t", fill=T, quote="")
#' get_lineage_ids(taxids, tax_df=tax_df)
#' }
#' @export
#'
get_lineage_ids <- function(taxids, tax_df){
  
  # taxids is a vector of taxids; there can be duplicated values
  lineages <- as.data.frame(taxids)
  colnames(lineages) <- c("tax_id")
  
  i <- 1 # i is the number of itaration. µIt should stop, if all lineages arrived to the root
  while(i < 100){
    # use i as name instead of tax_id
    new_colname <- as.character(i)
    # add parent_tax_id and rename columns
    lineages <- left_join(lineages, tax_df, by="tax_id")%>%
      select(-rank, -name_txt, -taxlevel) %>%
      # !! = interpret the variable
      rename(!!new_colname :=tax_id, "tax_id"=parent_tax_id)
    
    i <- i+1
    # stop if all lines has the same value (usually 1)
    tid_list <- unique(lineages$tax_id)
    if(length(tid_list) == 1 && tid_list[1] ==1){
      break
    }
  }
  # delete the last column, where all values are 1
  lineages <- lineages %>%
    select(-tax_id)
  # reverse order of columns
  lineages <- lineages[, ncol(lineages):1]
  # Apply the function to each row of the lineages data frame: 
  # delete all 1, shift the remaining elements of each row to the beginning, 
  # and replace missing values at the end of the row by NA
  lineages <- as.data.frame(t(apply(lineages, 1, delete_1_by_row)))
  # add as a first column the taxid, so they can be easily accessed
  lineages <- cbind(taxids, lineages)
  
  return(lineages)
}

#' delete_1_by_row
#' 
#' Delete all 1 from the beginning of a raw, shift the other values and 
#' replace de missing ones at the end by NA
#'  
#' @param row Vector of taxids.
#' @returns
#' Vector of taxids.
#' @examples
#' \dontrun{
#' delete_1_by_row(row)
#' }
#' @export
#'
delete_1_by_row <- function(row) {
  
  n <- length(row) 
  # Remove all occurrences of 1
  row <- row[row != 1]
  
  # Create a new row with NA at the end
  new_row <- c(row, rep(NA, n - length(row)))
  
  return(new_row)
}

#' Make LTG
#' 
#' Determine the Lowest Taxonomic Group (LTG) that contains 
#' phit percentage of the input taxids.
#'  
#' @param taxids Vector of taxIDs (taxonomic IDs). 
#' There can be duplicated values.
#' @param lineages Data frame: taxids in the first column followed by 
#' their lineages represented by taxids (starting at the lowest resolution).
#' @param phit Integer between 0 and 100: Percentage of taxids that should 
#' be included to the LTG.
#' @returns Numerical taxid of the LTG, or NA if LTG cannot be defined.
#' @examples
#' \dontrun{
#' taxids <- c(9593,9606,9606,9606,9606,9606,9606,9606,9606)
#' make_ltg(taxids, lineages=lineages, phit=80)
#' }
#' @export
#'
make_ltg <- function(taxids, lineages, phit=70){
  # taxids is a vector of taxids; there can be duplicated values
  # make a data frame from the vector
  lin <- as.data.frame(taxids)
  colnames(lin) <- c("staxid")
  
  # add lineage to each taxid
  lin <- left_join(lin, lineages, by=c("staxid" = "taxids")) %>%
    select(-where(~all(is.na(.)))) # delete columns if all elements are NA
  
  ltg <- NA
  if(length(unique(lin$staxid)) == 1){ # only one taxid among hits; avoid loops
    ltg <-lin$staxid[1]
  }else{
    for(i in 2:ncol(lin)){ # start from low resolution
      tmp <- as.data.frame(lin[,i])
      colnames(tmp) <- c("tax_id")
      # get unique taxids, and their numbers in the i-th column
      tmp <- tmp %>%
        group_by(tax_id) %>%
        summarize(nhit=length(tax_id)) %>%
        arrange(desc(nhit)) %>%
        ungroup()
      
      # stop, if the taxid with the highest number of sequences does not contain 
      # at least phit percent of the hits
      max_hitn <- as.numeric(tmp[1,"nhit"])
      if(is.na(tmp[1,"tax_id"])){ # the most frequent "taxid" is NA
        break
      }
      if(max_hitn/sum(tmp[,"nhit"]) < phit/100){# the most frequent taxid has less than phit%
        break
      }
      ltg <- as.numeric(tmp[1,"tax_id"])
      
    }
  }
  return(ltg)
}

#' Get Ranked Lineages
#' 
#' For each input taxid, get the major taxonomic ranks of their lineage.
#'  
#' @param taxids Vector of taxIDs.
#' @param tax_df Data frame with the following columns: 
#' tax_id, parent_tax_id, rank, name_txt, taxlevel 
#' (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 
#' 2: kingdom, 1: superkingdom, 0: root).
#' @param fill_lineage Boolean. If TRUE, fill in missing higher-level taxa 
#' in the lineage using the name of the next known lower-level taxon, prefixed 
#' by the current taxonomic level 
#' (e.g., No_kingdom_Chrysophyceae if kingdom is missing but class is known).
#' @returns Data frame with the ranked lineages of taxids. 
#' Columns: ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,
#' superkingdom_taxid,superkingdom,kingdom_taxid,kingdom,phylum_taxid,phylum,
#' class_taxid,class,order_taxid,order,family_taxid,family,
#' genus_taxid,genus,species_taxid,species).
#' @examples
#' \dontrun{
#' taxids <- c(9593,9606)
#' tax_df <- read.delim(taxonomy="xxxxxxxx", header=T, sep="\t", fill=T, quote="")
#' get_ranked_lineages(taxids, tax_df)
#' }
#' @export
#'

get_ranked_lineages <- function(taxids, tax_df, fill_lineage=TRUE){
  
  # taxids is a vector of taxids; there can be duplicated values
  ranked_lineages <- as.data.frame(taxids)%>%
    filter(!is.na(taxids)) %>%
    rename(tax_id=taxids)
  # make tmp data frame to keep a list of taxids
  tmp <- ranked_lineages
  # make tmp_lin data frame to keep a list of taxids and the lineage of each taxid 
  # (including, names, taxid, taxlevel)
  tmp_lin <- ranked_lineages
  
  # define first colums, with taxid, name, taxlevel
  ranked_lineages <- left_join(ranked_lineages, tax_df, by="tax_id")%>%
    select(-parent_tax_id)%>%
    rename(ltg_taxid=tax_id, ltg_name=name_txt, ltg_rank=rank, ltg_rank_index=taxlevel) %>%
    select(ltg_taxid, ltg_name, ltg_rank, ltg_rank_index)
  # add columns for each major taxlevel (taxid and name)
  now_cols <- c(
    "domain_taxid", "domain",
    "kingdom_taxid", "kingdom",
    "phylum_taxid", "phylum",
    "class_taxid", "class",
    "order_taxid", "order",
    "family_taxid", "family",
    "genus_taxid", "genus",
    "species_taxid", "species"
  )
  ranked_lineages[now_cols] <- NA
  
  i <- 1
  # get linegaes of each taxid
  while(i < 100){
    # get the tax name, and tax rank for each taxid
    tmp <- left_join(tmp, tax_df, by="tax_id")
    # info in tmp_lin
    tmp_lin <- cbind(tmp_lin, tmp$tax_id, tmp$name_txt, tmp$rank )
    # re-initilize tmp
    tmp <- tmp %>%
      select(parent_tax_id)%>%
      rename(tax_id=parent_tax_id)
    # stop if all linage ends with root
    if(all(tmp$tax_id ==1)){
      break
    }
    i<- i+1
  }
  
  # select only major taxonomic ranks from each line of tmp_lin; 
  # keep the results in ranked_lineages
  for (c in seq(from=6, to=ncol(ranked_lineages), by=2)) {# go though all major taxlevel
    taxrank <- colnames(ranked_lineages[c])
    for (i in 1:nrow(tmp_lin)) {
      row <- tmp_lin[i, ]  # Extract the current row
      col_index <- which(row == taxrank)  # Find the column index containing "species"
      
      if (length(col_index) > 0) {
        # Add taxon name and taxid to ranked_lineages
        ranked_lineages[i,c-1] <- tmp_lin[i,col_index-2]
        ranked_lineages[i,c] <- tmp_lin[i,col_index-1]
      }
    }
  }
  
  # if NA in a high level taxon and non NA in lower level taxon, 
  # replace NA by taxlevel_lower_level_taxon
  if(fill_lineage){
    ranked_lineages <- fill_NA_in_lineage(ranked_lineages)
  }
  
  return(ranked_lineages)
}

#' Fill in missing higher-level taxa
#' 
#' Fill in missing higher-level taxa in the lineage using the name of the 
#' next known lower-level taxon, prefixed by the current taxonomic level 
#'  
#' @param df A data frame containing taxonomic lineages. Major taxonomic levels 
#' are located in every second column, starting from column 6 (e.g. 6: domain, 
#' 8: kingdom, 10: phylum, 12: class, 14: order, 16:family, 18: genus, 20 species).
#' @returns Input data frame with missing taxon names replaced by the 
#' next known lower-level taxon, prefixed by the current taxonomic level
#' @examples
#' \dontrun{
#' df <- data.frame(
#' id = 1:3,
#' name = c("OTU1", "OTU2", "OTU3"),
#' sample1 = c(10, 5, 0),
#' sample2 = c(3, 7, 2),
#' domain_taxid = c(2, 2759, 2759),
#' domain = c("Bacteria", "Eukaryota", "Eukaryota"),
#' phylum_taxid = c(1224, 4762, NA),
#' phylum = c("Pseudomonadota", "Oomycota", NA),
#' class_taxid = c(1236, NA, 2825),
#' class = c("Gammaproteobacteria", NA, "Chrysophyceae"))
#' df1 <- fill_NA_in_lineage(df)
#' }
#' @export
#'
fill_NA_in_lineage <- function(df) {
  
  taxonomic_levels <- colnames(df)
  # Loop through the taxonomic levels from highest to second-lowest
  indices <- seq(from=6, to=(ncol(df)-2), by=2)
  
  for (i in indices) {
    current_level <- taxonomic_levels[i]
    for (j in seq(from=i+2, to=length(taxonomic_levels), by=2)) {
      lower_level <- taxonomic_levels[j]
      # Replace NA in current_level with 
      # paste0(current_level, "_", lower_level_value), if lower_level is not NA
      missing <- is.na(df[[current_level]]) & !is.na(df[[lower_level]])
      df[[current_level]][missing] <- paste0("No_", current_level, "_", df[[lower_level]][missing])
    }
  }
  return(df)
}


#' Adjust the resolution of LTG
#' 
#' If the LTG has a higher resolution than the ltgres parameter,
#' adjust the LTG and stop lineage at ltgres level.
#'  
#' @param taxres_df Data frame with the following columns: asv_id,ltg_taxid,
#' ltg_name,ltg_rank,ltg_rank_index,superkingdom_taxid,
#' superkingdom,kingdom_taxid,kingdom,phylum_taxid,phylum,class_taxid,class,
#' order_taxid,order,family_taxid,family,genus_taxid,genus,species_taxid,species,pid,
#' pcov,phit,taxn,seqn,refres,ltgres,asv.
#' @param tax_df Data frame with the following columns: 
#' tax_id, parent_tax_id, rank, name_txt, taxlevel 
#' (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 
#' 2: kingdom, 1: superkingdom, 0: root).
#' @returns The modified input data frame with lower resolution of LTGs if necessary.
#' @examples
#' \dontrun{
#' tax_df <- read.delim(taxonomy="xxxxxxxx", header=T, sep="\t", fill=T, quote="")
#' adjust_ltgres(taxres_df, tax_df)
#' }
#' @export
#'
adjust_ltgres <- function(taxres_df, tax_df){
  
  # link taxlevel index and tax rank
  taxlevel_index = data.frame(taxlevel_index=c(1,2,3,4,5,6,7,8),
                              taxrank=c("superkingdom","kingdom","phylum",
                                        "class","order","family","genus","species")
  )
  
  # add the name of the tax rank equivalent to the index in ltgref
  taxres_df <- left_join(taxres_df, taxlevel_index, by=c("ltgres" = "taxlevel_index"))
  
  for(i in 1:nrow(taxres_df)){ # all rows
    
    if(!is.na(taxres_df[i,"ltg_taxid"]) & 
       taxres_df[i,"ltg_rank_index"] > taxres_df[i,"ltgres"])
      { # if resolution of ltg is higher then ltgres
      # get the taxrank that corresponds to ltgres 
      tl <- taxres_df[i,"taxrank"]
      # get the index of the column that corresponds to the ltgres
      col_index <- which(colnames(taxres_df) == tl)
      # make a data frame with taxid, and get taxinfo from tax_df
      new_taxid <- as.data.frame(taxres_df[i, col_index-1]) 
      colnames(new_taxid) <- c("tax_id")
      new_taxid <- left_join(new_taxid, tax_df, by="tax_id") %>%
        select(tax_id, name_txt, rank, taxlevel)
      
      # replace ltg taxid and associated info
      taxres_df[i, 2:5] <- new_taxid[1,]
      # replace tax lineage over the ltgref by NA
      taxres_df[i, (col_index+1):(ncol(taxres_df)-9)] <- NA
    }# end if
  }# end for i
  
  taxres_df <- taxres_df %>%
    select(-taxrank)
  
  return(taxres_df)
}

#' Write ASV table
#' 
#' Write csv file with samples in columns, ASVs in lines, read_counts in cells.
#'  
#' @param read_count_samples_df Data frame with the following variables: 
#' asv_id, sample, read_count, asv.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param asv_tax Data frame or CSV file with taxonomic assignments with the following columns:
#' asv_id,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,
#' superkingdom_taxid,superkingdom,kingdom_taxid,kingdom,
#' phylum_taxid,phylum,class_taxid,class,order_taxid,order,
#' family_taxid,family,genus_taxid,genus,species_taxid,species.
#' pid,pcov,phit,taxn,seqn,refres,ltgres,asv. 
#' If given, the output is completed with taxonomic assignment of each ASV.
#' @param sortedinfo Data frame or csv file with columns: sample, sample_type.
#' Only necessary if add_empty_samples==T or add_expected_asv==T.
#' @param add_empty_samples logical: add a column for each samples 
#' in the original data set, even if they do not have reads after filtering.
#' @param add_sums_by_sample logical: add a line with the total number of reads
#'  in each sample, and another with the number of ASVs in each sample.
#' @param add_sums_by_asv logical: add a column with the total number of reads 
#' for each ASV, and another with the number of samples, where the ASV is present.
#' @param add_expected_asv logical: add a column for each mock sample in which 
#' keep and tolerate ASVs are flagged.
#' @param mock_composition Data frame or CSV file with the following columns: 
#' sample,action,asv. Action can take the following values: keep/tolerate.
#' Only necessary if add_expected_asv==T.
#' @param sep Field separator character in input and output csv files.
#' @returns Data frame corresponding to the output file. 
#' Samples in columns, ASVs in lines, read_counts in cells, plus
#' other information according to the input parameters.
#' @examples
#' \dontrun{
#' WriteASVtable(read_count_samples_df, 
#'     outfile="out/asv_table.csv", 
#'     asv_tax=asv_tax, 
#'     sortedinfo=sortedinfo_df,
#'     add_empty_samples=T, 
#'     add_sums_by_sample=T,
#'     add_sums_by_asv=T, 
#'     add_expected_asv=T,
#'     mock_composition="data/mock_compostion.csv"
#'     )
#' }
#' @export
#'
WriteASVtable <- function(read_count_samples_df, 
                          outfile="", 
                          asv_tax=NULL, 
                          sortedinfo="", 
                          add_empty_samples=F, 
                          add_sums_by_sample=F, 
                          add_sums_by_asv=F, 
                          add_expected_asv=F,
                          mock_composition="", 
                          sep=","
                          ){
  
  # make a wide data frame with samples in columns, ASVs in lines
  wide_read_count_df <- as.data.frame(pivot_wider(
    read_count_samples_df, 
    names_from = c(sample), 
    values_from = read_count, 
    values_fill=0, 
    names_sep = ".", 
    names_sort=T)
    )
  # put the asv column at the end

  # read the sortedinfo to a data frame 
  if(add_empty_samples | add_expected_asv){
    if(is.character(sortedinfo)){
      sortedinfo_df <- read.csv(sortedinfo, header=T, sep=sep)
    }else{
      sortedinfo_df <- sortedinfo
    }
  }
  
  
  if(add_empty_samples){
    # make vector with samples already in the data frame 
    # (asv_id and asv is also on the list, but it is not a pb)
    samples <- colnames(wide_read_count_df)
    # number of ASVs
    n <- nrow(wide_read_count_df)
    
    # make a vector with all unique samples in the sortedinfo
    all_samples <-unique(sortedinfo_df$sample)
    
    # add a column for each samples that are not yet in the data frame, 
    # with 0 read counts for all variants
    for(sample in all_samples){
      if(!(sample %in% samples)){
        wide_read_count_df[[sample]] <- rep(0, n)
      }
    }
  }
  
  # add a line with the total number of reads of each sample and 
  # another with the number of ASVs in the sample
  if(add_sums_by_sample){
    
    # make a data frame with same columns as wide_read_count_df
    sum_rc <- data.frame(matrix(0, nrow=2, ncol= ncol(wide_read_count_df)))
    colnames(sum_rc) <- colnames(wide_read_count_df)
    #  and total number of reads in line 1 
    sum_rc[1,1] <- "NA" # asv_id col
    sum_rc[1,2] <- NA # asv col
    # total number of reads for each sample (ignore cols 1 and 2, since it is asv_id ans asv)
    sum_rc[1,-c(1,2)] <- colSums(wide_read_count_df[,-c(1,2)])
    # Number of ASVs in each sample in line 2
    sum_rc[2,1] <- "NA"
    sum_rc[2,2] <- NA # asv col
    sum_rc[2,-c(1,2)] <- colSums(wide_read_count_df[,-c(1,2)] != 0)
    wide_read_count_df <- rbind(sum_rc, wide_read_count_df)
  }
  
  # add sum of read count and the number of occurrences for each asv
  if(add_sums_by_asv){
    # count the total number of reads for each asv and add a column
    wide_read_count_df$Total_number_of_reads <- rowSums(wide_read_count_df[,-c(1,2)])
    # This is the sum of number of ASVs per sample, Does not make much sens
    wide_read_count_df$Total_number_of_reads[2] <- NA 
    
    # count the number of samples where the ASV is present
    tmp <- read_count_samples_df %>%
      group_by(asv) %>%
      summarize(Number_of_samples=length(sample)) %>%
      ungroup()
    # add sample count to wide_read_count_df
    wide_read_count_df <- full_join(wide_read_count_df, tmp, by="asv")
  }
  
  if(add_expected_asv){
    
    # keep only mock samples in sortedinfo_df
    sortedinfo_df <- sortedinfo_df %>%
      filter(sample_type=="mock")
    # make a vector with all unique samples in the sortedinfo
    mock_samples <-unique(sortedinfo_df$sample)
    
    
    if(is.character(mock_composition)){
      mock_asv <-  read.csv(mock_composition, header=T, sep=sep)
      CheckFileinfo(file=mock_composition, file_type="mock_composition", sep=sep, quiet=TRUE)
    }else{
      mock_asv <- mock_composition
    }
    # keep only keep and tolerate action, in case the file contains other lines 
    mock_asv <- mock_asv%>%
      filter(action=="keep" | action=="tolerate")


    # add a column for each mock samples with keep or tolerate if relevant for each ASV 
    for(mock in mock_samples){
      # make a df containing only data for a given mock sample
      df <- mock_asv %>%
        filter(sample==mock) %>%
        select(action,asv)
      
      # add action to wide_read_count_df
      new_colname <- paste("action", mock,  sep=".")
      wide_read_count_df <- left_join(wide_read_count_df, df, by="asv") %>%
        rename_with(~new_colname, action)
    }
  }
  
  
  if(!is.null(asv_tax)){ #  taxonomic assignation is given
    if(is.character(asv_tax)  && asv_tax != ""){ # as a file
      asv_tax <- read.csv(asv_tax, header=T, sep=sep)
    }
      asv_tax$asv_id <- as.character(asv_tax$asv_id)
      wide_read_count_df$asv_id <- as.character(wide_read_count_df$asv_id)
      wide_read_count_df <- left_join(wide_read_count_df, asv_tax, by=c("asv_id", "asv"))
  }
  
  # put the asv column at the end
  wide_read_count_df <- wide_read_count_df %>%
    select(-asv, everything(), asv)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(wide_read_count_df, file=outfile, row.names = F, sep=sep)
  }
  return(wide_read_count_df)
}


#' OptimizePCRerror
#' 
#' Prepare a data frame that lists pairs of expected and unexpected ASVs 
#' in mock samples with maximum `max_mismatch` difference between them.
#' The `pcr_error_var_prop` parameter should be above the highest 
#' `pcr_error_var_prop` (unexpected_read_count/expected_read_count) 
#' value in the table.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param mock_composition Data frame or csv file with columns: 
#' sample, action (keep/tolerate), asv.
#' @param vsearch_path Character string: path to vsearch executables.
#' @param sep Field separator character in input and output csv files.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param max_mismatch Positive integer: maximum number of mismatches 
#' allowed between two asv to be compared.
#' @param min_read_count Positive integer: occurrences with fewer 
#' reads are ignored. 
#' @returns Data frame with the following columns: sample,expected_read_count,
#' unexpected_read_count,pcr_error_var_prop,expected_asv_id,unexpected_asv_id,
#' expected_asv,unexpected_asv
#' @examples
#' \dontrun{
#' OptimizePCRerror(read_count=read_count_df, 
#'     mock_composition="data/mock_composition.csv", 
#'     vsearch_path=vsearch_path, 
#'     max_mismatch=2, 
#'     min_read_count=5
#'     )
#' }
#' @export
#'

OptimizePCRerror <- function(read_count, 
                             mock_composition="", 
                             vsearch_path= "vsearch", 
                             sep=",", 
                             outfile="", 
                             max_mismatch=1, 
                             min_read_count=2
                             ){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  # can accept df or file as an input
  if(is.character(mock_composition)){
    # read known occurrences
    mock_composition_df <- read.csv(mock_composition, header=T, sep=sep)
  }else{
    mock_composition_df <- mock_composition
  }
  CheckFileinfo(file=mock_composition_df, file_type="mock_composition", sep=sep, quiet=TRUE)
  
  # read the mock composition file and keep only lines with keep and tolerate
  mock_composition_df <- mock_composition_df %>%
    filter(action=="keep" | action=="tolerate")
  unique_mock_list <- unique(mock_composition_df$sample)
  
  # sum read_counts over replicates 
  df <- read_count_df %>%
    group_by(sample, asv, asv_id) %>%
    summarize(read_count_sample=sum(read_count), .groups="drop_last") %>%
    filter(read_count_sample >=min_read_count) %>%
    filter(sample %in% unique_mock_list) %>%
    ungroup()
  
  # define an empty data frame for the output
  asv_pairs <- data.frame(
    sample= character(),
    expected_read_count= numeric(),
    unexpected_read_count= numeric(),
    pcr_error_var_prop= numeric(),
    expected_asv_id= numeric(),
    unexpected_asv_id= numeric(),
    expected_asv= character(),
    unexpected_asv= character())
  ###
  # loop over all mock samples
  ###
  for(mock in unique_mock_list){
    outdir_tmp <- paste('tmp_OptimizePCRError_', 
                        trunc(as.numeric(Sys.time())), 
                        sample(1:100, 1), 
                        sep=''
                        )
    outdir_tmp <- file.path(tempdir(), outdir_tmp)
    check_dir(outdir_tmp)
    # get the list of keep ASV in the given mock sample from mock_composition
    tmp_mock <- mock_composition_df %>%
      filter(sample==mock) %>%
      filter(action=="keep")
    asv_list_keep <- unique(tmp_mock$asv)
    # make fasta file with unique mock variants; use sequences as ids
    fas_keep <- file.path(outdir_tmp, paste(mock, "keepASV.fas", sep="_"))
    write_fasta(asv_list_keep, fas_keep, seq_as_id=T)
    
    # get the list of tolerate ASV in the given mock sample from mock_composition
    tmp_mock <- mock_composition_df %>%
      filter(sample==mock) %>%
      filter(action=="tolerate")
    asv_list_tolerate <- unique(tmp_mock$asv)
    
    # get list of ASVs present in the mock sample in read_count_df 
    # that are neither keep nor tolerate 
    tmp <- df %>%
      filter(sample==mock) %>%
      filter(!(asv %in% asv_list_keep)) %>%
      filter(!(asv %in% asv_list_tolerate))
    asv_list_delete <- unique(tmp$asv)   
    # make fasta file with unique variants that are neither keep nor tolerate in mock; 
    # use sequences as ids
    fas_delete <- paste(mock, "deleteASV.fas", sep="_")
    fas_delete <- file.path(outdir_tmp, fas_delete)
    write_fasta(asv_list_delete, fas_delete, seq_as_id=T)
    
    if(length(asv_list_delete)>0 && length(asv_list_keep)>0){ # sequences in both files
      # vsearch --usearch_global to find highly similar sequence pairs
      vsearch_out <- paste(mock, 'vsearch_out.out', sep="_")
      vsearch_out <- file.path(outdir_tmp, vsearch_out)
      vsearch <- paste(vsearch_path, " --usearch_global ", 
                       fas_delete, " --db ", fas_keep, 
                       ' --quiet --iddef 1 --self --id 0.90 --maxaccepts 0 --maxrejects 0 ',
                       ' --userfields "query+target+ids+aln"',
                       " --userout ", vsearch_out, 
                       sep="")
      #https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/system
      system(vsearch)

      if(file.exists(vsearch_out) && file.size(vsearch_out) > 0){
        # read vsearch results
        results_vsearch<- read.csv(vsearch_out, header = FALSE, sep="\t")
        colnames(results_vsearch) <- c("query","target","nb_ids","aln")
        # none of the values easily outputted by vsearch take into the external gaps as a diff 
        # => correct this, based on the alnlen and the number of identities
        results_vsearch$nb_diff <- nchar(results_vsearch$aln) - results_vsearch$nb_ids
        # keep only pairs with 1 difference 
        results_vsearch <- results_vsearch %>%
          filter(nb_diff <= max_mismatch)
        # delete unnecessary columns and add plate, marker, sample
        results_vsearch <- select(results_vsearch, -c(nb_ids, aln, nb_diff))
        if(nrow(results_vsearch) == 0){
          break
        }
        
        results_vsearch$sample <- rep(mock, nrow(results_vsearch))
        # add read_counts to results_vsearch
        results_vsearch <- rename(results_vsearch, asv = target)
        results_vsearch <- left_join(results_vsearch, df, by=c("sample", "asv")) 
        results_vsearch <- results_vsearch %>%
          select(sample, 
                 expected_read_count = read_count_sample, 
                 expected_asv = asv, 
                 expected_asv_id = asv_id, 
                 query
                 )
        
        results_vsearch <- rename(results_vsearch, asv = query)
        results_vsearch <- left_join(results_vsearch, df, by=c("sample", "asv"))
        results_vsearch <- results_vsearch %>%
          select(sample, 
                 expected_read_count, 
                 unexpected_read_count = read_count_sample, 
                 expected_asv_id,  
                 unexpected_asv_id = asv_id, 
                 expected_asv, 
                 unexpected_asv = asv
                 )      
        
        # delete row if the expected variant is non in sample
        results_vsearch <- results_vsearch %>%
          filter(!is.na(expected_read_count))
        
        results_vsearch$pcr_error_var_prop <- results_vsearch$unexpected_read_count / 
          results_vsearch$expected_read_count
        results_vsearch <- results_vsearch %>%
          arrange(desc(pcr_error_var_prop)) %>%
          select(sample, 
                 expected_read_count, 
                 unexpected_read_count, 
                 pcr_error_var_prop,
                 expected_asv_id, 
                 unexpected_asv_id, 
                 expected_asv, 
                 unexpected_asv
                 )
        # append results to existing asv_pairs
        asv_pairs <- rbind(asv_pairs, results_vsearch)
      }
    }
    # Delete the temp directory
    unlink(outdir_tmp, recursive = TRUE)
  }
  ###
  # end loop 
  ### 
  
  asv_pairs <- asv_pairs %>%
    arrange(desc(pcr_error_var_prop))
  
  # Delete the temp directory
  unlink(outdir_tmp, recursive = TRUE)
  
  if(outfile != "")
  {
    check_dir(outfile, is_file=TRUE)
    write.table(asv_pairs, file=outfile, sep=sep, row.names = F)
  }
  return(asv_pairs)
}

#' OptimizeLFNsampleReplicate
#' 
#' Prepare a data frame that lists all expected occurrences in all mock 
#' sample-replicates their read_counts and the proportion of 
#' their read_counts to the total number of reads in the sample-replicate. 
#' The cutoff parameter in LFNsampleReplicate should be bellow the smallest 
#' proportion in order to keep all expected ASVs in the data set.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, plate, marker, sample, replicate, read_count, asv.
#' @param mock_composition Data frame or csv file with columns: 
#' sample, action (keep/tolerate), asv.
#' @param sep Field separator character in input and output csv files.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @returns Data frame with the following columns: sample, replicate, action, 
#' asv_id, read_count, read_count_sample_replicate, lfn_sample_replicate_cutoff, 
#' asv
#' @examples
#' \dontrun{
#' OptimizeLFNsampleReplicate(read_count_df, mock_composition="data/mock_composition.csv")
#' }
#' @export
#'
OptimizeLFNsampleReplicate <- function(read_count, mock_composition="", sep=",", outfile=""){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  if(is.character(mock_composition)){
    # read known occurrences
    mock_composition_df <- read.csv(mock_composition, header=T, sep=sep)
  }else{
    mock_composition_df <- mock_composition
  }
  CheckFileinfo(file=mock_composition_df, file_type="mock_composition", sep=sep, quiet=TRUE)
  
  # read the mock composition file and keep only lines with keep
  mock_composition_df <- mock_composition_df %>%
    filter(action=="keep")
  # there is a asv_id column in mock_composition => delete it
  if("asv_id" %in% colnames(mock_composition_df )){ 
    mock_composition_df <- mock_composition_df %>%
    select(-asv_id)
  }

  
  # get a complete and unique list of sample, replicate
  sample_replicate_list <- read_count_df %>%
    ungroup() %>%
    select(sample, replicate) %>%
    unique
  
  # add replicate to mock_composition
  mock_composition_df <- left_join(mock_composition_df, 
                                   sample_replicate_list, 
                                   by=c("sample"), 
                                   relationship = "many-to-many"
                                   )
  unique_asv_keep <-  unique(mock_composition_df$asv)
  
  # get the total number of reads for each sample replicate for the mocks
  sample_replicate_rc <- read_count_df %>%
    group_by(sample, replicate) %>%
    summarize(read_count_sample_replicate= sum(read_count), .groups="drop_last") %>%
    filter(sample %in% mock_composition_df$sample) %>%
    ungroup()
  
  # sum read_counts over replicates 
  asv_keep_df <- left_join(mock_composition_df, 
                           read_count_df, 
                           by=c("sample", "replicate", "asv")
                           )
  asv_keep_df <- left_join(asv_keep_df, 
                           sample_replicate_rc, 
                           by=c("sample", "replicate")
                           )
  asv_keep_df$lfn_sample_replicate_cutoff <- 
    asv_keep_df$read_count/asv_keep_df$read_count_sample_replicate
  asv_keep_df$lfn_sample_replicate_cutoff <- 
    round(asv_keep_df$lfn_sample_replicate_cutoff-0.00005, digits=4)
  
  asv_keep_df <- asv_keep_df %>%
    arrange(lfn_sample_replicate_cutoff) %>%
    select(sample, 
           replicate, 
           action, 
           asv_id, 
           read_count,
           read_count_sample_replicate, 
           lfn_sample_replicate_cutoff,
           asv, 
           everything()
           )
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(asv_keep_df, file=outfile, sep=sep, row.names = F)
  }
  return(asv_keep_df)
}


#' MakeKnownOccurrences
#' 
#' Prepare three data frames:
#' * known_occurrences: Lists all occurrences that are expected in mock 
#' (True Positives) or 
#' False positives (unexpected variants in mocks, all variants in negative 
#' control samples, variants present in a wrong habitat).
#' * missing_occurrences: Lists all False Negatives (expected occurrences 
#' in mock samples that are missing).
#' * performance_metrics: Number of TP (expected in mock and present), 
#' FP (unexpected but present), 
#' FN (expected but absent), Precision (TP/(TP+FP)) and Sensitivity (TP/(TP+FN)).
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, plate, marker, sample, replicate (optional), read_count, asv.
#' @param sortedinfo Data frame or csv file with columns: 
#' sample, sample_type(mock/negative/real), habitat, replicate, (optional: file).
#' @param mock_composition Data frame or csv file with columns: 
#' sample, action (keep/tolerate), asv.
#' @param sep Field separator character in input and output csv files.
#' @param known_occurrences Character string: csv file containing 
#' known occurrences (expected occurrences in mock, and FP).
#' If empty, no file is written.
#' @param missing_occurrences Character string: csv file containing
#' the missing occurrences (FN). If empty, no file is written.
#' @param performance_metrics Character string: csv file containing performance 
#' metrics. If empty, no file is written.
#' @param habitat_proportion Real (between 0-1): for each asv, if the proportion 
#' of reads in a habitat is below this cutoff, 
#' it is considered as an artifact in all samples of the habitat.
#' @returns List of data frames:
#' * known_occurrences_df: sample, action, asv_id, asv
#' * missing_occurrences_df: sample, action, asv, asv_id
#' * performance_metrics_df: TP, FP, FN, Precision, Sensitivity
#' @examples
#' \dontrun{
#' results <- MakeKnownOccurrences(read_count_samples_df, 
#' sortedinfo=sortedinfo_df, 
#' mock_composition="data/mock_composition.csv", 
#' habitat_proportion=0.7
#' )
#' known_occurrences_df <- results[[1]]
#' missing_occurrences_df <- results[[2]]
#' performance_metrics_df <- results[[3]]
#' }
#' @export
#' 

MakeKnownOccurrences <- function(read_count, 
                                 sortedinfo, 
                                 mock_composition, 
                                 sep=",", 
                                 known_occurrences="", 
                                 missing_occurrences="", 
                                 performance_metrics="", 
                                 habitat_proportion=0.5){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # PoolReplicates before counting occurrences in samples
  if("replicate" %in% colnames(read_count_df)){  
    CheckFileinfo(file=read_count_df, 
                  file_type="read_count", 
                  quiet=TRUE)
    read_count_samples_df <- PoolReplicates(read_count_df)
  }else{
    CheckFileinfo(file=read_count_df, 
                  file_type="read_count_sample", 
                  quiet=TRUE)
    read_count_samples_df <- read_count_df
  }
 
  
  if(is.character(sortedinfo)){
    sortedinfo_df <- read.csv(sortedinfo, header=T, sep=sep)
  }else{
    sortedinfo_df <- sortedinfo
  }

  if(is.character(mock_composition)){
    mock_composition_df <- read.csv(mock_composition, header=T, sep=sep)
  }else{
    mock_composition_df <- mock_composition
  }
  CheckFileinfo(file=mock_composition_df, 
                file_type="mock_composition", 
                quiet=TRUE)
  
  # read info on samples types and keep only relevant info
  sortedinfo_df <- sortedinfo_df %>%
    select(sample, sample_type, habitat)
  # get unique lines to avoid replicates
  sortedinfo_df <- unique(sortedinfo_df)
  
  # define data frame for known occurrences
  occurrence_df <- read_count_samples_df
  # add habitat and sample_type to occurrence_df
  occurrence_df <- left_join(occurrence_df, sortedinfo_df, by="sample")
  # add action column
  occurrence_df$action <- rep(NA, nrow(occurrence_df))
  
  # flag occurrences in negative control samples as delete
  occurrence_df$action[which(occurrence_df$sample_type=="negative")] <- "delete"
  # flag all expected occurrences in mock samples as "keep", NA for tolerate, and delete for all others
  occurrence_df <- flag_from_mock(occurrence_df, mock_composition_df, sep=sep)
  # flag occurrences as delete with low read count in habitat, compared to the other habitats
  occurrence_df <- flag_from_habitat(occurrence_df, 
                                     habitat_proportion=habitat_proportion
                                     ) 
  
  # keep only relevant columns and lines, sort data
  occurrence_df <- occurrence_df %>%
    select(sample,action,asv_id,asv) %>%
    filter(!is.na(action)) %>%
    arrange(sample, action)

  
  # count the number of FP and expected TP
  FP <- nrow(occurrence_df %>%
               filter(action=="delete"))
  TP <- nrow(occurrence_df %>%
               filter(action=="keep"))
  
  # count the number of FN and write missing_occurrences, if filename is defined
  missing_occurrence_df <-make_missing_occurrences(read_count_samples=read_count_samples_df,
                                                   mock_composition=mock_composition_df, 
                                                   sep=sep, 
                                                   out=missing_occurrences
                                                   )
  FN <- nrow(missing_occurrence_df %>%
               filter(action=="keep"))
  # real TP is the expected occurrences - FN
  TP <- TP - FN
  Precision <- TP/(TP+FP)
  Sensitivity <- TP/(TP+FN)
  count_df <- data.frame("TP" = c(TP),
                         "FP" = c(FP),
                         "FN" = c(FN),
                         "Precision" = c(Precision),
                         "Sensitivity"= c(Sensitivity)
  )
  
  
  # write to outfiles (missing is written by function make_missing_occurrences)
  if(known_occurrences != ""){
    check_dir(known_occurrences, is_file=TRUE)
    write.table(occurrence_df, file=known_occurrences, row.names = F, sep=sep)
  }
  if(performance_metrics != ""){
    check_dir(performance_metrics, is_file=TRUE)
    write.table(count_df, file=performance_metrics, row.names = F, sep=sep)
  }
  df_list <- list(occurrence_df, missing_occurrence_df, count_df)
  return(df_list)
}

#' Flag occurrences in mock samples
#' 
#' Flag all occurrences in mock samples. 
#' Expected variants as 'keep', unexpected ASVs as 'delete', tolerate ASVs as NA. 
#' Tolerate ASVs are ASVs that can be in the mock, but the filtering 
#' should not be optimized to keep them 
#' (e.g. badly amplified species present in the mock).
#'  
#' @param occurrence_df Data frame with the following variables: 
#' asv_id, sample, read_count, asv, sample_type, habitat, action. 
#' @param mock_composition Data frame or csv file with columns: 
#' sample, action (keep/tolerate), asv.
#' @param sep Field separator character in input and output csv files.
#' @returns Data fram with the following columns: 
#' asv_id, sample, read_count, asv, sample_type, habitat, action.
#' @examples
#' \dontrun{
#' flag_from_mock(occurrence_df=occurrence_df, mock_composition="data/mock_composition.csv")
#' }
#' @export
#' 
flag_from_mock <- function(occurrence_df, mock_composition, sep=","){
  # can accept df or file as an input
  if(is.character(mock_composition)){
    # read known occurrences
    mock_composition_df <- read.csv(mock_composition, header=T, sep=sep)
  }else{
    mock_composition_df <- mock_composition
  }
  mock_composition_df <- mock_composition_df %>%
    rename(action_mock=action)
  
  if("asv_id" %in% colnames(mock_composition_df)){
    mock_composition_df <- mock_composition_df %>%
      select(-asv_id)
  }

   
  # add action_mock to occurrence_df; 
  # use full join, so expected ASV missing from occurrence_df will be added
  occurrence_df <- full_join(occurrence_df, mock_composition_df, by=c("sample", "asv"))
  # if expected ASV was missing from occurrence_df, complete the sample_type as mock
  occurrence_df$sample_type[
    which(is.na(occurrence_df$sample_type) & occurrence_df$action_mock=="keep")
    ] <- "mock"
  # set the action to delete, keep or tolerate in function of the mock composition
  occurrence_df$action[
    which((is.na(occurrence_df$action)) & occurrence_df$sample_type == "mock")
    ] <- "delete"
  occurrence_df$action[which(occurrence_df$action_mock =="keep")] <- "keep"
  occurrence_df$action[which(occurrence_df$action_mock =="tolerate")] <- NA
  # select original columns
  occurrence_df <- occurrence_df %>%
    select(asv_id, sample, read_count, asv, sample_type, habitat, action)
  
  return(occurrence_df)
}

#' Flag occurrences based on habitat
#' 
#' Flag False Positive occurrences based on the habitat of the samples.
#' All ASVs present in more than one habitat are checked. 
#' For each of these ASVs, if the proportion of reads in a habitat is below 
#' a cutoff (`habitat_proportion`), 
#' it is considered as an artifact in all samples of the habitat.
#'  
#' @param occurrence_df Data frame with the following variables: 
#' asv, plate, marker, sample, read_count, sample_type, habitat, action.
#' @param habitat_proportion For each of ASVs, if the proportion 
#' of reads in a habitat is below this cutoff, it is considered as an 
#' artifact in all samples of the habitat.
#' @returns The updated input data frame, with some occurrences flagged as 
#' delete in the action column.
#' @examples
#' \dontrun{
#' flag_from_habitat(occurrence_df, habitat_proportion=0.7)
#' }
#' @export
#' 
flag_from_habitat <- function(occurrence_df, habitat_proportion=0.5){
  
  # group by asv and habitat and count the total number of reads for 
  # each habitat-asv combination
  tmp <- occurrence_df %>%
    group_by(habitat, asv) %>%
    summarize(habitat_read_count=sum(read_count), .groups="drop_last") %>%
    filter(!is.na(habitat)) %>%
    ungroup()
  
  # count the number of habitats for each asv 
  # and keep only the ones present in at least two different habitats
  tmp2 <- tmp %>%
    group_by(asv) %>%
    summarize(nb_habitat=length(asv)) %>%
    filter(nb_habitat>1) %>%
    ungroup()
  # keep only selected asvs in tmp
  tmp <- tmp[tmp$asv %in% tmp2$asv, ]
  # get the total readcount for each asv in tmp
  tmp3 <- tmp %>%
    group_by(asv) %>%
    summarize(sum_read_count = sum(habitat_read_count)) %>%
    ungroup()
  # add total readcount of asv to tmp 
  # keep only lines where the asv in a given habitat is less frequent than in the other habitats
  tmp <- left_join(tmp, tmp3, by="asv")
  tmp <- tmp[tmp$habitat_read_count/tmp$sum_read_count < habitat_proportion, ]
  # keep only pertinent columns in tmp and add hab_action column with "delete"
  tmp <- tmp %>%
    select(habitat, asv)
  tmp$hab_action <- rep("delete", nrow(tmp))
  
  occurrence_df <- left_join(occurrence_df, tmp, by=c("habitat", "asv"))
  occurrence_df$action[which(occurrence_df$hab_action=="delete")] <- "delete"
  
  occurrence_df <- occurrence_df %>%
    select(-hab_action)
  
  return(occurrence_df)
}

#' Make Missing Occurrences
#' 
#' Prepare a data fram that list all expected occurrences that are missing 
#' (False negatives).
#'  
#' @param read_count_samples Data frame or csv file with the following variables: 
#' asv, plate, marker, sample, read_count.
#' @param mock_composition Data frame or csv file with columns: 
#' sample, action (keep/tolerate), asv.
#' @param sep Field separator character in input and output csv files.
#' @param out Character string: output file. If empty, no file is written. 
#' @returns Data frame with the following columns: 
#' sample, action, asv, asv_id
#' @examples
#' \dontrun{
#' make_missing_occurrences(read_count_samples=read_count_samples_df, 
#'     mock_composition="data/mock_composition.csv"
#'     )
#' }
#' @export
#'
make_missing_occurrences <- function(read_count_samples, mock_composition, sep=",", out=""){
  
  # can accept df or file as an input
  if(is.character(mock_composition)){
    # read known occurrences
    mock_comp <- read.csv(mock_composition, header=T, sep=sep)
  }else{
    mock_comp <- mock_composition
  }
  
  # read mock composition to a df
  mock_comp <- mock_comp %>%
    filter(action=="keep")
    if("asv_id" %in% colnames(mock_comp)){
      mock_comp <- mock_comp %>%
        select(-asv_id)
    }
  
  # can accept df or file as an input
  if(is.character(read_count_samples)){
    # read known occurrences
    read_count_samples_df <- read.csv(read_count_samples, header=T, sep=sep)
  }else{
    read_count_samples_df <- read_count_samples
  }
  
  # add read_count to df from read_count_samples_df, and keep only if value is NA
  df <- left_join(mock_comp, read_count_samples_df,  by=c("sample", "asv")) %>%
    filter(is.na(read_count)) %>%
    select(-"read_count")
#    select(sample, action, asv_id, asv)
  
  
  # write to outfile
  if(out != ""){
    check_dir(out, is_file=TRUE)
    write.table(df, file=out, row.names = F, sep=sep)
  }
  
#  FN <- nrow(df %>%
#               filter(action=="keep"))
  return(df)
}


#' OptimizeLFNreadCountLFNvariant
#' 
#' Suggest optimal parameters for `lfn_read_count_cutoff` and `lnf_variant_cutoff` 
#' functions.
#' The the `LFNreadCount` and `LFNvariant` is run for a series of parameter value 
#' combinations followed by `FilterMinReplicate`. 
#' For each parameter combination, the number of FN, TP, and FP is reported. 
#' Users should chose the parameter setting that minimizes, FN and FP.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param known_occurrences Data frame or file produced by 
#' `MakeKnownOccurrences` function, with known False Positives and True Positives.
#' @param sep Field separator character in input and output csv files.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param min_lfn_read_count_cutoff Positive integer: the lowest cutoff value for 
#' `LFNreadCount` function. 
#' @param max_lfn_read_count_cutoff Positive integer: the highest cutoff value 
#' for `LFNreadCount` function. 
#' @param increment_lfn_read_count_cutoff Positive integer: values from 
#' `min_lfn_read_count_cutoff` to `max_lfn_read_count_cutoff` 
#' are tested by `increment_lfn_read_count_cutof` of increment. 
#' @param min_lnf_variant_cutoff Real (0-1): the lowest cutoff value for 
#' `LFNvariant` function. 
#' @param max_lnf_variant_cutoff  Real (0-1): the highest value for `LFNvariant`
#'  function.
#' @param increment_lnf_variant_cutoff  Real (0-1): values from 
#' `min_lnf_variant_cutoff` to `max_lnf_variant_cutoff` are tested by 
#' `increment_lnf_variant_cutoff` increment. 
#' @param by_replicate logical: : compare read count of the occurrence to the 
#' read counts of the ASV-replicate in `LFNvariant` function.
#' @param min_replicate_number Positive integer: minimum number of replicates 
#' for `FilterMinReplicate`.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns data frame with the following columns: 
#' lfn_read_count_cutoff, lnf_variant_cutoff, FN, TP, FP.
#' @examples
#' \dontrun{
#' OptimizeLFNreadCountLFNvariant(read_count_df, 
#'     known_occurrences="data/known_occurrences.csv", 
#'     min_lfn_read_count_cutoff=10,
#'      max_lfn_read_count_cutoff=50, 
#'      increment_lfn_read_count_cutoff=10, 
#'      min_lnf_variant_cutoff=0.001, 
#'      max_lnf_variant_cutoff=0.005, 
#'      increment_lnf_variant_cutoff=0.001
#'      )
#' }
#' @export
#'

OptimizeLFNreadCountLFNvariant <- function(read_count, 
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
                                           ){
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
  CheckFileinfo(file=known_occurrences_df, 
                file_type="known_occurrences", 
                sep=sep, 
                quiet=TRUE
                )

  # make a series of cutoff values for LFNreadCount
  rc_cutoff_list <- seq(from=min_lfn_read_count_cutoff, 
                        to=max_lfn_read_count_cutoff, 
                        by=increment_lfn_read_count_cutoff
                        )
  # make a series of cutoff values for LFNreadCount
  var_cutoff_list <- seq(from=min_lnf_variant_cutoff, 
                         to=max_lnf_variant_cutoff, 
                         by=increment_lnf_variant_cutoff
                         )
  
  out_df <- data.frame(
    lfn_read_count_cutoff=numeric(),
    lnf_variant_cutoff=numeric(),
    FN=numeric(),
    TP=numeric(),
    FP=numeric()
  )
  # go through all parameter combination and count the number of TP and FN
  
  for(rc_cutoff in rc_cutoff_list){
    df_tmp <- read_count_df
    #LFNreadCount
    df_tmp <- LFNreadCount(df_tmp, rc_cutoff)
    for(var_cutoff in var_cutoff_list){
      # LFNvariant
      df_tmp <- LFNvariant(df_tmp, var_cutoff, by_replicate=by_replicate)
      # FilterMinReplicate
      df_tmp <- FilterMinReplicate(df_tmp, min_replicate_number)
      # PoolReplicates
      df_tmp_sample <- PoolReplicates(df_tmp, digits=0)
      # pool readcount info and known occurrences info
      ko <- full_join(df_tmp_sample, known_occurrences_df, by=c("sample", "asv")) %>%
        filter(!is.na(action)) %>% # keep only lines mentioned in the known occurrences
        # delete lines if asv is not present (read_count==NA) and the action is delete
        filter(!(is.na(read_count) & action=="delete")) 
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
      new_line <- data.frame(lfn_read_count_cutoff=rc_cutoff, 
                             lnf_variant_cutoff=var_cutoff ,
                             FN=FN_count, 
                             TP=TP_count, 
                             FP=FP_count
                             )
      if(!quiet){
        print(new_line)
      }
      out_df <- bind_rows(out_df, new_line )
    }
  }
  
  out_df <- out_df %>%
    arrange(FN, FP, lnf_variant_cutoff, lfn_read_count_cutoff)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(out_df, file=outfile, sep=sep, row.names = F)
  }
  return(out_df)
}


#' Pool Datasets
#' 
#' Take several input files, each in long format containing 
#' asv_id, sample, read_count and asv columns.
#' Pool the different data sets, if all have the same marker.
#' If more than one markers, ASVs identical on their overlapping 
#' regions are pooled into groups, and different ASVs of the same group 
#' are pooled under the centroid (longest ASV of the group).
#' Pooling can take the mean of the read counts of the ASV (default) or their sum.
#' 
#'  
#' @param files Data frame with the following variables: file (name of input files), marker.
#' Input files must have asv, sample and read_count columns.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param asv_with_centroids Character string: csv file name of the output file 
#' containing the same information as the concatenated input files, 
#' completed by centroid_id and centroid columns.
#' @param sep Field separator character in input and output csv files.
#' @param mean_over_markers logical: If TRUE, the mean read count is calculated 
#' over different ASVs of each cluster. Sum otherwise.
#' @param vsearch_path Character string: path to vsearch executables. 
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with asv_id, sample, read_count, asv columns.
#' @examples
#' \dontrun{
#' files <- data.frame(file=c("vtamR_test/out_mfzr/14_PoolReplicates.csv", 
#'     "vtamR_test/out_zfzr/14_PoolReplicates.csv"),
#'     marker=c("MFZR", "ZFZR"))
#' PoolDatasets(files, vsearch_path=vsearch_path)
#' }
#' @export
#'
PoolDatasets <- function(files, 
                         outfile="", 
                         asv_with_centroids="", 
                         sep=",", 
                         mean_over_markers=T, 
                         vsearch_path="vsearch", 
                         quiet=T
                         ){
  
  tmp_dir <-paste('tmp_pool_datasets_', 
                  trunc(as.numeric(Sys.time())), 
                  sample(1:100, 1), 
                  sep='')
  tmp_dir <- file.path(tempdir(), tmp_dir)
  check_dir(tmp_dir)
  
  ###
  # pool all data into one data frame (df), 
  # check if the all marker.sample combinations are unique among different data sets
  ###
  df <- data.frame("asv_id" = list(), 
                   "sample" = list(), 
                   "read_count" = list(), 
                   "asv"=list(), 
                   "marker"== list()
                   )
  samples <- c()
  for(i in 1:nrow(files)){
    file <- files[i, "file"]
    marker <- files[i, "marker"]
    #    print(file)
    tmp <- read.csv(file, sep=sep) %>%
      select(asv_id, sample, read_count, asv)
    tmp$marker <- rep(marker, nrow(tmp)) # add maker
    
    # make a list of marker.sample of the data set that just have been read
    local_samples <- unique(paste(tmp$marker, tmp$sample, sep="."))
    # see if they match earlier read marker.sample combinations
    shared_samples <- intersect(local_samples, samples) 
    if(length(shared_samples) > 0){
      print("The following samples are in at least 2 different data sets of the same marker. 
            Their read_counts will be summed. Use unique names if you want to keep them separate:")
      print(shared_samples)
    }
    samples <- c(samples, local_samples) 
    
    # add data set to df
    df <- rbind(df, tmp)
  }
  
  ###
  # Pool ASVs identical on their overlapping region, if more than one marker
  ###
  marker_list <- unique(df$marker)
  # more than one marker => pool sequences identical in their corresponding region
  if(length(marker_list) > 1){ 
    # get full list of ASVs
    asvs <- df %>%
      group_by(asv_id, asv) %>%
      summarize("rc" = sum(read_count), .groups="drop_last")  %>%
      ungroup()
    
    # arrange ASVs by decreasing sequence length and then by read_count
    asvs$length <- as.numeric(nchar(asvs$asv))
    asvs <- asvs %>%
      arrange(desc(length), desc(rc))
    
    # make a fasta file
    fasta <- file.path(tmp_dir, "vsearch_input.fasta")
    writeLines(paste(">", asvs$asv_id, "\n", asvs$asv, sep="" ), fasta)
    
    # cluster using cluster_smallmem and 1 as identity limit
    centroids_file <- file.path(tmp_dir, "consout.txt")
    #query sequences are shorter than subjects => centroids are in the subjects column
    blastout_file <- file.path(tmp_dir, "blastout.tsv")  
    vsearch_cmd <- paste(vsearch_path, " --cluster_smallmem ", 
                         fasta, " --consout ",centroids_file,
                         " --blast6out ", blastout_file,
                         " --id 1", 
                         sep=""
                         )
    if(!quiet){
      print(vsearch_cmd)
    }
    system(vsearch_cmd)
    
    ###
    # Make cent data frame with a complete list of ASVs and the centroïd for each of them.
    ###
    # read the ids of centoids, and get the list of centroids
    cent <- read.table(centroids_file)
    colnames(cent) <- c("centroid_id")
    cent <- cent %>%
      filter(grepl(">centroid=", centroid_id))
    cent$centroid_id <- gsub(">centroid=", "", cent$centroid_id)
    cent$nbseq <-   gsub(".+;seqs=", "", cent$centroid_id)
    cent$centroid_id <- gsub(";.+", "", cent$centroid_id)
    cent$centroid_id <- as.numeric(cent$centroid_id)
    cent$nbseq <- as.numeric(cent$nbseq)
    
    # add to centroid the asv_id that are in the same cluster
    blastout <- read.table(blastout_file) %>%
      select(1,2)
    colnames(blastout) <- c("asv_id", "centroid_id")
    cent <- left_join(cent, blastout, by= c("centroid_id"))
    # add centroid_id to asv_id column for singletons
    cent <- cent %>%
      mutate(asv_id = ifelse(is.na(asv_id), centroid_id, asv_id))
    # add a line for each non-singleton centroid, 
    # with centroid id in both the centroid and in query columns
    added_lines <- cent %>%
      filter(nbseq>1) %>%
      mutate(asv_id=centroid_id) %>%
      unique # add just one line per centroid, not several if many sequences in cluster
    cent<- rbind(cent, added_lines) %>%
      arrange(centroid_id)
    
    ###
    # Pool ASVs of the same cluster
    ###
    # add the centroid_id to each asv if df
    df <- left_join(df, cent, by=c("asv_id")) %>%
      select(-nbseq)
    # add the centroid to each centroid_id in df
    df <- left_join(df, asvs, by=c("centroid_id"="asv_id")) %>%
      select(-length, -rc) %>%
      rename("asv"=asv.x, "centroid"=asv.y) %>%
      select(centroid_id,asv_id,marker,sample,read_count,asv,centroid) %>%
      arrange(centroid_id, marker)
    
    if(mean_over_markers){
      df_pool <- df %>%
        group_by(centroid_id, sample) %>%
        summarize("read_count"=round(mean(read_count), digits=0), .groups =  "keep") %>%
        ungroup()
    }else{
        df_pool <- df %>%
          group_by(centroid_id, sample) %>%
          summarize("read_count"=sum(read_count), .groups =  "keep" )  %>%
          ungroup()
    }
    # add asv column and select columns
    # df_pool is a simple output with the format identical to the read_count_sample dfs,
    # but no info on the as that has been pooled together
    df_pool <- left_join(df_pool, asvs, by=c("centroid_id" = "asv_id")) %>%
      select("asv_id"=centroid_id, sample, read_count, asv)
    
    if(asv_with_centroids != ""){
      write.table(df, file=asv_with_centroids, sep=sep, row.names = F)
    }
  }else{# one marker
    df_pool <- df %>%
      select(asv_id, sample, read_count, asv)
  }
  
  unlink(tmp_dir, recursive = TRUE)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(df_pool, file=outfile, sep=sep, row.names = F)
  }

  return(df_pool)
}

#' History By
#' 
#' Filters a feature (asv_id/asv/sample/replicate), 
#' in all output files of intermediate filtering steps to retain only 
#' lines corresponding a value.
#'  
#' @param dir Character string: directory containing the output files 
#' of intermediate filtering steps. File names should start with a number 
#' followed by an underscore (e.g. 5_LFNsampleReplicate.csv).
#' @param feature Character string with the following possible values: 
#' "asv_id", "asv", "sample", "replicate", "read_count".
#' @param value Character string: values of feature that should selected. 
#' The output data frame contains all lines where value is present 
#' in feature in the input files.
#' @param sep Field separator character in input and output csv files.
#' @returns Data frame with the selected lines of the input files.
#' @examples
#' \dontrun{
#' HistoryBy(dir="out", feature="asv_id", value=1)
#' HistoryBy(dir="out", feature="sample", value="tpos1")
#' }
#' @export
#'
HistoryBy <- function(dir, feature, value, sep=","){

  check_dir(dir)
  files <- list.files(path=dir, pattern="^[0-9]+", full.names=FALSE)
  
  # get filenames to df and arrange the according to the number at the beginning of the filename
  df <- data.frame("files"= files)
  df$order <- gsub("_.*$", "", df$files)
  df$order <- as.numeric(df$order)
  df <- df %>%
    arrange(order)
  
  selected_lines <- data.frame(
    file= as.character(), 
    asv_id =as.numeric(),
    sample = as.character(), 
    replicate= as.character(), 
    read_count=as.numeric(),
    asv= as.character())
  
  for(i in 1:length(files)){
    file <- file.path(dir, df$files[i])
    
    data <- read.csv(file, sep=sep)
    data$file <- rep(df$files[i], nrow(data)) # add filename
    # add columns that can be missing eventually
    if(!"replicate" %in% colnames(data)){ 
      data$replicate <- rep(NA, nrow(data))
    }
    if(!"asv_id" %in% colnames(data)){
      data$asv_id <- rep(NA, nrow(data))
    }
    # check if the feature is in among the columns names of the data frame
    if(feature %in% colnames(data)){
      tmp <- data %>%
        filter(!!sym(feature)==value) %>% # filter using a symbol from feature
        select(file, asv_id, sample, replicate, read_count, asv)
      
      selected_lines <- rbind(selected_lines, tmp)
    }else{
      stop("ERROR: feature is not in file")
    }
  }
  return(selected_lines)
}

#' SummarizeBy
#' 
#' Summarizes the output of each intermediate filtering steps. 
#' Takes all files in the directory, that start with a number.
#' For each file, groups the lines by `grouped_by` variable, 
#' and gets the number of distinct values of `feature` for 
#' each group, or if `feature` is `read_count`, gets the total 
#' number of reads of each group.
#'  
#' @param dir Character string: directory containing the output files 
#' of intermediate filtering steps. File names should start with a number 
#' followed by an underscore (e.g. 5_LFNsampleReplicate.csv).
#' @param feature Character string with the following possible values: 
#' "asv_id", "asv", "sample", "replicate", "read_count".
#' @param grouped_by Character string with the following possible values: 
#' "asv_id", "asv", "sample", "replicate". Group data by this variable.
#' @param sep Field separator character in input and output csv files.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @returns Data frame with input file names in columns, `grouped_by` values in 
#' lines, and the number of features or the read counts in the cells.
#' @examples
#' \dontrun{
#' SummarizeBy(dir="vtamR_test/out_mfzr", feature="asv", grouped_by="sample")
#' SummarizeBy(dir="vtamR_test/out_mfzr", feature="read_count", grouped_by="sample")
#' }
#' @export
#'
SummarizeBy <- function(dir, feature, grouped_by, outfile="", sep=","){
  
  # read file names in dir
  check_dir(dir)
  files <- list.files(path=dir, pattern="^[0-9]+", full.names=FALSE)
  
  # get filenames to file_df and arrange the according to the number 
  # at the beginning of the file name
  file_df <- data.frame("files"= files)
  file_df$order <- gsub("_.*$", "", file_df$files)
  file_df$order <- as.numeric(file_df$order)
  file_df <- file_df %>%
    arrange(order)
  
  # define empty data frame
  df <- data.frame(
    "grouped_by"=character(),
    "count"=character(),
    "step"=character()
  )
  
  for(i in 1:length(files)){ # for each file
    # read file
    file <- file.path(dir, file_df$files[i])
    filename <- file_df$files[i]
    filename <- gsub("\\..*$", "", filename)
    tmp <- read.csv(file, sep=sep)
    
    if(grouped_by %in% colnames(tmp)){ # grouping variable is present in the file
      # feature variable is in the file
      if(feature %in% colnames(tmp)){
        if(feature == "read_count"){ # if feature is read_count, 
          # it should be summed instead of get the number of distinct values
          tmp <- tmp %>%
            group_by(!!sym(grouped_by)) %>%
            summarize(count = sum(!!sym(feature)))
        }else{ #  get the number of distinct values
          tmp <- tmp %>%
            group_by(!!sym(grouped_by)) %>%
            summarize(count = n_distinct(!!sym(feature)))
        }
      } else{ # feature variable is NOT in the file
        tmp <- tmp %>%
          select(!!sym(grouped_by)) %>%
          distinct()
        tmp$count <- NA
      }
      tmp$step <- rep(filename, nrow(tmp)) # add filename
      df <- rbind(df, tmp)
    }else{ # grouping variable is not present => go to next file
      msg <- paste("WARNING:",grouped_by, "variable is not present in", file, sep=" ")
      print(msg)
      next()
    }
  }
  # mkae wide format
  wide_df <- as.data.frame(pivot_wider(df, 
                                       names_from = c(step), 
                                       values_from = count, 
                                       values_fill=0, 
                                       names_sep = ".", 
                                       names_sort=F
                                       )
                           )
  # print outfile
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(wide_df, file=outfile, row.names = F, sep=sep)
  }
  return(wide_df)
}

#' Write data frame to fasta
#' 
#' Writes input data frame to a fasta file. Can produce gz compressed file 
#' or uncompressed files.
#' The output file name's extension is corrected according to the 
#' compression if necessary.
#'  
#' @param df input data frame with the following columns: header, sequence.
#' @param out Character string: name of the output file.
#' @param compress logical: Compress output files to gzip format.
#' @returns Character string: name of the output file updated according to the compression.
#' @examples
#' \dontrun{
#' df <- data.frame(header=c("seq1", "seq2"),
#'     sequence=c("AACTTGTTGTCACTGTAAACTGATGTA", "AACTTGTTGTCACTGTTTGACTGATGTA")
#'     )
#' write_df_to_fasta(df, out="out/test.fasta", compress=T)
#' }
#' @export
#'
write_df_to_fasta <- function(df, out, compress=F){
  
  if(compress){
    if(!endsWith(out, ".gz")){
      out <- paste(out, ".gz", sep="")
    }
    file_connection <- gzfile(out, "w")
  }else{
    if(endsWith(out, ".gz")){
      out <- sub(".gz", "", out)
    }
    file_connection <- file(out, "w")
  }
  
  df$header <- paste('>', df$header, sep="")
  writeLines(paste(df$header, df$sequence, sep="\n"), con=file_connection, sep="\n")
  close(file_connection)
  
  return(out)
}


#' Select sequences randomly
#' 
#' Selects `n` random sequences from the input fasta file.
#'  
#' @param file Character string: name of the input fasta file.
#' Can be  uncompressed or compressed in gz format (zip files are not supported).
#' @param n Positive integer: the number of sequences to be selected.
#' @param randseed Positive integer: seed for random sampling.
#'  0 by default means to use a pseudo-random seed. 
#'  A given non-zero seed produces always the same result.
#' @returns Data frame with two columns: headers and sequences.
#' @examples
#' \dontrun{
#' select_sequences(file="data/test.fasta", n=100, randseed=563)
#' }
#' @export
#'
select_sequences <- function(file, n=100, randseed=0){
  
  # read file to a df, with headers in one column and sequences in another
  seq_df <- read_fasta_to_df(file, dereplicate=FALSE)
  seq_n <- nrow(seq_df)
  if(seq_n > n){ # enough sequences
    if(randseed == 0){
      set.seed(Sys.time())
    }else{
      set.seed(randseed)
    }
    random_integers <- as.data.frame(sample(1:seq_n, size = n, replace = FALSE))
    colnames(random_integers) <- c("random_ind")
    seq_df$ind <- as.numeric(rownames(seq_df))
    # select lines that correspond to random numbers
    seq_df <- left_join(random_integers, seq_df, by=c("random_ind"= "ind")) %>%
      select(-random_ind)
    
  }else{
    msg <- paste("Only", seq_n, "sequences in ",file, 
                 "All sequences will used in subsequent steps", sep=" ")
    print(msg)
  }
  
  return(seq_df)
}


#' Read fasta file to data frame
#' 
#' Read a fasta file to a data frame 
#'  
#' @param file Character string: name of the input fasta file.
#' Can be  uncompressed or compressed in gz format (zip files are not supported).
#' @param dereplicate logical: If TRUE returns a data frame with asv and 
#' read_count columns. If FALSE data frame with header and sequence columns.
#' @returns Data frame with two columns: 
#' * headers and sequences if dereplicate==FALSE
#' * asv and read_count if dereplicate==TRUE
#' @examples
#' \dontrun{
#' read_fasta_to_df(file="data/test.fasta", dereplicate=F)
#' read_fasta_to_df(file="data/test.fasta", dereplicate=T)
#' }
#' @export
#'
read_fasta_to_df <- function(file, dereplicate=F){
  
  ### can deal with uncompressed files and gz compressed files. 
  # Zip files should be decompressed previously
  if(endsWith(file, ".gz")){
    file_connection <- gzfile(file, "rb") 
  }else{
    file_connection <- file(file, "r")
  }
  
  # read file to a vector. Ech element is a line
  file_contents <- readLines(file_connection, warn = FALSE)
  close(file_connection)
  
  # Identify lines starting with '>'
  header_indices <- grepl("^>", file_contents)
  # Use cumulative sum to create groups for each header
  group_indices <- cumsum(header_indices)
  # Split the file_contents into groups based on header indices
  grouped_lines <- split(file_contents, group_indices)
  rm(group_indices)
  
  # Create a data frame (columns: header, sequence)
  df <- data.frame(
    header = gsub("^>", "", file_contents[header_indices]), # Remove '>' from headers
    sequence = sapply(
      grouped_lines, function(x) if(length(x) > 1) paste(x[-1], collapse = "") else ""
      ),
    stringsAsFactors = FALSE
  )
  
  if(dereplicate){ #(columns: asv, read_count)
    df <- df %>%
      group_by(sequence) %>%
      summarize(read_count=n()) %>%
      select("asv"=sequence, read_count)
  }
  
  return(df)
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


#' Count reads
#' 
#' Count the number of sequences in a fasta or fastq file
#' or the number of lines in other files.
#' 
#' Efficient in linux-like systems, but slow on Windows for large files.
#' Can handle gz compressed and uncompressed files, but not zip files.
#'  
#' @param file Character string: name of the input file (including path).
#' @param file_type Character string with the following values: "fasta", "fastq".
#' For all other values, the number of lines is returned.
#' @returns Integer: The number of sequences for fasta or fastq files, 
#' number of lines for other files.
#' @examples
#' \dontrun{
#' count_reads_file(file="data/test.fasta", file_type="fasta")
#' count_reads_file(file="data/test.fasta", file_type="")
#' }
#' @export
#' 

count_reads_file <- function(file, file_type=""){
  
  if(endsWith(file, ".zip")){
    msg <- "File compression type is not supported"
    print(msg)
    return(0)
  }
  
  if(is_linux()){
    # compressed files
    if(endsWith(file, ".gz") || endsWith(file, ".bz") || endsWith(file, ".gz2")){
      if(file_type == "fastq"){
        cmd <- paste("zcat ", file, "| wc -l ", sep=" ")
        seq_count <- as.integer(system(cmd, intern=TRUE))
        seq_count <- seq_count/4
      }else if(file_type == "fasta"){
        cmd <- paste("zcat ", file, "| grep '^>' -P | wc -l", sep=" ")
        seq_count <- as.integer(system(cmd, intern=TRUE))
      }else{
        msg <- paste(file_type, "is neither fasta nor fastq. 
                     The number of liens in file will be returned for", file)
        print(msg)
        cmd <- paste("zcat ", file, "| wc -l ", sep=" ")
        seq_count <- as.integer(system(cmd, intern=TRUE))
      }
    }else{
      #uncompressed files
      if(file_type == "fastq"){
        cmd <- paste("wc -l", line, sep=" ")
        seq_count <- as.integer(system(cmd, intern=TRUE))
        seq_count <- seq_count/4
      }else if(file_type == "fasta"){
        cmd <- paste("grep '^>' -P", file, "| wc -l", sep=" ")
        seq_count <- as.integer(system(cmd, intern=TRUE))
      }else{
        msg <- paste(file_type, "is neither fasta nor fastq. 
                     The number of lines in file will be returned for", file)
        print(msg)
        cmd <- paste("wc -l", line, sep=" ")
        seq_count <- as.integer(system(cmd, intern=TRUE))
      }
    }
    return(seq_count)
  }else{
    print("WARNING: This command on non linux-like systems is slow 
          and might not work with very large files.")
    
    if(file_type == "fasta"){ # can deal with compressed and uncompressed files
      df <- read_fasta_to_df(file, dereplicate=F)
      seq_count <- nrow(df)
    }else { # fastq and others
      if(endsWith(file, ".gz") || endsWith(file, ".bz") || endsWith(file, ".gz2")){
        file_connection <- gzfile(file, "rb")
      }else{
        file_connection <- file(file, "r")
      }
      data <- readLines(file_connection, n = -1)
      close(file_connection)
      seq_count <- length(data)
      if(file_type == "fastq"){
        seq_count <- seq_count / 4
      }else{
        msg <- paste(file_type, "is neither fasta nor fastq. 
                     The number of lines in file will be returned for", file)
        print(msg)
      }
      
    }
    return(seq_count)
  } # end non-linux-like
  return(0)
}

#' Count reads in directory
#' 
#' Counts the number of sequences in fasta or fastq files, or the number of 
#' lines in other files. Reads all files in the input directory with 
#' `pattren` in their names.
#' 
#' Efficient in linux-like systems, but slow on Windows for large files.
#' Can handle gz compressed and uncompressed files, but not zip files.
#'  
#' @param dir Character string: name of the input directory.
#' @param pattern Regular expression: pattern in the name of the files in the 
#' input directory. Read only files with pattern in their names.
#' @param file_type Character string with the following values: "fasta", "fastq".
#' For all other values, the number of lines is returned.
#' @param sep Field separator character in input and output csv files.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with 2 columns: filename, read_count.
#' @examples
#' \dontrun{
#' CountReadsDir(dir="out", pattern=".fastq", file_type="fastq")
#' CountReadsDir(dir="out", pattern="^mfzr", file_type="fasta")
#' }
#' @export
#' 
CountReadsDir<- function(dir, 
                         pattern=".", 
                         file_type="fasta", 
                         outfile="", 
                         sep=",", 
                         quiet=T
                         ){
  
  check_dir(dir)
  files <- list.files(path = dir, pattern=pattern)
  df <- data.frame(
    "filename"=files,
    "read_count"=rep(NA, length(files))
  )
  
  for(i in 1:length(files)){
    file_p <- file.path(dir, files[i])
    if(!quiet){
      print(file_p)
    }
    n <- count_reads_file(file_p, file_type=file_type)
    df[i, "read_count"] <- n
  }
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(df, file=outfile, sep=sep, row.names = F)
  }
  return(df)
}

#' Check File formats
#' 
#' Check the format and coherence of different file types:
#' * Check if all obligatory columns are present (all file_type)
#' * Check if all sample names are alpha-numerical
#' * Check if asv, tag_fw, tag_rv primer_fw, primer_rv have only IUPAC codes
#' * Check if read_count is numerical
#' * Check sample type, habitat homogeneity across replicates 
#' (fastqinfo, fastainfo, sortedinfo)
#' * Check if fastq file pairs are coherent (e.g. 1 to 1 relation; fastqinfo)
#  * Check if files in fastq_fw, fastq_rv, fasta columns exist (fastqinfo, fastainfo, sortedinfo)
#' * Check if tag combinations are unique within a file(pair)s (fastqinfo, fastainfo)
#' * Check action (mock_composition, known_occurrences)
#' * Check if 1 to 1 relation between asv_id ad asv (read_count,read_count_sample,asv_list)
#' 
#'  
#' @param file Character string: name of the input file to be checked or 
#' data frame with the content of the input file.
#' @param dir Character string: Name of the directory containing the files 
#' in `fastq_fw`, `fastq_rv` and `fasta` columns in the input file.
#' @param file_type Character string with the following possible values: 
#' "fastqinfo", "fastainfo", "sortedinfo", "mock_composition", "known_occurrences", 
#' "read_count", "read_count_sample", "asv_list".
#' @param sep Field separator character in input and output csv files.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Error message if problem with the files and stop the script.
#' @examples
#' \dontrun{
#' CheckFileinfo(file="input/sortedinfo.csv", dir="fasta", file_type="sortedinfo")
#' CheckFileinfo(file=sortedinfo_df, dir="fasta", file_type="read_count_sample")
#' }
#' @export
#' 
CheckFileinfo <- function(file, dir="", file_type="fastqinfo", sep=",", quiet=FALSE){
  
  if(is.character(file)){
    # read known occurrences
    df <- read.csv(file, header=T, sep=sep)
  }else{
    df <- file
  }
  
  # define expected columns
  if(file_type == "fastqinfo"){
    column_heading <- c("tag_fw","primer_fw","tag_rv","primer_rv",
                        "sample","sample_type","habitat","replicate","fastq_fw","fastq_rv")
  }else if(file_type == "fastainfo"){
    column_heading <- c("tag_fw","primer_fw","tag_rv","primer_rv",
                        "sample","sample_type","habitat","replicate","fasta")
  }else if(file_type == "sortedinfo"){
    column_heading <- c("sample","sample_type","habitat","replicate","fasta")
  }else if(file_type == "mock_composition"){
    column_heading <- c("sample","action","asv")
  }else if(file_type == "known_occurrences"){
    column_heading <- c("sample","action","asv")
  }else if(file_type == "read_count"){
    column_heading <- c("asv","asv_id","sample","replicate","read_count")
  }else if(file_type == "read_count_sample"){
    column_heading <- c("asv","asv_id","sample","read_count")
  }else if(file_type == "asv_list"){
    column_heading <- c("asv","asv_id")
  }
  
  # check if all essential columns are present; 
  # File is used only to print its name in case of pb it the columns
  check_heading(column_heading, colnames(df), file=file)
  
  # check if all sample names are alphanumerical
  if("sample" %in% colnames(df)){
    tmp <- grep("[^A-z0-9_]", df$sample, perl = TRUE)
    if(length(tmp) > 0){
      msg <- paste("The following sample names contain non-alphanumerical characters:", paste(df$sample[tmp], collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Sample names are alphanumerical : OK")
      print(msg)
    }
  }
  
  # check if all characters in asv correspond to IUPAC nucleotide code 
  if("asv" %in% colnames(df)){
    tmp <- grep("[^ACGTRYSWKMBDHVN]", df$asv, perl = TRUE, ignore.case = TRUE)
    if(length(tmp) > 0){
      msg <- paste("The following ASV contain non-IUPAC characters:", paste(df$asv[tmp], collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Only IUPAC characters in ASV: OK")
      print(msg)
    }
  }
  
  # check if all characters in tag_fw correspond to IUPAC nucleotide code 
  if("tag_fw" %in% colnames(df)){
    tmp <- grep("[^ACGTRYSWKMBDHVN]", df$tag_fw, perl = TRUE, ignore.case = TRUE)
    if(length(tmp) > 0){
      msg <- paste("The following tag_fw contain non-IUPAC characters:", paste(df$tag_fw[tmp], collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Only IUPAC characters in tag_fw: OK")
      print(msg)
    }
  }
  
  # check if all characters in tag_rv correspond to IUPAC nucleotide code 
  if("tag_rv" %in% colnames(df)){
    tmp <- grep("[^ACGTRYSWKMBDHVN]", df$tag_rv, perl = TRUE, ignore.case = TRUE)
    if(length(tmp) > 0){
      msg <- paste("The following tag_rv contain non-IUPAC characters:", paste(df$tag_rv[tmp], collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Only IUPAC characters in tag_rv: OK")
      print(msg)
    }
  }
  
  # check if all characters in primer_fw correspond to IUPAC nucleotide code 
  if("primer_fw" %in% colnames(df)){
    tmp <- grep("[^ACGTRYSWKMBDHVN]", df$primer_fw, perl = TRUE, ignore.case = TRUE)
    if(length(tmp) > 0){
      msg <- paste("The following primer_fw contain non-IUPAC characters:", paste(df$primer_fw[tmp], collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Only IUPAC characters in primer_fw: OK")
      print(msg)
    }
  }
  
  # check if all characters in primer_rv correspond to IUPAC nucleotide code 
  if("primer_rv" %in% colnames(df)){
    tmp <- grep("[^ACGTRYSWKMBDHVN]", df$primer_rv, perl = TRUE, ignore.case = TRUE)
    if(length(tmp) > 0){
      msg <- paste("The following primer_rv contain non-IUPAC characters:", paste(df$primer_rv[tmp], collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Only IUPAC characters in primer_rv: OK")
      print(msg)
    }
  }
  
  # check if all read_count are numerical values 
  if("read_count" %in% colnames(df)){
    
    num <- is.numeric(df$read_count)
    na <- all(is.na(df$read_count))
    
    if(!num && !na){
      msg <- "Values should be numerical in read_count"
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
        msg <- paste("read_count is numerical: OK")
        print(msg)
    }
  }
  
  # Check sample type, habitat homogeneity across replicates
  if(file_type == "fastqinfo" || file_type == "fastainfo" || file_type == "sortedinfo" ){
    #sample_type
    tmp <- df %>%
      select("sample","sample_type","replicate") %>%
      group_by(sample) %>%
      summarise(same_sample_type = n_distinct(sample_type)) %>%
      filter(same_sample_type > 1)
    
    if(nrow(tmp) > 0){
      msg <- paste("Samples with inconsistent sample_type:", paste(tmp$sample, collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Coherence between samples and sample_type: OK")
      print(msg)
    }
    
    #habitat
    tmp <- df %>%
      select("sample","habitat","replicate") %>%
      group_by(sample) %>%
      summarise(same_sample_type = n_distinct(habitat)) %>%
      filter(same_sample_type > 1)
    if(nrow(tmp) > 0){
      msg <- paste("Samples with inconsistent habitat:", paste(tmp$sample, collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Coherence between samples and habitat: OK")
      print(msg)
    }
    
    # check sample_type
    sample_type_unique <- c("negative", "mock", "real")
    tmp <- unique(df$sample_type)
    incorrect_sample_type <- tmp[!tmp %in% sample_type_unique]
    if(length(incorrect_sample_type) > 0) {
      msg <- paste("The following sample types are not accepted:", 
                   paste(incorrect_sample_type, collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("sample_type: OK")
      print(msg)
    }
    
    # unique sample-replicate
    tmp <- df %>%
      select("sample","replicate") %>%
      group_by(sample, replicate) %>%
      summarize("n"=n(), .groups="drop_last") %>%
      filter(n>1)
    if(nrow(tmp) > 0){
      msg <- paste("Sample-replicate combinations must be unique:", 
                   paste(tmp$sample, collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Coherence between samples and replicates: OK")
      print(msg)
    }
  }
  
  # check if fastq file pairs are coherent (e.g. 1 to 1 relation)
  # check file extension. accept only .fastq or .fastq.gz
  if(file_type == "fastqinfo"){
    
    unique_files <- unique(df$fastq_fw)
    unique_files <- append(unique_files, unique(df$fastq_rv))
    bool <- TRUE
    for(file in unique_files){
      if( !( endsWith(file, ".fastq.gz") || endsWith(file, ".fastq") ) ){
        bool <- FALSE
        msg <- paste("Only fastq or fastq.gz formats are accepted in the 
                     fastq_fw and fastq_rv columns")
        tryCatch(stop(msg), error = function(e) message(msg))
      }
    }
    if (bool & !quiet){
      msg <- paste("File extension: OK")
      print(msg)
    }

    
    # check if fastq filepairs are coherent (e.g. 1 to 1 relation)
    tmp_rv <- df %>%
      select("fastq_fw","fastq_rv") %>%
      group_by(fastq_fw) %>%
      summarize("rv_count"=n_distinct(fastq_rv)) %>%
      filter(rv_count > 1)
    
    tmp_fw <- df %>%
      select("fastq_fw","fastq_rv") %>%
      group_by(fastq_rv) %>%
      summarize("fw_count"=n_distinct(fastq_fw)) %>%
      filter(fw_count > 1)
    
    if(nrow(tmp_rv)>0 || nrow(tmp_fw)>0) {
      msg <- paste("The following fastq files have more than one pairs:", 
                   paste(tmp_fw$fastq_rv, tmp_rv$fastq_fw, collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Coherence between fw and rv fastq filename: OK")
      print(msg)
    }
  }
  
  # check file extension. accept only .fasta .fas  .fasta.gz .fas.gz
  if(file_type == "fastainfo" || file_type == "sortedinfo"){
    bool <- TRUE
    unique_files <- unique(df$fasta)
    for(file in unique_files){
      if( !( endsWith(file, ".fasta.gz") || 
             endsWith(file, ".fasta") || 
             endsWith(file, ".fas")  || 
             endsWith(file, ".fas.gz")  
             ) 
          ){
        bool <- FALSE
        msg <- paste("Only fas, fasta, fas.gz or fasta.gz file extentions are accepted in ", 
                     file, sep="")
        tryCatch(stop(msg), error = function(e) message(msg))
      }
    }
    if(bool & !quiet){
      msg <- paste("File extension: OK")
      print(msg)
    }
  }

  
  # check if files exist
  if(file_type == "fastqinfo" || file_type == "fastainfo" || file_type == "sortedinfo"){
    
    if(file_type == "fastqinfo"){
      file_list_fw <- unique(df$fastq_fw)
      file_list_rv <- unique(df$fastq_rv)
      file_list <- c(file_list_fw, file_list_rv)
    }
    if(file_type == "fastainfo" || file_type == "sortedinfo"){
      file_list <- unique(df$fasta)
    }
    check_file_exists(dir=dir, file_list=file_list)
  }
  
  # check if tag combinations are unique within a file(pair)
  if(file_type == "fastqinfo" || file_type == "fastainfo"){
    if(file_type == "fastqinfo"){
      tmp <- df %>%
        select("tag_fw", "tag_rv", "file"=fastq_fw)
    }else{
      tmp <- df %>%
        select("tag_fw", "tag_rv", "file"=fasta)
    }
    tmp <- tmp %>%
      group_by(file, tag_fw, tag_rv) %>%
      summarize(count = n(), .groups="drop_last") %>%
      filter(count>1)
    tmp$res <- paste(tmp$tag_fw, tmp$tag_rv, tmp$file, sep=" ")
    
    if(nrow(tmp)>0) {
      msg <- paste("The following  within file tag combinations are not unique:", 
                   paste(tmp$res, collapse = "\n"))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Unique tag combinations : OK")
      print(msg)
    }
  }
  
  # check action
  if(file_type == "mock_composition" || file_type == "known_occurrences"){
    action_type <- c("keep", "delete", "tolerate")
    tmp <- unique(df$action)
    incorrect_action_type <- tmp[!tmp %in% action_type]
    if(length(incorrect_action_type) > 0) {
      msg <- paste("The following actions types are not accepted:", 
                   paste(incorrect_action_type, collapse = ", "))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("Action types : OK")
      print(msg)
    }
  }
  
  # check if 1 to 1 relation between asv_id ad asv
  if(file_type == "read_count" || 
     file_type == "read_count_sample" || 
     file_type == "asv_list" 
     ){
    
    tmp_asv_id <- df %>%
      select("asv_id","asv") %>%
      group_by(asv_id) %>%
      summarize("asv_count"=n_distinct(asv)) %>%
      filter(asv_count > 1)
    
    tmp_asv <- df %>%
      select("asv_id","asv") %>%
      group_by(asv) %>%
      summarize("asv_id_count"=n_distinct(asv_id)) %>%
      filter(asv_id_count > 1)
    
    if(nrow(tmp_asv_id)>0 || nrow(tmp_asv)>0) {
      msg <- paste("The following ASVs or asv_ids are not unique:", 
                   paste(tmp_asv_id$asv_id, tmp_asv$asv, collapse = "\n"))
      tryCatch(stop(msg), error = function(e) message(msg))
    }else if(!quiet){
      msg <- paste("1 to 1 relation between asv_id and asv : OK")
      print(msg)
    }
  }
  
}

#' Check existence of files
#' 
#' Check if all files in the input vector exist.
#'  
#' @param dir Character string: name of the directory containing the input files
#' @param file_list Vector of file names.
#' @returns Error message if some files do not exist and stop program.
#' @examples
#' \dontrun{
#' file_list <- c("14ben01-1.fasta", "14ben01-2.fasta")
#' check_file_exists(file_list=file_list, dir="vtamR_test/out_mfzr/sorted")
#' }
#' @export
#' 
check_file_exists <- function(file_list, dir=""){
  
  check_dir(dir)
  missing <- c()
  for(i in file_list){
    file_p <- file.path(dir, i)
    if(!file.exists(file_p)){
      missing <- append(missing, file_p)
    }
  }
  if(length(missing)>0 ) {
    msg <- paste("The following files do not exist :", paste(missing, collapse = ", "))
    tryCatch(stop(msg), error = function(e) message(msg))
  }
}

#' Check heading of a files
#' 
#' Compare the heading of files to an input vector. 
#' 
#'  
#' @param list1,list2 Vectors: `list1` is a heading of file, 
#' `list2` is the expected heading for the file.
#' @param file Character string: name of the file.
#' @returns Error message with elements of `list1` that are not on `list2`.
#' @examples
#' \dontrun{
#' list1 <- c("tag_fw","primer_fw","tag_rv","primer_rv","sample",
#' "sample_type","habitat","replicate","fasta")
#' list2 <- c("tag_fw","primer_fw","tag_rv","primer_rv","sample",
#' "sample_type","habitat","replicate","fastq_fw","fastq_rv")
#' check_heading(list1, list2)
#' }
#' @export
#' 
check_heading <- function(list1, list2, file="") {
  bool <- TRUE
  col <- character(0)
  for (i in list1) {
    if (!(i %in% list2)) {
      col <- append(col, i)
      bool <- FALSE
    }
  }
  col <- paste(col, collapse = ", ")
  if (!bool) {
    msg <- paste("The following column(s) are missing from", file, ":", col, sep = " ")
    tryCatch(stop(msg), error = function(e) message(msg))
  }
}

#' ClusterSize
#' 
#' Cluster ASV in input data frame with the cluster_size command of vsearch. 
#' Replace ASV by their centroids and regroup lines by centroid and sample.
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate (optional), read_count, asv.
#' @param id Real; Minimum identity between asv and centroid.
#' @param vsearch_path Character string: path to vsearch executables.
#' @param by_sample logical: run clustering separately for each sample.
#' @param num_threads Positive integer: Number of CPUs.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns read_count data frame, with ASV of the same cluster, sample, (replicate)
#' grouped to the same line
#' @examples
#' \dontrun{
#' clustered_df <- ClusterSize(read_count_df, id=0.97)
#' }
#' @export
#' 
ClusterSize <- function(read_count, 
                        id=0.97, 
                        vsearch_path="vsearch", 
                        by_sample = FALSE,
                        num_threads = 1,
                        outfile="", 
                        sep=",",
                        quiet=TRUE
                        ) {
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  if(by_sample){
    # make an empty output df, with the same columns and variable types as read_count_df
    out_df <- read_count_df %>%
      filter(asv=="")
    
    # get list of samples 
    sample_list <- unique(read_count_df$sample)
    
    # run cluster_size for each sample
    for(s in sample_list){
      if(!quiet){
        print(s)
      }
      
      # select occurrences for sample
      df_sample <- read_count_df %>%
        filter(sample==s)
      # run clustering
      df_sample <- run_clustersize(df_sample, 
                                   id = id,
                                   vsearch_path=vsearch_path, 
                                   num_threads=num_threads, 
                                   quiet=quiet)
      # add output of the sample to the total data frame
      out_df <- rbind(out_df, df_sample)
    }
  }else{ # run clustering for all samples together
    out_df <- run_clustersize(read_count_df, 
                              id = id,
                              vsearch_path=vsearch_path, 
                              num_threads=num_threads, 
                              quiet=quiet)
  }
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(out_df, file = outfile,  row.names = F, sep=sep)
  }
  return(out_df)
}

#' run_clustersize
#' 
#' Run vsearch clustersize on all ASV of the input data frame and pool
#' ASV of the same cluster under the asv_id of the centroid.
#'  
#' @param read_count_df Data frame with the following variables: 
#' asv_id, sample, replicate (optional), read_count, asv.
#' @param id Real; Minimum identity between asv and centroid.
#' @param vsearch_path Character string: path to vsearch executables.
#' @param num_threads Positive integer: Number of CPUs.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns read_count data frame, with ASV of the same cluster, sample, (replicate)
#' grouped to the same line
#' @examples
#' \dontrun{
#' clustered_df <- run_clustersize(read_count_df, id=0.97)
#' }
#' @export
run_clustersize <- function(read_count_df, 
                            id = 0.97,
                            vsearch_path="vsearch", 
                            num_threads=1, 
                            quiet=TRUE){
  
  ### make tmp dir
  outdir_tmp <- paste('tmp_cluster_size_', 
                      trunc(as.numeric(Sys.time())), 
                      sample(1:100, 1), 
                      sep='')
  outdir_tmp <- file.path(tempdir(), outdir_tmp)
  check_dir(outdir_tmp)
  
  ### define temporary file names
  fasta <- file.path(outdir_tmp, 'asv.fasta')
  blast6_file <- file.path(outdir_tmp, 'blast6.tsv')
  
  ### get unique list of asv with read_count and asv_id
  asv_rc <- read_count_df %>%
    group_by(asv) %>%
    summarize(
      read_count = sum(read_count), 
      asv_id = first(asv_id) 
    ) %>%
    arrange(desc(read_count))
  
  # make fasta file with abundances
  write_fasta_rc(asv_rc, fasta)
  
  # run vsearch cluster_size
  cmd <- paste(vsearch_path, 
               ' --cluster_size ', fasta,
               ' --blast6out ', blast6_file,
               ' --id ', id, 
               sep="")
  if(!quiet){
    print(cmd)
  }
  system(cmd)
  
  file_info <- file.info(blast6_file)
  
  if(file_info$size == 0){ # No output of clustering
      # Delete the temp directory
      unlink(outdir_tmp, recursive = TRUE)
      return(read_count_df)
  } else{
    # read clustering results to df_centroids
    # merged_id (include to a cluster), centroid_id (most abundant asv_id of the cluster)
    df_centroids <- read_blast6out(blast6_file)
    # add centroids to all asv_ids 
    read_count_df_tmp <- left_join(read_count_df, 
                                   df_centroids, 
                                   by=c("asv_id" = "merged_id")
    ) 
    
    read_count_df_tmp <- read_count_df_tmp %>%
      # if sequence is not in a cluster or if it is a centroid, asv not in blast6_file
      # => centroid_id is NA => change to asv_id
      mutate(centroid_id = if_else(is.na(centroid_id), asv_id, centroid_id)) %>%
      select(-asv, -asv_id) %>% # remore outdated columns
      rename(asv_id = centroid_id) # replace original asv_id by centroid_id
    
    
    # list of unique asv (before clustering)
    asv_rc <- asv_rc %>%
      select(-read_count) # colnames(asv_rc):asv, asv_id
    
    # add sequence of centroid to read_count_df
    read_count_df_tmp <- read_count_df_tmp %>%
      left_join(asv_rc, by=c("asv_id")) %>%
      select(asv_id, everything()) # put asv_id as a first column
    
    # regroup by asv (replaced by centroid) and sample and (replicate)
    if("replicate" %in% colnames(read_count_df_tmp)){
      read_count_df_tmp <- read_count_df_tmp %>%
        group_by(asv_id, sample, replicate) %>% 
        summarize(read_count = sum(read_count), asv = first(asv), .groups="drop_last")
    } else{ # not replicates
      read_count_df_tmp <- read_count_df_tmp %>%
        group_by(asv_id, sample) %>% 
        summarize(read_count = sum(read_count),asv = first(asv), .groups="drop_last")
    }
    
    # Delete the temp directory
    unlink(outdir_tmp, recursive = TRUE)
    return(read_count_df_tmp)
    }
    
}
#' 
#' 
#' 
#' 
#' read_blast6out
#' 
#' read output of vsearch --cluster_size to df
#'  
#' @param filename Tab separeted file with merged and centroid asv_ids in the first 2 columns
#' @returns data frame with merged_id and centroid_id
#' @export
#' 

read_blast6out <- function(filename) {
  
  df <- read.table(filename, header=FALSE, sep="\t")
  df <- df[,1:2]
  colnames(df) <- c("merged_id","centroid_id")
  df$merged_id <- gsub(';size=[0-9]+', '', df$merged_id)
  df$centroid_id <- gsub(';size=[0-9]+', '', df$centroid_id)
  
  df$merged_id <- as.integer(df$merged_id)
  df$centroid_id <- as.integer(df$centroid_id)
  
  return(df)
}

#' write_fasta_rc
#' 
#' Write a fasta file with definition lines in >label;size=### format
#'  
#' @param df Data frame with asv_id, read_count and asv columns
#' @param outfile haracter string: output fasta file name
#' @returns outfile fasta file with definition lines in >label;size=### format
#' @export
#' 
write_fasta_rc <- function(df, outfile) {
  # Open the file for writing
  file <- file(outfile, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(df$asv)) {
    seq <- df$asv[i]
    count <- df$read_count[i]
    seqid <- df$asv_id[i]
    header <- paste0(">", seqid, ';size=', count)
    writeLines(c(header, seq, ""), file)
  }
  # Close the file
  close(file)
}

#' write_fasta_df
#' 
#' Write a fasta file with definition lines in >label;size=### format
#'  
#' @param df Data frame with asv_id, asv and read_count (optional) columns
#' @param outfile character string: output fasta file name
#' @param read_count Boolean; If TRUE  definition line is >label;size=### format
#' @returns outfile fasta file with definition lines in >label;size=### format or >label format
#' @export
#' 
write_fasta_df <- function(df, outfile, read_count=FALSE) {
  # Open the file for writing
  
  if(read_count){
    df <- df %>%
      group_by(asv_id) %>%
      summarize(
        read_count = sum(read_count),
        asv = first(asv)
        ) 

  }else{
    df <- df %>%
      select(asv_id, asv) %>%
      distinct()
  }
  
  file <- file(outfile, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(df$asv)) {
    seq <- df$asv[i]
    seqid <- df$asv_id[i]
    if(read_count){
      count <- df$read_count[i]
      header <- paste0(">", seqid, ';size=', count)
    }else{
      header <- paste0(">", seqid)
    }
    writeLines(c(header, seq, ""), file)
  }
  # Close the file
  close(file)
}