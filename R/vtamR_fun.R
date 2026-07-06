#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first if_else
#' @importFrom dplyr slice_head
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym :=
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL


#' Run a command using system2
#' 
#' Execute a system command using `system2`, with optional control over 
#' the verbosity of the output.
#'  
#' @param path Character string specifying the path to the executable.
#' @param args Character vector of arguments passed to the command.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @examples 
#' \dontrun{
#' run_system2(path = "ls", args = c("-all"), quiet = FALSE)
#' }
#' 
#' @export

run_system2 <- function(path, args, quiet = FALSE) {
  
  # system2 cannot use ~ as a home
  path <- path.expand(path)
  
  if (!quiet) {
    # Show the full command that will be run
    cat("Running command:\n")
    cat(path, paste(shQuote(args), collapse = " "), "\n")
    
    system2(
      command = path,
      args = args,
      stdout = "",
      stderr = ""
    )
    
  } else {
    output <- suppressWarnings(system2(
      command = path,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    ))
    
    # Extract only error/warning/fail lines
    errors_only <- grep("error|warning|fail", output, ignore.case = TRUE, value = TRUE)
    
    # Remove lines containing "Mean expected error"
    errors_only <- grep("Mean expected error|Mean observed errors|Pairs that failed|expected error too high|failed merging", 
                        errors_only, ignore.case = TRUE, value = TRUE, 
                        invert = TRUE)
    
    if (length(errors_only) > 0) {
      cat(errors_only, sep = "\n")
    }
  }
}


#' Check directory
#' 
#' Check whether a directory exists and create it if necessary. If a file path 
#' is provided, the directory portion is extracted and created if needed.
#'  
#' @param path Character string specifying a directory or a file path.
#' @param is_file Logical. If `TRUE`, `path` is treated as a file path; 
#'   otherwise, it is treated as a directory path.
#' 
#' @return `NULL`
#' 
#' @examples 
#' \dontrun{
#' check_dir(path = "data")
#' }
#' 
#' @export

check_dir <- function(path, is_file=FALSE){
  
  if(is_file){
    dir_to_create <- dirname(path)
  }else{
    path <- sub("[/\\\\]+$", "", path) # remove / or \ at the end of dir name
    dir_to_create <- path
  }
  
  if(!dir.exists(dir_to_create)){
      dir.create(dir_to_create, recursive =TRUE)
  }
  return(invisible(path))
}


#' Compress or uncompress a file
#'
#' Compress or uncompress a file using `pigz` (if available, for faster 
#' performance) or fall back to `R.utils`.
#'
#' @param file Character string specifying the input file.
#' @param outfile Character string specifying the output file. If NULL 
#'   automatically derived from file.
#' @param remove Logical. If `TRUE`, remove the input file after a successful 
#'   operation.
#' @param method Character. Compression method: `pigz`, `gzip`, or `R`.  
#'   `pigz` requires `pigz` to be installed and available in the system 
#'   PATH (or specified via `pigz_path`).  
#'   `gzip` is available on Linux systems.  
#'   `R` uses `R.utils`, which is cross-platform but slower.  
#'   Relative speed: `R.utils` < `gzip` < `pigz`.
#' @param pigz_path Character string specifying the path to the `pigz`
#'   executable. Only required if `method = "pigz"` and it is not available 
#'   in the system PATH.
#' @param num_threads Positive integer specifying the number of CPU threads to 
#'   use. If `0`, all available CPUs are used.
#' @param quiet Logical. If `TRUE`, suppress informational messages.
#' @param compress Logical. If `TRUE`, compress the file; if 
#'   `FALSE`, decompress it.
#'
#' @return Invisibly returns a character string: the path to the output file.
#' 
#' @examples
#' \dontrun{
#' # Decompress a .gz file
#' smart_gzip("data.fastq.gz", compress = FALSE)
#'
#' # Compress a file
#' smart_gzip("data.fastq", compress = TRUE)
#' }
#' 
#' @export
smart_gzip <- function(file,
                       outfile = NULL,
                       remove = FALSE,
                       pigz_path = "pigz",
                       method = "R",
                       num_threads = 0,
                       quiet = TRUE,
                       compress = FALSE) {
  # Check input
  if (!file.exists(file)) stop("File not found: ", file)
  
  # Detect operation type and define flags
  if (compress) {
    mode_flag <- "-c"        # compress
    file_check <- !grepl("\\.gz$", file, ignore.case = TRUE)
  } else {
    mode_flag <- "-d"        # decompress
    file_check <- grepl("\\.gz$", file, ignore.case = TRUE)
  }
  
  if (!file_check)
    stop("File extension does not match expected mode: ",
         if (compress) "expected uncompressed input" else "expected .gz input")
  
  # Determine output file
  if (is.null(outfile)) {
    outfile <- if (compress) paste0(file, ".gz") else sub("\\.gz$", "", file)
  }else{ # check if output filename is coherent
    if (compress) {
      outfile_check <- grepl("\\.gz$", outfile, ignore.case = TRUE)
    } else {
      outfile_check <- !grepl("\\.gz$", outfile, ignore.case = TRUE)
    }
    if (!outfile_check)
      stop("File extension does not match expected mode: ",
           if (compress) "expected compressed (.gz ) output" else "expected uncompessed output")
  }
  outfile <- path.expand(outfile)
  
  # Detect number of threads
  if (num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  
  # Use pigz if provided
  if (method=="pigz") {
    if (!quiet) {
      message("Using pigz for ", if (compress) "compression" else "decompression")
    }
    # Build argument list
    args <- c(mode_flag, "-p", num_threads, "-c", file)
    
    # Run pigz with system2 (safe cross-platform)
    status <- system2(pigz_path, args = args, stdout = outfile)
    
    if (status != 0) stop("pigz failed with exit code: ", status)
    if (remove) unlink(file)
  } else if (method=="gzip"){
    if (!quiet) {
      message("Using gzip for ", if (compress) "compression" else "decompression")    
      
      if (compress) {
        status <- system2("gzip", args = c("-c", file), stdout = outfile)
        # If gzip ran successfully, remove the input file
        if (status == 0 && remove) unlink(file)
        
      } else {
        status <- system2("gunzip", args = c("-c", file), stdout = outfile)
        if (status == 0 && remove) unlink(file)
      }
    }
    
  } else {
    if (!quiet) {
      message("Using R.utils for ", if (compress) "compression" else "decompression")
    }
    
    if (compress) {
      R.utils::gzip(file, destname = outfile, remove = remove, overwrite = TRUE)
    } else {
      R.utils::gunzip(file, destname = outfile, remove = remove, overwrite = TRUE)
    }
  }
  
  return(invisible(outfile))
}


#' get_stat
#' 
#' Compute summary statistics (reads, variants, samples, and replicates) and 
#' append them to a statistics data frame.
#'  
#' @param read_count Data frame or path to a CSV file containing the following 
#'   variables: `asv_id`, `sample``, `replicate, 
#'   `read_count`, `asv`.
#' @param stat_df A data frame containing the following variables:
#' 
#'     `parameters` Parameters used for the analysis
#'     `asv_count` Number of ASVs detected
#'     `read_count` Total number of reads
#'     `sample_count` Number of samples analyzed
#'     `sample_replicate_count` Number of replicates per sample
#'   If provided, the function appends a new row to `stat_df`. If not provided,
#'   a new data frame with a single row is initialized and returned.
#' @param stage Character string specifying the name of the filtering step. It 
#'   is used as the row name in `stat_df`.
#' @param params Character string containing concatenated parameter values used 
#'   for the filtering step.
#' @param outfile Character string specifying the name of a CSV file to write 
#'   the updated data frame. If NULL, no file is written.
#' 
#' @return Data frame. The updated `stat_df` with an additional row.
#' 
#' @examples
#' \dontrun{
#' get_stat(
#'   read_count_df,
#'   stat_df,
#'   stage = "filter_occurrence_variant",
#'   params = "0.002;by_replicate=TRUE"
#' )
#' 
#' get_stat(
#'   read_count_df,
#'   stat_df = NULL,
#'   stage = "filter_indel",
#'   params = "0.002;by_replicate=TRUE",
#'   outfile = "out/ReadCount_stat.csv"
#' )
#' }
#' 
#' @export
#' 
get_stat <- function(read_count, stat_df=NULL, stage="", params=NA, outfile=NULL){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  if(is.null(stat_df)){
    stat_df <- data.frame(
      parameters = character(),
      asv_count = integer(),
      read_count = integer(),
      sample_count = integer(),
      sample_replicate_count = integer()
    )
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
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(stat_df, file = outfile,  row.names = F, sep=sep)
  }
  return(stat_df)
}

#' Merge forward and reverse reads
#' 
#' Merge paired-end FASTQ reads (forward and reverse) and convert the resulting 
#' sequences to FASTA format. Output FASTA files can optionally be compressed.
#' The generated `fastainfo.csv` file is similar to the input 
#' `fastqinfo` file, but FASTQ file columns are replaced by a 
#' `fasta` column containing the names of the output files.
#'   
#' @param fastqinfo Data frame or path to a CSV file with the following columns: 
#'   `tag_fw`, `primer_fw`, `tag_rv`, `primer_rv`, 
#'   `sample`, `sample_type` (mock/negative/real), 
#'   `habitat` (optional), `replicate`, `fastq_fw`, 
#'   `fastq_rv`.
#' @param fastq_dir Character string specifying the directory containing input 
#'   FASTQ files (listed in `fastqinfo$fastq_fw` and 
#'   `fastqinfo$fastq_rv`).
#' @param outdir Character string specifying the output directory.
#' @param vsearch_path Character string specifying the path to the 
#'   `vsearch` executable.
#' @param compress_method Character or logical. Compression method: 
#'   `"pigz"`, `"gzip"`, or `"R"`.  
#'   `"pigz"` requires `pigz` to be installed and available in the 
#'   system PATH (or specified via `pigz_path`).  
#'   `"gzip"` is available on Linux systems.  
#'   `"R"` uses `R.utils`, which is cross-platform but slower.  
#'   Relative speed: `R.utils` < `gzip` < `pigz`.
#' @param pigz_path Character string specifying the path to the `pigz` 
#'   executable. Only required if `compress_method = "pigz"` and it is not 
#'   available in the system PATH.
#' @param num_threads Positive integer specifying the number of CPU threads to 
#'   use. If `0`, all available CPUs are used.
#' @param fastq_ascii Integer (33 or 64) specifying the ASCII offset used for 
#'   FASTQ quality scores.
#' @param fastq_maxdiffs Positive integer specifying the maximum number of 
#'   mismatches allowed in the overlapping region.
#' @param fastq_maxee Positive integer specifying the maximum number of expected 
#'   errors allowed per sequence.
#' @param fastq_minlen Positive integer specifying the minimum allowed sequence 
#'   length.
#' @param fastq_maxlen Positive integer specifying the maximum allowed sequence 
#'   length.
#' @param fastq_minmergelen Positive integer specifying the minimum length of 
#'   merged sequences.
#' @param fastq_maxmergelen Positive integer specifying the maximum length of 
#'   merged sequences.
#' @param fastq_maxns Positive integer specifying the maximum number of 
#'   ambiguous bases (`N`) allowed per sequence.
#' @param fastq_truncqual Positive integer specifying the quality score threshold 
#'   for truncating sequences.
#' @param fastq_minovlen Positive integer specifying the minimum overlap length 
#'   required for merging reads.
#' @param fastq_allowmergestagger Logical. If `TRUE`, allow merging of 
#'   staggered read pairs (where the reverse read extends beyond the forward 
#'   read).
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param compress Logical. If `TRUE`, compress output files using gzip.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Data frame corresponding to the generated `fastainfo.csv` file.
#' 
#' @examples
#' \dontrun{
#' merge_fastq_pairs(
#'   fastqinfo_df,
#'   fastq_dir = "data/fastqdir",
#'   vsearch_path = "C:/Users/Public/vsearch",
#'   outdir = "data/fastadir",
#'   compress = TRUE,
#'   quiet = FALSE
#' )
#' 
#' merge_fastq_pairs(
#'   fastqinfo_df,
#'   fastq_dir = "data/fastqdir",
#'   outdir = "data/fastadir",
#'   fastq_maxdiffs = 5,
#'   fastq_maxee = 2,
#'   fastq_minlen = 60,
#'   fastq_maxlen = 100,
#'   fastq_minmergelen = 80,
#'   fastq_maxmergelen = 100,
#'   fastq_maxns = 1,
#'   fastq_truncqual = 20,
#'   fastq_minovlen = 20,
#'   fastq_allowmergestagger = TRUE
#' )
#' }
#' 
#' @export
#'
merge_fastq_pairs <- function(fastqinfo, 
                  fastq_dir, 
                  outdir, 
                  vsearch_path="vsearch",
                  compress_method="R",
                  pigz_path="pigz",
                  num_threads=0,
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
                  fastq_allowmergestagger=FALSE, 
                  sep=",", 
                  compress=FALSE, 
                  quiet=T){
  
  fastq_dir = check_dir(fastq_dir)
  outdir = check_dir(outdir)
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # can accept df or file as an input
  if(is.character(fastqinfo)){
    # read known occurrences
    fastqinfo_df <- read.csv(fastqinfo, header=T, sep=sep)
  }else{
    fastqinfo_df <- fastqinfo
  }
  check_file_info(file=fastqinfo_df, dir=fastq_dir, file_type="fastqinfo", sep=sep, quiet=TRUE)
  
  #  get unique list of fw_rv filenames and outfile names
  tmp <- fastqinfo_df %>%
    # only one sample-replicate combination per fw-rv file combination. Use sample-replicate.fasta as outfile
    mutate(fasta = paste(sample, "-", replicate, ".fasta", sep=""),
           read_count = NA) %>% 
    select(fasta, fastq_fw, fastq_rv)
  
  nb_sample_repl = length(unique(tmp$fasta))
  nb_fw_rv = nrow(tmp %>%
                    select(fastq_fw, fastq_rv)%>%
                    distinct()
  )
  if(nb_sample_repl != nb_fw_rv){# Same fw_rv file combination, have more than one sample-replicate
    # it should be demultiplexed latter. Use the fw filename for output
    tmp <- tmp %>%
      mutate(fasta = sub("\\..*$", ".fasta", fastq_fw)) %>%
      distinct()
  } 
  
  for(i in 1:nrow(tmp)){# for each file pairs
    
    # use the name of the fw fastq file and replace extension by fasta (uncompressed)
#    outfile <- sub("\\..*", ".fasta", tmp[i,1])
#    tmp$fasta[i] <- outfile
    outfile <- file.path(outdir, tmp[i,"fasta"])
    # add path to input filenames
    fw_fastq <- file.path(fastq_dir, tmp[i,"fastq_fw"])
    rv_fastq <- file.path(fastq_dir, tmp[i,"fastq_rv"])
    
    #Decompress input files, since they are cannot be treated directly by vsearch on the OS
    if(!is_linux() && endsWith(fw_fastq, ".gz")){
      fw_fastq <- smart_gzip(fw_fastq, 
                            remove = FALSE,
                            method=compress_method,
                            pigz_path=pigz_path,
                            num_threads = num_threads,
                            quiet = quiet,
                            compress = FALSE)
      rv_fastq <- smart_gzip(rv_fastq, 
                             remove = FALSE,
                             method=compress_method,
                             pigz_path=pigz_path,
                             num_threads = num_threads,
                             quiet = quiet,
                             compress = FALSE)
      
        }
    
    if(fw_fastq == outfile){ # stop the run if input and output files have the same name
      stop("ERROR: Input and output directories for fastq and fasta files are indentical. 
           Please, give a different output directory")
    }
    
    ##### run vsearch
    # Build argument vector
    args <- c(
      "--fastq_mergepairs", fw_fastq,
      "--reverse", rv_fastq ,
      "--fastaout", outfile,
      "--quiet",
      "--fastq_ascii", fastq_ascii,
      "--fastq_maxdiffs", fastq_maxdiffs, 
      "--fastq_maxee", fastq_maxee, 
      "--fastq_minlen", fastq_minlen, 
      "--fastq_maxlen",fastq_maxlen, 
      "--fastq_minmergelen",fastq_minmergelen,
      "--fastq_maxmergelen",fastq_maxmergelen,
      "--fastq_maxns", fastq_maxns, 
      "--fastq_truncqual", fastq_truncqual, 
      "--fastq_minovlen", fastq_minovlen
    )
    if(num_threads > 0){
      args <- append(args, c("--threads", num_threads))
    }
    if(fastq_allowmergestagger){ # if reads are longer than the amplicon
      args <- append(args, c("--fastq_allowmergestagger"))
    }
    run_system2(vsearch_path, args, quiet=quiet)
    
    seq_n <- count_seq(outfile)
    tmp$read_count[i] <- seq_n

    # vsearch produces uncompressed files even if input is compressed => compress output file
    if(compress){
        out <- smart_gzip(outfile, 
                          remove = TRUE,
                          method=compress_method,
                          pigz_path=pigz_path,
                          num_threads = num_threads,
                          quiet = quiet,
                          compress = TRUE)
        # correct output filename in fastainfo if necessary
        if(!endsWith(tmp$fasta[i], ".gz")){ 
          tmp$fasta[i] <- paste(tmp$fasta[i], ".gz", sep="")
        }
      }
      
    original_fw_fastq <- file.path(fastq_dir, tmp[i,"fastq_fw"])
    if( original_fw_fastq != fw_fastq){# the input fastq has been unzipped for vsearch => rm unzipped file to free space
      file.remove(fw_fastq)
      file.remove(rv_fastq)
    }
  } # end loop over files
  # make fastainfo file
  # rename read_count if present in input fastqinfo
  if("read_count" %in% colnames(fastqinfo_df)){
    fastqinfo_df <- fastqinfo_df %>%
      rename(read_count_input = read_count)
  }
  fastainfo_df <- left_join(fastqinfo_df, tmp, by=c("fastq_fw", "fastq_rv")) %>%
    select(-fastq_fw, -fastq_rv)
  write.table(fastainfo_df, file = file.path(outdir, "fastainfo.csv"),  row.names = F, sep=sep)
  
  return(fastainfo_df)
  
}

#' Test if OS is Linux-like
#' 
#' Determine whether the operating system is Unix/Linux-like.
#'  
#' @return Logical. `TRUE` if the `"sysname"` returned by 
#'   `Sys.info()` starts with one of: `linux`, `sunos`, 
#'   `darwin`, `gnu`, or `unix`; otherwise `FALSE`.
#' 
#' @examples
#' \dontrun{
#' os_linux <- is_linux()
#' }
#' 
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

#' Count sequences in FASTA
#' 
#' Count the number of sequences in a FASTA file.
#'  
#' The input file can be uncompressed or gzip-compressed; other compression 
#' formats are not supported. On Unix/Linux-like systems, the function uses 
#' the `grep` and `wc` shell commands for fast processing. On other 
#' operating systems, it falls back to `count.fields`, which is slower.
#'  
#' @param file Character string specifying the input FASTA file (including path).
#' 
#' @return Integer. The number of sequences in the input file.
#' 
#' @examples
#' \dontrun{
#' n <- count_seq(file = "data/test.fasta")
#' }
#' 
#' @export
#' 
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
#' Trim primer sequences from an input FASTA file.
#'  
#' The input FASTA file can be uncompressed or gzip-compressed; other 
#' compression formats are not supported. The compression of the output file 
#' is determined by the output filename.
#'   
#' @param fasta Character string specifying the input FASTA file (including path). 
#' @param outfile Character string specifying the output FASTA file (including path).
#' @param primer_fw Character string specifying the forward primer (IUPAC ambiguity codes are accepted).
#' @param primer_rv Character string specifying the reverse primer (IUPAC ambiguity codes are accepted).
#' @param vsearch_path Character string specifying the path to the 
#'   `vsearch` executable. 
#' @param cutadapt_path Character string specifying the path to the 
#'   `cutadapt` executable. 
#' @param num_threads Positive integer specifying the number of CPU threads to 
#'   use. If `0`, all available CPUs are used.
#' @param compress_method Character or logical. Compression method: 
#'   `"pigz"`, `"gzip"`, or `"R"`.  
#'   `"pigz"` requires `pigz` to be installed and available in the 
#'   system PATH (or specified via `pigz_path`).  
#'   `"gzip"` is available on Linux systems.  
#'   `"R"` uses `R.utils`, which is cross-platform but slower.  
#'   Relative speed: `R.utils` < `gzip` < `pigz`.  
#'   Only used if `check_reverse = TRUE`.
#' @param pigz_path Character string specifying the path to the `pigz` 
#'   executable. Only required if `compress_method = "pigz"` and it is not 
#'   available in the system PATH.
#' @param check_reverse Logical. If `TRUE`, also check reverse-complemented 
#'   sequences from the input FASTA file.
#' @param primer_to_end Logical. If `TRUE`, primers are assumed to be 
#'   directly adjacent to tags (i.e., no heterogeneity spacer).
#' @param cutadapt_error_rate Numeric value between 0 and 1 specifying the 
#'   maximum allowed error rate between primers and reads (exact match is 
#'   required for tags).
#' @param cutadapt_minimum_length Positive integer specifying the minimum length 
#'   of trimmed sequences.
#' @param cutadapt_maximum_length Positive integer specifying the maximum length 
#'   of trimmed sequences.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return `NULL`. Produces an output FASTA file.
#' 
#' @examples
#' \dontrun{
#' trim_primers_file(
#'   fasta = "data/test.fasta",
#'   outfile = "out/test_trimmed.fasta",
#'   primer_fw = "TCCACTAATCACAARGATATTGGTAC",
#'   primer_rv = "WACTAATCAATTWCCAAATCCTCC",
#'   check_reverse = TRUE,
#'   primer_to_end = TRUE,
#'   cutadapt_minimum_length = 300,
#'   cutadapt_maximum_length = 400
#' )
#' }
#' 
#' @export
#' 
trim_primers_file <- function(fasta, 
                               outfile, 
                               primer_fw, 
                               primer_rv, 
                               cutadapt_path="cutadapt", 
                               vsearch_path="vsearch", 
                               compress_method="R",
                               pigz_path="pigz",
                               num_threads = 0,
                               check_reverse=F, 
                               primer_to_end=T, 
                               cutadapt_error_rate=0.1,
                               cutadapt_minimum_length=50,
                               cutadapt_maximum_length=500, 
                               sep = ",",
                               quiet=T
                               ){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
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
    g <- paste("^", primer_fw, "...", primer_rv_rc, "$", sep="")
  } else{
    g <- paste(primer_fw, ";min_overlap=", nchar(primer_fw), "...", primer_rv_rc, ";min_overlap=", nchar(primer_rv_rc), sep="")
  }
  ##### run cutadapt
  args <- c(
    "--cores", num_threads,
    "-e", cutadapt_error_rate,
    "--no-indels",
    "--trimmed-only",
    "--minimum-length", cutadapt_minimum_length,
    "--maximum-length", cutadapt_maximum_length, 
    "-g ", shQuote(g),
    "--output", outfile,
    fasta
  )
  if(quiet){
    args <- append(args, c("--quiet"))
  }
  run_system2(cutadapt_path, args, quiet=quiet)
  
  if(check_reverse){
    # exchange fw and rv primers
    primer_rv_rc  <- reverse_complement(primer_fw)
    primer_fw <- primer_rv
    # change output filename
    out_rv <- sub("\\.", "_rv.", outfile)
    
    if(primer_to_end){
      g <- paste("^", primer_fw, "...", primer_rv_rc, "$", sep="")
    } else{
      g <- paste(primer_fw, ";min_overlap=", nchar(primer_fw), "...", primer_rv_rc, ";min_overlap=", nchar(primer_rv_rc), sep="")
    }
    ##### run cutadapt
    args <- c(
      "--cores", num_threads,
      "-e", cutadapt_error_rate,
      "--no-indels",
      "--trimmed-only",
      "--minimum-length", cutadapt_minimum_length,
      "--maximum-length", cutadapt_maximum_length, 
      "-g ", shQuote(g),
      "--output", out_rv,
      fasta
    )
    if(quiet){
        args <- append(args, c("--quiet"))
    }
    run_system2(cutadapt_path, args, quiet=quiet)
    
    # reverse complement rv file and append it to the outfile
    if(file.size(out_rv) > 0){ # there are sequences in the reverse trimmed file
      # reverse complement sequences in out_rv file
      out_rv_rc <- sub("\\.", "_rc.", out_rv)

      ##### run vsearch
      args <- c(
        "--fastx_revcomp", out_rv, 
        "--fastaout", out_rv_rc, 
        "--quiet"
      )
      if(num_threads > 0){
        args <- append(args, c("--threads", num_threads))
      }
      run_system2(vsearch_path, args, quiet=quiet)
      
      # append content of minus_rc to plus file
      file.append(outfile, out_rv_rc)
      unlink(out_rv_rc)
    }
    unlink(out_rv)
    
    if(outfile != original_output ){ # Output should be compressed
      outfile <- smart_gzip(file=outfile,
                            remove = TRUE,
                            method = compress_method,
                            pigz_path = pigz_path,
                            num_threads = num_threads,
                            quiet = quiet,
                            compress = TRUE)
    }
  }
}

#' Trim primers from FASTA files
#' 
#' Trim primer sequences from each FASTA file listed in the input 
#' `fastainfo` data frame or CSV file. Only trimmed reads are retained. 
#' Both orientations can be checked if required. The number of reads in the 
#' output files is also computed.
#'   
#' @param fastainfo Data frame or path to a CSV file with the following columns: 
#'   `tag_fw`, `primer_fw`, `tag_rv`, `primer_rv`, 
#'   `sample`, `sample_type`, `habitat`, `replicate`, 
#'   `fasta`, `read_count`.
#' @param fasta_dir Character string specifying the directory containing the 
#'   input FASTA files.
#' @param outdir Character string specifying the output directory for trimmed 
#'   FASTA files.
#' @param cutadapt_path Character string specifying the path to the 
#'   `cutadapt` executable. 
#' @param vsearch_path Character string specifying the path to the 
#'   `vsearch` executable.
#' @param check_reverse Logical. If `TRUE`, also check reverse-complemented 
#'   sequences from the input FASTA files.
#' @param primer_to_end Logical. If `TRUE`, primers are assumed to be 
#'   directly adjacent to tags (i.e., no heterogeneity spacer).
#' @param cutadapt_error_rate Numeric value between 0 and 1 specifying the 
#'   maximum allowed error rate between primers and reads (exact match is 
#'   required for tags).
#' @param cutadapt_minimum_length Positive integer specifying the minimum length 
#'   of trimmed sequences.
#' @param cutadapt_maximum_length Positive integer specifying the maximum length 
#'   of trimmed sequences.
#' @param compress Logical. If `TRUE`, compress output files using gzip.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param compress_method Character or logical. Compression method: 
#'   `"pigz"`, `"gzip"`, or `"R"`.  
#'   `"pigz"` requires `pigz` to be installed and available in the 
#'   system PATH (or specified via `pigz_path`).  
#'   `"gzip"` is available on Linux systems.  
#'   `"R"` uses `R.utils`, which is cross-platform but slower.  
#'   Relative speed: `R.utils` < `gzip` < `pigz`. 
#' @param pigz_path Character string specifying the path to the `pigz` 
#'   executable. Only required if `compress_method = "pigz"` and it is not 
#'   available in the system PATH.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Data frame. The updated `fastainfo` data frame with modified 
#'   file names and sequence counts.
#' 
#' @examples
#' \dontrun{
#' fastainfo_df <- trim_primers(
#'   fastainfo,
#'   fasta_dir,
#'   outdir,
#'   compress = TRUE,
#'   check_reverse = TRUE,
#'   primer_to_end = FALSE,
#'   cutadapt_error_rate = 0.1,
#'   cutadapt_minimum_length = 300,
#'   cutadapt_maximum_length = 350,
#'   quiet = TRUE
#' )
#' }
#' 
#' @export
#' 
trim_primers <- function(fastainfo, 
                       fasta_dir, 
                       outdir, 
                       compress=F, 
                       compress_method="R",
                       pigz_path="pigz",
                       cutadapt_path="cutadapt", 
                       vsearch_path="vsearch", 
                       check_reverse=F, 
                       primer_to_end=T, 
                       cutadapt_error_rate=0.1, 
                       cutadapt_minimum_length=50, 
                       cutadapt_maximum_length=500, 
                       sep = ",",
                       quiet=T
                       ){
  
  fasta_dir = check_dir(fasta_dir)
  outdir = check_dir(outdir)
  
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  check_file_info(file=fastainfo_df, dir=fasta_dir, file_type="fastainfo", sep=sep, quiet=TRUE)
  
  # upper case for all primers and tags
  fastainfo_df$primer_fw <- toupper(fastainfo_df$primer_fw)
  fastainfo_df$primer_rv <- toupper(fastainfo_df$primer_rv)
  # make a column for output filenames
  fastainfo_df$filename <- NA
  
  # check dirs
  outdir = check_dir(outdir)
  fasta_dir = check_dir(fasta_dir)
  
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
    trim_primers_file(input, 
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
                       compress_method=compress_method,
                       pigz_path=pigz_path,
                       quiet=quiet
                       )
    # count reads
    seq_n <- count_seq(output)
    fastainfo_df$read_count[i] <- seq_n
  }
  fastainfo_df <- fastainfo_df %>%
    select(sample, sample_type, habitat, replicate, fasta=filename, read_count)
  
  write.table(fastainfo_df, 
              file = file.path(outdir, "sampleinfo.csv"),  
              row.names = F, 
              sep=sep
              )
  return(fastainfo_df)
}

#' Demultiplex and trim tags and primers
#' 
#' Demultiplex each input FASTA file using tag combinations located at the 
#' ends of merged reads, then trim primer sequences from the resulting reads.
#'  
#' The output `sampleinfo.csv` file is similar to the input `fastainfo`, 
#' but without the tag and primer columns.
#'  
#' @param fastainfo Data frame or path to a CSV file with the following columns: 
#'   `tag_fw`, `primer_fw`, `tag_rv`, `primer_rv`, 
#'   `sample`, `sample_type` (mock/negative/real), 
#'   `habitat` (optional), `replicate`, `fasta`.
#' @param fasta_dir Character string specifying the directory containing input 
#'   FASTA files (listed in the `fasta` column of `fastainfo`).
#' @param vsearch_path Character string specifying the path to the 
#'   `vsearch` executable. 
#' @param cutadapt_path Character string specifying the path to the 
#'   `cutadapt` executable. 
#' @param compress_method Character or logical. Compression method: 
#'   `"pigz"`, `"gzip"`, or `"R"`.  
#'   `"pigz"` requires `pigz` to be installed and available in the 
#'   system PATH (or specified via `pigz_path`).  
#'   `"gzip"` is available on Linux systems.  
#'   `"R"` uses `R.utils`, which is cross-platform but slower.  
#'   Relative speed: `R.utils` < `gzip` < `pigz`.  
#'   Only needed if `check_reverse = TRUE`.
#' @param pigz_path Character string specifying the path to the `pigz` 
#'   executable. Only required if `compress_method = "pigz"` and it is not 
#'   available in the system PATH.
#' @param num_threads Positive integer specifying the number of CPU threads to 
#'   use. If `0`, all available CPUs are used.
#' @param outdir Character string specifying the output directory.
#' @param check_reverse Logical. If `TRUE`, also check reverse-complemented 
#'   sequences from the input FASTA files.
#' @param tag_to_end Logical. If `TRUE`, tags are expected to be located 
#'   at the extremities of reads (starting at the first base).
#' @param primer_to_end Logical. If `TRUE`, primers are assumed to follow 
#'   directly after tags (i.e., no heterogeneity spacer).
#' @param cutadapt_error_rate Numeric value between 0 and 1 specifying the 
#'   maximum allowed error rate between primers and reads (exact match is 
#'   required for tags).
#' @param cutadapt_minimum_length Positive integer specifying the minimum length 
#'   of trimmed sequences.
#' @param cutadapt_maximum_length Positive integer specifying the maximum length 
#'   of trimmed sequences.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param compress Logical. If `TRUE`, compress output files using gzip.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Data frame corresponding to the output `sampleinfo.csv` file 
#'   and one FASTA file per tag combination for each input FASTA file, 
#'   containing trimmed reads.
#' 
#' @examples
#' \dontrun{
#' fastainfo_df <- demultiplex_and_trim(
#'   fastainfo = fastainfo_df,
#'   fasta_dir = "data/fasta",
#'   outdir = "data/sorted",
#'   check_reverse = FALSE,
#'   tag_to_end = TRUE,
#'   primer_to_end = TRUE,
#'   cutadapt_minimum_length = 300,
#'   cutadapt_maximum_length = 350,
#'   sep = ","
#' )
#' }
#' 
#' @export
#' 
demultiplex_and_trim <- function(fastainfo, 
                      fasta_dir, 
                      outdir, 
                      cutadapt_path="cutadapt",
                      vsearch_path="vsearch", 
                      compress_method="R",
                      pigz_path="pigz",
                      num_threads=0,
                      check_reverse=FALSE, 
                      tag_to_end=TRUE, 
                      primer_to_end=TRUE, 
                      cutadapt_error_rate=0.1,
                      cutadapt_minimum_length=50,
                      cutadapt_maximum_length=500,
                      sep=",",
                      compress=FALSE,
                      quiet=T
                      ){
  
  fasta_dir = check_dir(fasta_dir)
  outdir = check_dir(outdir)
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  
  check_file_info(file=fastainfo_df, dir=fasta_dir, file_type="fastainfo", sep=sep, quiet=TRUE)
  
  #########
  # demultiplex_and_trim_strand_plus does the whole demultilexing, trimming and compress on the + strand
  # If sequences are not oriented, the -strand should be checked => 
  # run demultiplex_and_trim_strand_plus of plus strand and on - strand after 
  # swapping fw and rev tags and primers,
  # take the reverse complement of the -strand results (vsearch)
  # pool the results of the 2 strands
  # compress if necessary
  
  # run on strand +
  if(check_reverse){
    #### use +strand, output to sorted_dir, uncompressed
    sampleinfo_df <- demultiplex_and_trim_strand_plus(fastainfo_df, 
                                          fasta_dir=fasta_dir, 
                                          outdir=outdir, 
                                          cutadapt_path=cutadapt_path, 
                                          num_threads=num_threads,
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
    outdir = check_dir(outdir)
    rc_dir <- paste('rc_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    rc_dir <- file.path(tempdir(), rc_dir)
    rc_dir = check_dir(rc_dir)
    # run demultiplex_and_trim on for reverse strand
    sampleinfo_df <- demultiplex_and_trim_strand_plus(fastainfo_df_tmp, 
                                          fasta_dir=fasta_dir, 
                                          outdir=rc_dir, 
                                          cutadapt_path=cutadapt_path, 
                                          num_threads = num_threads,
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
        args <- c(
          "--fastx_revcomp", minus, 
          "--fastaout", minus_rc
        )
        if(quiet){
          args <- append(args, c("--quiet"))
        }
        run_system2(vsearch_path, args, quiet=quiet)
        
        # append content of minus_rc to plus file
        file.append(plus, minus_rc)
      }
    }
    
    # delete temporary reverse_comp dir
    unlink(rc_dir, recursive = TRUE)
    
    ### compress
    if(compress){
      
      for(i in 1:nrow(sampleinfo_df)){
        
        file <- sampleinfo_df$fasta[i]
        sampleinfo_df$fasta[i] <- paste(file, ".gz", sep="") # correct output filename
        file <- file.path(outdir, file) # add path
        file_gz <- smart_gzip(file,
                              remove = TRUE,
                              method = compress_method,,
                              pigz_path = pigz_path,
                              num_threads = num_threads,
                              quiet = quiet,
                              compress = TRUE)
      }
    }
  }
  else{
    # check only + strand
    sampleinfo_df <- demultiplex_and_trim_strand_plus(fastainfo_df, 
                                          fasta_dir=fasta_dir,
                                          outdir=outdir, 
                                          cutadapt_path=cutadapt_path, 
                                          num_threads = num_threads,
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
  
  sampleinfo_df <- add_read_counts(sampleinfo_df, dir=outdir)
  write.table(sampleinfo_df, file = file.path(outdir, "sampleinfo.csv"),  row.names = F, sep=sep)
  
  return(sampleinfo_df)
}

#' Add read_count to df
#' 
#' Count the number of reads in all FASTA files listed in the `fasta` column 
#' of a data frame, and add a `read_count` column to it.
#'  
#' @param df Data frame with a `fasta` column containing the names of FASTA files.
#' @param dir Character string specifying the directory containing the FASTA files.
#' 
#' @return Data frame with an additional `read_count` column.
#' 
#' @examples
#' \dontrun{
#' df <- add_read_counts(df, sorted_dir)
#' }
#' 
#' @export

add_read_counts <- function(df, dir){
  
  df$read_count <- NA
  
  fastas <- unique(df$fasta)
  
  for(file in fastas){
    file_path <- file.path(dir, file)
    read_n <- count_seq(file_path)
    df$read_count[which(df$fasta==file)] <- read_n
  }
return(df)
}

#' Demultiplex and trim tags and primers (no reverse strand check)
#' 
#' Same as `demultiplex_and_trim`, but without checking the reverse-complement 
#' of the sequences. Demultiplex each input FASTA file using tag combinations 
#' located at the extremities of merged reads, then trim primer sequences.
#'
#' Input files can be compressed or uncompressed. Output compression is 
#' controlled by `compress`.
#' 
#' The output `sampleinfo.csv` file is similar to the input `fastainfo` file, 
#' but without tag and primer columns.
#'  
#' @param fastainfo Data frame or path to a CSV file with the following columns: 
#'   `tag_fw`, `primer_fw`, `tag_rv`, `primer_rv`, 
#'   `sample`, `sample_type` (mock/negative/real), 
#'   `habitat` (optional), `replicate`, `fasta`.
#' @param fasta_dir Character string specifying the directory containing input 
#'   FASTA files (listed in the `fasta` column of `fastainfo`).
#' @param cutadapt_path Character string specifying the path to the 
#'   `cutadapt` executable.
#' @param num_threads Positive integer specifying the number of CPU threads to 
#'   use. If `0`, all available CPUs are used.
#' @param outdir Character string specifying the output directory.
#' @param tag_to_end Logical. If `TRUE`, tags are assumed to be located at the 
#'   extremities of reads (starting at the first base).
#' @param primer_to_end Logical. If `TRUE`, primers are assumed to follow 
#'   directly after tags (i.e., no heterogeneity spacer).
#' @param cutadapt_error_rate Numeric value between 0 and 1 specifying the 
#'   maximum allowed error rate between primers and reads (exact match is 
#'   required for tags).
#' @param cutadapt_minimum_length Positive integer specifying the minimum length 
#'   of trimmed sequences.
#' @param cutadapt_maximum_length Positive integer specifying the maximum length 
#'   of trimmed sequences.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param compress Logical. If `TRUE`, compress output files using gzip.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Data frame corresponding to the output `sampleinfo.csv` file and one 
#'   FASTA file per tag combination for each input FASTA file, containing 
#'   trimmed reads.
#' 
#' @examples
#' \dontrun{
#' fastainfo_df <- demultiplex_and_trim_strand_plus(
#'   fastainfo = fastainfo_df,
#'   fasta_dir = "data/fasta",
#'   outdir = "data/sorted",
#'   tag_to_end = TRUE,
#'   primer_to_end = TRUE,
#'   cutadapt_minimum_length = 300,
#'   cutadapt_maximum_length = 350,
#'   sep = ","
#' )
#' }
#' 
#' @export
#' 
demultiplex_and_trim_strand_plus <- function(fastainfo, 
                                 fasta_dir, 
                                 outdir, 
                                 cutadapt_path="cutadapt", 
                                 num_threads=0,
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
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  fasta_dir = check_dir(fasta_dir)
  outdir = check_dir(outdir)
  
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
    tag_file <- write_cutadapt_adapter_fasta(fastainfo_df, 
                                   fasta_file=fasta_file, 
                                   tag_to_end=tag_to_end, 
                                   outdir=tmp_dir_fasta
                                   )

     # add path
    fasta_file <- file.path(fasta_dir, fasta_file)
    # demultiplex fasta, write output to tmp file

    
    ##### run cmd
    g <- paste("file:", tag_file, sep="")
    out <- file.path(tmp_dir_fasta, "tagtrimmed-{name}.fasta") 
    args <- c(
      "-e", "0",
      "--no-indels",
      "--trimmed-only",
      "-g",  shQuote(g),
      "-o",  shQuote(out),
      fasta_file
    )
    if(num_threads > 0){
      args <- append(args, c("--cores", num_threads), after=2)
    }
    if(quiet){
      args <- append(args, c("--quiet"), after=2)
    }
    run_system2(cutadapt_path, args, quiet=quiet)
    
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
          g <- paste("^", primer_fwl, "...", primer_rvl_rc, "$", sep="")
        }
        else{
          g <- paste(primer_fwl, ";min_overlap=", nchar(primer_fwl),"...", primer_rvl_rc,";min_overlap=",nchar(primer_rvl_rc),sep="")
        }

        ##### run cmd
        args <- c(
          "-e", cutadapt_error_rate,
          "--no-indels",
          "--trimmed-only",
          "--minimum-length", cutadapt_minimum_length,
          "--maximum-length", cutadapt_maximum_length, 
          "-g", shQuote(g),
          "--output", primer_trimmed_file, 
          tag_trimmed_file
        )
        if(num_threads > 0){
          args <- append(args, c("--cores", num_threads), after=2)
        }
        if(quiet){
          args <- append(args, c("--quiet"), after=2)
        }
        run_system2(cutadapt_path, args, quiet=quiet)

      } # end tag-trimmed 
    # delete the tmp dir with the tag-trimmed files
    unlink(tmp_dir_fasta, recursive = TRUE)
  }# end fasta
  
  # make sampleinfo file
  fastainfo_df <- fastainfo_df %>%
    select(-fasta) %>%
    select(sample, sample_type, habitat, replicate, "fasta" = filename)
    
  return(fastainfo_df)
}

#' Make a FASTA file with adapters
#' 
#' Create a FASTA file containing tag combinations formatted for `cutadapt`. 
#' This file is used by `demultiplex_and_trim` to demultiplex input FASTA files.
#' 
#' @param fastainfo_df Data frame with columns: `tag_fw`, `tag_rv`, `fasta`.
#' @param fasta_file Character string specifying the FASTA file to be demultiplexed 
#'   (must be present in the `fasta` column of `fastainfo_df`).
#' @param outdir Character string specifying the output directory.
#' @param tag_to_end Logical. If `TRUE`, tags are assumed to be located at the 
#'   extremities of reads.
#' 
#' @return Path to the generated `tags.fasta` file. Returns `NA` if all tags 
#'   are `NA` in `fastainfo_df` for the given `fasta_file`.
#' 
#' @examples 
#' \dontrun{
#' write_cutadapt_adapter_fasta(
#'   fastainfo_df = fastainfo_df, 
#'   fasta_file = "test.fasta", 
#'   tag_to_end = FALSE, 
#'   outdir = "data/out"
#' )
#' }
#' 
#' @export
#' 
write_cutadapt_adapter_fasta <- function(fastainfo_df, fasta_file, outdir, tag_to_end=T){
  
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
  outdir = check_dir(outdir)
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


#' Reverse complement a sequence
#' 
#' Compute the reverse complement of a DNA sequence. IUPAC ambiguity codes are accepted.
#' 
#' @param sequence Character string specifying a DNA sequence.
#' 
#' @return Character string. The reverse complement of the input sequence.
#' 
#' @examples
#' reverse_complement(sequence = "AAATGCRC")
#' 
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

#' Read all FASTA files to a data frame and dereplicate
#' 
#' Read all FASTA files listed in the `fasta` column of a `sampleinfo` data 
#' frame (or CSV file), then dereplicate reads into ASVs. The number of reads 
#' per ASV is counted for each input file, and a unique `asv_id` is assigned 
#' to each ASV.
#'  
#' If `input_asv_list` is provided (containing previously observed ASV and 
#' `asv_id` pairs), existing `asv_id` values are reused when possible and new 
#' unique `asv_id` values are assigned to novel ASVs. If `output_asv_list` is 
#' provided, an updated file containing all `asv`–`asv_id` pairs is written.
#' 
#' @param sampleinfo Data frame or path to a CSV file with columns: 
#'   `sample`, `replicate`, `fasta`, and optionally `sample_type`, `habitat`. 
#'   The `fasta` column contains the names of FASTA files to be dereplicated.
#' @param dir Character string specifying the directory containing input FASTA files.
#' @param outfile Character string specifying the CSV output file containing the 
#'   resulting data frame (`asv_id`, `sample`, `replicate`, `read_count`). If 
#'   empty, no file is written.
#' @param input_asv_list Data frame or path to a CSV file containing previously 
#'   observed `asv`–`asv_id` pairs. Optional; used to harmonize `asv_id` 
#'   across datasets.
#' @param output_asv_list Character string specifying the output file containing 
#'   the updated `asv`–`asv_id` pairs. Optional.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Data frame with columns: `asv_id`, `sample`, `replicate`, 
#'   `read_count`, `asv`.
#' 
#' @examples
#' \dontrun{
#' dereplicate(sampleinfo = sampleinfo, dir = "data/sorted")
#' }
#' 
#' @export
#' 
dereplicate <- function(sampleinfo, 
                        dir, 
                        outfile=NULL, 
                        input_asv_list=NULL, 
                        output_asv_list=NULL, 
                        sep=",", 
                        quiet=T
                        ){
  # can accept df or file as an input
  if(is.character(sampleinfo)){
    # read known occurrences
    sampleinfo_df <- read.csv(sampleinfo, header=T, sep=sep)
  }else{
    sampleinfo_df <- sampleinfo
  }
  check_file_info(file=sampleinfo_df, 
                dir=dir, 
                file_type="sampleinfo", 
                sep=sep, 
                quiet=TRUE
                )
  
  # read all fasta files in sampleinfo to a read_count_df
  if(nchar(dir)>0){
    dir = check_dir(dir)
  }
  # define empty read_count_df to pool the results of variables
  read_count_df <- data.frame(asv=character(),
                              read_count=integer(),
                              sample=character(),
                              replicate=character())
  # read all fasta files in sampleinfo and count the reads
  for(i in 1:length(sampleinfo_df$fasta)){
    fas <- file.path(dir, sampleinfo_df$fasta[i])
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
    read_count_df_tmp$sample <- sampleinfo_df[i,"sample"]
    read_count_df_tmp$replicate <- sampleinfo_df[i,"replicate"]
    read_count_df <- rbind(read_count_df, read_count_df_tmp)
  }
  rm(read_count_df_tmp)
  # reorder columns
  read_count_df <- read_count_df[, c("asv", "sample", "replicate", "read_count")]
  
  # add input_asv_id and write output_asv_list if filename is given
  read_count_df <- add_ids(read_count_df, 
                           input_asv_list=input_asv_list, 
                           sep=sep,  
                           output_asv_list=output_asv_list
                           )
  
  # write read_count table
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Add numerical identifiers to ASVs
#' 
#' Add `asv_id` values to a data frame or CSV file containing an `asv` column. 
#' Previously existing `asv`–`asv_id` pairs (from earlier datasets) can be taken 
#' into account via `input_asv_list`.
#'  
#' If `output_asv_list` is provided, the `input_asv_list` is updated with newly 
#' observed ASVs and their corresponding `asv_id` values.
#' 
#' @param read_count Data frame or path to a CSV file with columns: 
#'   `asv`, `sample`, `replicate`, `read_count`.
#' @param input_asv_list Data frame or path to a CSV file containing `asv`–`asv_id` 
#'   pairs. Optional; used to harmonize `asv_id` values across datasets.
#' @param output_asv_list Character string specifying the output file containing 
#'   the updated `asv`–`asv_id` pairs. Optional.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Data frame with an added `asv_id` column.
#' 
#' @examples
#' \dontrun{
#' add_ids(read_count_df, input_asv_list = asv_list)
#' }
#' 
#' @export
#' 
add_ids <- function(read_count, 
                    input_asv_list=NULL, 
                    output_asv_list=NULL, 
                    sep=",", 
                    quiet=T
                    ){
  
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  ### read earlier asv_id - asv pairs
  if(is.character(input_asv_list)){
    if(input_asv_list == ""){ # no filename
      asv_df <- data.frame("asv_id"=integer(),
                           "asv"=as.character())
    }else{  # read known asv
      asv_df <- read.csv(input_asv_list, header=T, sep=sep)
    }
  }else if (is.null(input_asv_list)){
    asv_df <- data.frame("asv_id"=integer(),
                         "asv"=as.character())
  }else{
    asv_df <- input_asv_list
  }
  t <- check_one_to_one(asv_df) # stop execution, if FALSE
  

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
  
  # if input_asv_list should be updated, write it to a new file
  if(!is.null(output_asv_list)){
    write.table(asv_df, file=output_asv_list, row.names = FALSE, sep=sep)
  }
  
  return(read_count_df)
}


#' Read sequences from a FASTA file. 
#'  
#' Input FASTA files can be gzip-compressed or uncompressed.
#' 
#' @param filename Character string specifying the input FASTA file (including full path).
#' @param dereplicate Logical. If `TRUE`, return ASVs with read counts instead of individual reads.
#' 
#' @return Data frame. If `dereplicate = FALSE`, one read per row is returned. If 
#'   `dereplicate = TRUE`, a data frame with `asv` and `read_count` columns is returned.
#' 
#' @examples
#' \dontrun{
#' read_df <- read_fasta_seq(filename = "data/test.fasta", dereplicate = FALSE)
#' asv_df <- read_fasta_seq(filename = "data/test.fasta", dereplicate = TRUE)
#' }
#' 
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

#' Check one-to-one relationship
#' 
#' Check whether there is a one-to-one relationship between unique ASVs and 
#' unique `asv_id` values in the input data frame.
#'  
#' The same `asv`–`asv_id` combination may appear multiple times in the data frame.
#' 
#' @param df Data frame containing at least the columns `asv_id` and `asv` 
#' 
#' @return Logical. Returns `TRUE` if a one-to-one relationship exists between 
#'   `asv` and `asv_id`. Otherwise, the function stops execution with an error.
#' 
#' @examples
#' \dontrun{
#' check_one_to_one(df = read_count_df)
#' }
#' 
#' @export
#' 
check_one_to_one <- function(df){
  
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
    stop("Some of the the asv have to multile asv_id")
  }
  return(TRUE)
}

#' Update ASV list
#' 
#' Merge unique `asv`–`asv_id` pairs from two input data frames or CSV files. 
#' The function checks for consistency and ensures there are no conflicts 
#' within or between inputs.
#' 
#' The result is a complete, non-redundant list of `asv`–`asv_id` pairs.
#' 
#' @param asv_list1 Data frame or path to a CSV file with columns `asv_id`, `asv`; 
#' same format as `asv_list2`.
#' @param asv_list2 Data frame or path to a CSV file with columns `asv_id`, `asv`.
#' @param outfile Character string specifying the output CSV file name.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param return_df Logical. If `TRUE`, the function returns a data frame.
#' 
#' @return Data frame or CSV file containing all unique `asv_id`–`asv` pairs 
#'   from the inputs. If any conflict is detected within or between inputs, the 
#'   function stops with an error.
#' 
#' @examples
#' \dontrun{
#' update_asv_list(
#'   asv_list1 = read_count_df, 
#'   asv_list2 = "data/asv_list.csv", 
#'   outfile = "out/updated_asv_list.csv"
#' )
#' }
#' 
#' @export
#' 
update_asv_list <- function(asv_list1, asv_list2, outfile, sep=",", return_df=FALSE){
  
  if(is.character(asv_list1)){
    df1 <- read.csv(asv_list1, header=T, sep=sep)
  }else{
    df1 <- asv_list1
  }
  # select columns and make unique
  df1 <- df1 %>%
    select(asv_id, asv) %>%
    distinct()
  t <- check_one_to_one(df1)
  
  if(is.character(asv_list2)){
    df2 <- read.csv(asv_list2, header=T, sep=sep)
  }else{
    df2 <- asv_list2
  }  
  # select columns and make unique
  df2 <- df2 %>%
    select(asv_id, asv) %>%
    distinct()
  t <- check_one_to_one(df2)
  
  # pool 
  df1 <- rbind(df1, df2) %>%
    distinct() %>%
    arrange(asv_id)
  t <- check_one_to_one(df1)
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(df1, file=outfile, row.names = FALSE, sep=sep)
  }
  if(return_df){
    return(df1)
  }
}

#' filter_contaminant
#' 
#' Remove ASVs that are likely contaminants, defined as ASVs that have a higher 
#' read count in at least one negative control than in any biological sample.
#' 
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param sampleinfo Data frame or path to a CSV file with at least the following columns: 
#'   `sample`, `sample_type` (negative/mock/real).
#' @param outfile Character string specifying the CSV file to write the filtered 
#'   data frame. If NULL, no file is written.
#' @param conta_file Character string specifying a CSV file to store the filtered-out 
#'   contaminant ASVs. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count` data frame with contaminant ASVs removed.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_contaminant(read_count_df, sampleinfo = sampleinfo)
#' }
#' 
#' @export
#' 
filter_contaminant <- function (read_count, sampleinfo, outfile=NULL, 
                                       conta_file=NULL,sep=",") {
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  if(is.character(sampleinfo)){
    # read known occurrences
    sampleinfo <- read.csv(sampleinfo, header=T, sep=sep)
  }
  
  sampleinfo <- sampleinfo %>%
    filter(sample_type == "negative")
  
  negative_control_samples <- unique(sampleinfo$sample)
  
  # get for each asv the sample that have the highest read_count
  df <- read_count_df %>%
    group_by(asv_id) %>%
    arrange(desc(read_count)) %>%
    summarize(sample=first(sample)) %>%
    ungroup() %>%
    filter(sample %in% negative_control_samples)
  
  asv_conta <- unique(df$asv_id)
  
  # make a file with contaminants
  if(!is.null(conta_file)){
    check_dir(conta_file, is_file=TRUE)
    
    contaminant_df <- read_count_df %>%
      filter(asv_id %in% asv_conta) %>%
      group_by(asv_id) %>%
      arrange(asv_id, desc(read_count))
    
    write.table(contaminant_df, file = conta_file,  row.names = F, sep=sep)
  }

  # delete potential contaminants
  read_count_df <- read_count_df %>%
    filter(!(asv_id %in% asv_conta))

  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Filter out ASVs with Low Total Read Count
#' 
#' Remove ASVs with fewer than a given number of reads across the entire dataset.
#' 
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param cutoff Positive integer specifying the minimum total number of reads 
#'   required for an ASV to be retained. ASVs below this threshold are removed.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count` data frame.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_asv_global(read_count_df, cutoff = 2)
#' }
#' 
#' @export
#' 
filter_asv_global <- function (read_count, cutoff=10, outfile=NULL, sep=",") {
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
    filter(read_count_all >= cutoff) %>%
    ungroup()
  read_count_df <- filter(read_count_df, (asv %in% df$asv))
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}


#' Filter low-abundance occurrences using a read count cutoff
#' 
#' Remove occurrences (presence of an ASV in a sample-replicate) with fewer 
#' than a given number of reads.
#' 
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param cutoff Positive integer specifying the minimum number of reads required 
#'   for an occurrence to be retained. Occurrences below this threshold are removed.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count` data frame.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count <- filter_occurrence_read_count(read_count_df, cutoff = 20)
#' }
#' 
#' @export
#' 
filter_occurrence_read_count <- function (read_count, cutoff=10, outfile=NULL, sep=",") {
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  read_count_df <- filter(read_count_df,  (read_count >= cutoff))
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Filter low-abundance occurrences relative to total sample-replicate read count
#' 
#' 
#' Remove occurrences (presence of an ASV in a sample-replicate) where the 
#' relative abundance within the sample-replicate is below a given threshold.
#' Specifically, occurrences are removed when:
#' `read_count / sum(read_count in sample-replicate) < cutoff`.
#' 
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param cutoff Numeric value between 0 and 1 specifying the minimum proportion 
#'   of reads required for an occurrence to be retained. Occurrences below this 
#'   threshold are removed.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count` data frame.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_occurrence_sample(read_count_df, cutoff = 0.005)
#' }
#' 
#' @export
#' 
filter_occurrence_sample <- function (read_count, cutoff=0.001, outfile=NULL, sep=",") {
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
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Generate ASV-specific cutoff values
#' 
#' In some datasets, certain ASVs occur frequently across many samples with 
#' high read counts. In such cases, tag-jump or inter-sample contamination 
#' may persist even after applying a fixed threshold in 
#' `filter_occurrence_variant`. This function computes ASV-specific cutoff 
#' values to address such cases, focusing on ASVs with known false-positive 
#' occurrences (as identified by `classify_control_occurrences`).
#'    
#' For each ASV, the function:
#' - identifies all false-positive occurrences,
#' - takes the maximum read count among these false positives,
#' - divides this value by the total read count of the ASV across the dataset,  
#'   or within each replicate if `by_replicate = TRUE`.
#'   
#' Because some false positives may not be due to tag-jump contamination, the 
#' resulting cutoff values can be overly stringent. Therefore:
#' - filter the dataset as thoroughly as possible before computing cutoffs,
#' - set a reasonable upper limit using the `max_cutoff` parameter.
#'
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv`, `sample`, `replicate`, `read_count`.
#' @param max_cutoff Numeric value specifying the maximum allowed cutoff value.
#' @param mock_composition Data frame or path to a CSV file with columns: 
#'   `sample`, `action` (keep/tolerate), `asv`.
#' @param habitat_proportion Numeric value between 0 and 1. For each ASV, if the 
#'   proportion of reads in a habitat is below this threshold, it is considered 
#'   a false positive in all samples of that habitat.
#' @param by_replicate Logical. If `TRUE`, compute cutoffs separately for each replicate.
#' @param outfile Character string specifying the output file. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Data frame with columns: `asv_id`, `replicate` (if `by_replicate = TRUE`), `cutoff`.  
#' Only ASVs with known false-positive occurrences are included.
#' 
#' @seealso `filter_occurrence_variant()`
#' 
#' @examples
#' \dontrun{
#' cutoffs <- compute_asv_specific_cutoff(
#'   read_count = read_count_df, 
#'   mock_composition = "data/mock_composition.csv"
#' )
#' }
#' 
#' @export
#'
compute_asv_specific_cutoff <- function(read_count, 
                              max_cutoff=0.05,
                              mock_composition=NULL,
                              habitat_proportion=0.5,
                              by_replicate=FALSE, 
                              outfile=NULL, 
                              sep=",")  {
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  ### classify_control_occurrences to make known_occurrences_df
  results <- classify_control_occurrences(read_count_df, 
                                  sampleinfo=sampleinfo, 
                                  mock_composition=mock_composition,
                                  habitat_proportion=habitat_proportion,
                                  quiet=TRUE)
  
  known_occurrences_df <- results[[1]]
  
  # total number of read by asv, or asv.replicate
  if(by_replicate){
    asv_total_rc <- read_count_df %>%
      group_by(asv_id, replicate) %>%
      summarize(total_rc = sum(read_count), .groups="drop")
  }else{
    asv_total_rc <- read_count_df %>%
      group_by(asv_id) %>%
      summarize(total_rc = sum(read_count))%>%
      ungroup()
  }
  
  # list of asv.sample FP
  delete_occurrences_df <- known_occurrences_df %>%
    filter(action=="delete") %>%
    select(sample,asv_id) %>%
    distinct() %>%
    mutate(asv_sample = paste(asv_id, sample, sep="."))
  
  # 
  if(by_replicate){
    asv_spec_cutoff_df <- read_count_df %>%
      mutate(asv_sample = paste(asv_id, sample, sep=".")) %>%
      filter(asv_sample %in% delete_occurrences_df$asv_sample) %>%
      group_by(asv_id, replicate) %>%
      filter(read_count==max(read_count))%>%
      ungroup() %>%
      left_join(asv_total_rc, by=c("asv_id", "replicate")) %>%
      mutate(cutoff_asv_spec = round((read_count/total_rc)+0.00005, digits=4)) %>%
      select(asv_id, replicate, cutoff_asv_spec)
    
  }else{
    asv_spec_cutoff_df <- read_count_df %>%
      mutate(asv_sample = paste(asv_id, sample, sep=".")) %>%
      filter(asv_sample %in% delete_occurrences_df$asv_sample) %>%
      group_by(asv_id) %>%
      filter(read_count==max(read_count))%>%
      ungroup() %>%
      left_join(asv_total_rc, by=c("asv_id")) %>%
      mutate(cutoff_asv_spec = round((read_count/total_rc)+0.00005, digits=4)) %>%
      select(asv_id, cutoff_asv_spec)
  }
  
  # adjust too high values to max_cutoff
  asv_spec_cutoff_df <- asv_spec_cutoff_df %>%
    mutate(cutoff_asv_spec = if_else(
      cutoff_asv_spec > max_cutoff, max_cutoff, cutoff_asv_spec))
  
  # write to outfile
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(asv_spec_cutoff_df, file=outfile, row.names = F, sep=sep)
  }
  
  return(asv_spec_cutoff_df)
}

#' Filter low-abundance occurrences relative to total variant read count
#' 
#' Filter potential false positives arising from tag-jumps or low-level 
#' inter-sample contamination.
#'   
#' When `by_replicate = FALSE`, occurrences are removed when the ratio 
#' `read_count / total read count of the ASV in the dataset` is lower than `cutoff`.
#'   
#' When `by_replicate = TRUE`, occurrences are removed when the ratio 
#' `read_count / total read count of the ASV within the same replicate` is lower than `cutoff`.
#'   
#' A warning is issued if the total read count of an ASV after filtering falls 
#' below `min_read_count_prop` of its original abundance, as this may indicate 
#' that the cutoff value is too stringent.
#'   
#' By default, a single global cutoff value (`cutoff`) is applied to all ASVs. 
#' Alternatively, ASV-specific cutoff values can be provided via 
#' `asv_specific_cutoffs` (data frame or CSV file). When both are provided, the 
#' global `cutoff` is used when no ASV-specific cutoff is available or when the 
#' ASV-specific cutoff is lower than the global value.
#'
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param cutoff Numeric value between 0 and 1 specifying the minimum proportion 
#'   of reads required for an occurrence to be retained. Occurrences below this 
#'   threshold are removed.
#' @param asv_specific_cutoffs Data frame or path to a CSV file specifying ASV-specific 
#'   cutoff values. May include a `replicate` column if `by_replicate = TRUE`.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param lost_asv_file Character string specifying the CSV file used to store 
#'   information about ASVs whose retained read counts fall below 
#'   `min_read_count_prop` of their original abundance.
#' @param by_replicate Logical. If `TRUE`, compare read counts to ASV abundance 
#'   within each replicate.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param min_read_count_prop Numeric value specifying the minimum proportion of 
#'   total reads that must be retained per ASV after filtering. ASVs below this 
#'   threshold trigger a warning.
#' 
#' @return Filtered `read_count` data frame.
#' 
#' @examples
#' \dontrun{
#' read_count_df <- filter_occurrence_variant(
#'   read_count = "read_counts.csv",
#'   cutoff = 0.001,
#'   asv_specific_cutoffs = "asv_cutoffs.csv",
#'   by_replicate = FALSE
#' )
#' }
#' 
#' @export
#' 
filter_occurrence_variant <- function(read_count, 
                       cutoff=NULL, 
                       asv_specific_cutoffs = NULL,
                       outfile=NULL,
                       lost_asv_file =NULL,
                       by_replicate=FALSE, 
                       sep=",", 
                       min_read_count_prop=0.7){
  
  #### get read_count_df
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  ##### check coherence of parameters
  if(is.null(asv_specific_cutoffs)){
    if(is.null(cutoff)){
      stop("ERROR: cutoff and asv_specific_cutoffs are both NULL, Please, specify at least one of them.")
    }
  }else{
    # make asv_specific_cutoffs_df
    if (is.character(asv_specific_cutoffs)){
      # read known occurrences
      asv_specific_cutoffs_df <- read.csv(asv_specific_cutoffs, header=T, sep=sep)
    }else{
      asv_specific_cutoffs_df <- asv_specific_cutoffs
    }
    
    if(!("replicate" %in% colnames(asv_specific_cutoffs_df)) & by_replicate==TRUE){
      stop("ERROR: When by_replicate is TRUE, asv_specific_cutoffs should have a replicate column.")
    }
    if("replicate" %in% colnames(asv_specific_cutoffs_df) & by_replicate==FALSE){
      stop("ERROR: When by_replicate id FALSE, asv_specific_cutoffs should not have a replicate column.")
    }
  }
  
  
  #### input read count and sample count for later comparison
  asvs <- read_count_df %>%
    group_by(asv_id) %>%
    summarize("sample_count_input" = length(unique(sample)), "read_count_input"=sum(read_count)) %>%
    filter(read_count_input > 10) %>%
    ungroup()
  
  #### make df with asv total read count 
  if(by_replicate){
    sum_by_asv <- read_count_df %>%
      group_by(asv_id,replicate) %>%
      summarize(asv_sum = sum(read_count), .groups="drop")
  } else{
    sum_by_asv <- read_count_df %>%
      group_by(asv_id) %>%
      summarize(asv_sum = sum(read_count)) %>%
      ungroup()
  }
  
  #### Simple case of fixed cutoff
  if(is.null(asv_specific_cutoffs)){
    if(by_replicate){
      read_count_df <- left_join(read_count_df, sum_by_asv, by=c("asv_id", "replicate")) %>%
        filter(read_count/asv_sum >= cutoff) %>%
        select(-asv_sum)
    } else{
      read_count_df <- left_join(read_count_df, sum_by_asv, by=c("asv_id")) %>%
        filter(read_count/asv_sum >= cutoff)%>%
        select(-asv_sum)
    }
  }
  
  #### ASV specific cutoff
  # add  asv_specific_cutoffs_df to sum_by_asv
  if(!is.null(asv_specific_cutoffs)){
    
    # add cutoff_asv_spec from input asv_specific_cutoffs
    if(by_replicate){
      sum_by_asv <- left_join(sum_by_asv, asv_specific_cutoffs_df, by=c("asv_id", "replicate"))
    }else{
      sum_by_asv <- left_join(sum_by_asv, asv_specific_cutoffs_df, by=c("asv_id"))
    }
    # Change NA to 0 if no fixed cutoff
    if(is.null(cutoff)){ # no fix cutoff
      sum_by_asv <- sum_by_asv %>%
        mutate(cutoff_asv_spec = if_else(is.na(cutoff_asv_spec), 0, cutoff_asv_spec))
    }
    else{  # Change NA to cutoff if fixed cutoff 
      # and modify cutoff_asv_spec to fixed cutoff, it fixed cutoff if higher
      ### add fixed cutoff, when not specified in the input asv_specific_cutoffs_df
      sum_by_asv <- sum_by_asv %>%
        mutate(cutoff_asv_spec = if_else(is.na(cutoff_asv_spec), cutoff, cutoff_asv_spec)) %>%
        mutate(cutoff_asv_spec = if_else(cutoff_asv_spec<cutoff, cutoff, cutoff_asv_spec))
    }
    
    ### filter
    if(by_replicate){
      read_count_df <- left_join(read_count_df, sum_by_asv, by = c("asv_id", "replicate"), relationship = "many-to-many")
    }else{
      read_count_df <- left_join(read_count_df, sum_by_asv, by =c("asv_id"), relationship = "many-to-many")
    }
    read_count_df <- read_count_df %>%
      filter(read_count/asv_sum >= cutoff_asv_spec) %>%
      select(-cutoff_asv_spec, -asv_sum)
  }
  
  ###
  # Check if filter did not eliminate occurrences with relatively high read_count 
  ###
  asvs_output <- read_count_df %>%
    group_by(asv_id) %>%
    summarize("sample_count_output" = length(unique(sample)), 
              "read_count_output"=sum(read_count)) %>%
    ungroup()
  
  # join sample and read counts before and after filtering
  asvs <- left_join(asvs, asvs_output, by="asv_id") %>%
    mutate(sample_count_output = if_else(is.na(sample_count_output), 0, sample_count_output)) %>%
    mutate(read_count_output = if_else(is.na(read_count_output), 0, read_count_output)) %>%
    mutate(sample_prop = sample_count_output / sample_count_input) %>%
    mutate(read_count_prop = read_count_output / read_count_input) %>%
    filter(read_count_prop<min_read_count_prop) %>%
    arrange(read_count_prop, sample_prop) %>%
    select("asv_id", "read_count_input", "read_count_output", 
           "sample_count_input", "sample_count_output")
  
  if(nrow(asvs) > 0){
    cat("WARNING: The following ASVs have lost a high proportion of their 
          reads during this filtering step. 
          The cutoff value of filter_occurrence_variant function might need to be reduced.")
    print(asvs)
  }
  
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  
  if(!is.null(lost_asv_file)){
    check_dir(lost_asv_file, is_file=TRUE)
    write.table(asvs, file = lost_asv_file,  row.names = F, sep=sep)
  }
  
  return(read_count_df)
}


#' Retain occurrences shared across multiple filtered datasets
#' 
#' Pool multiple `read_count_df` data frames and retain only occurrences 
#' present in all filters.
#' 
#' @param ... Data frames with the following columns: 
#'   `asv_id`, `sample`, `replicate` (optional), `read_count`, `asv`.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count_df` data frame containing only shared occurrences.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- pool_filters(read_count_df1, read_count_df2)
#' }
#' 
#' @export
#' 
pool_filters <- function(... , outfile=NULL, sep=","){
  df_list <- list(...)
  merged <-  df_list[[1]]
  for(i in 2:length(df_list)){
    suppressMessages( merged <- inner_join(merged, df_list[[i]]) )
  }

  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(merged, file = outfile,  row.names = F, sep=sep)
  }
  return(merged)
}

#' Filter occurrences by minimum number of replicates per sample
#' 
#' Remove occurrences where an ASV is not present in at least `cutoff` 
#' replicates within a sample.
#'  
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param cutoff Positive integer specifying the minimum number of replicates 
#'   in which an ASV must be detected within a sample to be retained.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count` data frame.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_min_replicate(read_count_df, cutoff = 3)
#' }
#' 
#' @export
#'
filter_min_replicate <- function(read_count, cutoff=2, outfile=NULL, sep=","){
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
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Filter ASVs by frame consistency (indel detection)
#' 
#' Remove ASVs whose length is not compatible with the dominant reading frame 
#' of the dataset. Specifically, ASVs are filtered out if their length modulo 3 
#' differs from that of the majority of ASVs, which can indicate indel errors.
#'  
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count` data frame.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_indel(read_count_df)
#' }
#' 
#' @export
#' 
filter_indel <- function(read_count, outfile=NULL, sep=","){
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
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}


#' Codon stop set for a genetic code
#' 
#' Return the set of stop codons corresponding to a given genetic code number 
#' as defined by NCBI.
#'  
#' Genetic codes follow the NCBI translation tables:
#' https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
#'  
#' @param genetic_code Integer specifying the NCBI genetic code number.
#' 
#' @return Character vector of stop codons.
#' 
#' @examples
#' \dontrun{
#' get_stop_codons(genetic_code = 1)
#' }
#' 
#' @export
#' 
get_stop_codons <- function(genetic_code=5){
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

#' Filter ASVs containing stop codons
#' 
#' Remove ASVs that contain stop codons in all three reading frames of the 
#' forward (direct) strand, which typically indicates non-coding sequences or 
#' sequencing errors.
#'  
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate` (optional), `read_count`, `asv`.
#' @param genetic_code Positive integer specifying the NCBI genetic code number.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count` data frame.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_stop_codon(read_count_df, genetic_code = 5)
#' }
#' 
#' @export
#' 
filter_stop_codon <- function(read_count, outfile=NULL, genetic_code=5, sep=","){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  codon_stops <- get_stop_codons(genetic_code=genetic_code)
  if(length(codon_stops) == 0){
    print("WARNING: The Genetic Code Number provided does not correspond to any code in 
          https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes\n
          The filter_stop_codon step is skipped")
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
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Write FASTA file
#' 
#' Write a vector of sequences to a FASTA file, either using the sequences 
#' themselves as identifiers (`seq_as_id = TRUE`) or generating arbitrary 
#' sequence identifiers (`seq_as_id = FALSE`).
#'  
#' @param sequences Character vector of DNA sequences.
#' @param filename Character string specifying the output FASTA file name.
#' @param seq_as_id Logical. If `TRUE`, sequences are used as FASTA identifiers; 
#'   otherwise, default identifiers are assigned.
#' 
#' @return `NULL`. Writes a FASTA file to disk.
#' 
#' @examples
#' \dontrun{
#' write_fasta_vector(
#'   sequences = read_count_df$asv, 
#'   filename = "out/seq.fasta", 
#'   seq_as_id = TRUE
#' )
#' }
#' 
#' @export
#' 
write_fasta_vector <- function(sequences, filename, seq_as_id=F) {
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

#' Flag PCR errors using vsearch
#' 
#' Identify potential PCR errors: ASVs that are highly similar (within 
#' `max_mismatch`) to a more abundant ASV, based on abundance ratios in the 
#' input dataset.
#' 
#' A sequence is flagged as a potential PCR error when it is closely related to 
#' a more abundant ASV and its relative abundance is below `pcr_error_var_prop`.
#' A new column is added to the input data frame: `PCRerror` (1 = potential PCR 
#' error, 0 = otherwise).
#'  
#' @param unique_asv_df Data frame with the following columns: `asv`, `read_count`. 
#'   ASVs must be unique.
#' @param pcr_error_var_prop Numeric value between 0 and 1 specifying the 
#'   maximum allowed abundance ratio between similar ASVs. Less abundant ASVs 
#'   below this threshold are flagged as PCR errors.
#' @param max_mismatch Positive integer specifying the maximum number of 
#'   mismatches (including gaps) used to define similarity between ASVs.
#' @param vsearch_path Character string specifying the path to the 
#'   `vsearch` executable.
#' @param num_threads Positive integer specifying the number of CPU threads to 
#'   use. If `0`, all available CPUs are used.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Input data frame with an added `PCRerror` column 
#'   (1 = potential PCR error, 0 = otherwise).
#' 
#' @examples
#' \dontrun{
#' unique_asv_df <- read_count_df %>%
#'   dplyr::group_by(asv) %>%
#'   dplyr::summarise(read_count = sum(read_count))
#' 
#' unique_asv_df_flagged <- flag_pcr_error(
#'   unique_asv_df,
#'   vsearch_path = vsearch_path,
#'   pcr_error_var_prop = 0.2,
#'   max_mismatch = 2
#' )
#' }
#' 
#' @export
#' 
flag_pcr_error <- function(unique_asv_df,
                                 vsearch_path="vsearch", 
                                 num_threads=0,
                                 pcr_error_var_prop=0.1, 
                                 max_mismatch=1,
                                 quiet=TRUE
                                 ){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # no ASV in the unique_asv_df => return a dataframe with 0 for all ASVs in PCRerror column
  if(length(unique_asv_df$asv) == 0){ 
    unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
    return(unique_asv_df)
  }
  
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste('tmp_PCRerror_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  outdir_tmp <- file.path(tempdir(), outdir_tmp)
  outdir_tmp = check_dir(outdir_tmp)
  
  # make fasta file with unique reads; use sequences as ids
  fas <- file.path(outdir_tmp, 'unique.fas')
  write_fasta_vector(unique_asv_df$asv, fas, seq_as_id=T)
  # vsearch --usearch_global to find highly similar sequence pairs
  vsearch_out <- file.path(outdir_tmp, 'unique_vsearch_out.out')

  ##### run cmd
  args <- c(
    "--usearch_global", fas, 
    "--db", fas, 
    "--userout ",  vsearch_out,
    "--iddef", "1",
    "--self",
    "--id", "0.90",
    "--maxaccepts", 0,
    "--maxrejects", 0,
    "--userfields", shQuote("query+target+ids+aln") 
  )
  if(num_threads > 0){
    args <- append(args, c("--threads", num_threads))
  }
  if(quiet){
    args <- append(args, c("--quiet"))
  }
  run_system2(vsearch_path, args, quiet=quiet)
  
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

#' Filter PCR errors
#' 
#' Remove ASVs flagged as potential PCR errors based on sequence similarity 
#' (`max_mismatch`) and relative abundance (`pcr_error_var_prop`).
#'  
#' ASVs are considered PCR errors when they are highly similar to a more abundant 
#' ASV and their abundance ratio is below or equal to `pcr_error_var_prop`.
#'  
#' The analysis can be performed across the full dataset (`by_sample = FALSE`) 
#' or independently within each sample (`by_sample = TRUE`).
#'  
#' When `by_sample = TRUE`, an ASV is removed only if it is flagged as a PCR error 
#' in at least a proportion `sample_prop` of samples.
#'  
#' @param read_count Data frame or path to a CSV file with the following columns:  
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param pcr_error_var_prop Numeric value between 0 and 1 specifying the maximum 
#'   abundance ratio between similar ASVs for PCR error assignment. Less abundant 
#'   ASVs at or below this threshold are flagged as PCR errors.
#' @param max_mismatch Positive integer specifying the maximum number of mismatches 
#'   (including gaps) used to define sequence similarity.
#' @param by_sample Logical. If `TRUE`, PCR error detection is performed separately 
#'   within each sample.
#' @param sample_prop Numeric value between 0 and 1 specifying the minimum proportion 
#'   of samples in which an ASV must be flagged as a PCR error (when `by_sample = TRUE`) 
#'   to be removed.
#' @param min_read_count Positive integer specifying the minimum read count threshold; 
#'   occurrences below this value are ignored, to spead up the analyses.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param vsearch_path Character string specifying the path to the `vsearch` executable.
#' @param num_threads Positive integer specifying the number of CPU threads to use. 
#'   If `0`, all available CPUs are used.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param quiet Logical. If `TRUE`, suppress informational messages and only show 
#'   warnings or errors.
#' 
#' @return Filtered `read_count` data frame with PCR error ASVs removed.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_pcr_error(
#'   read_count_df,
#'   vsearch_path = vsearch_path,
#'   pcr_error_var_prop = 0.2,
#'   max_mismatch = 2,
#'   by_sample = TRUE,
#'   sample_prop = 0.8
#' )
#' 
#' filtered_read_count_df <- filter_pcr_error(
#'   read_count_df,
#'   vsearch_path = vsearch_path,
#'   pcr_error_var_prop = 0.2,
#'   max_mismatch = 2,
#'   by_sample = FALSE
#' )
#' }
#' 
#' @export
#' 
filter_pcr_error <- function(read_count,
                           outfile=NULL, 
                           vsearch_path="vsearch", 
                           num_threads=0,
                           pcr_error_var_prop=0.1,
                           max_mismatch=1, 
                           by_sample=TRUE, 
                           min_read_count=10,
                           sample_prop=0.8, 
                           sep=",",
                           quiet=TRUE
                           ){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  
  if(pcr_error_var_prop >= 1){
    stop("pcr_error_var_prop must be between 0-1.")
  }
  if(pcr_error_var_prop > 0.5){
    warning("pcr_error_var_prop above 0.5 is unusually high and may be unrealistic.")
  }
  
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
    filter(read_count >= min_read_count) %>%
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
        filter(read_count >= min_read_count) %>%
        arrange(desc(read_count)) %>%
        ungroup()
      
      # flag PCR errors; 
      # add one column to unique_asv_df for each sample with 1 if ASV is flagged in the sample, 
      # 0 otherwise
      unique_asv_df_sample <- flag_pcr_error(unique_asv_df_sample, 
                                                   vsearch_path=vsearch_path, 
                                                   num_threads=num_threads,
                                                   pcr_error_var_prop=pcr_error_var_prop, 
                                                   max_mismatch=max_mismatch,
                                                   quiet=quiet
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
    unique_asv_df <- flag_pcr_error(unique_asv_df, 
                                          vsearch_path=vsearch_path, 
                                          num_threads=num_threads,
                                          pcr_error_var_prop=pcr_error_var_prop, 
                                          max_mismatch=max_mismatch
                                          )
  }
  
  # count the number of times each ASV has been flagged and when it has not. 
  # Ignore NA, when the ASV is not present in the sample
  unique_asv_df$yes <- rowSums(unique_asv_df[3:ncol(unique_asv_df)] == 1, na.rm = TRUE)
  unique_asv_df$no <- rowSums(unique_asv_df[3:(ncol(unique_asv_df)-1)] == 0, na.rm = TRUE)
  # keep only ASVs, 
  # that are flagged in sample_prop proportion of the samples where they are present  
  unique_asv_df <- unique_asv_df %>%
    filter(yes/(yes+no) >= sample_prop)
  
  # eliminate potential PCRerrors from read_count_df
  read_count_df <- read_count_df %>%
    filter(!asv %in% unique_asv_df$asv)
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Flag chimeric sequences
#' 
#' Identify potential chimeric ASVs in a dataset of unique sequences.
#' A new column is added to the input data frame indicating whether each ASV 
#' is a putative chimera.
#'  
#' Chimeras are detected using abundance- and similarity-based criteria: an ASV 
#' is flagged when it can be explained as a combination of more abundant parent 
#' sequences under the specified `abskew` threshold.
#'  
#' @param unique_asv_df Data frame with the following columns: `asv`, `read_count`. 
#'   ASVs must be unique.
#' @param abskew Positive integer specifying the minimum abundance ratio required 
#'   between parent sequences and a chimera candidate.
#' @param vsearch_path Character string specifying the path to the `vsearch` executable.
#' @param quiet Logical. If `TRUE`, suppress informational messages and only show 
#'   warnings or errors.
#' @param num_threads Positive integer specifying the number of CPU threads to use. 
#'   If `0`, all available CPUs are used.
#' 
#' @return Input data frame with an added `chimera` column 
#'   (1 = potential chimera, 0 = otherwise).
#' 
#' @examples
#' \dontrun{
#' unique_asv_df <- read_count_df %>%
#'   dplyr::group_by(asv) %>%
#'   dplyr::summarise(read_count = sum(read_count))
#' 
#' flag_chimera(unique_asv_df, vsearch_path = vsearch_path, abskew = 2)
#' }
#' 
#' @export
#'
flag_chimera <- function(unique_asv_df, vsearch_path="vsearch", abskew=2, 
                        quiet=TRUE, num_threads=0){

  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # no ASV in the unique_asv_df => return a data frame with 0 for all ASVs in Chimera column
  if(length(unique_asv_df$asv) == 0){ 
    unique_asv_df$chimera <- rep(0, length(unique_asv_df$asv))
    return(unique_asv_df)
  }
  
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste('tmp_filter_chimera_', 
                      trunc(as.numeric(Sys.time())), 
                      sample(1:100, 1), 
                      sep=''
                      )
  outdir_tmp <- file.path(tempdir(), outdir_tmp)
  outdir_tmp = check_dir(outdir_tmp)
  
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
  
  vsearch_out <- file.path(outdir_tmp, "uchime3_denovo_out.tsv")
  ##### run uchime3_denovo
  # uchime3_denovo in vsearch v2.14.0 does not support multithreading 
  args <- c(
    "--uchime3_denovo", fas, 
    "--abskew", abskew,
    "--uchimeout", vsearch_out
  )
  if(quiet){
    args <- append(args, c("--quiet"))
  }
  run_system2(vsearch_path, args, quiet=quiet)
  
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

#' Filter chimeric sequences
#' 
#' Remove ASVs identified as chimeras based on abundance- and similarity-based 
#' detection using `vsearch`.
#'  
#' Chimeras are detected using the `abskew` parameter, which defines the minimum 
#' abundance ratio required between parental sequences and a chimera candidate.
#'  
#' Detection can be performed across the full dataset (`by_sample = FALSE`) or 
#' independently within each sample (`by_sample = TRUE`).
#'  
#' When `by_sample = TRUE`, an ASV is removed if it is flagged as a chimera in 
#' at least a proportion `sample_prop` of the samples in which it is present.
#'  
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param abskew Positive integer specifying the minimum abundance ratio used to 
#'   identify chimeric sequences.
#' @param by_sample Logical. If `TRUE`, chimera detection is performed separately 
#'   within each sample.
#' @param sample_prop Numeric value between 0 and 1 specifying the minimum proportion 
#'   of samples in which an ASV must be flagged as a chimera (when `by_sample = TRUE`) 
#'   to be removed.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param vsearch_path Character string specifying the path to the `vsearch` executable.
#' @param num_threads Positive integer specifying the number of CPU threads to use. 
#'   If `0`, all available CPUs are used.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param quiet Logical. If `TRUE`, suppress informational messages and only show 
#'   warnings or errors.
#' 
#' @return Filtered `read_count` data frame with chimeric ASVs removed.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_chimera(
#'   read_count_df,
#'   vsearch_path = vsearch_path,
#'   by_sample = TRUE,
#'   sample_prop = 0.7,
#'   abskew = 4
#' )
#' 
#' filtered_read_count_df <- filter_chimera(
#'   read_count_df,
#'   vsearch_path = vsearch_path,
#'   by_sample = FALSE,
#'   abskew = 4
#' )
#' }
#' 
#' @export
#' 
filter_chimera <- function(read_count, 
                          outfile=NULL, 
                          vsearch_path="vsearch",
                          num_threads=0,
                          by_sample=T, 
                          sample_prop=0.8, 
                          abskew=2, 
                          sep=",",
                          quiet=TRUE
                          ){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
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
      unique_asv_df_sample <- flag_chimera(unique_asv_df_sample, 
                                          vsearch_path=vsearch_path,
                                          num_threads = num_threads,  
                                          abskew=abskew,
                                          quiet=quiet
                                          )
      
      # remove read_count column
      unique_asv_df_sample <- select(unique_asv_df_sample, -c("read_count"))
      # add a column for each sample to unique_asv_df, with 1 if ASV is flagged in the sample,
      # 0 otherwise
      unique_asv_df <- left_join(unique_asv_df, unique_asv_df_sample, by = "asv")
    }
  }else{ # whole dataset
    # add a column to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
    unique_asv_df <- flag_chimera(unique_asv_df, 
                                 vsearch_path=vsearch_path, 
                                 num_threads = num_threads,  
                                 abskew=abskew,
                                 quiet=quiet)
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
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Calculate Renkonen distance
#' 
#' Compute the Renkonen distance between two ASV abundance profiles.
#'  
#' The Renkonen distance is based on the absolute differences in relative 
#' abundances of shared ASVs and ranges from 0 (identical compositions) to 1 
#' (completely different compositions).
#'  
#' @param df1,df2 Data frames with `asv` and `read_count` columns.
#' 
#' @return Numeric value between 0 and 1 representing the Renkonen distance 
#' between the two samples.
#' 
#' @examples
#' \dontrun{
#' df1 <- read_count_df %>%
#'   dplyr::filter(sample == "tpos1") %>%
#'   dplyr::group_by(asv) %>%
#'   dplyr::summarise(read_count = sum(read_count))
#' 
#' df2 <- read_count_df %>%
#'   dplyr::filter(sample == "tnegtag") %>%
#'   dplyr::group_by(asv) %>%
#'   dplyr::summarise(read_count = sum(read_count))
#' 
#' renkonen_dist(df1, df2)
#' }
#' 
#' @export
#'
renkonen_dist <- function(df1, df2){
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
  return(rdist)
}

#' Compute Renkonen distances
#' 
#' Calculate pairwise Renkonen distances between sample-replicates.
#' 
#' Distances can be computed either between all sample-replicate pairs 
#' (`compare_all = TRUE`) or only between replicates within the same sample 
#' (`compare_all = FALSE`).
#' 
#' @param read_count_df Data frame with `asv`, `sample`, `replicate`, and 
#'   `read_count` columns.
#' @param compare_all Logical. If `TRUE`, compute Renkonen distances for all 
#'   pairs of sample-replicates. If `FALSE`, compute distances only between 
#'   replicates of the same sample.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' 
#' @return Data frame with columns: `sample1`, `sample2`, `replicate1`, 
#'   `replicate2`, `renkonen_d`, `sample_comp` (indicating `"within"` if 
#'   sample1 equals sample2, `"between"` otherwise).
#' 
#' @examples
#' \dontrun{
#' compute_renkonen_distances(read_count_df, compare_all = FALSE)
#' }
#' 
#' @export
#' 
compute_renkonen_distances <- function(read_count_df, 
                                  compare_all=FALSE,
                                  outfile=NULL){
  
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
        rdist <- renkonen_dist(dfi, dfj)
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
          rdist <- renkonen_dist(dfi, dfj)
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
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(renkonen_df, file = outfile,  row.names = F, sep=sep)
  }
  return(renkonen_df)
}

#' Filter replicates based on Renkonen distance
#' 
#' Remove replicates that are inconsistent with other replicates of the same sample 
#' based on Renkonen distance. Replicates that show high dissimilarity to most other 
#' replicates are considered outliers and are filtered out.
#'  
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param cutoff Numeric value between 0 and 1 specifying the maximum acceptable 
#'   Renkonen distance for retaining replicates. Replicates exceeding this threshold 
#'   relative to most others in the same sample are removed.
#' @param renkonen_distance_quantile Numeric value between 0 and 1 specifying a 
#'   quantile-based cutoff. If `cutoff` is not provided, this quantile of the 
#'   Renkonen distance distribution is used to define the threshold (e.g., 0.9 
#'   corresponds to the 90th percentile).
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Filtered `read_count` data frame with inconsistent replicates removed.
#' 
#' @examples
#' \dontrun{
#' filtered_read_count_df <- filter_replicate(read_count_df, cutoff = 0.6)
#' filtered_read_count_df <- filter_replicate(
#'   read_count_df, 
#'   renkonen_distance_quantile = 0.9
#' )
#' }
#' 
#' @export
#' 
filter_replicate <- function(read_count, 
                           outfile=NULL,
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
  renkonen_df <- compute_renkonen_distances(read_count_df, compare_all=FALSE) %>%
    select("sample" = sample1, replicate1, replicate2, renkonen_d) %>%
    arrange(renkonen_d)
  
  # determine the cut off renkonen distance; values > cutoff are considered as high
  if(is.na(cutoff)){  
    last_row <- floor(length(renkonen_df$renkonen_d) * renkonen_distance_quantile)
    cutoff <- renkonen_df$renkonen_d[last_row]
  }
  msg <- paste("The cutoff value for Renkonen distances is ", cutoff)
  print(msg)
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
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)  
}

#' Pool replicates by sample
#' 
#' Aggregate replicates within each sample by combining read counts for each ASV.
#' For a given sample-ASV pair, read counts can be summarized using different 
#' strategies such as mean, sum, maximum, or minimum across replicates.
#'  
#' @param read_count Data frame or path to a CSV file with the following columns: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`, `cluster_id` (optional).
#' @param method Character string specifying the aggregation method used to pool 
#'   replicates. Must be one of `"mean"`, `"max"`, `"sum"`, or `"min"`.
#' @param digits Positive integer specifying the number of decimal places used 
#'   when rounding mean read counts.
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' 
#' @return Data frame with columns: `asv`, `sample`, `read_count` (aggregated 
#'   across replicates), and optional `cluster_id`.
#' 
#' @examples
#' \dontrun{
#' pool_replicates(read_count_df)
#' }
#' 
#' @export
#'
pool_replicates <- function(read_count, method="mean", digits=0, outfile=NULL, sep=","){
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  t <- check_one_to_one(read_count_df)
  
  # method
    method <- match.arg(method, c("mean", "max", "sum", "min"))
    fun <- switch(method,
                  mean = function(x) mean(x, na.rm = TRUE),
                  max  = function(x) max(x, na.rm = TRUE),
                  sum  = function(x) sum(x, na.rm = TRUE),
                  min  = function(x) min(x, na.rm = TRUE))
    
  
  # make df with unique asv_id cluster_id combinations
  if("cluster_id" %in% colnames(read_count_df)){
    cluster_df <- read_count_df %>%
      select(asv_id, cluster_id) %>%
      distinct()
    
    n <- length(unique(cluster_df$asv_id))
    if(n != nrow(cluster_df)){
      cat("WARNING: Some of the the asv have multile cluster_id. 
             This can happen if ASVs have been clustered sample by sample.
             cluster_id will not be written to the output.")
      read_count_df <- read_count_df %>%
        select(-cluster_id)
    }
  }
  
  read_count_samples_df <- read_count_df %>%
    group_by(asv_id,sample,asv) %>%
    summarize(read_count = fun(read_count), .groups="drop") %>%
    select(asv_id, sample, read_count, asv)
  
  read_count_samples_df$read_count <- round(read_count_samples_df$read_count, digits =digits)
  
  if("cluster_id" %in% colnames(read_count_df)){
    read_count_samples_df <- left_join(read_count_samples_df, cluster_df, by="asv_id")
  }
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_samples_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_samples_df)
}

#' Assign taxonomy using BLAST-based Lowest Common Ancestor (LTG) method
#' 
#' Assign the Lowest Taxonomic Group (LTG) for each ASV using a BLAST-based 
#' Lowest Common Ancestor approach.
#'  
#' This function implements the mkLTG algorithm described in 
#' Meglécz (2024) https://rdcu.be/dxABF and 
#' https://github.com/meglecz/mkLTG. It evaluates BLAST hits across multiple 
#' identity thresholds, each associated with specific filtering parameters.
#'  
#' For each identity threshold (`pid`), the following parameters are used:
#' * `pcov`: percentage of query coverage
#' * `phit`: proportion of validated hits required for LTG assignment
#' * `taxn`: minimum number of taxa among validated hits
#' * `seqn`: minimum number of sequences among validated hits
#' * `refres`: minimum resolution of validated hits
#' * `ltgres`: maximum resolution of the resulting LTG
#'  
#' @param asv Data frame or path to a CSV file containing at least `asv` and `asv_id` columns.
#' @param taxonomy TSV file containing the following columns: 
#'   `tax_id`, `parent_tax_id`, `rank`, `name_txt`, `old_tax_id` (merged tax IDs), 
#'   `taxlevel` (8: species ... 0: root).
#' @param blast_db Character string specifying the BLAST database name.
#' @param blast_path Character string specifying the path to the BLAST executable.
#' @param ltg_params Data frame or path to a CSV file defining identity thresholds 
#'   (`pid`) and associated parameters (`pcov`, `phit`, `taxn`, `seqn`, `refres`, `ltgres`).
#' @param outfile Character string specifying the CSV file to write the output 
#'   data frame. If NULL, no file is written.
#' @param fill_lineage Logical. If `TRUE`, missing higher taxonomic levels are 
#'   filled using the nearest known lower-level taxon with a prefix indicating 
#'   the missing rank.
#' @param num_threads Positive integer specifying the number of CPU threads to use. 
#'   If `0`, all available CPUs are used.
#' @param tax_sep Character string specifying the field separator used in the taxonomy file.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param quiet Logical. If `TRUE`, suppress informational messages and only show 
#'   warnings or errors.
#' 
#' @return Data frame with columns: 
#' `asv_id`, `ltg_taxid`, `ltg_name`, `ltg_rank`, `ltg_rank_index`,
#' `domain_taxid`, `domain`, `kingdom_taxid`, `kingdom`,
#' `phylum_taxid`, `phylum`, `class_taxid`, `class`, `order_taxid`, `order`,
#' `family_taxid`, `family`, `genus_taxid`, `genus`, `species_taxid`, `species`,
#' `pid`, `pcov`, `phit`, `taxn`, `seqn`, `refres`, `ltgres`, `asv`
#' 
#' @examples
#' \dontrun{
#' assign_taxonomy_ltg(
#'   asv = read_count_df,
#'   taxonomy = "xxxxxx",
#'   blast_db = "xxxxxxxxx",
#'   num_threads = 4
#' )
#' }
#' 
#' @export
#'
assign_taxonomy_ltg <- function(
    asv, 
    taxonomy, 
    blast_db, 
    blast_path="blastn", 
    ltg_params=NULL, 
    outfile=NULL, 
    fill_lineage=TRUE,
    num_threads=0, 
    tax_sep="\t", 
    sep=",",
    quiet=TRUE
    ){

  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
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
t <- check_one_to_one(asv_df)

if(is.null(ltg_params)){
  ltg_params_df = data.frame( pid=c(100,97,95,90,85,80),
                              pcov=c(70,70,70,70,70,70),
                              phit=c(70,70,70,70,70,70),
                              taxn=c(1,1,2,3,4,4),
                              seqn=c(1,1,2,3,4,4),
                              refres=c(8,8,8,7,6,6),
                              ltgres=c(8,8,8,8,7,7)
  )
} else if (is.character(ltg_params)){ 
  
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
outdir_tmp = check_dir(outdir_tmp)

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
lineages <- get_taxonomic_lineage(unique(blast_res$staxids), tax_df)

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
      ltg <- find_ltg(df_intern$staxids, lineages, phit = phitl)
      # fill out line with the ltg and the parameters that were used to get it
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
taxres_df <- adjust_ltg_resolution(taxres_df, tax_df)
# taxres_df data frame with the following columns: 
# asv_id,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,domain_taxid,
# domain,kingdom_taxid,kingdom,phylum_taxid,phylum,class_taxid,class,
# order_taxid,order,family_taxid,family,genus_taxid,genus,species_taxid,species,pid,
# pcov,phit,taxn,seqn,refres,ltgres,asv

# delete temporary  dir
unlink(outdir_tmp, recursive = TRUE)

if(!is.null(outfile)){
  check_dir(outfile, is_file=TRUE)
  write.table(taxres_df, file = outfile,  row.names = FALSE, sep=sep)
}

return(taxres_df)
}

#' Run BLAST
#' 
#' Perform BLAST searches using ASV sequences as query inputs.
#'  
#' This function runs BLAST and returns a table of matching hits for downstream 
#' taxonomic assignment or filtering.
#'  
#' @param df Data frame containing `asv` and `asv_id` columns.
#' @param blast_db Character string specifying the BLAST database (including path if needed).
#' @param outdir Character string specifying the output directory.
#' @param blast_path Character string specifying the path to the BLAST executable.
#' @param qcov_hsp_perc Real number between 0 and 100 specifying the minimum query coverage.
#' @param perc_identity Real number between 0 and 100 specifying the minimum percentage identity.
#' @param num_threads Positive integer specifying the number of CPU threads to use. 
#'   If `0`, all available CPUs are used.
#' @param quiet Logical. If `TRUE`, suppress informational messages and only show 
#'   warnings or errors.
#' 
#' @return Data frame containing BLAST results with columns: `qseqid`, `pident`, 
#'   `qcovhsp`, `staxids`.
#' 
#' @examples
#' \dontrun{
#' run_blast(
#'   df = read_count_df,
#'   blast_db = "xxxxxxxx",
#'   blast_path = "blastn",
#'   qcov_hsp_perc = 80,
#'   perc_identity = 90,
#'   num_threads = 4
#' )
#' }
#' 
#' @export
#'
run_blast <- function(df, 
                      blast_db, 
                      outdir, 
                      blast_path="blastn", 
                      qcov_hsp_perc=70, 
                      perc_identity=70, 
                      num_threads=0, 
                      quiet=T
                      ){

  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  outdir = check_dir(outdir)
  
  # make fasta file with unique reads; use numbers as ids  
#  seqs <- unique(df$asv)
  fas <- file.path(outdir, 'unique.fas')
  write_fasta_with_counts(df, outfile=fas, read_count=FALSE)
#  write_fasta_vector(seqs, fas, seq_as_id=T)
  # define the name of the output file
  blast_out <- file.path(outdir, 'blast.out')
  
  task = "megablast"
  e = 1e-20
  dust = "yes"
  max_target_seqs=500
  
  # Build argument vector
  args <- c(
    "-task", task,
    "-db", blast_db ,
    "-query", fas,
    "-evalue", e,
    "-out", blast_out,
    "-outfmt", shQuote("6 qseqid pident qcovhsp staxids"),
    "-dust", dust,
    "-qcov_hsp_perc", qcov_hsp_perc,
    "-perc_identity", perc_identity,
    "-max_target_seqs", max_target_seqs
  )
  if(num_threads > 0){
    args <- append(args, c("-num_threads", num_threads))
  }
  
  run_system2(blast_path, args, quiet=quiet)
  
  # read BLAST results; 
  # take care of lines where there is several taxids in the staxids column 
  # This can happen in ncbi nt
  blast_res <- read_blast_results(file=blast_out)
  return(blast_res)
}

#' Read BLAST results into a data frame
#' 
#' Import BLAST output into a structured data frame. If a hit contains multiple 
#' taxonomy IDs (`staxids`), each taxid is split into a separate row.
#'  
#' @param file Character string specifying the BLAST output file. The file must 
#'   be tab-separated and contain the columns: `qseqid`, `pident`, `qcovhsp`, `staxids`.
#' @return Data frame with columns: `qseqid`, `pident`, `qcovhsp`, `staxids`.
#' @examples
#' \dontrun{
#' read_blast_results(file = "blastout.txt")
#' }
#' 
#' @export
#'
read_blast_results <- function(file){
  
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

#' Expand rows by splitting taxid field
#' 
#' Split entries containing multiple taxonomy IDs into separate rows, ensuring 
#' a one-to-one mapping between hits and taxids.
#'  
#' @param row A data frame row containing the columns: `qseqid`, `pident`, 
#'   `qcovhsp`, `staxids`.
#' @return Data frame with one row per taxid, with columns: `qseqid`, `pident`, 
#'   `qcovhsp`, `staxids`.
#' @examples
#' \dontrun{
#' expand_rows(row)
#' }
#' 
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

#' Update obsolete taxids
#' 
#' Replace outdated NCBI taxonomy IDs with their current valid equivalents.
#' This step ensures that merged or deprecated taxids are mapped to updated identifiers.
#'  
#' @param df Data frame with the following columns: `qseqid`, `pident`, `qcovhsp`, `staxids`.
#' @param old_taxid Data frame with columns `tax_id` and `old_tax_id` used to map 
#'   deprecated taxids to valid ones.
#' @return Data frame with updated `staxids`, where obsolete taxids have been replaced 
#'   by current valid taxids.
#' @examples
#' \dontrun{
#' update_taxids(df = blastout_df, old_taxid = old_taxid_df)
#' }
#' 
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

#' Retrieve complete taxonomic lineages from taxids
#' 
#' Construct the full taxonomic lineage for each input taxid using a reference 
#' taxonomy table. Lineages are returned from the highest taxonomic level to 
#' the queried taxid.
#'  
#' @param taxids Vector of taxids (taxonomic identifiers).
#' @param tax_df Data frame with the following columns: `tax_id`, `parent_tax_id`, 
#'   `rank`, `name_txt`, `taxlevel` 
#'   (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 
#'   2: kingdom, 1: domain, 0: root).
#' @return Data frame where each row corresponds to an input taxid and contains 
#'   its full lineage as a sequence of taxids from the root to the queried taxid.
#' @examples
#' \dontrun{
#' taxids <- c(197147, 43823)
#' tax_df <- read.delim("xxxxxxxxxxxxxxx", header = TRUE, sep = "\t", fill = TRUE, quote = "")
#' get_taxonomic_lineage(taxids, tax_df = tax_df)
#' }
#' 
#' @export
#'
get_taxonomic_lineage <- function(taxids, tax_df){
  
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
  lineages <- as.data.frame(t(apply(lineages, 1, remove_leading_ones)))
  # add as a first column the taxid, so they can be easily accessed
  lineages <- cbind(taxids, lineages)
  
  return(lineages)
}

#' Remove leading ones from a vector
#' 
#' Remove all leading `1` values from a vector, shift remaining values to the left, 
#' and pad the end with `NA` to preserve the original length.
#'  
#' @param row Vector of taxids.
#' @return Vector of taxids with leading `1` values removed and replaced by `NA` at the end.
#' @examples
#' \dontrun{
#' remove_leading_ones(row)
#' }
#' 
#' @export
#'
remove_leading_ones <- function(row) {
  
  n <- length(row) 
  # Remove all occurrences of 1
  row <- row[row != 1]
  
  # Create a new row with NA at the end
  new_row <- c(row, rep(NA, n - length(row)))
  
  return(new_row)
}

#' Determine the Lowest Taxonomic Group (LTG)
#' 
#' Identify the Lowest Taxonomic Group (LTG) that contains at least `phit` percent 
#' of the input taxids.
#'  
#' @param taxids Vector of taxids (taxonomic identifiers), possibly containing duplicates.
#' @param lineages Data frame with taxids in the first column followed by their 
#'   taxonomic lineages (represented as taxids, starting from the lowest resolution).
#' @param phit Integer between 0 and 100 specifying the minimum percentage of taxids 
#'   that must be included in the LTG.
#' @return Numeric taxid corresponding to the LTG, or `NA` if no LTG can be determined.
#' @examples
#' \dontrun{
#' taxids <- c(189839, 1077837)
#' tax_df <- read.delim(xxxxxxxxxxxx, package = "vtamR"), header = TRUE, sep = "\t", fill = TRUE, quote = "")
#' lineages <- get_taxonomic_lineage(taxids, tax_df = tax_df)
#' find_ltg(taxids, lineages = lineages, phit = 80)
#' }
#' 
#' @export
#'
find_ltg <- function(taxids, lineages, phit=70){
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

#' Get ranked taxonomic lineages
#' 
#' Retrieve major taxonomic ranks for each input taxid based on a reference 
#' taxonomy table.
#'  
#' @param taxids Vector of taxids (taxonomic identifiers).
#' @param tax_df Data frame with the following columns: `tax_id`, `parent_tax_id`, 
#'   `rank`, `name_txt`, `taxlevel` 
#'   (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 
#'   2: kingdom, 1: domain, 0: root).
#' @param fill_lineage Logical. If `TRUE`, fill missing higher-level taxa using the 
#'   name of the next known lower-level taxon, prefixed by the corresponding rank 
#'   (e.g., `No_kingdom_Chrysophyceae` if kingdom is missing but class is known).
#' @return Data frame containing ranked lineages with columns: 
#'   `ltg_taxid`, `ltg_name`, `ltg_rank`, `ltg_rank_index`,
#'   `domain_taxid`, `domain`, `kingdom_taxid`, `kingdom`,
#'   `phylum_taxid`, `phylum`, `class_taxid`, `class`,
#'   `order_taxid`, `order`, `family_taxid`, `family`,
#'   `genus_taxid`, `genus`, `species_taxid`, `species`.
#' @examples
#' \dontrun{
#' taxids <- c(9593, 9606)
#' tax_df <- read.delim("xxxxxxxx", header = TRUE, sep = "\t", fill = TRUE, quote = "")
#' get_ranked_lineages(taxids, tax_df)
#' }
#' 
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
    ranked_lineages <- fill_missing_taxa(ranked_lineages)
  }
  
  return(ranked_lineages)
}

#' Fill missing higher-level taxa in lineages
#' 
#' Replace missing higher-level taxon names using the name of the next available 
#' lower-level taxon, prefixed by the corresponding taxonomic rank.
#'  
#' @param df Data frame containing taxonomic lineages. Major taxonomic levels are 
#'   located in every second column, starting from column 6 
#'   (e.g. 6: domain, 8: kingdom, 10: phylum, 12: class, 14: order, 
#'   16: family, 18: genus, 20: species).
#' @return Input data frame with missing taxon names filled using lower-level taxa, 
#'   prefixed by their corresponding rank.
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   id = 1:3,
#'   name = c("OTU1", "OTU2", "OTU3"),
#'   sample1 = c(10, 5, 0),
#'   sample2 = c(3, 7, 2),
#'   domain_taxid = c(2, 2759, 2759),
#'   domain = c("Bacteria", "Eukaryota", "Eukaryota"),
#'   phylum_taxid = c(1224, 4762, NA),
#'   phylum = c("Pseudomonadota", "Oomycota", NA),
#'   class_taxid = c(1236, NA, 2825),
#'   class = c("Gammaproteobacteria", NA, "Chrysophyceae")
#' )
#' df1 <- fill_missing_taxa(df)
#' }
#' 
#' @export
#'
fill_missing_taxa <- function(df) {
  
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

#' Adjust LTG resolution
#' 
#' Reduce the resolution of the Lowest Taxonomic Group (LTG) when it exceeds 
#' the level specified by `ltgres`, truncating the lineage at the desired rank.
#'  
#' @param taxres_df Data frame with the following columns: `asv_id`, `ltg_taxid`,
#'   `ltg_name`, `ltg_rank`, `ltg_rank_index`, `domain_taxid`,
#'   `domain`, `kingdom_taxid`, `kingdom`, `phylum_taxid`, `phylum`, 
#'   `class_taxid`, `class`, `order_taxid`, `order`, `family_taxid`, `family`, 
#'   `genus_taxid`, `genus`, `species_taxid`, `species`, `pid`,
#'   `pcov`, `phit`, `taxn`, `seqn`, `refres`, `ltgres`, `asv`.
#' @param tax_df Data frame with the following columns: `tax_id`, `parent_tax_id`, 
#'   `rank`, `name_txt`, `taxlevel` 
#'   (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 
#'   2: kingdom, 1: domain, 0: root).
#' @return Data frame with LTG resolution adjusted where necessary.
#' @examples
#' \dontrun{
#' tax_df <- read.delim("xxxxxxxx", header = TRUE, sep = "\t", fill = TRUE, quote = "")
#' adjust_ltg_resolution(taxres_df, tax_df)
#' }
#' 
#' @export
#'
adjust_ltg_resolution <- function(taxres_df, tax_df){
  
  # link taxlevel index and tax rank
  taxlevel_index = data.frame(taxlevel_index=c(0,1,2,3,4,5,6,7,8),
                              taxrank=c("root","domain","kingdom","phylum",
                                        "class","order","family","genus","species")
  )
  
  # add the name of the tax rank equivalent to the index in ltgref
  taxres_df <- left_join(taxres_df, taxlevel_index, by=c("ltgres" = "taxlevel_index"))
  
  for(i in 1:nrow(taxres_df)){ # all rows
    if(!is.na(taxres_df[i,"ltg_taxid"]) & 
       taxres_df[i,"ltg_rank_index"] > taxres_df[i,"ltgres"])
      { # if resolution of ltg is higher then ltgres
      # get the taxrank (name) that corresponds to ltgres 
      tl <- as.character(taxres_df[i,"taxrank"])
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
#' Generate and optionally write an ASV abundance table with samples as columns, 
#' ASVs as rows, and read counts as values.
#'  
#' @param read_count Data frame or CSV file with the following variables: 
#'   `asv_id`, `sample`, `replicate` (optional), `read_count`, `asv`, `cluster_id` (optional).
#' @param outfile Character string specifying the output CSV file. If NULL, no file is written.
#' @param asv_tax Data frame or CSV file containing taxonomic assignments. 
#'   Must include at least the columns `asv_id` and `asv`. Additional columns 
#'   describe the taxonomy of each ASV (e.g. ranks, identifiers, 
#'   or other annotations).
#'   If provided, taxonomic annotations are appended to the output.
#' @param sampleinfo Data frame or CSV file with columns: `sample`, `sample_type`.
#'   Required if `add_empty_samples = TRUE` or `add_expected_asv = TRUE`.
#' @param pool_replicates Logical. If `TRUE`, aggregate read counts across replicates 
#'   (see `method`). If `FALSE`, keep replicates separate (columns are `sample.replicate`).
#' @param method Character string specifying how replicate read counts are aggregated. 
#'   Must be one of `"mean"`, `"max"`, `"sum"`, or `"min"`.
#' @param add_empty_samples Logical. If `TRUE`, include all samples from the original dataset, 
#'   even if they have zero reads after filtering.
#' @param add_sums_by_sample Logical. If `TRUE`, add rows with total read counts and 
#'   number of ASVs per sample.
#' @param add_sums_by_asv Logical. If `TRUE`, add columns with total read counts per ASV 
#'   and the number of samples in which each ASV is present.
#' @param add_expected_asv Logical. If `TRUE`, add columns indicating expected ASVs 
#'   (e.g. in mock samples).
#' @param mock_composition Data frame or CSV file with columns: `sample`, `action`, `asv`. 
#'   `action` can be `keep` or `tolerate`. Required if `add_expected_asv = TRUE`.
#' @param sep Field separator character used in input and output CSV files.
#' 
#' @return Invisible data frame corresponding to the ASV table. Columns represent samples, 
#'   rows represent ASVs, and cells contain read counts, optionally extended with 
#'   taxonomic and summary information.
#' 
#' @examples
#' \dontrun{
#' write_asv_table(
#'   read_count_samples_df,
#'   outfile = "out/asv_table.csv",
#'   asv_tax = asv_tax,
#'   sampleinfo = sampleinfo_df,
#'   add_empty_samples = TRUE,
#'   add_sums_by_sample = TRUE,
#'   add_sums_by_asv = TRUE,
#'   add_expected_asv = TRUE,
#'   mock_composition = "data/mock_compostion.csv"
#' )
#' }
#' 
#' @export
#'
write_asv_table <- function(read_count, 
                          outfile=NULL, 
                          asv_tax=NULL, 
                          sampleinfo=NULL, 
                          pool_replicates=FALSE,
                          method="mean",
                          add_empty_samples=FALSE, 
                          add_sums_by_sample=FALSE, 
                          add_sums_by_asv=FALSE, 
                          add_expected_asv=FALSE,
                          mock_composition=NULL, 
                          sep=","
                          ){
  
  
  if(is.character(read_count)){
    read_count_samples_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_samples_df <- read_count
  }
  
  # check asv_id - asv
  t <- check_one_to_one(read_count_samples_df)
  

  #### Check if cluster info is available, and make df with unique asv_id an cluster_id, to add them to output
  add_cluster_id <- FALSE
  # make df with unique asv_id cluster_id combinations
  if("cluster_id" %in% colnames(read_count_samples_df)){
    cluster_df <- read_count_samples_df %>%
      select(asv_id, cluster_id) %>%
      distinct()
    add_cluster_id <- TRUE
    
    n <- length(unique(cluster_df$asv_id))
    if(n != nrow(cluster_df)){
      cat("WARNING: Some of the the asv have multile cluster_id. 
             This can happen if ASVs have been clustered sample by sample.
             cluster_id will not be written to the output.")
      add_cluster_id <- FALSE
    }
    # delete cluster_id
    read_count_samples_df <- read_count_samples_df %>%
      select(-cluster_id)
  }
  
  
  # read the sampleinfo to a data frame 
  if(add_empty_samples | add_expected_asv){
    if(!is.null(sampleinfo)){
      if(is.character(sampleinfo)){
        sampleinfo_df <- read.csv(sampleinfo, header=T, sep=sep)
      }else{
        sampleinfo_df <- sampleinfo
      }
    }
  }
  
  ### deal with replicates
  adjust_sample <- FALSE
  total_rc <- "sum_rc"
  if("replicate" %in% colnames(read_count_samples_df)){
    if(pool_replicates){ # take the mean read count of the replicates of the same sample
      read_count_samples_df <- pool_replicates(read_count_samples_df, method=method)
      total_rc <- paste("sum", method, "rc", sep="_")
    }else{ # make sample column and replace sample by sample.replicate 
      read_count_samples_df <- read_count_samples_df %>%
        mutate(sample = paste(sample, replicate, sep=".")) %>%
        select(-replicate)
      ### modify sample column in sampleinfo as well
      if(add_empty_samples | add_expected_asv){
        sampleinfo_df <- sampleinfo_df %>%
          mutate(sample = paste(sample, replicate, sep=".")) %>%
          select(-replicate)
      }
      adjust_sample <- TRUE
    }
  }
  
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


  if(add_empty_samples){
    # make vector with samples already in the data frame 
    # (asv_id and asv is also on the list, but it is not a pb)
    samples <- colnames(wide_read_count_df)
    # number of ASVs
    n <- nrow(wide_read_count_df)
    
    # make a vector with all unique samples in the sampleinfo
    all_samples <-unique(sampleinfo_df$sample)
    
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
    sum_rc[1,1] <- NA # asv_id col
    sum_rc[1,2] <- NA # asv col
    # total number of reads for each sample (ignore cols 1 and 2, since it is asv_id ans asv)
    sum_rc[1,-c(1,2)] <- colSums(wide_read_count_df[,-c(1,2)])
    # Number of ASVs in each sample in line 2
    sum_rc[2,1] <- NA
    sum_rc[2,2] <- NA # asv col
    sum_rc[2,-c(1,2)] <- colSums(wide_read_count_df[,-c(1,2)] != 0)
    wide_read_count_df <- rbind(sum_rc, wide_read_count_df)
  }
  
  # add sum of read count and the number of occurrences for each asv
  if(add_sums_by_asv){
    # count the number of reads for each ASV
    asvs <- read_count_samples_df %>%
      group_by(asv) %>%
      summarize(!!total_rc :=sum(read_count)) %>%
      ungroup()
    # add sample count to wide_read_count_df
    wide_read_count_df <- full_join(wide_read_count_df, asvs, by="asv")
    
    
    
    # count the number of samples where the ASV is present
    tmp <- read_count_samples_df %>%
      group_by(asv) %>%
      summarize(nb_samples=length(sample)) %>%
      ungroup()
    # add sample count to wide_read_count_df
    wide_read_count_df <- full_join(wide_read_count_df, tmp, by="asv")
  }
  
  if(add_expected_asv){
    
    # keep only mock samples in sampleinfo_df
    sampleinfo_df <- sampleinfo_df %>%
      filter(sample_type=="mock")
    # make a vector with all unique samples in the sampleinfo
    mock_samples <-unique(sampleinfo_df$sample)
    
    if(adjust_sample){ # samples are sample.replicate in wide_read_count_df
      mock_samples <- sub("\\..+$", "", mock_samples)
      mock_samples <- unique(mock_samples)
    }
    
    if(is.null(mock_composition)){
      stop("When add_expected_asv is TRUE, mock_composition must be provided")
    }
    
      if(is.character(mock_composition)){
        mock_asv <-  read.csv(mock_composition, header=T, sep=sep)
        check_file_info(file=mock_composition, file_type="mock_composition", sep=sep, quiet=TRUE)
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
  # add cluster_id
  if(add_cluster_id){
    wide_read_count_df <- left_join(wide_read_count_df, cluster_df, by="asv_id")
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
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(wide_read_count_df, file=outfile, row.names = F, sep=sep)
  }
  return(invisible(wide_read_count_df))
}

#' Suggest PCR error cutoff
#' 
#' Identify pairs of expected and unexpected ASVs in mock samples that differ 
#' by at most `max_mismatch`, and compute their read count ratios to help 
#' determine an appropriate `pcr_error_var_prop` threshold.
#' 
#' The suggested cutoff should be set above the maximum observed ratio 
#' (`unexpected_read_count / expected_read_count`) in the output table.
#' 
#' @param read_count Data frame or CSV file with the following variables: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param mock_composition Data frame or CSV file with columns: 
#'   `sample`, `action` (`keep` or `tolerate`), `asv`.
#' @param vsearch_path Character string specifying the path to vsearch executables.
#' @param num_threads Positive integer specifying the number of CPU threads to use. 
#'   If `0`, all available CPUs are used.
#' @param sep Field separator character used in input and output CSV files.
#' @param outfile Character string specifying the output CSV file. If NULL, 
#'   no file is written.
#' @param max_mismatch Positive integer specifying the maximum number of mismatches 
#'   allowed between ASVs to be compared.
#' @param min_read_count Positive integer specifying the minimum read count threshold; 
#'   occurrences below this value are ignored.
#' @param quiet Logical. If `TRUE`, suppress informational messages and only 
#'   show warnings or errors.
#' @return Data frame with the following columns: `sample`, `expected_read_count`,
#'   `unexpected_read_count`, `pcr_error_var_prop`, `expected_asv_id`, 
#'   `unexpected_asv_id`, `expected_asv`, `unexpected_asv`.
#' @examples
#' \dontrun{
#' suggest_pcr_error_cutoff(
#'   read_count = read_count_df,
#'   mock_composition = "data/mock_composition.csv",
#'   vsearch_path = vsearch_path,
#'   max_mismatch = 2,
#'   min_read_count = 5
#' )
#' }
#' 
#' @export
#'
suggest_pcr_error_cutoff <- function(read_count, 
                             mock_composition, 
                             vsearch_path= "vsearch", 
                             num_threads=0,
                             sep=",", 
                             outfile=NULL, 
                             max_mismatch=1, 
                             min_read_count=10,
                             quiet=TRUE
                             ){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  
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
  check_file_info(file=mock_composition_df, file_type="mock_composition", sep=sep, quiet=TRUE)
  
  # read the mock composition file and keep only lines with keep and tolerate
  mock_composition_df <- mock_composition_df %>%
    filter(action=="keep" | action=="tolerate")
  unique_mock_list <- unique(mock_composition_df$sample)
  
  
  #### Test if all expected (keep) variants are present in read_count_df
  # Extract unique ASVs from the mock composition
  unique_mock_asv <- mock_composition_df %>%
    select(asv) %>%
    distinct()
  
  # Keep ASVs that are present in the input data frame
  present <- read_count_df %>%
    select(asv_id, asv) %>%
    filter(asv %in% unique_mock_asv$asv) %>%
    distinct()
  
  # Identify ASVs from the mock that are missing in the input data
  missing <- unique_mock_asv %>%
    filter(!(asv %in% present$asv))
  
  # If any ASVs are missing, issue a warning (not an error)
  if (nrow(missing) > 0) {
    warning(
      "The following expected ASVs are not present in the read_count dataframe:\n",
      paste(missing$asv, collapse = ", "),
      "\nPlease check whether these sequences are correct."
    )
  }
  
  
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
    outdir_tmp <- paste('tmp_suggest_pcr_error_cutoff_', 
                        trunc(as.numeric(Sys.time())), 
                        sample(1:100, 1), 
                        sep=''
                        )
    outdir_tmp <- file.path(tempdir(), outdir_tmp)
    outdir_tmp = check_dir(outdir_tmp)
    # get the list of keep ASV in the given mock sample from mock_composition
    tmp_mock <- mock_composition_df %>%
      filter(sample==mock) %>%
      filter(action=="keep")
    asv_list_keep <- unique(tmp_mock$asv)
    # make fasta file with unique mock variants; use sequences as ids
    fas_keep <- file.path(outdir_tmp, paste(mock, "keepASV.fas", sep="_"))
    write_fasta_vector(asv_list_keep, fas_keep, seq_as_id=T)
    
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
    write_fasta_vector(asv_list_delete, fas_delete, seq_as_id=T)
    
    if(length(asv_list_delete)>0 && length(asv_list_keep)>0){ # sequences in both files
      # vsearch --usearch_global to find highly similar sequence pairs
      vsearch_out <- paste(mock, 'vsearch_out.out', sep="_")
      vsearch_out <- file.path(outdir_tmp, vsearch_out)

      ##### run cmd
      args <- c(
        "--usearch_global", fas_delete,
        "--db", fas_keep, 
        "--iddef", 1,
        "--self",
        "--id", 0.90,
        "--maxaccepts", 0,
        "--maxrejects", 0,
        "--userfields", shQuote("query+target+ids+aln"),
        "--userout", vsearch_out
      )
      if(num_threads > 0){
        args <- append(args, c("--threads", num_threads))
      }
      if(quiet){
        args <- append(args, c("--quiet"))
      }
      run_system2(vsearch_path, args, quiet=quiet)
      

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
        # delete unnecessary columns and add sample
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
  
  if(!is.null(outfile))
  {
    check_dir(outfile, is_file=TRUE)
    write.table(asv_pairs, file=outfile, sep=sep, row.names = F)
  }
  return(asv_pairs)
}

#' Suggest cutoff for `filter_occurrence_sample`
#' 
#' Compute read count proportions of expected ASV occurrences in mock 
#' sample-replicates to help determine an appropriate cutoff for 
#' `filter_occurrence_sample`.
#' 
#' The cutoff should be set below the smallest observed proportion to ensure 
#' that all expected ASVs are retained in the dataset.
#'  
#' @param read_count Data frame or CSV file with the following variables: 
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param mock_composition Data frame or CSV file with columns: 
#'   `sample`, `action` (`keep` or `tolerate`), `asv`.
#' @param sep Field separator character used in input and output CSV files.
#' @param outfile Character string specifying the output CSV file. If NULL, 
#'   no file is written.
#' @return Data frame with the following columns: `sample`, `replicate`, `action`, 
#'   `asv_id`, `read_count`, `read_count_sample_replicate`, `sample_cutoff`, `asv`.
#' @examples
#' \dontrun{
#' suggest_sample_cutoff(
#'   read_count_df,
#'   mock_composition = "data/mock_composition.csv"
#' )
#' }
#' 
#' @export
#'
suggest_sample_cutoff <- function(read_count, mock_composition, sep=",", outfile=NULL){
  
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
  check_file_info(file=mock_composition_df, file_type="mock_composition", sep=sep, quiet=TRUE)
  
  
  #### Test if all expected (keep) variants are present in read_count_df
  # Extract unique ASVs from the mock composition
  unique_mock_asv <- mock_composition_df %>%
    select(asv) %>%
    distinct()
  
  # Keep ASVs that are present in the input data frame
  present <- read_count_df %>%
    select(asv_id, asv) %>%
    filter(asv %in% unique_mock_asv$asv) %>%
    distinct()
  
  # Identify ASVs from the mock that are missing in the input data
  missing <- unique_mock_asv %>%
    filter(!(asv %in% present$asv))
  
  # If any ASVs are missing, issue a warning (not an error)
  if (nrow(missing) > 0) {
    warning(
      "The following expected ASVs are not present in the read_count dataframe:\n",
      paste(missing$asv, collapse = ", "),
      "\nPlease check whether these sequences are correct."
    )
  }
  
  
  #########
  # read the mock composition file and keep only lines with keep
  mock_composition_df <- mock_composition_df %>%
    filter(action=="keep")
  # there is a asv_id column in mock_composition => delete it
  if("asv_id" %in% colnames(mock_composition_df )){ 
    mock_composition_df <- mock_composition_df %>%
    select(-asv_id)
  }
  
  read_count_mock_df <- read_count_df %>%
    filter(sample %in% mock_composition_df$sample)

  # get the unique list of asv and asvid (it will be use to add asv_id at the end,even if asv not in mock)
  asv_ids <- read_count_mock_df %>%
    select(asv_id, asv) %>%
    distinct()
  
  # get a complete and unique list of sample, replicate of mocks
  sample_replicate_list <- read_count_mock_df %>%
    select(sample, replicate) %>%
    distinct()
  
  # add replicate to mock_composition
  mock_composition_df <- left_join(mock_composition_df, 
                                   sample_replicate_list, 
                                   by=c("sample"), 
                                   relationship = "many-to-many"
                                   )
  
  # get the total number of reads for each sample replicate for the mocks
  sample_replicate_rc <- read_count_mock_df %>%
    group_by(sample, replicate) %>%
    summarize(read_count_sample_replicate= sum(read_count), .groups="drop_last") %>%
    ungroup()
  
  # add read_count
  asv_keep_df <- left_join(mock_composition_df, 
                           read_count_mock_df, 
                           by=c("sample", "replicate", "asv")
                           )
  # add sum read_counts over replicates 
  asv_keep_df <- left_join(asv_keep_df, 
                           sample_replicate_rc, 
                           by=c("sample", "replicate")
                           )
  # rm asv_id and add it from asv_ids: if ASV is not prensent in mock, but present
  # in another, the asv_id will not be NA
  asv_keep_df <- asv_keep_df %>%
    select(-asv_id) %>%
    left_join(asv_ids, by= "asv") %>%
    mutate(read_count = if_else(is.na(read_count), 0, read_count))
  
  asv_keep_df$sample_cutoff <- 
    asv_keep_df$read_count/asv_keep_df$read_count_sample_replicate
  asv_keep_df$sample_cutoff <- 
    round(asv_keep_df$sample_cutoff-0.00005, digits=4)
  
  asv_keep_df <- asv_keep_df %>%
    arrange(sample_cutoff) %>%
    select(sample, 
           replicate, 
           action, 
           asv_id, 
           sample_cutoff,
           read_count,
           read_count_sample_replicate, 
           asv, 
           everything()
           )
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(asv_keep_df, file=outfile, sep=sep, row.names = F)
  }
  return(asv_keep_df)
}

#' Classify control occurrences and compute performance metrics
#' 
#' Identify expected and unexpected occurrences in control samples and summarize 
#' pipeline performance using true positives, false positives, and false negatives.
#' 
#' This function also computes key evaluation metrics including precision and sensitivity.
#'  
#' @param read_count Data frame or CSV file with the following variables: 
#'   `asv_id`, `sample`, `replicate` (optional), `read_count`, `asv`.
#' @param sampleinfo Data frame or CSV file with columns: 
#'   `sample`, `sample_type` (`mock`/`negative`/`real`), `habitat`, `replicate` (optional), 
#'   and optionally `file`.
#' @param mock_composition Data frame or CSV file with columns: 
#'   `sample`, `action` (`keep`/`tolerate`), `asv`.
#' @param sep Field separator character used in input and output CSV files.
#' @param known_occurrences Character string specifying output CSV file for known occurrences 
#'   (expected occurrences in mocks and false positives). If NULL, no file is written.
#' @param false_negatives Character string specifying output CSV file for false negatives. 
#'   If NULL, no file is written.
#' @param performance_metrics Character string specifying output CSV file for performance metrics. 
#'   If NULL, no file is written.
#' @param habitat_proportion Numeric between 0 and 1. For each ASV, if the proportion of reads 
#'   within a habitat is below this threshold, it is considered an artifact in all samples 
#'   of that habitat.
#' @param quiet Logical. If `TRUE`, suppress informational messages and only show warnings or errors.
#' 
#' @return A list containing:
#'   * `known_occurrences_df`: sample, action, asv_id, asv
#'   * `false_negatives_df`: sample, action, asv, asv_id
#'   * `performance_metrics_df`: TP, FP, FN, Precision, Sensitivity
#' 
#' @examples
#' \dontrun{
#' results <- classify_control_occurrences(
#'   read_count_samples_df,
#'   sampleinfo = sampleinfo_df,
#'   mock_composition = "data/mock_composition.csv",
#'   habitat_proportion = 0.7
#' )
#' known_occurrences_df <- results[[1]]
#' false_negatives_df <- results[[2]]
#' performance_metrics_df <- results[[3]]
#' }
#' 
#' @export
#' 

classify_control_occurrences <- function(read_count, 
                                 sampleinfo, 
                                 mock_composition, 
                                 sep=",", 
                                 known_occurrences=NULL, 
                                 false_negatives=NULL, 
                                 performance_metrics=NULL, 
                                 habitat_proportion=0.5,
                                 quiet=TRUE){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # pool_replicates before counting occurrences in samples
  if("replicate" %in% colnames(read_count_df)){  
    check_file_info(file=read_count_df, 
                  file_type="read_count", 
                  quiet=TRUE)
    read_count_samples_df <- pool_replicates(read_count_df)
  }else{
    check_file_info(file=read_count_df, 
                  file_type="read_count_sample", 
                  quiet=TRUE)
    read_count_samples_df <- read_count_df
  }
 
  
  if(is.character(sampleinfo)){
    sampleinfo_df <- read.csv(sampleinfo, header=T, sep=sep)
  }else{
    sampleinfo_df <- sampleinfo
  }

  if(is.character(mock_composition)){
    mock_composition_df <- read.csv(mock_composition, header=T, sep=sep)
  }else{
    mock_composition_df <- mock_composition
  }
  check_file_info(file=mock_composition_df, 
                file_type="mock_composition", 
                quiet=TRUE)
  
  # read info on samples types and keep only relevant info
  sampleinfo_df <- sampleinfo_df %>%
    select(sample, sample_type, habitat)
  # get unique lines to avoid replicates
  sampleinfo_df <- unique(sampleinfo_df)
  
  # define data frame for known occurrences
  occurrence_df <- read_count_samples_df
  # add habitat and sample_type to occurrence_df
  occurrence_df <- left_join(occurrence_df, sampleinfo_df, by="sample")
  # add action column
  occurrence_df$action <- rep(NA, nrow(occurrence_df))
  
  # flag occurrences in negative control samples as delete
  occurrence_df$action[which(occurrence_df$sample_type=="negative")] <- "delete"
  # flag all expected occurrences in mock samples as "keep", NA for tolerate, and delete for all others
  occurrence_df <- flag_mock_asv(occurrence_df, mock_composition_df, sep=sep)
  # flag occurrences as delete with low read count in habitat, compared to the other habitats
  occurrence_df <- flag_by_habitat(occurrence_df, 
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
  
  # count the number of FN and write false_negatives, if filename is defined
  missing_occurrence_df <-detect_false_negatives(read_count_samples=read_count_samples_df,
                                                   mock_composition=mock_composition_df, 
                                                   sep=sep, 
                                                   out=false_negatives,
                                                   quiet=quiet
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
  
  
  # write to outfiles (missing is written by function detect_false_negatives)
  if(!is.null(known_occurrences)){
    check_dir(known_occurrences, is_file=TRUE)
    write.table(occurrence_df, file=known_occurrences, row.names = F, sep=sep)
  }
  if(!is.null(performance_metrics)){
    check_dir(performance_metrics, is_file=TRUE)
    write.table(count_df, file=performance_metrics, row.names = F, sep=sep)
  }
  df_list <- list(occurrence_df, missing_occurrence_df, count_df)
  return(df_list)
}


#' Flag occurrences in mock samples
#' 
#' Classify all occurrences observed in mock samples as expected or unexpected
#' based on a reference mock composition table.
#' Expected variants are flagged as `"keep"`, unexpected ASVs as `"delete"`,
#' and `"tolerate"` ASVs are left as `NA`.
#' 
#' Tolerate ASVs correspond to taxa that may be present in mock samples but
#' should not influence optimization of filtering thresholds
#' (e.g. poorly amplified taxa included in the mock design).
#' 
#' @param occurrence_df Data frame with columns `asv_id`, `sample`, `read_count`,
#'   `asv`, `sample_type`, `habitat`, `action`.
#' @param mock_composition Data frame or CSV file with columns `sample`,
#'   `action` (keep/tolerate), and `asv`.
#' @param sep Field separator character used in input and output CSV files.
#' @return Data frame with columns `asv_id`, `sample`, `read_count`, `asv`,
#'   `sample_type`, `habitat`, `action`.
#' @examples
#' \dontrun{
#' flag_mock_asv(
#'   occurrence_df = occurrence_df,
#'   mock_composition = "data/mock_composition.csv"
#' )
#' }
#' @export
#' 
flag_mock_asv <- function(occurrence_df, mock_composition, sep=","){
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
#' Identify and flag false-positive occurrences based on the distribution of
#' ASVs across habitats. 
#'  
#' ASVs present in multiple habitats are evaluated: 
#' For each of these ASVs, the proportion of reads within each habitat is computed.
#' If this proportion is below `habitat_proportion`, the occurrence is
#' considered a likely artifact and flagged as `"delete"` in all samples of
#' that habitat.
#' 
#' @param occurrence_df Data frame with columns `asv`, `sample`, `read_count`,
#'   `sample_type`, `habitat`, `action`.
#' @param habitat_proportion Numeric value between 0 and 1. Minimum proportion
#'   of reads required for an ASV to be considered valid within a habitat.
#' @return Updated input data frame with occurrences flagged as `"delete"` in
#'   the `action` column when they fail the habitat criterion.
#' @examples
#' \dontrun{
#' flag_by_habitat(occurrence_df, habitat_proportion = 0.7)
#' }
#' @export
#' 
flag_by_habitat <- function(occurrence_df, habitat_proportion=0.5){
  
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

#' Detect false negatives in mock
#' 
#' Build a data frame of expected occurrences that are absent from the data
#' (false negatives) based on a reference mock composition.
#' 
#' @param read_count_samples Data frame or CSV file with columns
#'   `asv`, `sample`, `read_count`.
#' @param mock_composition Data frame or CSV file with columns
#'   `sample`, `action` (keep/tolerate), and `asv`.
#' @param sep Field separator character used in input and output CSV files.
#' @param out Character string naming the output file. If NULL, no file is written.
#' @param quiet logical; if TRUE, suppress informational messages and show only
#'   warnings or errors.
#' @return Data frame with columns `sample`, `action`, `asv`, `asv_id`.
#' @examples
#' \dontrun{
#' detect_false_negatives(
#'   read_count_samples = read_count_samples_df,
#'   mock_composition = "data/mock_composition.csv"
#' )
#' }
#' @export
#'
detect_false_negatives <- function(read_count_samples, mock_composition, sep=",", out=NULL, quiet=TRUE){
  
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
  
  # make asv_id asv pairs for identifing the asv_id of the missing occurrences
  asvs <- read_count_samples_df %>%
    select(asv_id, asv) %>%
    distinct()
  
  # add read_count to df from read_count_samples_df, and keep only if value is NA
  df <- left_join(mock_comp, read_count_samples_df,  by=c("sample", "asv")) %>%
    filter(is.na(read_count)) %>%
    select(-read_count, -asv_id) %>% # delete asv_id, since it is NA
    left_join(asvs, by=c("asv")) # add asv_id if exists
  
  if (!quiet & nrow(df) > 0) {
    missing_asv_text <- paste(capture.output(print(df)), collapse = "\n")
    warning(
      paste0(
        "\n  Some expected ASVs are missing from the mock samples.\n",
        "----------------------------------------------------------\n",
        missing_asv_text,
        "\n----------------------------------------------------------\n"
      ),
      call. = FALSE
    )
  }

  
  
  # write to outfile
  if(!is.null(out)){
    check_dir(out, is_file=TRUE)
    write.table(df, file=out, row.names = F, sep=sep)
  }
  
#  FN <- nrow(df %>%
#               filter(action=="keep"))
  return(df)
}

#' Suggest cutoff values for `filter_occurrence_variant` and `filter_occurrence_read_count`
#' 
#' Explore combinations of parameters for `filter_occurrence_read_count` and
#' `filter_occurrence_variant` (followed by `filter_min_replicate`) in order
#' to identify optimal filtering thresholds.
#' 
#' For each parameter combination, the number of false negatives (FN),
#' true positives (TP), and false positives (FP) is computed. Users are
#' encouraged to select parameter settings that minimize both FN and FP.
#' 
#' If `known_occurrences` is not provided, it is constructed from
#' `read_count`, `mock_composition`, `sampleinfo`, and `habitat_proportion`.
#' 
#' @param read_count Data frame or CSV file with columns
#'   `asv_id`, `sample`, `replicate`, `read_count`, `asv`.
#' @param outdir Character string naming the output directory.
#' @param known_occurrences Data frame or file produced by
#'   `classify_control_occurrences()`, containing known true positives and
#'   false positives. Optional; if NULL, it is computed from
#'   `mock_composition` and `sampleinfo`.
#' @param mock_composition Data frame or CSV file with columns
#'   `sample`, `action` (keep/tolerate), and `asv`.
#' @param sampleinfo Data frame or CSV file with columns
#'   `sample`, `sample_type` (mock/negative/real), `habitat`, `replicate`
#'   (optional), and optional file metadata.
#' @param habitat_proportion Numeric value between 0 and 1. Used to infer
#'   habitat-based artifacts when constructing known occurrences.
#' @param sep Field separator character used in input and output CSV files.
#' @param min_read_count_cutoff Positive integer specifying the minimum value
#'   tested for `filter_occurrence_read_count()`.
#' @param max_read_count_cutoff Positive integer specifying the maximum value
#'   tested for `filter_occurrence_read_count()`.
#' @param increment_read_count_cutoff Positive integer defining the step size
#'   between tested read count cutoff values.
#' @param min_variant_cutoff Numeric value between 0 and 1 specifying the
#'   minimum value tested for `filter_occurrence_variant()`.
#' @param max_variant_cutoff Numeric value between 0 and 1 specifying the
#'   maximum value tested for `filter_occurrence_variant()`.
#' @param increment_variant_cutoff Numeric value between 0 and 1 defining the
#'   step size between tested variant cutoff values.
#' @param by_replicate logical; argument passed to filter_occurrence_variant
#' @param `min_replicate_number()` Positive integer specifying the minimum number of
#'   replicates for `filter_min_replicate()`.
#' @param quiet logical; if TRUE, suppress informational messages and show only
#'   warnings or errors.
#' @return Data frame with columns `read_count_cutoff`, `variant_cutoff`,
#'   `FN`, `TP`, and `FP`.
#' @examples
#' \dontrun{
#' suggest_variant_readcount_cutoffs(
#'   read_count_df,
#'   known_occurrences = "data/known_occurrences.csv",
#'   min_read_count_cutoff = 10,
#'   max_read_count_cutoff = 50,
#'   increment_read_count_cutoff = 10,
#'   min_variant_cutoff = 0.001,
#'   max_variant_cutoff = 0.005,
#'   increment_variant_cutoff = 0.001
#' )
#' }
#' @export
#'
suggest_variant_readcount_cutoffs <- function(read_count, 
                                           outdir, 
                                           known_occurrences = NULL, 
                                           mock_composition = NULL,
                                           sampleinfo = NULL,
                                           habitat_proportion = 0.5,
                                           sep=",",
                                           min_read_count_cutoff=10, 
                                           max_read_count_cutoff=100, 
                                           increment_read_count_cutoff=5, 
                                           min_variant_cutoff=0.001, 
                                           max_variant_cutoff=0.01, 
                                           increment_variant_cutoff=0.001, 
                                           by_replicate=FALSE, 
                                           min_replicate_number=1, 
                                           quiet=T
){

  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  outdir <- path.expand(outdir)
  # can accept df or file as an input
  if(is.null(known_occurrences)){ # Calculate known_occurrences from data
    
    known_occurrences <- file.path(outdir, "known_occurrences.csv")
    false_negatives <- file.path(outdir, "false_negatives.csv")
    performance_metrics <- file.path(outdir, "performance_metrics.csv")
    
    results <- classify_control_occurrences(
      read_count_df,
      sampleinfo = sampleinfo,
      mock_composition = mock_composition,
      known_occurrences = known_occurrences,
      false_negatives = false_negatives,
      performance_metrics = performance_metrics,
      sep = ",",
      habitat_proportion = habitat_proportion,
      quiet=TRUE
    )
    
    known_occurrences_df <- results[[1]]
    false_negatives_df <- results[[2]]
    if(nrow(false_negatives_df)>0){
      warning(
        paste0(
          "\n  Some expected ASVs are missing from the mock samples.\n",
          "  Check the ", false_negatives, " file! \n"
        ),
        call. = FALSE
      )
    }
  }else{# known_occurrences is given
    if(is.character(known_occurrences)){
      # read known occurrences
      known_occurrences_df <- read.csv(known_occurrences, header=T, sep=sep)
    }else{
      known_occurrences_df <- known_occurrences
    }
    check_file_info(file=known_occurrences_df, 
                  file_type="known_occurrences", 
                  sep=sep, 
                  quiet=TRUE
    )
  }
  
  # make a series of cutoff values for filter_occurrence_read_count
  rc_cutoff_list <- seq(from=min_read_count_cutoff, 
                        to=max_read_count_cutoff, 
                        by=increment_read_count_cutoff
  )
  # make a series of cutoff values for filter_occurrence_read_count
  var_cutoff_list <- seq(from=min_variant_cutoff, 
                         to=max_variant_cutoff, 
                         by=increment_variant_cutoff
  )
  
  out_df <- data.frame(
    read_count_cutoff=numeric(),
    variant_cutoff=numeric(),
    FN=numeric(),
    TP=numeric(),
    FP=numeric()
  )
  # go through all parameter combination and count the number of TP and FN
  
  for(rc_cutoff in rc_cutoff_list){
    df_tmp <- read_count_df
    #filter_occurrence_read_count
    df_tmp <- filter_occurrence_read_count(df_tmp, rc_cutoff)
    for(var_cutoff in var_cutoff_list){
      # filter_occurrence_variant
      df_tmp <- filter_occurrence_variant(df_tmp, var_cutoff, by_replicate=by_replicate)
      # filter_min_replicate
      df_tmp <- filter_min_replicate(df_tmp, min_replicate_number)
      # pool_replicates
      df_tmp_sample <- pool_replicates(df_tmp, method="max",digits=0) # the method does really not matter here
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
      new_line <- data.frame(read_count_cutoff=rc_cutoff, 
                             variant_cutoff=var_cutoff ,
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
    arrange(FN, FP, variant_cutoff, read_count_cutoff)
  
  outfile <- file.path(outdir, "suggest_variant_readcount_cutoffs.csv")
  check_dir(outfile, is_file=TRUE)
  write.table(out_df, file=outfile, sep=sep, row.names = F)
  
  return(out_df)
}

#' Pool multiple datasets
#' 
#' Combine multiple datasets generated from the same genetic marker (i.e.,
#' the same genomic region and primer pair).
#' 
#' When identical sample–replicate combinations are present across datasets,
#' their read counts are aggregated using the selected method: sum, mean,
#' minimum, or maximum.
#' 
#' Input files must be in long format and include the columns:
#' `asv_id`, `sample`, `read_count`, and `asv`, with an optional `replicate`
#' column.
#' 
#' Consistency between `asv_id` and `asv` across datasets is checked prior
#' to pooling.
#' 
#' @param files Character vector of file paths to the datasets to pool.
#' Each file must follow the same format, containing `asv_id`, `sample`,
#' `read_count`, `asv`, and optionally `replicate`.
#' @param outfile Character string specifying the output CSV file name.
#' If NULL, no file is written.
#' @param method Character string specifying how read counts from identical
#' sample–replicates are aggregated. Must be one of `"mean"`, `"max"`,
#' `"sum"`, or `"min"`.
#' @param sep Character string specifying the field separator used in input
#' and output CSV files.
#' @param quiet Logical; if TRUE, suppress informational messages and show only
#' warnings or errors.
#' @return A data frame with columns `asv_id`, `sample`, `replicate`
#' (optional), `read_count`, and `asv`.
#' @examples
#' \dontrun{
#' files <- c("vtamR_test/run1/1_Input.csv",
#'            "vtamR_test/run2/1_Input.csv")
#' df <- pool_datasets(files, method = "sum")
#' }
#' @export
#' 
pool_datasets <- function(files, 
                         outfile=NULL, 
                         method="mean",
                         sep=",", 
                         quiet=T
                         ){
  
  # method
  method <- match.arg(method, c("mean", "max", "sum", "min"))
  fun <- switch(method,
                mean = function(x) mean(x, na.rm = TRUE),
                max  = function(x) max(x, na.rm = TRUE),
                sum  = function(x) sum(x, na.rm = TRUE),
                min  = function(x) min(x, na.rm = TRUE))
  
  
  # read the first file
  df_pool <- read.table(files[1], sep=sep, header=TRUE)
  cols_sorted <- sort(colnames(df_pool))
  cols <- colnames(df_pool)

  for(i in 2:length(files)){
    
      df <- read.table(files[i], sep=sep, header=TRUE)
      
      # test if same columns
      cols_sorted_tmp <-  sort(colnames(df))
      if(!identical(cols_sorted_tmp, cols_sorted)){
        stop("All input files must have the same columns")
      }
      
      # use the same column order as in the first file
      df <- df %>%
        select(!!cols)
      # concatenation 
      df_pool <- rbind(df_pool, df)
      # check if coherence among asv and asv_id
      if(!check_one_to_one(df_pool)){
        stop("Incoherence between asv and asv_id among different datasets")
      }
    }
  
  if("replicate" %in% cols){
    df_pool <- df_pool %>%
      group_by(asv_id, sample, replicate, asv) %>%
      summarise(read_count = fun(read_count), .groups="drop")
  }else{
    df_pool <- df_pool %>%
      group_by(asv_id, sample, asv) %>%
      summarise(read_count = fun(read_count), .groups="drop") %>%
      mutate(read_count = round(read_count, digits=0))
  }
 
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(df_pool, file=outfile, sep=sep, row.names = F)
  }

  return(df_pool)
}

#' Pool markers
#' 
#' Combine multiple datasets from the same samples generated using different
#' markers targeting overlapping genomic regions.
#' 
#' Input files must be in long format and contain the columns:
#' `asv_id`, `sample`, `read_count`, `asv`, and an optional `replicate`.
#' 
#' ASVs sharing identical sequences across overlapping regions are grouped
#' together. Within each group, a centroid sequence is defined as the ASV
#' with the highest total read count. Read counts are then aggregated across
#' ASVs within each group for each sample–replicate combination.
#' 
#' To avoid ambiguity between identical `asv_id`s originating from different
#' markers, optional marker identifiers can be appended to `asv_id`.
#' 
#' @param files Character vector of file paths to the datasets to pool.
#' Each file must follow the same format, containing `asv_id`, `sample`,
#' `read_count`, `asv`, and optionally `replicate`.
#' @param marker_ids Optional character or integer vector specifying marker IDs.
#' Must be positive non-zero integers and provided in the same order as
#' `files`. These IDs are appended to `asv_id` to ensure uniqueness across
#' markers. If not provided, `asv_id`s must already be unique across datasets.
#' @param outfile Character string specifying the output CSV file name.
#' If NULL, no file is written.
#' @param asv_with_centroids Character string specifying an optional output CSV
#' file containing the merged dataset annotated with `centroid_id` and
#' `centroid` columns.
#' @param method Character string specifying how read counts are aggregated
#' within groups. Must be one of `"mean"`, `"max"`, `"sum"`, or `"min"`.
#' @param vsearch_path Character string specifying the path to the `vsearch`
#' executable.
#' @param num_threads Positive integer specifying the number of CPUs to use.
#' If 0, all available CPUs are used.
#' @param sep Character string specifying the field separator used in input
#' and output CSV files.
#' @param quiet Logical; if TRUE, suppress informational messages and show only
#' warnings or errors.
#' 
#' @return A data frame with columns `asv_id`, `sample`, `read_count`,
#' `asv`, and optional `replicate`, where ASVs belonging to the same group
#' within each sample–replicate are merged into a single row. The `read_count`
#' reflects the selected aggregation method (mean, max, min, or sum).
#' 
#' @examples
#' \dontrun{
#' files <- c(
#'   "~/vtamR_demo_out_zfzr/filter/7_filter_chimera.csv",
#'   "~/vtamR_demo/filter/7_filter_chimera.csv"
#' )
#' marker_ids <- c(1, 2)
#' pool_markers(files, marker_ids = marker_ids, method = "mean")
#' }
#' @export
#' 
pool_markers <- function(files, 
                         marker_ids = NULL,
                         outfile=NULL, 
                         asv_with_centroids=NULL,
                         method="mean", 
                         vsearch_path="vsearch", 
                         num_threads=0,
                         sep=",", 
                         quiet=T){
  
  # method
  method <- match.arg(method, c("mean", "max", "sum", "min"))
  fun <- switch(method,
                mean = function(x) mean(x, na.rm = TRUE),
                max  = function(x) max(x, na.rm = TRUE),
                sum  = function(x) sum(x, na.rm = TRUE),
                min  = function(x) min(x, na.rm = TRUE))
  
  
  #########################
  # concatenate all data in csv
  # read the first file
  df_pool <- read.table(files[1], sep=sep, header=TRUE)
  cols_sorted <- sort(colnames(df_pool))
  cols <- colnames(df_pool)
  if(!is.null(marker_ids)){
    df_pool <- df_pool %>%
      mutate(asv_id = paste(asv_id, marker_ids[1], sep=".")) %>%
      mutate(asv_id = as.numeric(asv_id))
    
  }

  for(i in 2:length(files)){
    
    df <- read.table(files[i], sep=sep, header=TRUE)
    # test if same columns
    cols_sorted_tmp <-  sort(colnames(df))
    if(!identical(cols_sorted_tmp, cols_sorted)){
      stop("All input files must have the same columns")
    }
    # use the same column order as in the first file
    df <- df %>%
      select(!!cols)
    # prefix asv_id
    if(!is.null(marker_ids)){
      df <- df %>%
        mutate(asv_id = paste(asv_id, marker_ids[i], sep=".")) %>%
        mutate(asv_id = as.numeric(asv_id))
    }
    # concatenation 
    df_pool <- rbind(df_pool, df)
    # check if coherence among asv and asv_id
    if(!check_one_to_one(df_pool)){
      stop("Incoherence between asv and asv_id among different datasets")
    }
  }
  ############################
  # Make asv_with_centroids df: concatenated input + centroid_id + centroid
  df_pool <- cluster_asv(read_count=df_pool,
                     group = FALSE,
                     by_sample = FALSE,
                     method = "vsearch",
                     path = vsearch_path,
                     identity = 1,
                     quiet = TRUE
                     )
  
  asv_with_centroids_df <- df_pool %>%
    rename(centroid_id = cluster_id) %>%
    # add centroid seq
    group_by(centroid_id) %>%
    mutate(centroid = first(asv[centroid_id==asv_id])) %>%
    ungroup()
  
  ############################
  ### group lines of the same cluster
  if("replicate" %in% cols){
    df_pool <- asv_with_centroids_df %>%
      select(-asv_id, -asv) %>%
      rename(asv_id = centroid_id, asv=centroid) %>%
      group_by(asv_id, sample, replicate, asv) %>%
      summarize(read_count = fun(read_count), .groups = "drop") %>%
      mutate(read_count = round(read_count, digits=0)) %>%
      select(asv_id, sample, replicate, read_count, asv)
  }else{
    df_pool <- asv_with_centroids_df %>%
      select(-asv_id, -asv) %>%
      rename(asv_id = centroid_id, asv=centroid) %>%
      group_by(asv_id, sample, asv) %>%
      summarize(read_count = fun(read_count), .groups = "drop") %>%
      mutate(read_count = round(read_count, digits=0)) %>%
      select(asv_id, sample, read_count, asv)
  }
    
  if(!is.null(asv_with_centroids)){
    check_dir(asv_with_centroids, is_file=TRUE)
    write.table(asv_with_centroids_df, file=asv_with_centroids, sep=sep, row.names = F)
  }
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(df_pool, file=outfile, sep=sep, row.names = F)
  }
  
  return(df_pool)
}


#' Filter processing history by feature
#' 
#' Filter a feature (`asv_id`, `asv`, `sample`, `replicate`, or `read_count`)
#' across all intermediate filtering output files to retain only rows matching
#' a specified value.
#' 
#' This function scans all output files from intermediate filtering steps in a
#' specified directory and extracts rows where the selected feature matches
#' any of the values provided in a vector.
#' 
#' By default, input filenames must start by a number (e.g. `5_filter_occurrence_sample.csv`).
#' See `pattern` to change this behavior.
#' 
#' @param dir Character string specifying the directory containing intermediate
#' filtering output files.
#' @param pattern A regular expression used to select filenames to be scanned.
#' @param feature Character string specifying the feature to filter by.
#' Must be one of `"asv_id"`, `"asv"`, `"sample"`, `"replicate"`, or
#' `"read_count"`.
#' @param values Numerical or Character vector with values to match in the selected feature.
#' Only rows containing this value are retained.
#' @param sep Field separator character used in input and output CSV files.
#' @return An invisible data frame containing all rows matching the selected
#' feature value across all intermediate files.
#' @examples
#' \dontrun{
#' history_by(dir = "out", feature = "asv_id", value = 1)
#' history_by(dir = "out", feature = "sample", value = "tpos1")
#' }
#' @export
#
history_by <- function(dir, pattern="^\\d", feature, values, sep=","){

  dir = check_dir(dir)
  files <- list.files(path=dir, pattern=pattern, full.names=FALSE)
  
  # get filenames to df and arrange the according to the number at the beginning of the filename
  df <- data.frame("files"= files) %>%
    arrange(files)
  
  df$order <- gsub("[^0-9].*$", "", df$files)
  if(!any(df$order == "")){
    df$order <- as.numeric(df$order)
    df <- df %>%
      arrange(order, files)
  }
  
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
#        filter(!!sym(feature)==value) %>% # filter using a symbol from feature
        filter(!!sym(feature) %in% values) %>% # filter using a symbol from feature
        select(file, asv_id, sample, replicate, read_count, asv)
      
      selected_lines <- rbind(selected_lines, tmp)
    }else{
      stop("ERROR: feature is not in file")
    }
  }
  return(invisible(selected_lines))
}

#' Summarize intermediate filtering steps by feature
#' 
#' Summarize the outputs of intermediate filtering steps across all files in a
#' directory.
#' By default, input filenames must start by a number (e.g. `5_filter_occurrence_sample.csv`).
#' See `pattern` to change this behavior.
#' 
#' For each file, rows are grouped by `grouped_by`, and either:
#' - the number of distinct values of `feature` is computed, or
#' - if `feature = "read_count"`, the total read counts are summed per group.
#' 
#' @param dir Character string specifying the directory containing intermediate
#' filtering output files. Files must start with a numeric prefix followed by
#' an underscore (e.g. `5_filter_occurrence_sample.csv`).
#' @param pattern A regular expression used to select filenames to be scanned.
#' @param feature Character string specifying the feature to summarize. Must be
#' one of `"asv_id"`, `"asv"`, `"sample"`, `"replicate"`, or `"read_count"`.
#' @param grouped_by Character string specifying the grouping variable. Must be
#' one of `"asv_id"`, `"asv"`, `"sample"`, or `"replicate"`.
#' @param sep Field separator character used in input and output CSV files.
#' @param outfile Character string specifying the output CSV file name.
#' If NULL, no file is written.
#' @return An invisible data frame where columns correspond to input files,
#' rows correspond to `grouped_by` values, and cell values represent either
#' counts of `feature` occurrences or summed read counts.
#' @examples
#' \dontrun{
#' summarize_by(dir = "vtamR_test/out_mfzr", feature = "asv", grouped_by = "sample")
#' summarize_by(dir = "vtamR_test/out_mfzr", feature = "read_count", grouped_by = "sample")
#' }
#' @export
#'
summarize_by <- function(dir, pattern = "^\\d", feature, grouped_by, outfile=NULL, sep=","){
  
  # read file names in dir
  dir = check_dir(dir)
  files <- list.files(path=dir, pattern=pattern, full.names=FALSE)
  
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
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(wide_df, file=outfile, row.names = F, sep=sep)
  }
  return(invisible(wide_df))
}

#' Write data frame to FASTA file
#' 
#' Export a data frame to a FASTA file. Output can be written as plain text
#' or gzip-compressed. The file extension is automatically adjusted based on
#' the compression setting if necessary.
#' 
#' @param df Input data frame with columns `header` and `sequence`.
#' @param out Character string specifying the output file name.
#' @param compress logical; if TRUE, the output file is compressed using gzip.
#' @return Character string giving the final output file name, updated to
#' reflect compression if applicable.
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   header = c("seq1", "seq2"),
#'   sequence = c(
#'     "AACTTGTTGTCACTGTAAACTGATGTA",
#'     "AACTTGTTGTCACTGTTTGACTGATGTA"
#'   )
#' )
#' write_fasta_from_df(df, out = "out/test.fasta", compress = TRUE)
#' }
#' @export
#'
write_fasta_from_df <- function(df, out, compress=F){
  
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

#' Read FASTA file to data frame
#' 
#' Import a FASTA file into a data frame. The file can be plain text or
#' gzip-compressed (ZIP files are not supported).
#' 
#' @param file Character string specifying the input FASTA file path.
#' @param dereplicate logical; if TRUE, returns a data frame of unique ASVs
#' with read counts. If FALSE, returns a data frame of individual sequences.
#' @return A data frame with:
#'   - `header` and `sequence` columns if `dereplicate = FALSE`
#'   - `asv` and `read_count` columns if `dereplicate = TRUE`
#' @examples
#' \dontrun{
#' read_fasta_to_df(file = "data/test.fasta", dereplicate = FALSE)
#' read_fasta_to_df(file = "data/test.fasta", dereplicate = TRUE)
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
  
  # read file to a vector. Each element is a line
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

#' Count reads in sequence or text files
#' 
#' Count the number of sequences in FASTA or FASTQ files, or the number of lines 
#' in other file types.
#' 
#' This function is optimized for Linux-like systems and may be slower on Windows 
#' for large files.
#' It supports both gz-compressed and uncompressed files (ZIP archives are not supported).
#'  
#' @param file Character string: path to the input file.
#' @param file_type Character string specifying the file type: `"fasta"` or `"fastq"`.
#'   For any other value, the function returns the number of lines in the file.
#' @return Integer: number of sequences for FASTA/FASTQ files, or number of 
#' lines for other file types.
#' @examples
#' \dontrun{
#' count_reads(file = "data/test.fasta", file_type = "fasta")
#' }
#' @export
#' 
count_reads <- function(file, file_type="fastq"){
  
  if (endsWith(file, ".zip")) {
    stop("File compression type is not supported.")
  }
  
  if(is_linux()){
    # compressed files
    if(endsWith(file, ".gz") || endsWith(file, ".bz") || endsWith(file, ".gz2")){
      if(file_type == "fastq"){
        cmd <- paste("zcat ", file, "| wc -l ", sep=" ")
        seq_count <- scan(text = system(cmd, intern = TRUE), what = integer(), nmax = 1, quiet = TRUE)
#        seq_count <- as.integer(system(cmd, intern=TRUE))
        seq_count <- seq_count/4
      }else if(file_type == "fasta"){
        cmd <- paste("zcat ", file, "| grep '^>' -P | wc -l", sep=" ")
        seq_count <- scan(text = system(cmd, intern = TRUE), what = integer(), nmax = 1, quiet = TRUE)
      }else{
        msg <- paste(file_type, "is neither fasta nor fastq. 
                     The number of liens in file will be returned for", file)
        print(msg)
        cmd <- paste("zcat ", file, "| wc -l ", sep=" ")
        seq_count <- scan(text = system(cmd, intern = TRUE), what = integer(), nmax = 1, quiet = TRUE)
      }
    }else{
      #uncompressed files
      if(file_type == "fastq"){
        cmd <- paste("wc", file, "-l", sep=" ")
        seq_count <- scan(text = system(cmd, intern = TRUE), what = integer(), nmax = 1, quiet = TRUE)
        seq_count <- seq_count/4
      }else if(file_type == "fasta"){
        cmd <- paste("grep '^>' -P", file, "| wc -l", sep=" ")
        seq_count <- scan(text = system(cmd, intern = TRUE), what = integer(), nmax = 1, quiet = TRUE)
      }else{
        msg <- paste(file_type, "is neither fasta nor fastq. 
                     The number of lines in file will be returned for", file)
        print(msg)
        cmd <- paste("wc", file, "-l", sep=" ")
        seq_count <- scan(text = system(cmd, intern = TRUE), what = integer(), nmax = 1, quiet = TRUE)
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
}

#' Count reads in files within a directory
#' 
#' Count the number of sequences in FASTA or FASTQ files, or the number of lines 
#' in other file types,
#' for all files in a directory matching a given pattern.
#' 
#' This function is optimized for Linux-like systems and may be slower on Windows 
#' for large files.
#' It supports both gz-compressed and uncompressed files (ZIP archives are not supported).
#'  
#' @param dir Character string: path to the input directory.
#' @param pattern Regular expression: pattern used to select files in the directory. Only files whose names match the pattern are processed.
#' @param file_type Character string specifying file type: `"fasta"` or `"fastq"`.
#'   For any other value, the function returns the number of lines in each file.
#' @param sep Field separator character in input and output CSV files.
#' @param outfile Character string: output CSV file name. If NULL, no file is written.
#' @param quiet Logical: if TRUE, suppress informational messages and show only warnings or errors.
#' @return Data frame with two columns: `filename`, `read_count`.
#' @examples
#' \dontrun{
#' count_reads_in_dir(dir = "out", pattern = "\\.fastq", file_type = "fastq")
#' count_reads_in_dir(dir = "out", pattern = "^mfzr", file_type = "fasta")
#' }
#' @export
#' 
count_reads_in_dir<- function(dir, 
                         pattern=".", 
                         file_type="fasta", 
                         outfile=NULL, 
                         sep=",", 
                         quiet=T
                         ){
  
  dir = check_dir(dir, is_file=FALSE)
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
    n <- count_reads(file_p, file_type=file_type)
    df[i, "read_count"] <- n
  }
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(df, file=outfile, sep=sep, row.names = F)
  }
  return(df)
}

#' Validate input file formats and consistency
#' 
#' Performs a series of checks on input files or data frames to ensure format validity
#' and internal consistency across different pipeline inputs.
#' 
#' The function performs the following checks depending on `file_type`:
#' * Presence of all required columns for each file type
#' * Validity of sample names (alphanumeric format)
#' * Validity of sequence-related fields (e.g. `asv`, `tag_fw`, `tag_rv`, `primer_fw`, `primer_rv`)
#'   restricted to IUPAC nucleotide codes
#' * Numeric format of `read_count` values
#' * Consistency of sample type and habitat across replicates
#'   (`fastqinfo`, `fastainfo`, `sampleinfo`)
#' * One-to-one pairing of FASTQ files (e.g. forward/reverse reads; `fastqinfo`)
#' * Existence of referenced files in `fastq_fw`, `fastq_rv`, and `fasta` columns
#'   (`fastqinfo`, `fastainfo`, `sampleinfo`)
#' * Uniqueness of tag combinations within file pairs (`fastqinfo`, `fastainfo`)
#' * Validity of action values (`mock_composition`, `known_occurrences`)
#' * One-to-one mapping between `asv_id` and `asv`
#'   (`read_count`, `read_count_sample`, `asv_list`)
#' 
#' @param file Character string or data frame: input file path or already-loaded data frame.
#' @param dir Character string: directory containing files referenced in `fastq_fw`, `fastq_rv`,
#'   or `fasta` columns.
#' @param file_type Character string specifying the type of input file. Must be one of:
#'   `"fastqinfo"`, `"fastainfo"`, `"sampleinfo"`, `"mock_composition"`,
#'   `"known_occurrences"`, `"read_count"`, `"read_count_sample"`, `"asv_list"`.
#' @param sep Field separator character used in CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and show only warnings or errors.
#' @return Return an error message and stops execution if inconsistencies are detected.
#' @examples
#' \dontrun{
#' check_file_info(file = "input/sampleinfo.csv", dir = "fasta", file_type = "sampleinfo")
#' check_file_info(file = sampleinfo_df, dir = "fasta", file_type = "read_count_sample")
#' }
#' @export
#' 
check_file_info <- function(file, dir, file_type="fastqinfo", sep=",", quiet=FALSE){
  
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
  }else if(file_type == "sampleinfo"){
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
  if(file_type == "fastqinfo" || file_type == "fastainfo" || file_type == "sampleinfo" ){
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
      msg <- paste("Sample-replicate combinations should be unique:", 
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
  if(file_type == "fastainfo" || file_type == "sampleinfo"){
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
  if(file_type == "fastqinfo" || file_type == "fastainfo" || file_type == "sampleinfo"){
    
    if(file_type == "fastqinfo"){
      file_list_fw <- unique(df$fastq_fw)
      file_list_rv <- unique(df$fastq_rv)
      file_list <- c(file_list_fw, file_list_rv)
    }
    if(file_type == "fastainfo" || file_type == "sampleinfo"){
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
#' Verifies that all files listed in `file_list` exist in the specified directory.
#' 
#' @param dir Character string: path to the directory containing the input files.
#' @param file_list Character vector of file names to check.
#' @return Throws an error and stops execution if one or more files are missing.
#' @examples
#' \dontrun{
#' file_list <- c("14ben01-1.fasta", "14ben01-2.fasta")
#' check_file_exists(file_list = file_list, dir = "vtamR_test/out_mfzr/sorted")
#' }
#' @export
#' 
check_file_exists <- function(file_list, dir){
  
  dir = check_dir(dir)
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


#' Check file column headings
#' 
#' Compares the column names (header) of a file or vector against an expected set of names.
#' 
#' @param list1 Character vector: actual column names found in the file.
#' @param list2 Character vector: expected column names to compare against.
#' @param file Character string: name of the file being checked (used for error reporting).
#' @return Throws an error listing elements in `list1` that are missing from `list2`.
#' @examples
#' \dontrun{
#' list1 <- c("tag_fw","primer_fw","tag_rv","primer_rv","sample",
#'            "sample_type","habitat","replicate","fasta")
#' list2 <- c("tag_fw","primer_fw","tag_rv","primer_rv","sample",
#'            "sample_type","habitat","replicate","fastq_fw","fastq_rv")
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

#' Read VSEARCH cluster size output (outmft6 format)
#' 
#' Reads the output of `vsearch --cluster_size` into a data frame.
#' 
#' @param filename Character string: path to the tab-separated output file 
#' produced by `vsearch --cluster_size`.
#'   The file must contain merged and centroid ASV identifiers in the first two columns.
#' @return Data frame with columns `merged_id` and `centroid_id`.
#' @export
#' 
read_vsearch_cluster_size_outmft6 <- function(filename) {
  
  df <- read.table(filename, header=FALSE, sep="\t")
  df <- df[,1:2]
  colnames(df) <- c("merged_id","centroid_id")
  df$merged_id <- gsub(';size=[0-9]+', '', df$merged_id)
  df$centroid_id <- gsub(';size=[0-9]+', '', df$centroid_id)
  
  df$merged_id <- as.numeric(df$merged_id)
  df$centroid_id <- as.numeric(df$centroid_id)
  
  return(df)
}

#' Write FASTA file with optional read counts
#' 
#' Writes a FASTA file from a data frame, optionally including read counts in the sequence headers.
#' 
#' If `read_count = TRUE`, definition lines are written in the format `>label;size=###`.
#' If `read_count = FALSE`, definition lines are written in the format `>label`.
#' 
#' @param df Data frame containing at least `asv_id` and `asv` columns, and optionally `read_count`.
#' @param outfile Character string: output FASTA file name.
#' @param read_count Logical: if TRUE, include read counts in FASTA headers using `;size=` format.
#' @return NULL; writes a FASTA file to disk.
#' @export
#' 
write_fasta_with_counts <- function(df, outfile, read_count=FALSE) {
  # Open the file for writing
  
  if(read_count){
    df <- df %>%
      group_by(asv_id, asv) %>%
      summarize(read_count = sum(read_count)) 

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

#' Match ASV variants to expected mock species
#'
#' Identifies expected variants in mock samples by matching ASVs to a
#' custom reference database.
#'
#' This function performs the following steps:
#'
#' - Builds a small BLAST database from sequences of species expected in the mock samples
#'   (or closely related taxa). Reference sequences should at least 
#'   partially cover the target amplicon
#'   (approximately 70% coverage) and may be shorter or longer than the ASVs.
#' - Assigns taxonomy to all ASVs detected in mock samples using this custom database.
#' - Selects the most abundant ASV for each taxon.
#' - Generates a `mock_composition` template file for user validation and manual curation.
#'
#' The output directory contains two files:
#'
#' - *mock_taxassign.csv*:
#'   Taxonomic assignment results, including total read counts per ASV in mock samples
#'   and the number of mock samples in which each ASV is detected.
#'
#' - *mock_composition_template_to_check.csv*:
#'   A template for constructing the `mock_composition` file.
#'   It includes the most abundant sequence for each `ltg_name` identified in the taxonomic assignment,
#'   repeated across mock samples.
#'
#'   When multiple mock samples have different expected compositions, users should remove lines
#'   corresponding to taxa not expected in each sample.
#'
#'   Note: this file does not include reference species that did not match any ASV.
#'   This may occur if the species was not amplified in the mock samples or if the reference
#'   sequence in the custom database is incorrect.
#'
#' @param read_count Data frame or CSV file containing columns `asv_id`, `sample`, `asv`, and `read_count`.
#' @param fas Character string: path to a FASTA file containing reference sequences from expected mock species
#'   or closely related taxa. FASTA headers must follow the format:
#'   `>HQ563207.1 taxID=1592914`, where `taxID` is a valid NCBI taxonomic identifier
#'   (https://www.ncbi.nlm.nih.gov/taxonomy).
#' @param sampleinfo Data frame or CSV file containing at least `sample` and `sample_type` columns.
#' @param taxonomy Character string: path to a TSV taxonomy file containing columns:
#'   `tax_id`, `parent_tax_id`, `rank`, `name_txt`, `old_tax_id` (merged into `tax_id`),
#'   and `taxlevel` (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum,
#'   2: kingdom, 1: domain, 0: root).
#'   A COInr taxonomy file can be used.
#' @param outdir Character string: output directory for generated files.
#' @param blast_path Character string: path to BLAST executable.
#' @param num_threads Positive integer: number of CPU threads to use. If 0, all available CPUs are used.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and show only warnings or errors.
#' @return Data frame with columns: `sample`, `action`, `asv`, `taxon`, `asv_id`.
#' @examples
#' \dontrun{
#' match_variants_to_mock_species(
#'   read_count = read_count_df,
#'   fas = "reference.fasta",
#'   taxonomy = "taxonomy.tsv",
#'   blast_path = blast_path,
#'   sampleinfo = sampleinfo_df
#' )
#' }
#' @export
#'
match_variants_to_mock_species <- function(
    read_count,
    fas,
    sampleinfo,
    taxonomy,
    outdir,
    blast_path = "blastn",
    num_threads=0,
    sep=",",
    quiet=TRUE
){
  
  ##### Make blast db from mock fasta
  
  ## make TSV with seqID and taxID
  taxids <- file.path(tempdir(), "taxid.tsv" )
  write_taxid_mapping(file=fas, outfile=taxids)
  
  blastdb_path <- sub("blastn$", "makeblastdb", blast_path)
  write_output <- TRUE
  if(is.null(outdir)){
    outdir = file.path(tempdir(), "mock")
    write_output <- FALSE
  }
  bdmock <- file.path(outdir, "db_mock", "db")
  check_dir(bdmock, is_file = TRUE)
  
  args = c(
    "-dbtype", "nucl",
    "-in", fas,
    "-taxid_map", taxids,
    "-out", bdmock,
    "-parse_seqids"
  )
  run_system2(blastdb_path, args, quiet=quiet)
  
  if(is.character(sampleinfo)){
    # read known occurrences
    sampleinfo_df <- read.csv(sampleinfo, header=T, sep=sep)
  }else{
    sampleinfo_df <- sampleinfo
  }
  
  ### Get list of mocks
  mocks <- sampleinfo_df %>%
    filter(sample_type == "mock") %>%
    select(sample) %>%
    distinct()
  
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  ### Select sequences in mock samples
  # select mock samples
  mock_df <- read_count_df %>%
    filter(sample %in% mocks$sample)
  
  #### taxassign
    ltg_params_df = data.frame( pid=c(100,97,95,90,80),
                                pcov=c(70,70,70,70,70),
                                phit=c(0,0,0,0,0),
                                taxn=c(1,1,1,1,1),
                                seqn=c(1,1,1,1,1),
                                refres=c(1,1,1,1,1),
                                ltgres=c(8,8,8,8,8))

  taxa <- assign_taxonomy_ltg(asv= mock_df, 
                    taxonomy = taxonomy, 
                    blast_db = bdmock, 
                    blast_path=blast_path, 
                    ltg_params=ltg_params_df, 
                    quiet=quiet, 
                    fill_lineage=TRUE)
  
  
  
  ### filter and organize taxa
  # For each ASV, get the sum of the read counts in the mocks, and the number of mocks, where the ASV is present.
  df <- mock_df %>%
    group_by(asv_id) %>%
    summarize(total_read_count_mock = sum(read_count), number_mock=n_distinct(sample))
  # complete taxa with read and sample count info and arrange
  taxa <- left_join(taxa, df, by="asv_id") %>%
    select(asv_id, total_read_count_mock, number_mock, ltg_name, species, genus, family, order, class, phylum, pid, asv) %>%
    arrange(ltg_name, desc(total_read_count_mock))
  
  # select columns and sort the lines by species and decreasing read counts
  taxa_select <- taxa %>%
    filter(!is.na(ltg_name)) %>%
    group_by(ltg_name) %>%
    slice_head(n = 1) %>%
    select(asv, taxon=ltg_name, asv_id) %>%
    ungroup()
  
  mock_composition_template <- merge(mocks, taxa_select, by = NULL)%>%
    mutate(action="keep") %>%
    arrange(sample, taxon) %>%
    select(sample, action, asv,taxon, asv_id)
  
  #  mocks <- rbind(mocks, data.frame(sample=c("tpos2")))
  
  
  ### make mock composition template
  
  if(write_output){
    out_comp <- file.path(outdir, "mock_taxassign.csv")
    write.table(taxa, file=out_comp, row.names = FALSE, sep=sep)
    out <- file.path(outdir, "mock_composition_template_to_check.csv")
    write.table(mock_composition_template, file=out, row.names = FALSE, sep=sep)
  }
  
  return(mock_composition_template)
}

#' Create a taxID mapping file from a FASTA file
#'
#' Extracts sequence identifiers and NCBI taxonomic identifiers from FASTA headers
#' and writes them as a two-column tab-separated file.
#'
#' The FASTA headers must contain a taxonomic identifier in the format:
#' `>SequenceName taxID=12345`
#'
#' @param file Character string: path to a FASTA file containing sequences with
#'   taxonomic identifiers in the header (format `taxID=<number>`).
#' @param outfile Character string: path to the output file to be created.
#'   The file contains two tab-separated columns: `seq_id` and `tax_id`,
#'   without header or quotes.
#' @return NULL; writes a tab-delimited file to `outfile`.
#' @details
#' - Uses the internal helper `read_fasta_to_df()` to parse the FASTA file.
#' - Only the first two space-separated elements of each FASTA header are used;
#'   any additional information is ignored.
#' - The `taxID=` prefix is removed before converting values to numeric.
#' @examples
#' \dontrun{
#' write_taxid_mapping("input_sequences.fasta", "taxid_mapping.tsv")
#' }
#' @export
#' 
write_taxid_mapping <- function(file, outfile){
  
  df <- read_fasta_to_df(file)
  
  df <- df %>%
    separate(header, into=c("seq_id", "tax_id"), extra="drop", sep=" ", remove=TRUE) %>%
    mutate(tax_id = sub("taxID=", "", tax_id)) %>%
    mutate(tax_id = as.numeric(tax_id)) %>%
    select(seq_id, tax_id)
  
  write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}


#' Randomly sample sequences from a FASTA file
#'
#' Randomly selects `n` sequences from an input FASTA file.
#'
#' The input FASTA file may be uncompressed or gzipped. Output compression is
#' inferred automatically from the output file extension (i.e., `.gz` produces a
#' gzipped file).
#'
#' This implementation is platform-independent but may be slower than
#' system-optimized alternatives.
#'
#' If the input file contains fewer than or exactly `n` sequences, it is copied
#' to the output file, preserving the requested output compression format.
#'
#' @param fasta Character string: path to the input FASTA file (optionally gzipped).
#' @param outfile Character string: path to the output FASTA file. If it ends in `.gz`,
#'   the output will be gzip-compressed.
#' @param n Positive integer: number of sequences to randomly sample.
#' @param randseed Integer or NULL: random seed for reproducibility. If NULL (default),
#'   sampling is non-deterministic.
#' @param quiet Logical: if TRUE, suppress informational messages; only warnings and errors are shown.
#' @return Invisibly returns the number of sequences written to the output file.
#' @examples
#' \dontrun{
#' random_sample_r(
#'   fasta = "all_sequences.fasta.gz",
#'   outfile = "subset_100.fasta.gz",
#'   n = 100,
#'   randseed = 123,
#'   quiet = FALSE
#' )
#' }
#' @export
#' 
random_sample_r <- function(fasta, outfile, n=1000000, randseed = NULL, quiet=TRUE) {
  if (!is.null(randseed)) set.seed(randseed)
  
  # --- Count sequences
  con <- if (grepl("\\.gz$", fasta)) gzfile(fasta, "rt") else file(fasta, "rt")
  total <- 0L
  if(!quiet){cat("Counting sequences.\n")}
  while (length(lines <- readLines(con, n = 100000)) > 0) { # read by chunks
    total <- total + sum(startsWith(lines, ">"))
  }
  close(con)
  
  if(!quiet){ cat("Total sequences:", total, "\n")}
  
  # --- Check if n >= total
  if (n >= total) {
    txt <- sprintf("%s contains %d sequences.\nThe input file is copied to output.",
                   fasta, total)
    warning(txt, call. = FALSE)
    
#    txt <- paste("\n", fasta, "contains", total, "sequences.\n", "The input file is copied to output\n")
#    warning(txt)
    
    # determine compression for output
    if (grepl("\\.gz$", outfile)) {
      out_con <- gzfile(outfile, "wt")
    } else {
      out_con <- file(outfile, "wt")
    }
    
    # read input (compressed or not) and write to output
    con <- if (grepl("\\.gz$", fasta)){
      gzfile(fasta, "rt")
    }else {
      file(fasta, "rt")
    }
    while (length(lines <- readLines(con, n = 100000)) > 0) {
      writeLines(lines, out_con)
    }
    close(con)
    close(out_con)
    
    return(invisible(total))
  }
  
  # --- Sample indices
  keep_idx <- sort(sample.int(total, n))
  
  # --- Extract sampled sequences
  if(!quiet){cat("Extracting sampled sequences.\n")}
  con <- if (grepl("\\.gz$", fasta)) gzfile(fasta, "rt") else file(fasta, "rt")
  out_con <- if (grepl("\\.gz$", outfile)) gzfile(outfile, "wt") else file(outfile, "wt")
  
  seq_idx <- 0L # integer
  keep_pointer <- 1L # integer; position in the vector of keep_idx
  current_seq <- character()
  write_seq <- FALSE
  
  while (length(line <- readLines(con, n = 1)) > 0) {
    if (startsWith(line, ">")) {
      # write previous sequence if needed
      if (length(current_seq) > 0 && write_seq) {
        writeLines(current_seq, out_con)
      }
      seq_idx <- seq_idx + 1L
      # (keep_pointer <= length(keep_idx) Make sure not to refer to an element in keep_idx, that do not exists
      write_seq <- (keep_pointer <= length(keep_idx) && seq_idx == keep_idx[keep_pointer])
      if (write_seq && keep_pointer <= length(keep_idx)) keep_pointer <- keep_pointer + 1L
      current_seq <- if (write_seq) line else character(0)
    } else if (write_seq) {
      current_seq <- c(current_seq, line)
    }
  }
  # write last sequence if needed
  if (length(current_seq) > 0 && write_seq) {
    writeLines(current_seq, out_con)
  }
  
  close(con)
  close(out_con)
  
  return(invisible(n))
}

#' Randomly sample sequences from a FASTA file (Linux, VSEARCH-based)
#'
#' Randomly selects `n` sequences from an input FASTA file using `fastx_subsample`
#' from VSEARCH.
#'
#' This function is Linux-specific. For a cross-platform alternative, use
#' `random_sample_r()`.
#'
#' The input FASTA file can be uncompressed or gzipped. Output compression is
#' inferred from the output file extension (i.e., `.gz` produces a gzipped file).
#'
#' If the input file contains fewer than or exactly `n` sequences, it is copied
#' to the output file, preserving the requested compression format.
#'
#' @param fasta Character string: path to the input FASTA file (optionally gzipped).
#' @param outfile Character string: path to the output FASTA file. If it ends in `.gz`,
#'   the output will be gzip-compressed.
#' @param n Positive integer: number of sequences to randomly sample.
#' @param vsearch_path Character string: path to the VSEARCH executable.
#' @param randseed Integer or NULL: random seed for sampling. If NULL or 0 (default),
#'   a pseudo-random seed is used. A non-zero seed ensures reproducible results.
#' @param compress_method Character or logical: compression method used for output files.
#'   Must be one of `"pigz"`, `"gzip"`, or `"R"`.
#'   - `"pigz"` requires `pigz` to be installed and available in the system PATH
#'     (or specified via `pigz_path`).
#'   - `"gzip"` is Linux-only.
#'   - `"R"` uses `R.utils` and is cross-platform but slower.
#'   Performance (fastest → slowest): `pigz` > `gzip` > `R`.
#' @param pigz_path Character string: path to the `pigz` executable. Only required if
#'   `pigz` is used and not available in the system PATH.
#' @param num_threads Positive integer: number of CPU threads to use. If 0, all available CPUs are used.
#' @param quiet Logical: if TRUE, suppress informational messages; only warnings and errors are shown.
#' @return Invisibly returns the number of sequences written to the output file.
#' @examples
#' \dontrun{
#' random_sample_linux(
#'   fasta = "all_sequences.fasta.gz",
#'   outfile = "subset_100.fasta.gz",
#'   vsearch_path = "vsearch",
#'   n = 100,
#'   randseed = 123,
#'   quiet = FALSE
#' )
#' }
#' @export
#' 
random_sample_linux <- function(fasta, 
                                   outfile, 
                                   n=1000000, 
                                   vsearch_path="vsearch", 
                                   compress_method= "R",
                                   pigz_path="pigz",
                                   randseed = NULL, 
                                   quiet=TRUE, 
                                   num_threads=0) {
  
  if(!is_linux()){
    stop("This parameter setting is suppored only on linux")
  }
  
  check_dir(outfile, is_file=TRUE)
  
  if(is.null(randseed)){
    randseed <- 0
  }
  
  # --- Count sequences
  if(endsWith(fasta, '.gz')){
    cmd <- paste("zcat", fasta, "| grep '>' | wc -l", sep=" ")
  }else{
    cmd <- paste("grep '>' ",fasta, "| wc -l", sep=" ")
  }
  total <- as.integer(system(cmd, intern=TRUE))
  
  if(!quiet){ cat("Total sequences:", total, "\n")}
  
  # --- Check if n >= total
  if (n >= total) {
    txt <- paste(fasta, "contains", total, "sequences.\n", "The input file is copied to output\n")
    warning(txt)
    
    # copy input to outfile, and compress/uncompress if necessary
    if (grepl("\\.gz$", fasta) && !grepl("\\.gz$", outfile)) { #input gz, output not
      outfile <- smart_gzip(file=fasta,
                             outfile = outfile,
                             remove = FALSE,
                             pigz_path = pigz_path,
                             method = compress_method,
                             num_threads = num_threads,
                             quiet = quiet,
                             compress = FALSE)
    } else if (!grepl("\\.gz$", fasta) && grepl("\\.gz$", outfile)) { #input not uncompressed - output gz
      outfile <- smart_gzip(file=fasta,
                            outfile = outfile,
                            remove = FALSE,
                            pigz_path = pigz_path,
                            method = compress_method,
                            num_threads = num_threads,
                            quiet = quiet,
                            compress = TRUE)
    } else{ # same compression, simply copy file
      file.copy(from = fasta, to = outfile, overwrite = TRUE)
    }
    return(invisible(total))
  }
  
  # --- total > n => random sample with vsearch
  # do not transform large numbers to scientific forms, since it would lead to an error in vsearch
  options(scipen=100)
  if(grepl("\\.gz$", outfile)){ # output should be compressed
    outfile_tmp <- gsub("\\.gz", "", outfile) # vsearch makes decompressed files
  }else{
    outfile_tmp <- outfile
  }
  
  ##### run vsearch
  # Build argument vector
  args <- c("--fastx_subsample", fasta,
            "--fastaout", outfile_tmp,
            "--sample_size",  n,
            "--randseed", randseed
  )
  if(num_threads > 0){
    args <- append(args, c("--threads", num_threads), after = 2)
  }
  run_system2(vsearch_path, args, quiet=quiet)
  options(scipen=0)
  
  if(grepl("\\.gz$", outfile)){ # output should be compressed
    outfile <- smart_gzip(file=outfile_tmp,
                         outfile = outfile,
                         remove = TRUE,
                         pigz_path = pigz_path,
                         method = compress_method,
                         num_threads = num_threads,
                         quiet = quiet,
                         compress = TRUE)
  }
  return(invisible(n))
}

#' Randomly subsample sequences from FASTA files
#' 
#' Randomly selects `n` sequences (without replacement) from each FASTA file
#' listed in the `fastainfo` data frame.
#' This function is useful for standardizing sequencing depth across libraries.
#'
#' It can be applied before or after demultiplexing (`demultiplex_and_trim`),
#' but is recommended after `merge_fastq_pairs` and before `demultiplex_and_trim`.
#'
#' Two implementations are available:
#'   - VSEARCH-based (`use_vsearch = TRUE`): fast, but only available on Linux-like systems.
#'   - R-based (`use_vsearch = FALSE`): slower, but fully cross-platform.
#'  
#' @param fastainfo Data frame or CSV file containing a `fasta` column with input file names.
#'   Files may be gzipped.
#' @param fasta_dir Character string: directory containing the input FASTA files.
#' @param n Positive integer: number of sequences to sample from each file.
#' @param outdir Character string: directory where output FASTA files are written.
#' @param use_vsearch Logical: if TRUE, uses VSEARCH for subsampling (Linux only).
#'   Otherwise uses a cross-platform R-based implementation.
#' @param vsearch_path Character string: path to the VSEARCH executable.
#' @param randseed Integer or NULL: random seed for reproducibility. If NULL or 0,
#'   a pseudo-random seed is used. A non-zero seed ensures reproducible results.
#' @param num_threads Positive integer: number of CPU threads to use. If 0, all available CPUs are used.
#' @param compress Logical: if TRUE, output FASTA files are gzip-compressed.
#' @param compress_method Character string: compression method (used only when `use_vsearch = TRUE`).
#'   Options are:
#'   - `"pigz"`: fastest; requires `pigz` installed (or `pigz_path` specified)
#'   - `"gzip"`: slower; available on Linux-like systems
#'   - `"R"`: cross-platform; uses `R.utils` but is slowest
#'   Performance order: `pigz` > `gzip` > `R`
#' @param pigz_path Character string: path to the `pigz` executable. Required only if `pigz`
#'   is used and not available in the system PATH.
#' @param quiet Logical: if TRUE, suppress informational messages; only warnings and errors are shown.
#' @param sep Character string: field separator used in input and output CSV files.
#' @return Updated input data frame with adjusted file names (if needed) and updated read counts.
#' @examples
#' \dontrun{
#' # Fast version (Linux)
#' fastainfo_df <- subsample_fasta(
#'   fastainfo_df, merged_dir, "random_seq",
#'   use_vsearch = TRUE, n = 10000, compress = TRUE
#' )
#'
#' # Cross-platform version
#' fastainfo_df <- subsample_fasta(
#'   fastainfo_df, merged_dir, "random_seq",
#'   use_vsearch = FALSE, n = 10000, compress = TRUE
#' )
#' }
#' @export
#' 
subsample_fasta <- function(fastainfo, 
                       n,
                       fasta_dir,
                       outdir, 
                       use_vsearch=FALSE,
                       vsearch_path="vsearch",
                       randseed=NULL, 
                       compress_method="R",
                       pigz_path="pigz",
                       num_threads=0,
                       compress=FALSE, 
                       sep=",",
                       quiet=TRUE){
  
  # can accept df or file as an input
  if(is.character(fastainfo)){
    # read known occurrences
    fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  }else{
    fastainfo_df <- fastainfo
  }
  
  fasta_dir = check_dir(fasta_dir)
  outdir = check_dir(outdir)
  
  unique_fasta <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(unique_fasta)){ # go through all fasta files
    input_fasta <- unique_fasta[i]
    input_fasta_p <- file.path(fasta_dir, input_fasta)
    
    # outfine name is the same as the input, but different folder
    # adjusted compression
    outfile <- input_fasta
    if(!compress && grepl("\\.gz$", input_fasta)){
      outfile <- sub("\\.gz$", "", outfile)
    }
    if(compress && !grepl("\\.gz$", input_fasta)){
      outfile <- paste(outfile, "gz", sep=".")
    }
    outfile_p <- file.path(outdir, outfile)
    
    ##### Change algo if windows and use_vsearch=FALSE to avoid stopping the run 
    if(!is_linux() && use_vsearch){
      warning("The fastx_subsample commande in VSEARCH is not available for Windows\n
              A slower, but cross-platform function (random_sample_r)\n
              will be used for random sampling sequences (random_sample_r).")
      use_vsearch <- FALSE
    }
    
    if(use_vsearch){

      seqn <- random_sample_linux(fasta=input_fasta_p,
                                     outfile=outfile_p,
                                     n=n,
                                     vsearch_path=vsearch_path,
                                     compress_method=compress_method,
                                     pigz_path=pigz_path,
                                     randseed = randseed,
                                     quiet=quiet,
                                     num_threads=num_threads)
    }else{
      seqn <- random_sample_r(fasta=input_fasta_p, 
                                 outfile = outfile_p,
                                 n=n, 
                                 randseed = randseed,
                                 quiet=quiet)
    }
    
    fastainfo_df$fasta[which(fastainfo_df$fasta == input_fasta)] <- outfile
    fastainfo_df$read_count[which(fastainfo_df$fasta == outfile)] <- seqn
  } # end for
  write.table(fastainfo_df, file = file.path(outdir, "fastainfo.csv"),  row.names = F, sep=sep)
  return(fastainfo_df)
}


#' Concatenate contents of files with identical names
#'
#' Reads files from multiple directories and concatenates the contents of files
#' sharing the same filename into a single output file.
#'
#' For each unique filename matching `pattern`, all corresponding files found in
#' `dirs` are read in order and merged into a file of the same name in `outdir`.
#' 
#' This function works both on uncompressed and gz compressed files.
#'
#' @param dirs A character vector of input directories to search for files.
#' @param outdir A character string specifying the output directory where pooled
#'   files will be written.
#' @param pattern A regular expression used to select filenames within the input
#'   directories.
#' @param quiet Logical: if TRUE, suppress informational messages; only warnings and errors are shown.
#'
#' @return Invisibly data frame of files and directories 
#'
#' @examples
#' \dontrun{
#' concatenate_files(
#'   dirs = c("run1", "run2", "run3"),
#'   outdir = "pooled",
#'   pattern = "\\.fastq$"
#' )
#' }
#'
#' @export
concatenate_files <- function(dirs, outdir, pattern = "\\.", quiet=TRUE) {
  
  outdir = check_dir(outdir, is_file=FALSE)
  
  df <- data.frame(
    file = character(),
    dir = character(),
    stringsAsFactors = FALSE
  )
  
  # make df with filenames and dir names as columns
  for (dir in dirs) {
    
    files <- list.files(
      path = dir,
      pattern = pattern
    )
    
    tmp <- data.frame(
      file = files,
      dir = dir,
      stringsAsFactors = FALSE
    )
    
    df <- rbind(df, tmp)
  }
  
  # loop over unique filenames
  for (f in unique(df$file)) {
    if(!quiet){
      print(f)
    }
    sub_df <- df %>%
      filter(file == f)
    
    out <- file.path(outdir, f)
    
    gz <- grepl("\\.gz$", f)
    if (gz) {
      con_out <- gzfile(out, open = "wt")
    } else {
      con_out <- file(out, open = "w")
    }
    
    for (i in 1:nrow(sub_df)) {
      
      path <- file.path(sub_df$dir[i], sub_df$file[i])
      
      if (gz) {
        con_in <- gzfile(path, open = "rt")
      } else {
        con_in <- file(path, open = "r")
      }
      lines <- readLines(con_in)
      writeLines(lines, con_out)
      close(con_in)
    }
    close(con_out)
  }
  invisible(df)
}

#' Demultiplex fastq file pairs and trim tags and primers (reverse strand not checked)
#' 
#' Same as `demultiplex_fastq_pairs`, but without checking the reverse-complement 
#' of the sequences.
#'  
#' FASTQ file pairs are first demultiplexed by requiring a perfect match between
#' the tag sequence and the 5' end of the read. The resulting files are then
#' trimmed to remove primer sequences using less stringent matching parameters
#' (controlled by `cutadapt_error_rate`). A match between the primer and the 5'
#' end of the read is required for the read to be trimmed and retained. Matching
#' of the 3' primer is optional
#'  
#' Input files can be compressed or uncompressed. Output compression is 
#' controlled by `compress`.
#'  
#' @param fastqinfo Data frame or path to a CSV file with the following columns: 
#'   `tag_fw`, `primer_fw`, `tag_rv`, `primer_rv`, 
#'   `sample`, `sample_type` (mock/negative/real), 
#'   `habitat` (optional), `replicate`, `fastq_fw`, 
#'   `fastq_rv`
#' @param fastq_dir Character string specifying the directory containing input 
#'   FASTQ files (listed in the `fastq_fw` and `fastq_rv` columns of `fastqinfo`).
#' @param cutadapt_path Character string specifying the path to the 
#'   `cutadapt` executable.
#' @param num_threads Positive integer specifying the number of CPU threads to 
#'   use. If `0`, all available CPUs are used.
#' @param outdir Character string specifying the output directory.
#' @param tag_to_end Logical. If `TRUE`, tags are assumed to be located at the 
#'   extremities of reads (starting at the first base).
#' @param primer_to_end Logical. If `TRUE`, primers are assumed to follow 
#'   directly after tags (i.e., no heterogeneity spacer).
#' @param cutadapt_error_rate Numeric value between 0 and 1 specifying the 
#'   maximum allowed error rate between primers and reads (exact match is 
#'   required for tags).
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param compress Logical. If `TRUE`, compress output files using gzip.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Data frame similar to the input `fastqinfo` file, 
#'   but contains the output fastq file names and read counts.
#' 
#' @examples
#' \dontrun{
#' fastqinfo_df <- demultiplex_fastq_pairs_strand_plus(
#'   fastqinfo = fastqinfo_df,
#'   fastq_dir = "data/fastq",
#'   outdir = "data/fastq_demultiplexed",
#'   tag_to_end = TRUE,
#'   primer_to_end = TRUE,
#'   sep = ","
#' )
#' }
#' 
#' @export
#' 

demultiplex_fastq_pairs_strand_plus <- function(fastqinfo, 
                                                fastq_dir, 
                                                outdir, 
                                                cutadapt_path="cutadapt", 
                                                num_threads=0,
                                                tag_to_end=T, 
                                                primer_to_end=T, 
                                                cutadapt_error_rate=0.1,
                                                sep=",",  
                                                compress=F, 
                                                quiet=T
){
  
  # do the complete job of demultiplexing and trimming of input file without checking the reverse sequences
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  fastq_dir = check_dir(fastq_dir)
  outdir = check_dir(outdir)
  
  # can accept df or file as an input
  if(is.character(fastqinfo)){
    # read known occurrences
    fastqinfo_df <- read.csv(fastqinfo, header=T, sep=sep)
  }else{
    fastqinfo_df <- fastqinfo
  }
  
  check_file_info(fastqinfo_df, fastq_dir, file_type="fastqinfo", sep=",", quiet=TRUE)
  
  # upper case for all primers and tags
  fastqinfo_df$tag_fw <- toupper(fastqinfo_df$tag_fw)
  fastqinfo_df$tag_rv <- toupper(fastqinfo_df$tag_rv)
  fastqinfo_df$primer_fw <- toupper(fastqinfo_df$primer_fw)
  fastqinfo_df$primer_rv <- toupper(fastqinfo_df$primer_rv)
  # make columns for output filenames
  fastqinfo_df$fastq_fw_demultiplexed <- NA
  fastqinfo_df$fastq_rv_demultiplexed <- NA
  
  # get unique list of input fastq file pairs
  fastqs <- fastqinfo_df %>%
    select(fastq_fw, fastq_rv) %>%
    distinct()
  
  for(i in 1:nrow(fastqs)){ # for each input fastq pair
    # select lines in fastqinfo_df that corresponds to a given input fasta file
    fastq_fw_local <- fastqs$fastq_fw[i]
    fastq_rv_local <- fastqs$fastq_rv[i]
    df <- fastqinfo_df %>%
      filter(fastq_fw==fastq_fw_local & fastq_rv==fastq_rv_local)
    
    # Make a tmp_dir_fastq in tempdir specific to a fastq file pair. It will contain the tagtrimmed files.
    # This can be deleted at the end and avoid reusing tagtrimmed files created for a previous fastq files
    tmp_fastq <- paste(fastq_fw_local, "_", trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    tmp_dir_fastq <- file.path(tempdir(), tmp_fastq)
    
    # Delete tmp_dir_fastq if exists (previous run crushed before deleting it) 
    if (dir.exists(tmp_dir_fastq)) {
      unlink(tmp_dir_fastq, recursive = TRUE)
    }
    # Create it 
    dir.create(tmp_dir_fastq)
    
    # make a tag_fw.fasta and tag_rv.fasta files with all tag of the fastq pairs to be demultiplexed
    tag_files <- write_cutadapt_adapter_fastq(fastqinfo_df, 
                                              fastq_file=fastq_fw_local, 
                                              tag_to_end=tag_to_end, 
                                              outdir=tmp_dir_fastq
    )
    
    # add path
    fastq_fw_local <- file.path(fastq_dir, fastq_fw_local)
    fastq_rv_local <- file.path(fastq_dir, fastq_rv_local)
    
    ##### run demultiplexing
    g <- paste("file:", tag_files[1], sep="")
    G <- paste("file:", tag_files[2], sep="")
    out_fw <- file.path(tmp_dir_fastq, "tagtrimmed-{name1}-{name2}_fw.fastq") 
    out_rv <- file.path(tmp_dir_fastq, "tagtrimmed-{name1}-{name2}_rv.fastq") 
    args <- c(
      "-e", "0",
      "--no-indels",
      "--trimmed-only",
      "-g",  shQuote(g),
      "-G",  shQuote(G),
      "-o",  shQuote(out_fw),
      "-p",  shQuote(out_rv),
      fastq_fw_local, fastq_rv_local
    )
    if(num_threads > 0){
      args <- append(args, c("--cores", num_threads), after=2)
    }
    if(quiet){
      args <- append(args, c("--quiet"), after=2)
    }
    run_system2(cutadapt_path, args, quiet=quiet)
    
    ############ primer trimming with less stingent conditions
    # for a given marker, there is only one primer combination
    primer_fwl <- df[1,"primer_fw"]
    primer_rvl <- df[1,"primer_rv"]
    
    for(f in 1:nrow(df)){# go through each de-multiplexed, tag-trimmed file and trim primers
      outfile_fw <- paste(df[f,"sample"], "-", df[f,"replicate"], "_fw", sep="")
      outfile_rv <- paste(df[f,"sample"],  "-", df[f,"replicate"], "_rv", sep="")
      outfile_fw <- paste(outfile_fw, ".fastq", sep="")
      outfile_rv <- paste(outfile_rv, ".fastq", sep="")
      if(compress){
        outfile_fw <- paste(outfile_fw, ".gz", sep="") # cutadapt detects from filename, if outfiles should be compressed
        outfile_rv <- paste(outfile_rv, ".gz", sep="")
      }
      # complete fastqinfo_df with output fastq names
      fastqinfo_df$fastq_fw_demultiplexed[
        which(fastqinfo_df$sample==df[f,"sample"] & 
                fastqinfo_df$replicate==df[f,"replicate"])
      ]<- outfile_fw
      fastqinfo_df$fastq_rv_demultiplexed[
        which(fastqinfo_df$sample==df[f,"sample"] & 
                fastqinfo_df$replicate==df[f,"replicate"])
      ]<- outfile_rv
      
      # add path to output file
      primer_trimmed_fw <- file.path(outdir, outfile_fw)
      primer_trimmed_rv <- file.path(outdir, outfile_rv)
      tag_trimmed_fw <- paste("tagtrimmed-", 
                              df[f,"tag_fw"], "-", 
                              df[f,"tag_rv"], 
                              "_fw.fastq", 
                              sep="")
      tag_trimmed_rv <- paste("tagtrimmed-", 
                              df[f,"tag_fw"], "-", 
                              df[f,"tag_rv"], 
                              "_rv.fastq", 
                              sep="")
      tag_trimmed_fw <- file.path(tmp_dir_fastq, tag_trimmed_fw)
      tag_trimmed_rv <- file.path(tmp_dir_fastq, tag_trimmed_rv)
      
      primer_rvl_rc <- reverse_complement(primer_rvl)
      primer_fwl_rc <- reverse_complement(primer_fwl)
      
      ##### run primer trimming
      if(primer_to_end){
#        g <- paste("^", primer_fwl, sep="")
#        G <- paste("^", primer_rvl, sep="")
        a <- paste("^", primer_fwl, "...", primer_rvl_rc, sep="")
        A <- paste("^", primer_rvl, "...", primer_fwl_rc, sep="")
      }
      else{
#        g <- paste(primer_fwl, ";min_overlap=", nchar(primer_fwl), sep="")
#        G <- paste(primer_rvl, ";min_overlap=", nchar(primer_rvl), sep="")
        a <- paste(primer_fwl, ";required;min_overlap=", nchar(primer_fwl), "...",  primer_rvl_rc, sep="")
        A <- paste(primer_rvl, ";required;min_overlap=", nchar(primer_rvl), "...",  primer_fwl_rc, sep="")
      }
      args <- c(
        "-e", cutadapt_error_rate,
        "--no-indels",
        "--trimmed-only",
#        "-g",  shQuote(g),
#        "-G",  shQuote(G),
        "-a",  shQuote(a),
        "-A",  shQuote(A),
        "-o",  primer_trimmed_fw,
        "-p",  primer_trimmed_rv,
        tag_trimmed_fw, tag_trimmed_rv
      )
      if(num_threads > 0){
        args <- append(args, c("--cores", num_threads), after=2)
      }
      if(quiet){
        args <- append(args, c("--quiet"), after=2)
      }
      run_system2(cutadapt_path, args, quiet=quiet)
      
    } # end tag-trimmed 
    # delete the tmp dir with the tag-trimmed files
    unlink(tmp_dir_fastq, recursive = TRUE)
  }# end fastq
  
  # make sampleinfo file
  fastqinfo_df <- fastqinfo_df %>%
    select(-fastq_fw, -fastq_rv) %>%
    rename(fastq_fw = fastq_fw_demultiplexed, fastq_rv = fastq_rv_demultiplexed)
  
  return(fastqinfo_df)
  
}

#' Make a FASTA file with forward or reverse adapters
#' 
#' Create two FASTA files containing forward or reverse tags
#' formatted for `cutadapt`. 
#' This file is used by `demultiplex_and_trim_fastq` to demultiplex input FASTQ files pairs.
#' 
#' @param fastqinfo_df Data frame with columns: `tag_fw`, `tag_rv`, `fastq_fw` and `fastq_rv`.
#' @param fastq Character string specifying the `fastq_fw` file to be demultiplexed 
#'   (must be present in the `fastq_fw` column of `fastqinfo_df`).
#' @param outdir Character string specifying the output directory.
#' @param tag_to_end Logical. If `TRUE`, tags are assumed to be located at the 
#'   extremities of reads.
#' 
#' @return Vector with the output fasta files, or `NA` if all tags 
#'   are `NA` in `fastqinfo_df` for the given `fastq` file.
#' 
#' @examples 
#' \dontrun{
#' write_cutadapt_adapter_fastq(
#'   fastqinfo_df = fastqinfo_df, 
#'   fastq_file = "fastq_file", 
#'   tag_to_end = FALSE, 
#'   outdir = "data/out"
#' )
#' }
#' 
#' @export

write_cutadapt_adapter_fastq <- function(
  fastqinfo_df, 
  fastq_file, 
  outdir, 
  tag_to_end=TRUE
){
  
  # select all tags for the fastq file
  tags <- fastqinfo_df %>%
    filter(fastq_fw==fastq_file) %>%
    select(tag_fw, tag_rv)
  
  
  tags_fw <- unique(toupper(tags$tag_fw))
  tags_rv <- unique(toupper(tags$tag_rv))
  
  # return NA if all tags are NA
  if(length(tags_fw) == 1){ # only one tag
    if(is.na(tags_fw[1])){ # no tags
      return(NA)
    }
  }
  if(length(tags_rv) == 1){ # only one tag
    if(is.na(tags_rv[1])){ # no tags
      return(NA)
    }
  }
  
  # Specify the file path
  outdir = check_dir(outdir)
  tag_file_fw <- file.path(outdir, "tags_fw.fasta")
  if (tag_to_end) {
    text <- as.vector(rbind(
      paste0(">", tags_fw),
      paste0("^", tags_fw)
    ))
  } else {
    text <- as.vector(rbind(
      paste0(">", tags_fw),
      paste0(tags_fw, ";min_overlap=", nchar(tags_fw))
    ))
  }
  writeLines(text, tag_file_fw)
  
  tag_file_rv <- file.path(outdir, "tags_rv.fasta")
  
  if (tag_to_end) {
    text <- as.vector(rbind(
      paste0(">", tags_rv),
      paste0("^", tags_rv)
    ))
  } else {
    text <- as.vector(rbind(
      paste0(">", tags_rv),
      paste0(tags_rv, ";min_overlap=", nchar(tags_rv))
    ))
  }
  writeLines(text, tag_file_rv)
  
  files <- c(tag_file_fw, tag_file_rv)
  return(files)
}

#' Demultiplex FASTQ file pairs and trim tags and primers
#'
#' FASTQ file pairs are first demultiplexed by requiring a perfect match between
#' the tag sequence and the 5' end of the read. The resulting files are then
#' trimmed to remove primer sequences using less stringent matching parameters
#' (controlled by `cutadapt_error_rate`). A match between the primer and the 5'
#' end of the read is required for the read to be trimmed and retained. Matching
#' of the 3' primer is optional.
#'  
#' If `check_reverse = TRUE`, forward and reverse reads are swapped, demultiplexed 
#' again, and trimmed for 5' primers.
#' 
#' When the same set of tags is used at both ends of the reads, demultiplexing 
#' both the original and swapped orientations may produce incorrect assignments.
#' These incorrectly assigned reads are subsequently removed during primer trimming.
#'
#' Input files can be compressed or uncompressed. Output compression is 
#' controlled by `compress`.
#'  
#' @param fastqinfo Data frame or path to a CSV file with the following columns: 
#'   `tag_fw`, `primer_fw`, `tag_rv`, `primer_rv`, 
#'   `sample`, `sample_type` (mock/negative/real), 
#'   `habitat` (optional), `replicate`, `fastq_fw`, 
#'   `fastq_rv`
#' @param fastq_dir Character string specifying the directory containing input 
#'   FASTQ files (listed in the `fastq_fw` and `fastq_rv` columns of `fastqinfo`).
#' @param cutadapt_path Character string specifying the path to the 
#'   `cutadapt` executable. 
#' @param compress Logical. If `TRUE`, compress output files using gzip.
#' @param num_threads Positive integer specifying the number of CPU threads to 
#'   use. If `0`, all available CPUs are used.
#' @param outdir Character string specifying the output directory.
#' @param check_reverse Logical. If `TRUE`, also check reverse-complemented 
#'   sequences from the input FASTA files.
#' @param tag_to_end Logical. If `TRUE`, tags are expected to be located 
#'   at the extremities of reads (starting at the first base).
#' @param primer_to_end Logical. If `TRUE`, primers are assumed to follow 
#'   directly after tags (i.e., no heterogeneity spacer).
#' @param cutadapt_error_rate Numeric value between 0 and 1 specifying the 
#'   maximum allowed error rate between primers and reads (exact match is 
#'   required for tags).
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#' @param quiet Logical. If `TRUE`, suppress informational messages and 
#'   only display warnings or errors.
#' 
#' @return Data frame similar to the input `fastqinfo` file, 
#'   but contains the output fastq file names and read counts.
#' 
#' @examples
#' \dontrun{
#' fastqinfo_df <- demultiplex_fastq_pairs(
#'   fastqinfo = fastqinfo_df,
#'   fastq_dir = "data/fastq",
#'   outdir = "data/fastq_demultiplexed",
#'   tag_to_end = TRUE,
#'   primer_to_end = TRUE,
#'   sep = ","
#' )
#' }
#' 
#' @export

demultiplex_fastq_pairs <- function(
  fastqinfo, 
  fastq_dir, 
  outdir, 
  cutadapt_path="cutadapt",
  check_reverse=FALSE, 
  num_threads=0,
  tag_to_end=TRUE, 
  primer_to_end=TRUE, 
  cutadapt_error_rate=0.1,
  sep=",",
  compress=FALSE,
  quiet=T
  ){
  
  fastq_dir <- check_dir(fastq_dir)
  outdir <- check_dir(outdir)
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # can accept df or file as an input
  if(is.character(fastqinfo)){
    # read known occurrences
    fastqinfo_df <- read.csv(fastqinfo, header=T, sep=sep)
  }else{
    fastqinfo_df <- fastqinfo
  }
  
  check_file_info(file=fastqinfo_df, dir=fastq_dir, file_type="fastqinfo", sep=sep, quiet=TRUE)
  
  #########
  # demultiplex_fastq_pairs_strand_plus does the whole demultilexing, trimming and compress on the + strand
  # If sequences are not oriented, the -strand should be checked => 
  # swap fw and reverse input files and run again demultiplex_fastq_pairs_strand_plus
  # pool the output
  
  
  # run demultiplex_and_trim_strand_plus of plus strand and on - strand after 
  # swapping fw and rev tags and primers,
  # take the reverse complement of the -strand results (vsearch)
  # pool the results of the 2 strands
  # compress if necessary
  
  # run on strand +
  if(check_reverse){
    
    # make temp dirs
    fw_tmp_dir <- paste('fw_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    fw_tmp_dir <- file.path(tempdir(), fw_tmp_dir)
    check_dir(fw_tmp_dir)
    rv_tmp_dir <- paste('rv_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    rv_tmp_dir <- file.path(tempdir(), rv_tmp_dir)
    check_dir(rv_tmp_dir)
    
    #### use +strand, output to sorted_dir, uncompressed
    fastqinfo_demultiplexed_fw <- demultiplex_fastq_pairs_strand_plus(fastqinfo_df, 
                                                                      fastq_dir=fastq_dir, 
                                                                      outdir=fw_tmp_dir, 
                                                                      cutadapt_path=cutadapt_path, 
                                                                      num_threads=num_threads,
                                                                      tag_to_end=tag_to_end, 
                                                                      primer_to_end=primer_to_end, 
                                                                      cutadapt_error_rate=cutadapt_error_rate, 
                                                                      sep=sep, 
                                                                      compress=compress,
                                                                      quiet=quiet
    )
    
    #### use - strand
    # swap fw and rv tags and primers
    fastqinfo_df_tmp <- fastqinfo_df %>%
      rename(fastq_fw = fastq_rv, fastq_rv = fastq_fw)
    
    
    # run demultiplex_and_trim on for reverse strand
    fastqinfo_demultiplexed_rv <- demultiplex_fastq_pairs_strand_plus(fastqinfo_df_tmp, 
                                                                      fastq_dir=fastq_dir, 
                                                                      outdir=rv_tmp_dir, 
                                                                      cutadapt_path=cutadapt_path, 
                                                                      num_threads = num_threads,
                                                                      tag_to_end=tag_to_end, 
                                                                      primer_to_end=primer_to_end, 
                                                                      cutadapt_error_rate=cutadapt_error_rate, 
                                                                      sep=sep,
                                                                      compress=compress, 
                                                                      quiet=quiet
    )
    
    # concatenate files
    for(i in 1:nrow(fastqinfo_demultiplexed_fw)){
      
      # fw reads from original and swapped demultiplexing
      filename_base <- fastqinfo_demultiplexed_fw$fastq_fw[i]
      fw_in <- file.path(fw_tmp_dir, filename_base)
      rv_in <- file.path(rv_tmp_dir, filename_base)
      out <- file.path(outdir, filename_base)
      concat_files(files= c(fw_in, rv_in), outfile=out)
      
      filename_base <- fastqinfo_demultiplexed_fw$fastq_rv[i]
      # rv reads from original and swapped demultiplexing
      fw_in <- file.path(fw_tmp_dir, filename_base)
      rv_in <- file.path(rv_tmp_dir, filename_base)
      out <- file.path(outdir, filename_base)
      concat_files(files= c(fw_in, rv_in), outfile=out)
    }
    # fastqinfo_demultiplexed_fw and fastqinfo_demultiplexed_rv are identical
    fastqinfo_demultiplexed <- fastqinfo_demultiplexed_fw
    
    
    # delete temporary dirs
    unlink(fw_tmp_dir, recursive = TRUE)
    unlink(rv_tmp_dir, recursive = TRUE)
    
  }
  else{
    # check only + strand
    fastqinfo_demultiplexed <- demultiplex_fastq_pairs_strand_plus(fastqinfo_df, 
                                                                   fastq_dir=fastq_dir,
                                                                   outdir=outdir, 
                                                                   cutadapt_path=cutadapt_path, 
                                                                   num_threads = num_threads,
                                                                   tag_to_end=tag_to_end, 
                                                                   primer_to_end=primer_to_end, 
                                                                   cutadapt_error_rate=cutadapt_error_rate, 
                                                                   sep=sep, 
                                                                   compress=compress, 
                                                                   quiet=quiet
    )
  }
  
  #fastqinfo_demultiplexed <- add_read_counts(fastqinfo_demultiplexed, dir=outdir)
  
  read_count <- count_reads_in_dir(
    dir=outdir, 
    pattern="_fw.fastq", 
    file_type="fastq"
  )
  fastqinfo_demultiplexed <- fastqinfo_demultiplexed %>%
    left_join(read_count, by=c("fastq_fw" = "filename"))
  
  write.table(fastqinfo_demultiplexed, file = file.path(outdir, "fastqinfo.csv"),  row.names = F, sep=sep)
  
  return(fastqinfo_demultiplexed)
}



#' Concatenate files in a portable, streaming way
#'
#' This function concatenates multiple files by copying their raw bytes
#' sequentially into a single output file. It works for both plain text
#' files and gzip-compressed files (`.gz`), and is fully cross-platform
#' (Windows, macOS, Linux) without relying on system commands.
#'
#' For gzip files, concatenation produces a valid gzip stream because the
#' format supports concatenated members.
#'
#' @param files Character vector of input file paths.
#' @param outfile Character string specifying the output file path.
#' @param chunk_size Integer. Number of bytes to read per iteration.
#'   Larger values are faster but use more memory. Default is 1 MB.
#'
#' @return The output file path (invisibly).
#'
#' @details
#' The function performs binary-safe streaming using `readBin` and
#' `writeBin`. It does not decompress or interpret file contents.
#'
#' All input files must exist. The function stops otherwise.
#'
#' Mixing compressed and uncompressed files is technically possible but
#' generally not meaningful.
#'
#' @examples
#' \dontrun{
#' concat_files(c("a.txt", "b.txt"), "out.txt")
#'
#' concat_files(c("sample1.fastq.gz", "sample2.fastq.gz"),
#'              "merged.fastq.gz")
#' }
#'
#' @export
concat_files <- function(files, outfile, chunk_size = 1024^2) {
  
  # Basic input validation
  stopifnot(is.character(files), length(files) > 0)
  stopifnot(is.character(outfile), length(outfile) == 1)
  stopifnot(all(file.exists(files)))
  
  # Open output connection in binary write mode
  out <- file(outfile, open = "wb")
  on.exit(close(out), add = TRUE)
  
  # Loop over input files
  for (f in files) {
    
    # Open input file in binary read mode
    in_con <- file(f, open = "rb")
    
    # Stream file content in chunks
    repeat {
      
      # Read a block of raw bytes
      chunk <- readBin(in_con, what = "raw", n = chunk_size)
      
      # Stop when end-of-file is reached
      if (!length(chunk)) break
      
      # Write bytes to output
      writeBin(chunk, out)
    }
    
    # Close input connection for this file
    close(in_con)
  }
  
  # Return output path invisibly
  invisible(normalizePath(outfile, mustWork = FALSE))
  
  
}
