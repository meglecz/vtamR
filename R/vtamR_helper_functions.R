#' @importFrom dplyr filter mutate group_by select summarize summarise arrange last
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first if_else
#' @importFrom dplyr add_row slice_head
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym :=
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL

#'
#'
#' Download from Zenodo
#'
#' Downloads a file from Zenodo and optionally extracts it if it is a
#' `.tar.gz` archive.
#'
#' @param filename Character string specifying the name of the file to download.
#' @param url Character string; URL of the Zenodo record (copy from browser address bar).
#' @param dest_dir Character string; path to the directory where the file should be downloaded.
#' @param untar Logical; if `TRUE`, the archive is extracted after download.
#' @param quiet Logical; if `TRUE`, suppress informational messages and show only warnings or errors.
#'
#' @return Invisible full path to the directory where the file was downloaded.
#'
#' @examples
#' \dontrun{
#' download_zenodo(
#'   filename = "COInr_2025_05_23.tar.gz",
#'   url = "https://zenodo.org/records/15515860",
#'   dest_dir = "~/vtamR",
#'   untar = TRUE,
#'   quiet = FALSE
#' )
#' }
#'
#' @export
#' 
download_zenodo <- function(
    filename,
    url,
    dest_dir = ".",
    untar = TRUE,
    quiet = FALSE
) {
  
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest_dir <- normalizePath(dest_dir, mustWork = FALSE)
  
  url <- paste0(url, "/files/", filename, "?download=1")
  destfile <- file.path(dest_dir, filename)
  
  if (!quiet) message("Downloading ", filename, " to ", dest_dir)
  utils::download.file(url, destfile)
  
  if (!quiet) message("Extracting archive...")
  utils::untar(destfile, exdir = dest_dir)
  
  unlink(destfile)
  
  if (!quiet) message("Extraction complete: ", dest_dir)
  return(invisible(dest_dir))
}

#' Download from OSF
#'
#' Downloads a file from the Open Science Framework (OSF).
#' If the file is a `.tar.gz` archive, it can optionally be extracted after download.
#'
#' @param filename Character string specifying the name of the file to download.
#' @param url Character string; direct download URL of the OSF file.
#' @param dest_dir Character string; path to the directory where the file should be downloaded.
#' @param untar Logical; if `TRUE`, the archive is extracted after download.
#' @param quiet Logical; if `TRUE`, suppress informational messages and show only warnings or errors.
#'
#' @return Invisible full path to the directory where the file was downloaded.
#'
#' @examples
#' \dontrun{
#' download_osf(
#'   filename = "COInr_for_vtam_2025_05_23_dbV5.tar.gz",
#'   url = "https://osf.io/download/jyhz6/",
#'   dest_dir = "~/vtamR/OSF",
#'   untar = TRUE,
#'   quiet = FALSE
#' )
#' }
#'
#' @export
#' 
download_osf <- function(
    filename,
    url,
    dest_dir = ".",
    untar = TRUE,
    quiet = FALSE
) {
  
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest_dir <- normalizePath(dest_dir, mustWork = FALSE)
  
  destfile <- file.path(dest_dir, filename)
  
  if (!quiet) message("Downloading ", filename, " to ", dest_dir)
  utils::download.file(url, destfile, method="auto", mode = "wb")
  
  if (!quiet) message("Extracting archive...")
  utils::untar(destfile, exdir = dest_dir)
  
  unlink(destfile)
  
  if (!quiet) message("Extraction complete: ", dest_dir)
  return(invisible(dest_dir))
}


#' Collect Information on Loaded R Packages
#'
#' Collects information on all packages that are currently attached to the
#' search path or loaded as namespaces in the current R session.
#'
#' For each package, the function reports whether it is attached or only
#' loaded, whether it is a base, recommended, or contributed package, and
#' its installed version. This information can be useful for documenting the
#' software environment used to run an analysis.
#' 
#' @param outfile Character string specifying the output file. If NULL 
#'   automatically derived from file.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#'   
#'
#' @return A data frame with one row per package and the following columns:
#' \describe{
#'   \item{Package}{Package name.}
#'   \item{Status}{Package status in the current R session:
#'     \code{"attached"} if the package is on the search path,
#'     or \code{"loaded only"} if only its namespace is loaded.}
#'   \item{Type}{Package type:
#'     \code{"base"}, \code{"recommended"}, or
#'     \code{"contributed"} (CRAN, Bioconductor, GitHub, or other
#'     user-installed packages).}
#'   \item{Version}{Installed package version.}
#' }
#'
#' @details
#' Attached packages are available on the R search path and their exported
#' objects can be used without qualification. Loaded-only packages are
#' available as namespaces, typically because they are imported by another
#' package, but are not attached to the search path.
#'
#' @examples
#' pkg_info <- collect_package_info()
#' head(pkg_info)
#'
#' @export
collect_package_info <- function(outfile=NULL, sep=",") {

    ## Attached packages
    attached <- sub("^package:", "", grep("^package:", search(), value = TRUE))
    
    ## Loaded namespaces
    loaded <- loadedNamespaces()
    
    ## All packages
    pkgs <- sort(unique(c(attached, loaded)))
    
    ## Status
    status <- ifelse(pkgs %in% attached,
                     "attached",
                     "loaded only")
    
    ## Installed package information
    ip <- installed.packages()
    
    ## Version
    version <- vapply(
      pkgs,
      function(pkg) as.character(packageVersion(pkg)),
      character(1)
    )
    
    
    type <- ifelse(is.na(ip[pkgs, "Priority"]),
                   "contributed",
                   ip[pkgs, "Priority"])
    
    df <- data.frame(
      Package  = pkgs,
      Status   = status,
      Type = type,
      Version  = version,
      stringsAsFactors = FALSE
    )
    
    if(!is.null(outfile)){
      write.table(df, file=outfile, sep=sep, row.names = FALSE)
    }
    invisible(df)
}
  

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
#' @keywords internal
#' @noRd

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
#' @keywords internal

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
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
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
get_stat <- function(read_count, stat_df=NULL, stage="", params=NA, outfile=NULL, sep=","){
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

#' Identify the version of a third-party program
#'
#' This function runs an external program with its version option and extracts
#' the version number from the returned message.
#'
#' @param name Character string giving the name of the argument used to define
#'   the program path. This is used to handle programs with different version
#'   options (e.g. \code{blast_path} uses \code{-version}).
#' @param path Character string giving the path to the third-party program
#'   executable.
#'
#' @return A character string containing the detected program version.
#'
#' @examples
#' \dontrun{
#' version <- get_program_version(
#'   name = "vsearch_path",
#'   path = "/usr/local/bin/vsearch"
#' )
#' }
#'
#' @export
#'
get_program_version <- function(name, path) {
  
  path <- path.expand(path)
  
  if(!endsWith(name, "_path")){
    return(NA_character_)
  }
  
  if (name == "blast_path") {
    out <- system2(
      path,
      "-version",
      stdout = TRUE,
      stderr = TRUE
    )
  } else {
    out <- system2(
      path,
      "--version",
      stdout = TRUE,
      stderr = TRUE
    )
  }
  
  out <- paste(out, collapse = " ")
  
  version <- regmatches(
    out,
    regexpr("[0-9]+\\.[0-9]+[0-9\\.\\+]*", out, perl = TRUE)
  )
  
  return(version)
}

#' Collect function call information
#'
#' Collects the name of the calling function, the names and values of all its
#' arguments (including default values), and the time at which the function was
#' called. The information is returned as a data frame and can subsequently be
#' written to a log file with `complete_log()`.
#'
#' @return A data frame containing the calling function name, argument names,
#'   argument values, and the timestamp. Returned invisibly.
#'
#'
#' @keywords internal
collect_log <- function(file = NULL) {
  
  if(is.na(file)){ # if user defines NA => do not write log
    return(invisible(NA))
  }
  ## Name of calling function
  fun_name <- as.character(sys.call(-1)[[1]])
  
  ## Calling environment
  env <- parent.frame()
  
  ## Function definition and call
  fun <- sys.function(-1)
  call <- match.call(
    definition = fun,
    call = sys.call(-1),
    expand.dots = FALSE
  )
  
  ## Formal arguments
  fmls <- names(formals(fun))
  
  ## Get argument values (except ...)
  arg_names <- setdiff(fmls, "...")
  log_values <- mget(arg_names, envir = env, inherits = FALSE)
  
  ## Default names
  log_names <- arg_names
  
  ## Replace explicitly supplied arguments by their expressions
  call_args <- as.list(call)[-1]
  
  for (x in intersect(names(call_args), arg_names)) {
    log_names[log_names == x] <- deparse1(call_args[[x]])
  }
  
  
  ## Deal with ...
  if ("..." %in% fmls && "..." %in% names(call_args)) {
    
    dot_values <- eval(
      quote(list(...)),
      envir = env
    )
    
    dot_exprs <- call_args[["..."]]
    
    if (length(dot_values) > 0) {
      
      log_values <- c(log_values, dot_values)
      
      log_names <- c(
        log_names,
        vapply(dot_exprs, deparse1, character(1))
      )
    }
  }
  
  
  time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  
  to_string <- function(x, obj_name) {
    
    if (is.null(x)) {
      "NULL"
      
    } else if (length(x) == 0) {
      "<empty>"
      
    } else if (is.atomic(x)) {
      paste(as.character(x), collapse = ", ")
      
    } else {
      sprintf("%s <%s>", obj_name, class(x)[1])
    }
  }
  
  
  ## Safety check
  stopifnot(
    length(log_values) == length(log_names)
  )
  
  
  log_local <- data.frame(
    function_name = rep(fun_name, length(log_values)),
    argument_name = names(log_values),
    value = mapply(
      to_string,
      log_values,
      log_names,
      SIMPLIFY = TRUE
    ),
    start_time = rep(time, length(log_values)),
    end_time = NA_character_,
    runtime = NA_character_,
    version = NA_character_,
    stringsAsFactors = FALSE
  )
  
  log_local$version <- mapply(
    get_program_version,
    log_local$argument_name,
    log_local$value
  )
  
  return(log_local)
}

#' Complete the log and write it to a CSV file
#'
#' Completes the log data frame returned by `collect_log()` by adding the end
#' time and execution time, then appends the log to a CSV file.
#'
#' @param log A data frame returned by `collect_log()`, containing the calling
#'   function name, argument names, argument values, and the start time.
#' @param file Character string specifying the path to the CSV log file. If
#'   `NULL`, no log file is written.
#' @param sep Character string specifying the field separator used in input and 
#'   output CSV files.
#'
#' @return The completed log data frame, returned invisibly.
#'
#' @keywords internal
 
write_log <- function(log, file = NULL, sep=",") {
  
  if(is.na(file)){ # user do not want log
    return(invisible(NA))
  }
  
#  if (is.null(file))
#    return(invisible(log))
  
  check_dir(file, is_file = TRUE)
  
  end_time <- Sys.time()
  
  log$runtime <- as.numeric(
    difftime(
      end_time,
      as.POSIXct(log$start_time, format = "%Y-%m-%d %H:%M:%S"),
      units = "secs"
    )
  )
  log$runtime <- round(log$runtime, 0)
  
  log$end_time <- format(end_time, "%Y-%m-%d %H:%M:%S")
  
  write.table(
    log,
    file = file,
    append = file.exists(file),
    sep = sep,
    row.names = FALSE,
    col.names = !file.exists(file),
    quote = TRUE
  )
  
  invisible(log)
}

#' Resolve the path to the log file
#'
#' Internal helper that determines which log file path should be used,
#' following the priority order: an explicitly supplied `file` argument,
#' then the `vtamR.log_file` package option, then a default
#' `"vtamR_log.csv"` in the current working directory. If the resolved
#' path is not `NA`, the function ensures that the parent directory of
#' the log file exists (creating it if necessary).
#'
#' @param file Character string specifying the path to the CSV log file,
#'   `NA` to disable logging, or `NULL` (the default) to resolve the path
#'   automatically from the `vtamR.log_file` option or the default
#'   `"vtamR_log.csv"`.
#'
#' @return A character string with the resolved log file path, or `NA` if
#'   logging is disabled.
#'
#' @keywords internal
get_log_file <- function(file = NULL){
  
  if(!is.null(file)){ # user defined something (NA => no log file, or direct path)
    file = file
  } else {
    package_log <- getOption("vtamR.log_file")
    if(!is.null(package_log)){ # package log defined
      file = package_log
    }else{
      file = "vtamR_log.csv"
    }
  }
  
  if(!is.na(file)){ # user did not deacivate log
    check_dir(file, is_file = TRUE) # make dir if do not exists
  }
  return(file)
}


#' Convert a table to a plain text representation
#'
#' Internal helper function that converts a data frame or table-like object into
#' a single character string. The table is printed without row names and the
#' resulting lines are collapsed using newline separators.
#'
#' This function is mainly used to embed tabular information into text-based
#' MIEM reporting fields.
#'
#' @param x A table-like object that can be printed, typically a data frame.
#'
#' @return A single character string containing the printed representation of
#'   the table.
#'
#' @keywords internal
#' @noRd
table_to_text <- function(x){
  txt <- paste(
    capture.output(base::print.data.frame(x, row.names = FALSE)),
    collapse = "\n"
  )
  txt
}

#' Update a MIEM field with information
#'
#' Internal helper function that updates the `Information` column of a MIEM
#' table for a specified reporting field. The function searches for the
#' requested field in the `Step` column and inserts the provided text into the
#' corresponding row.
#'
#' An error is returned if the requested MIEM field does not exist in the table.
#'
#' @param miem A data frame containing MIEM fields. It must contain `Step` and
#'   `Information` columns.
#' @param field A character string specifying the MIEM field to update. The
#'   value must match an entry in the `Step` column.
#' @param text A character string containing the information to add to the
#'   selected MIEM field.
#'
#' @return The updated MIEM data frame.
#'
#' @keywords internal
#' @noRd

set_miem_field <- function(miem, field, text) {
  
  i <- match(field, miem$Step)
  
  if (is.na(i))
    stop("Unknown MIEM field: ", field)
  
  miem$Information[i] <- text
  
  return(miem)
  
}


#' Generate a MIEM bioinformatics report
#'
#' Generate a MIEM (Minimum Information for an eDNA Metabarcoding Study)
#' reporting table from a vtamR log file.
#'
#' The function extracts software versions, pipeline parameters and
#' taxonomic assignment settings from a vtamR log and summarizes them
#' in a table suitable for inclusion in the Supplementary Material of
#' a manuscript.
#'
#' @param log_file Character string specifying the path to the CSV log file.
#'   The path is resolved with the following priority:
#'   \enumerate{
#'     \item the `log_file` argument, if explicitly provided by the user;
#'     \item the package-level option/variable storing a default log path
#'       (if set);
#'     \item `"vtamR_log.csv"` in the current working directory, used as a
#'       last resort if neither of the above is set.
#'   }
#' @param r_versions A data.frame containing installed R package versions
#'   or the path to a CSV file.
#' @param outfile Optional output CSV file.
#' @param sep Field separator used when reading/writing CSV files.
#'
#' @return
#' A data.frame with three columns:
#'
#' * Step
#' * Information
#' * Reporting Requirements
#'
#'
#' @export
#' 
miem_bioinformatics <- function(
  log_file = NULL,
  r_versions,
  outfile = "",
  sep = ","){
  
  log_file <- get_log_file(file = log_file)
  if(is.na(log_file)){
    msg <- paste0("ERROR: `log_file` cannot be NA. Please, provide the path by the `log_file` argument")
    stop(msg)
  }
  
  # read input 
  log_df <-read_input(log_file, sep = sep)
  r_df <- read_input(r_versions, sep = sep)
  
  # initialize miem
  miem <- initialize_miem_bioinformatics()
  
  # get software and package versions
  versions <- extract_versions(log_df, r_df) 
  vtamR_version <- unname(versions["vtamR"])
  blast_version <- unname(versions["blast_path"])
  cutadapt_version <- unname(versions["cutadapt_path"])
  pigz_version <- unname(versions["pigz_path"])
  swarm_version <- unname(versions["swarm_path"])
  vsearch_version <- unname(versions["vsearch_path"])
  rRDP_version <- unname(versions["rRDP"])
  rRDPData_version <- unname(versions["rRDPData"])

  # LTG --------------------------------------------------
  
  
  tmp <- get_log_entries(
    log_df,
    functions = c("assign_taxonomy_ltg"),
    arguments = NULL,
    remove_duplicates = TRUE) 
  
  if(nrow(tmp > 0)){
    
    db <- basename(tmp[tmp$argument_name == "blast_db", "value"])
    tax <- basename(tmp[tmp$argument_name == "taxonomy", "value"])
    ltg_par <- tmp[tmp$argument_name == "ltg_params", "value"]
    
    
    # Taxonomic assignment method
    msg <- paste0(
      "Taxonomic assignment was performed using the Lowest Taxonomic Group (LTG) ",
      "method implemented in the vtamR R package (v", vtamR_version, ") ",
      "[Meglécz, 2023](https://link.springer.com/article/10.1007/s42977-024-00201-x). ",
      "The LTG method is a BLAST-based Lowest Common Ancestor ",
      "(LCA) approach relying on NCBI-BLAST (v", blast_version, ") similarity searches."
    )
    miem <- set_miem_field(miem, field="Taxonomic assignment method", msg)
    
    # Taxonomic assignment parameters (and thresholds)
    if(ltg_par == "NULL"){
      msg <- paste0(
        "Default parameters of the assign_taxonomy_ltg function from the vtamR ",
        "package (v", vtamR_version, ") were used for taxonomic assignment."
      )
    }else{
      msg <- "USER INPUT REQUIRED: Report the custom ltg_params values used for taxonomic assignment."
    }
    miem <- set_miem_field(miem, field="Taxonomic assignment parameters (and thresholds)", msg)
    
    # Database creation: Source of sequences and steps to identify locus of interest
    msg <- paste0(
      "The ", db, " database was used with the ", tax,
      " taxonomy file. ",
      "USER INPUT REQUIRED: Report the database origin and any sequence curation ",
      "performed. If COInr was used, cite that it contains COI sequences compiled ",
      "from the NCBI nucleotide (nt) and BOLD databases ",
      "[Meglécz, 2023](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756)."
    )
    miem <- set_miem_field(miem, field="Database creation: Source of sequences and steps to identify locus of interest", msg)
    
    # Database creation: Link to database or repository
    msg <- paste0(
      "USER INPUT REQUIRED: Report the exact database version and corresponding DOI or URL. ",
      "If COInr was used, cite the appropriate Zenodo release ",
      "(https://zenodo.org/records/20020232) corresponding to the version downloaded ",
      "with the vtamR (v", vtamR_version, ") download_osf function."
    )
    miem <- set_miem_field(miem, field="Database creation: Link to database or repository", msg)
    
    # Database creation: Method for sequence curation
    msg <- paste0(
      "USER INPUT REQUIRED: Report the sequence curation procedure applied. ",
      "For COInr, note that sequence redundancy within taxa is reduced, but the ",
      "database is not curated for potential taxonomic mislabeling."
    )
    miem <- set_miem_field(miem, field="Database creation: Method for sequence curation", msg)
  }
  
  # RDP --------------------------------------------------
  
  tmp <- get_log_entries(
    log_df,
    functions = c("assign_taxonomy_rdp"),
    arguments = NULL,
    remove_duplicates = TRUE) 
  
  if(nrow(tmp > 0)){
    
    db <- tmp[tmp$argument_name == "dir", "value"]
    confidence <- tmp[tmp$argument_name == "confidence", "value"]
    chloroplast <- tmp[tmp$argument_name == "rm_chloroplast", "value"]
    
    # Taxonomic assignment method
    msg <- paste0(
      "Sequences were assigned to taxa using the RDP Naive Bayesian classifier ",
      "implemented in the rRDP package (v", rRDP_version,
      "; https://bioconductor.org/packages/3.21/bioc/html/rRDP.html)."
    )
    miem <- set_miem_field(miem, field="Taxonomic assignment method", msg)
    
    # Taxonomic assignment parameters (and thresholds)
    msg <- paste0(
      "Default parameters of the rRDP package were used, retaining only taxonomic ",
      "assignments with bootstrap support greater than or equal to ",
      confidence, "."
    )
    if (chloroplast) {
      msg <- paste0(
        msg,
        " Taxonomic assignments classified as Chloroplast were set to NA."
      )
    }
    miem <- set_miem_field(miem, field="Taxonomic assignment parameters (and thresholds)", msg)
    
    # Database creation: Source of sequences and steps to identify locus of interest
    if (db == "NULL") {
      msg <- paste0(
        "The reference database distributed with the rRDPData package (v",
        rRDPData_version,
        "; https://bioconductor.org/packages/3.21/data/experiment/html/rRDPData.html) ",
        "was used for taxonomic assignment."
      )
    } else {
      msg <- paste0(
        "USER INPUT REQUIRED: The reference database information cannot be ",
        "retrieved from the log file. Report the database name, version, ",
        "source, and any modifications performed before use."
      )
    }
    miem <- set_miem_field(miem, field="Database creation: Source of sequences and steps to identify locus of interest", msg)
    
    # Database creation: Link to database or repository
    if (db == "NULL") {
      msg <- paste0(
        "rRDPData (v", rRDPData_version,
        "; https://bioconductor.org/packages/3.21/data/experiment/html/rRDPData.html)."
      )
    } else {
      msg <- "USER INPUT REQUIRED: Provide the URL or DOI of the reference database used."
    }
    miem <- set_miem_field(miem, field="Database creation: Link to database or repository", msg)
  }
  
  # Primer removal (trimming): program, version, parameters #############################################
  trim_functions <- c("trim_primers", 
                      "demultiplex_and_trim_fasta",
                      "demultiplex_and_trim_fastq")
  params <- c("check_reverse",
              "primer_to_end",
              "cutadapt_error_rate",
              "cutadapt_minimum_length",
              "cutadapt_maximum_length",
              "tag_to_end")
  tmp <- get_log_entries(
    log_df,
    functions = trim_functions,
    arguments = params,
    remove_duplicates = TRUE) 
  
  if(nrow(tmp) > 0){
#    txt <- paste(capture.output(print(tmp,row.names = FALSE)),collapse = "\n")
    txt <- table_to_text(tmp)
    msg <- paste0(
      "Primer and tag removal was performed using Cutadapt (v", cutadapt_version,
      "), integrated within the following vtamR (v", vtamR_version,
      ") function(s). The parameters used were:",
      "\n",
      txt
    )
    miem <- set_miem_field(miem, field="Primer removal (trimming): program, version, parameters", msg)
  }
  
  # QC program: program, version, parameters #############################################
  params <- c("fastq_maxee",
              "fastq_minlen",
              "fastq_maxlen",
              "fastq_minmergelen",
              "fastq_maxmergelen",
              "fastq_maxns",
              "fastq_truncqual",
              "fastq_maxdiffs",
              "fastq_minovlen")
  
  tmp <- get_log_entries(
    log_df,
    functions = "merge_fastq_pairs",
    arguments = params,
    remove_duplicates = TRUE) 
  
  if(nrow(tmp) >0){
    txt <- table_to_text(tmp)
    msg <- paste0(
      "Quality filtering was performed using the fastq_mergepairs function from ",
      "VSEARCH (v", vsearch_version, ") through the merge_fastq_pairs function ",
      "implemented in vtamR (v", vtamR_version, "). The parameters used were:",
      "\n",
      txt
    )
    miem <- set_miem_field(miem, field="QC program: program, version, parameters", msg)
  }
  
  # Read pair merging: program, version, parameters #############################################
  params <- c("fastq_maxdiffs",
              "fastq_minovlen",
              "fastq_minlen",
              "fastq_maxlen",
              "fastq_minmergelen",
              "fastq_maxmergelen",
              "fastq_allowmergestagger")
  
  tmp <- get_log_entries(
    log_df,
    functions = "merge_fastq_pairs",
    arguments = params,
    remove_duplicates = TRUE) 
  if(nrow(tmp) >0){
    txt <- table_to_text(tmp)
    msg <- paste0(
      "Read pair merging was performed using the fastq_mergepairs function of ",
      "VSEARCH (v", vsearch_version, ") called through the merge_fastq_pairs ",
      "function implemented in vtamR (v", vtamR_version, "). The parameters used were:",
      "\n",
      txt
    )
    miem <- set_miem_field(miem, field="Read pair merging: program, version, parameters", msg)
  }
  
  # Chimera removal: program, version, parameters #############################################
  
  params <- c("abskew",
              "by_sample",
              "filter_occurrence",
              "sample_prop")
  
  tmp <- get_log_entries(
    log_df,
    functions = "filter_chimera",
    arguments = params,
    remove_duplicates = TRUE) 
  
  if(nrow(tmp) >0){
    txt <- table_to_text(tmp)
    
    msg <- paste0(
      "Chimera removal was performed using the uchime3_denovo function from ",
      "VSEARCH (v", vsearch_version, ") through the filter_chimera function ",
      "implemented in vtamR (v", vtamR_version, "). The parameters used were:",
      "\n",
      txt
    )
    miem <- set_miem_field(miem, field="Chimera removal: program, version, parameters", msg)
  }
  
  # Clustering: OTUs or ASVs (and thresholds) #############################################
  
  params <- c("method",
              "by_sample",
              "swarm_d",
              "fastidious",
              "identity")
  
  tmp <- get_log_entries(
    log_df,
    functions = "cluster_asv",
    arguments = params,
    remove_duplicates = TRUE) 
  
  if(nrow(tmp) >0){
    cluster_method <- tmp[tmp$argument_name == "method", "value"]
    
    if(cluster_method == "vsearch"){
      tmp <- tmp %>%
        filter(argument_name %in% c("by_sample", "identity"))
      txt <- table_to_text(tmp)
      
      msg <- paste0(
        "Clustering of ASVs into mOTUs was performed using the cluster_size function ",
        "of VSEARCH (v", vsearch_version, ") through the cluster_asv function ",
        "implemented in vtamR (v", vtamR_version, "). The parameters used were:",
        "\n",
        txt
      )
      miem <- set_miem_field(miem, field="Clustering: OTUs or ASVs (and thresholds)", msg)
    }else{
      tmp <- tmp %>%
        filter(argument_name %in% c("by_sample", "swarm_d", "fastidious"))
      txt <- table_to_text(tmp)
      
      msg <- paste0(
        "Clustering of ASVs into mOTUs was performed using swarm (v", 
        swarm_version, ") through the cluster_asv function implemented in vtamR ",
        "(v", vtamR_version, "). The parameters used were:",
        "\n",
        txt
      )
      miem <- set_miem_field(miem, field="Clustering: OTUs or ASVs (and thresholds)", msg)
    }
  }
  
  
  # Additional filtering: removal of singletons or other methods #############################################
  
  functions <- c("denoise_by_swarm",
                 "filter_asv_global",
                 "filter_stop_codon",
                 "filter_indel",
                 "filter_min_replicate",
                 "filter_replicate")
  
  delete_params <- c("read_count",
                     "swarm_path",
                     "vsearch_path",
                     "pigz_path",
                     "num_threads",
                     "outfile",
                     "sep",
                     "quiet",
                     "log_file",
                     "sampleinfo",
                     "conta_file",
                     "mock_composition",
                     "outdir",
                     "known_occurrences")
  
  tmp <- get_log_entries(
    log_df,
    functions = functions,
    arguments = delete_params,
    keep_arguments = FALSE,
    remove_duplicates = TRUE) 
  
  if(nrow(tmp) > 0){
    
    txt <- table_to_text(tmp)
    
    msg <- paste0(
      "Additional filtering was performed using the following vtamR functions ",
      "(v", vtamR_version, ") with the parameters reported below. ",
      "The denoise_by_swarm function uses swarm (v", swarm_version, "). ",
      "For detailed information on the order of filtering steps, refer to the ",
      "log file provided in the Supplementary Materials.",
      "\n",
      txt
    )
    miem <- set_miem_field(miem, field="Additional filtering: removal of singletons or other methods", msg)
  }
  
  
  # Additional filtering: decontamination using sequenced controls #############################################
  
  functions <- c("filter_contaminant",
                 "filter_pcr_error",
                 "filter_occurrence_read_count",
                 "filter_occurrence_sample",
                 "filter_occurrence_variant",
                 "suggest_pcr_error_cutoff",
                 "suggest_variant_readcount_cutoffs",
                 "suggest_sample_cutoff")
  
  delete_params <- c("read_count",
                     "swarm_path",
                     "vsearch_path",
                     "pigz_path",
                     "num_threads",
                     "outfile",
                     "sep",
                     "quiet",
                     "log_file",
                     "sampleinfo",
                     "conta_file",
                     "mock_composition",
                     "outdir",
                     "known_occurrences")
  
  tmp <- get_log_entries(
    log_df,
    functions = functions,
    arguments = delete_params,
    keep_arguments = FALSE,
    remove_duplicates = TRUE) 
  
  if(nrow(tmp) >0){
    
    txt <- table_to_text(tmp)
    
    msg <- paste0(
      "The parameters used for the following 'filter_xxx' functions were selected ",
      "based on the output of the corresponding 'suggest_xxx' functions. ",
      "These functions estimate optimal filtering parameters from the composition ",
      "of control samples in order to minimize false positives and false negatives. ",
      "All functions are implemented in vtamR (v", vtamR_version, "). ",
      "The filter_pcr_error function uses the usearch_global function from ",
      "VSEARCH (v", vsearch_version, "). ",
      "For detailed information on the order of filtering steps, refer to the ",
      "log file provided in the Supplementary Materials.",
      "\n",
      txt
    )
    miem <- set_miem_field(miem, field="Additional filtering: decontamination using sequenced controls", msg)
  }
  
  #### write csv
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(miem, file = outfile,  row.names = F, sep=sep)
  }
  invisible(miem)
}


#' Initialize the default MIEM bioinformatics reporting table
#'
#' Creates and returns a data frame containing the MIEM (Minimum Information
#' for an eDNA Metabarcoding study) bioinformatics reporting fields and their
#' corresponding reporting requirements. The `Information` column is initialized
#' with `NA` values and is intended to be populated later in the workflow.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{Step}{The MIEM bioinformatics reporting item.}
#'   \item{Information}{A placeholder (`NA_character_`) for user-supplied information.}
#'   \item{Reporting Requirements}{Whether the item must be reported or is only
#'   required if applicable.}
#' }
#'
#' @keywords internal
#' @noRd

initialize_miem_bioinformatics <- function(){
  miem_fields <- c(
    "Database creation: Source of sequences and steps to identify locus of interest",
    "Database creation: Method for sequence curation",
    "Database creation: Link to database or repository",
    "Primer removal (trimming): program, version, parameters",
    "QC program: program, version, parameters",
    "Read pair merging: program, version, parameters",
    "Chimera removal: program, version, parameters",
    "Clustering: OTUs or ASVs (and thresholds)",
    "Additional filtering: removal of singletons or other methods",
    "Additional filtering: decontamination using sequenced controls",
    "Taxonomic assignment method",
    "Taxonomic assignment parameters (and thresholds)",
    "Read normalization: methods"
  )
  
  mien_requirements <- c(
    rep("Report", 8),
    "If Applicable",
    "If Applicable",
    "Report",
    "Report",
    "If Applicable"
  )
  miem <- data.frame(
    "Step" = miem_fields,
    "Information" = NA_character_,
    "Reporting Requirements" = mien_requirements,
    stringsAsFactors = FALSE
    
  )
  return(miem)
  
}


#' Read input data from a file path or an existing data frame
#'
#' Internal helper function that accepts either a path to a delimited text file
#' or an already loaded data object. If `x` is a character string, the function
#' reads the file using [utils::read.table()] with the specified separator.
#' Otherwise, the input object is returned unchanged.
#'
#' @param x A character string giving the path to a CSV file, or an existing
#'   data object (typically a data frame).
#' @param sep A single character used to separate fields in the input file.
#'   Defaults to `","`.
#'
#' @return A data frame or data object containing the input data.
#'
#' @keywords internal
read_input <- function(x, sep = ",", header=TRUE) {
  
  if(is.character(x)){
    df <- read.table(x, header=header, sep=sep, stringsAsFactors = FALSE)
  }else{
    df <- x
  }
  return(df)
}

#' Extract software versions from a workflow log
#'
#' Internal helper function that extracts version information for software
#' recorded in a workflow log and combines it with versions of R packages.
#' Third-party software versions are identified from the log data, while R
#' package versions are obtained from the provided package information table.
#'
#' If multiple versions of the same third-party program are detected, a warning
#' is issued and only the last recorded version is retained.
#'
#' @param log_df A data frame containing workflow log information. It must
#'   contain `argument_name` and `version` columns for third-party software
#'   extraction.
#' @param r_df A data frame containing R package information. It must contain
#'   `Package` and `Version` columns.
#'
#' @return A named character vector containing software names as names and
#'   corresponding version numbers as values. Both third-party software and
#'   R package versions are included.
#'
#' @keywords internal
#' @noRd

extract_versions <- function(log_df, r_df) {
  
  # named vector
  r_versions <- setNames(
    r_df$Version,
    r_df$Package
  )
  
  third_party <- log_df %>%
    select(argument_name, version) %>%
    filter(!is.na(version)) %>%  # get only prorgam paths
    distinct()
  
  if (anyDuplicated(third_party$argument_name)) {
    warning(
      "Some third-party programs were used with multiple versions. ",
      "Only the last called version will be reported.",
      call. = FALSE
    )
    
    third_party <- third_party %>%
      group_by(argument_name) %>%
      summarise(
        version = last(version),
        .groups = "drop"
      )
  }
  
  tp_versions <-
    setNames(
      third_party$version,
      third_party$argument_name
    )
  
  versions <- c(tp_versions, r_versions)
  as.character(versions)
  
  return(versions)
}

#' Extract selected entries from a workflow log
#'
#' Internal helper function that extracts log entries corresponding to selected
#' functions and, optionally, selected arguments. The returned data frame
#' contains the function name, argument name, and associated value recorded in
#' the log.
#'
#' Arguments can either be retained (`keep_arguments = TRUE`) or excluded
#' (`keep_arguments = FALSE`). Duplicate entries can optionally be removed.
#'
#' @param log_df A data frame containing workflow log information. It must
#'   contain `function_name`, `argument_name`, and `value` columns.
#' @param functions A character vector of function names to extract from the
#'   log.
#' @param arguments Optional. A character vector of argument names to retain or
#'   exclude, depending on the value of `keep_arguments`.
#' @param keep_arguments Logical. If `TRUE` (default), keep only entries whose
#'   argument names are included in `arguments`. If `FALSE`, remove entries
#'   whose argument names are included in `arguments`.
#' @param remove_duplicates Logical. If `TRUE` (default), remove duplicated entries from
#'   the output.
#'
#' @return A data frame containing the extracted log entries with columns:
#'   `function_name`, `argument_name`, and `value`.
#'
#' @keywords internal

get_log_entries <- function(
  log_df,
  functions,
  arguments = NULL,
  keep_arguments = TRUE,
  remove_duplicates = TRUE) {
  
  x <- log_df %>%
    filter(function_name %in% functions)
  
  if (!is.null(arguments)){
    if(keep_arguments){
      x <- x %>%
        filter(argument_name %in% arguments)
    } else{
      x <- x %>%
        filter(!(argument_name %in% arguments))
    }
  }
  
  x <-  x %>%
    select(
      function_name,
      argument_name,
      value
    )
  
  if (remove_duplicates){
    x <- distinct(x)
  }
  
  return(x)
}

#' Generate summary files for the MIEM sequencing summary statistics
#'
#' Generates summary statistics and quality-control metrics from the
#' intermediate and final outputs of the vtmaR pipeline. The resulting
#' files can be used to complete the "Results - Sequencing Summary
#' Statistics" section of the MIEM checklist.
#'
#' Depending on the input files provided, the function summarizes
#' sequencing read counts, read counts per sample and sample type,
#' taxonomic assignments, and results from mock-community controls.
#' All generated results are written to files in `outdir`.
#'
#' @param info_files Character vector of information files
#'   (`fastqinfo`, `fastainfo`, or `sampleinfo`) containing information
#'   such as FASTQ/FASTA file names, sample names, and sample types
#'   (`real`, `mock`, or `negative`). These files are generated at
#'   different preprocessing steps. In this function they are used to 
#'   calculate the total number of reads after the preprocessing steps of vtamR.
#'   The last one is also used to determine sample information.
#' @param fastq_dir Character string specifying the directory containing
#'   FASTQ files. Used to calculate read counts when they are not already
#'   available in `info_files`.
#' @param read_count_files Character vector of read-count files generated
#'   during the pipeline. Each file contains `asv_id`, `sample`,
#'   `read_count`, and optionally `replicate` and `cluster_id`.
#'   These files are used to summarize the number of reads retained
#'   during preprocessing and filtering.
#'
#'   A read-count file from any filtering step can be provided. The two
#'   most informative files are usually the first file (typically the
#'   output of `dereplicate`) and the last file, as they provide the
#'   number of reads and ASVs after the preprocessing (merge, demultiplex and quality
#'   filter) and after vtamR filtering (e.g. chimera, low, frequency noise), respectively.
#'
#'   The last file in the vector is also used to count the number of
#'   ASVs or mOTUs assigned to each taxonomic rank and to identify
#'   false-positive and false-negative assignments.
#'
#'   If the last file is the output of `cluster_asv` and contains a
#'   `cluster_id` column, the number of mOTUs rather than ASVs is used
#'   for the taxonomic summaries.
#' @param taxa Data frame containing taxonomic assignments. It must contain
#'   `asv_id` (or `cluster_id`) as well as the taxonomic ranks 
#'   `domain`, `phylum`, `class`, `order`,
#'   `family`, `genus`, and `species` and may contain `ltg_rank_index`.
#'   It is used to summarize the number
#'   of ASVs or mOTUs assigned at each taxonomic rank.
#' @param outdir Character string specifying the directory where result
#'   files will be written. The directory is created if necessary.
#'   Defaults to the current working directory.
#' @param mock_composition Mock-community composition used to evaluate
#'   control samples and identify false-positive and false-negative
#'   assignments.
#' @param sep Character string used as the field separator for input
#'   and output files. Defaults to `","`.
#'
#' @details
#' The function generates different result files depending on the
#' arguments supplied.
#'
#' If `info_files` is provided, the total number of reads associated
#' with each information file is calculated using
#' [count_reads_from_info()]. Existing read counts are used when
#' available; otherwise, reads are counted directly from the FASTQ files.
#' The results are written to `preprocess_read_count.csv`.
#'
#' If both `info_files` and `read_count_files` are provided, read counts
#' are summarized for each read-count file and separately for real, mock,
#' and negative samples. The results are
#' written to `read_count_by_sample.csv`.
#'
#' If both `taxa` and `read_count_files` are provided, the number of ASVs
#' or mOTUs assigned at each major taxonomic rank is calculated using
#' [count_taxassing_by_rank()] from the last file in `read_count_files`. 
#' The results are written to
#' `ASV_or_mOTU_count_by_taxonomic_rank.csv`.
#'
#' If `info_files`, `read_count_files` and `mock_composition` are
#' provided, mock-community control results are evaluated using
#' [classify_control_occurrences()] from the last read_count_files. 
#' False-positive and false-negative
#' occurrences are combined and written to
#' `false_positives_and_negatives.csv`.
#'
#' @return Invisibly returns `NULL`. Results are written as output files
#'   to `outdir`.
#'
#' @examples
#' \dontrun{
#' miem_results(
#'   info_files = info_files,
#'   fastq_dir = "fastq",
#'   read_count_files = read_count_files,
#'   taxa = taxa,
#'   outdir = "miem",
#'   mock_composition = mock_composition
#' )
#' }
#'
#' @export

miem_results <- function(
  info_files = NULL, 
  fastq_dir = ".", 
  read_count_files = NULL, 
  taxa = NULL, 
  outdir = ".", 
  mock_composition = NULL, 
  sep=","){
  
  outdir <- check_dir(outdir)
  
  ##########################################################
  ### Total number of raw sequence reads produced
  ### Total number of reads assigned to indices
  ### Total number of reads that made it through bioinformatic filtering
  
  # for each info_files add the read_count if exists, otherwise count reads in fasta/fatsq files
  # There should be only one info file where the read_count in not already present (the first fastq used in the pipeline). For this file
  # use de fastq_dir to access the fastq files liste in the info_file
  
  if(!is.null(info_files)){
    count_reads_from_info(info_files = info_files, fastq_dir = fastq_dir, outdir = outdir, sep = sep)
  }
  
  ##########################################################
  ### Total number of reads used for final/ subsequent analyses
  ### Average number of reads per sample
  ### Minimum and maximum number of reads per sample
  # for each read_count_files count the 
  # - Total number of reads
  # - for each sample_type (real/mock/negative)
  #    - Min, max, mean, median, number of samples
  # use the last info file to get sample_types
  
  if(!is.null(info_files) && !is.null(read_count_files)){
    read_count_by_sample(info_files = info_files, read_count_files = read_count_files, outdir = outdir, sep = sep)
  }
  
  ### Total number of OTUs or ASVs assigned to taxa (and to what level of taxonomy)
  ### Number of OTUs or ASVs unassigned
  # from taxa, final_read_count, count the number of ASV assigned to each major taxonomic level
  
  if(!is.null(taxa) && !is.null(read_count_files)){
    
    final_read_count <- read_count_files[length(read_count_files)]
    
    taxa_df <- format_taxa(taxa, sep=sep)
    count_by_taxrank <- count_taxassing_by_rank(final_read_count, taxa_df, sep = sep)
    
    outfile = file.path(outdir, "ASV_or_mOTU_count_by_taxonomomic_rank.csv")
    write.table(count_by_taxrank, file = outfile, row.names= FALSE, sep = sep)
  }
  
  ### Results from Controls
  # Results are in performance_metrics.csv, and details in known_occurrences.csv and false_negatives.csv
  # if mock_composition is given, run classify_control_occurrences using final_read_count as read_count 
  # and last info file to get sample_types (sampleinfo)
  
  if(!is.null(info_files) && !is.null(read_count_files) && !is.null(mock_composition)){
    
    final_read_count <- read_count_files[length(read_count_files)]
    
    sampleinfo <- read_input(info_files[length(info_files)], sep = sep)
    results <- classify_control_occurrences(
      read_count = final_read_count, 
      sampleinfo = sampleinfo, 
      mock_composition = mock_composition, 
      sep = sep,
      log_file = NA)
    known_occurrences <- results[[1]]
    false_neagtives <- results[[2]]
    performance <- results[[3]]
    
    fp <- known_occurrences %>%
      filter(action == "delete") %>%
      mutate(occurrence_type = "FP") %>%
      select(occurrence_type, sample, asv_id, asv)
    
    if(!"asv_id" %in% colnames(false_neagtives)){
      false_neagtives <- false_neagtives %>%
        mutate(asv_id = NA)
    }
    false_neagtives <- false_neagtives %>%
      mutate(occurrence_type = "FN") %>%
      select(occurrence_type, sample, asv_id, asv)
    
    tmp <- rbind(fp, false_neagtives)
    
    outfile <- file.path(outdir, "false_positives_and_nevatives.csv")
    write.table(tmp, file = outfile, sep = sep, row.names = FALSE)
  }
  
  
}

#' Summarize read counts by sample type
#'
#' Computes per-file read-count and ASV-count summaries, including
#' minimum, maximum, mean, and median reads per sample for real, mock,
#' and negative-control samples.
#'
#' The sample type information is taken from the last file in
#' `info_files`. Results are written to `read_count_by_sample.csv`
#' in `outdir`.
#'
#' @param read_count_files Character vector of read-count files to summarize.
#' @param info_files Character vector of sample information files.
#'   The last file is used to determine sample types, thus it should contain
#'   `sample` and `sample_type` (real/mock/negative) columns.
#' @param outdir Character string specifying the output directory.
#'   Defaults to the current working directory.
#' @param sep Character string used as the field separator for input
#'   and output files. Defaults to `","`.
#'
#' @return A data frame containing read-count and ASV-count summaries
#'   for each input read-count file, and minimum, maximum, mean, median of the
#'   number of reads per sample for each sample type. 
#'    The same data are written to `read_count_by_sample.csv` in `outdir`.
#'
#' @details
#' For each read-count file, reads are first summed across ASVs for
#' each sample. These per-sample read counts are then summarized
#' separately for real, mock, and negative samples.
#'
#' The output also contains the total number of reads and the total
#' number of unique ASVs in each read-count file.
#'
#' @keywords internal
read_count_by_sample  <- function(read_count_files, info_files, outdir = ".", sep = ","){
  
  sampleinfo <- info_files[length(info_files)]
  sampleinfo_df <- read_input(sampleinfo, sep = sep) %>%
    select(sample_type, sample) %>%
    distinct()
  
  read_count_sample <- data.frame(
    read_count_file = as.character(),
    total_read_count = numeric(),
    total_asv_count = numeric(),
    # real
    min_reads_per_sample_real = numeric(),
    max_reads_per_sample_real = numeric(),
    mean_reads_per_sample_real = numeric(),
    median_reads_per_sample_real = numeric(),
    n_samples_real = numeric(),
    # mock
    min_reads_per_sample_mock = numeric(),
    max_reads_per_sample_mock = numeric(),
    mean_reads_per_sample_mock = numeric(),
    median_reads_per_sample_mock = numeric(),
    n_samples_mock = numeric(),
    # negative
    min_reads_per_sample_negative = numeric(),
    max_reads_per_sample_negative = numeric(),
    mean_reads_per_sample_negative = numeric(),
    median_reads_per_sample_negative = numeric(),
    n_samples_negative = numeric()
  )
  
  for(file in read_count_files){
    read_count_df <- read_input(file, sep = sep)
    total = sum(read_count_df$read_count)
    asv_count = length(unique(read_count_df$asv))
    
    read_count_df <- read_count_df %>%
      group_by(sample) %>%
      summarise(read_count = sum(read_count), .groups = "drop") %>%
      left_join(sampleinfo_df, by = "sample") %>%
      group_by(sample_type) %>%
      summarize(
        min_reads_per_sample = min(read_count),
        max_reads_per_sample = max(read_count),
        mean_reads_per_sample = round(mean(read_count), digits = 0),
        median_reads_per_sample = round(median(read_count), digits = 0),
        n_samples = n()
      )
    if(!"mock" %in% read_count_df$sample_type){
      read_count_df <- read_count_df %>%
        add_row(sample_type = "mock", min_reads_per_sample = NA, max_reads_per_sample=NA, mean_reads_per_sample=NA, median_reads_per_sample=NA,  n_samples = 0)
    }
    if(!"negative" %in% read_count_df$sample_type){
      read_count_df <- read_count_df %>%
        add_row(sample_type = "negative", min_reads_per_sample = NA, max_reads_per_sample=NA, mean_reads_per_sample=NA, median_reads_per_sample=NA,  n_samples = 0)
    }
    if(!"real" %in% read_count_df$sample_type){
      read_count_df <- read_count_df %>%
        add_row(sample_type = "real", min_reads_per_sample = NA, max_reads_per_sample=NA, mean_reads_per_sample=NA, median_reads_per_sample=NA,  n_samples = 0)
    }
    read_count_sample <- read_count_sample %>%
      add_row(
        read_count_file = file,
        total_read_count = total,
        total_asv_count = asv_count,
        min_reads_per_sample_real = read_count_df$min_reads_per_sample[read_count_df$sample_type == "real"],
        max_reads_per_sample_real = read_count_df$max_reads_per_sample[read_count_df$sample_type == "real"],
        mean_reads_per_sample_real = read_count_df$mean_reads_per_sample[read_count_df$sample_type == "real"],
        median_reads_per_sample_real = read_count_df$median_reads_per_sample[read_count_df$sample_type == "real"],
        n_samples_real = read_count_df$n_samples[read_count_df$sample_type == "real"],
        # mock
        min_reads_per_sample_mock = read_count_df$min_reads_per_sample[read_count_df$sample_type == "mock"],
        max_reads_per_sample_mock = read_count_df$max_reads_per_sample[read_count_df$sample_type == "mock"],
        mean_reads_per_sample_mock = read_count_df$mean_reads_per_sample[read_count_df$sample_type == "mock"],
        median_reads_per_sample_mock = read_count_df$median_reads_per_sample[read_count_df$sample_type == "mock"],
        n_samples_mock = read_count_df$n_samples[read_count_df$sample_type == "mock"],
        # negative
        min_reads_per_sample_negative = read_count_df$min_reads_per_sample[read_count_df$sample_type == "negative"],
        max_reads_per_sample_negative = read_count_df$max_reads_per_sample[read_count_df$sample_type == "negative"],
        mean_reads_per_sample_negative = read_count_df$mean_reads_per_sample[read_count_df$sample_type == "negative"],
        median_reads_per_sample_negative = read_count_df$median_reads_per_sample[read_count_df$sample_type == "negative"],
        n_samples_negative = read_count_df$n_samples[read_count_df$sample_type == "negative"]
      )
  }
  
  outfile = file.path(outdir, "read_count_by_sample.csv")
  write.table(read_count_sample, file = outfile, row.names= FALSE, sep = sep)
  invisible(read_count_sample)
}

#' Count reads from sample information files
#'
#' Calculates the total number of reads associated with each sample
#' information file. If an input file already contains a `read_count`
#' column, the existing counts are summed. Otherwise, reads are counted
#' directly from the FASTQ files listed in the `fastq_fw` column.
#'
#' The resulting read counts are written to
#' `preprocess_read_count.csv` in `outdir`.
#'
#' @param info_files Character vector of sample information files.
#' @param fastq_dir Character string specifying the directory containing
#'   FASTQ files. Defaults to the current working directory.
#' @param outdir Character string specifying the output directory.
#'   Defaults to the current working directory.
#' @param sep Character string used as the field separator when reading
#'   and writing files. Defaults to `","`.
#'
#' @return A data frame with one row per information file and two columns:
#'   `info_filename`, containing the input file name, and `read_count`,
#'   containing the total number of reads. The same data frame is written
#'   to `preprocess_read_count.csv` in `outdir`.
#'
#' @details
#' For information files containing a `read_count` column, the function
#' retains the last two columns, removes duplicate rows, and sums the
#' read counts.
#'
#' For information files without a `read_count` column, the function
#' extracts unique FASTQ file names from the `fastq_fw` column and counts
#' the reads in each FASTQ file using [count_reads()].
#'
#' @keywords internal

count_reads_from_info <- function(info_files = NULL, fastq_dir = ".", outdir = ".", sep = ","){
  preprocess_read_count_df = data.frame(
    info_filename = character(),
    read_count = numeric()
  )
  for(file in info_files){
    info_df <- read_input(file, sep = sep)
    if("read_count" %in% colnames(info_df)){
      info_df <- info_df %>%
        select((ncol(.) - 1):ncol(.)) %>% # keep the last filename column and the read_count
        distinct() # get unique list
      
      total_read_count <- sum(info_df$read_count)
      preprocess_read_count_df <- preprocess_read_count_df %>%
        add_row(info_filename = file, read_count = total_read_count)
    }else{ # reads_ should be counted from fastq files
      
      fastq_files <- read_input(file, sep = sep) %>%
        select(fastq_fw) %>%
        distinct()
      
      total_read_count = 0
      for (fastq in fastq_files$fastq_fw) {
        n <- count_reads(
          file.path(fastq_dir, fastq),
          file_type = "fastq"
        )
        total_read_count <- total_read_count + n
      }
      preprocess_read_count_df <- preprocess_read_count_df %>%
        add_row(info_filename = file, read_count = total_read_count)
    }
  }
  outfile = file.path(outdir, "preprocess_read_count.csv")
  write.table(preprocess_read_count_df, file = outfile, row.names= FALSE, sep = sep)
  invisible(preprocess_read_count_df)
}


#' Format taxonomic assignment results
#'
#' Reads and formats a taxonomic assignment file into a standardized data
#' frame containing ASV identifiers, taxonomic rank indices, and taxonomic
#' levels. Supports taxonomic assignments generated using either LTG
#' (`ltg_rank_index`) or standard taxonomic rank columns (`domain` through
#' `species`).
#'
#' @param taxa Character string giving the path to the taxonomic assignment
#'   file.
#' @param sep Character string used to separate fields in the input file.
#'   Defaults to `","`.
#'
#' @return A data frame containing the ASV identifier (`asv_id`), the
#'   corresponding taxonomic rank index (`rank_index`), and the taxonomic
#'   level (`taxonomic_level`). For LTG assignments, the rank index is
#'   derived from `ltg_rank_index`. For standard taxonomic assignments,
#'   the rank index is determined from the highest available taxonomic
#'   level.
#'
#' @keywords internal
#' 
format_taxa <- function(taxa, sep=","){
  
  taxa_df <- read_input(taxa, sep = sep)
  
  ########### Read and format taxonomy results tog get asv_id
  tax_ind <- data.frame(
    rank_index = c(1,2,3,4,5,6,7,8),
    taxonomic_level = c("root","domain","phylum","class","order","family","genus","species")
  )
  ### LTG
  if("ltg_rank_index" %in% colnames(taxa_df)){
    
    taxa_df <- read_input(taxa, sep = sep) %>%
      select(asv_id, "rank_index" = ltg_rank_index) %>%
      mutate(rank_index = if_else(is.na(rank_index), 1, floor(rank_index))) %>%
      left_join(tax_ind, by="rank_index")
    
  } else {
    
    taxa_df <- read_input(taxa, sep = sep) %>%
      select(-asv) %>%
      mutate(rank_index = 8- rowSums(is.na(select(., domain:species)))) %>%
      left_join(tax_ind, by="rank_index") %>%
      select(asv_id, rank_index, taxonomic_level)
  }
  return(taxa_df)
}

#####################################################################

#' Count ASVs or mOTUs by taxonomic rank
#'
#' Counts the number of distinct ASVs or mOTUs assigned to each taxonomic
#' level. If the input data contains a `cluster_id` column, it is used as
#' the identifier for mOTUs instead of `asv_id`.
#'
#' @param read_count Character string giving the path to the read-count
#'   file containing ASV or mOTU identifiers.
#' @param taxa_df A data frame containing `asv_id`, `rank_index`, and 
#' `taxonomic_level` columns.
#' @param sep Character string used to separate fields in the input files.
#'   Defaults to `","`.
#'
#' @return A data frame containing the number of ASVs or mOTUs assigned to
#'   each taxonomic level, ordered from the highest to the lowest rank.
#'
#' @keywords internal

count_taxassing_by_rank <- function(read_count, taxa_df, sep = ","){
  
  asv_by_rank <- data.frame(
    "taxonomic_level" = character(),
    "ASV_or_mOTU_number" = numeric()
  )
  
  read_count_df <- read_input(read_count, sep = sep)
  
  # if output of cluster is not grouped, replace asv_id column by cluster_id
  if("cluster_id" %in% colnames(read_count_df)){
    read_count_df <- read_count_df %>%
      select(-asv_id) %>%
      rename(asv_id = cluster_id)
  }
  
  asv_by_rank <- read_count_df %>%
    select(asv_id) %>%
    distinct() %>%
    left_join(taxa_df, by="asv_id") %>%
    group_by(taxonomic_level) %>%
    summarize(ASV_or_mOTU_number = n(), rank_index = first(rank_index)) %>%
    arrange(desc(rank_index)) %>%
    select(-rank_index)
  return(asv_by_rank)
}



#' Format OTU and taxonomy tables for vegan
#'
#' Converts OTU/ASV abundance and taxonomy input tables into the format
#' required for downstream analysis with the \pkg{vegan} package. The OTU
#' table is returned with samples as rows and ASVs (or clusters) as columns,
#' while the taxonomy table contains taxonomy information for the ASVs
#' present in the OTU table.
#'
#' When \code{rm_coltrol = TRUE}, only samples classified as \code{"real"} in
#' the \code{sample_type} input are retained, thus control samples are deleted.
#' If a \code{replicate} column is
#' present, sample and replicate identifiers are combined into a single
#' \code{sample} identifier using a hyphen.
#'
#' If a \code{cluster_id} column is present in the OTU table, reads are
#' aggregated by cluster and sample. Otherwise, reads are aggregated by
#' ASV and sample.
#'
#' @param otu Character string giving the path to the OTU/ASV abundance
#'   table, or a data frame.
#'   The table must contain \code{asv_id}, \code{sample}, and
#'   \code{read_count} columns. An optional \code{replicate} or
#'   \code{cluster_id} column may also be present.
#' @param tax Character string giving the path to the taxonomy table, or a data frame. 
#'   The table must contain an
#'   \code{asv_id} column and taxonomy columns corresponding to the expected
#'   ranks.
#' @param rm_coltrol Logical; if \code{TRUE}, remove control samples and
#'   retain only samples whose \code{sample_type} is \code{"real"}.
#'   Defaults to \code{TRUE}.
#' @param sample_type Character string giving the path to the sample_type
#'   table, or a data frame. Required when
#'   \code{rm_coltrol = TRUE}. The table must contain \code{sample} and
#'   \code{sample_type} columns.
#' @param outfile_motu Optional character string specifying the output file
#'   path for the formatted OTU table. If \code{NULL}, the table is not
#'   written to disk.
#' @param outfile_taxa Optional character string specifying the output file
#'   path for the formatted taxonomy table. If \code{NULL}, the table is not
#'   written to disk.
#' @param sep Character used to separate fields in input and output files.
#'   Defaults to \code{","}.
#' @param log_file Character string specifying the path to the CSV log file.
#'   The path is resolved with the following priority:
#'   \enumerate{
#'     \item the `log_file` argument, if explicitly provided by the user;
#'     \item the package-level option/variable storing a default log path
#'       (if set);
#'     \item `"vtamR_log.csv"` in the current working directory, used as a
#'       last resort if neither of the above is set.
#'   }
#'   To disable logging entirely, set `log_file = NA`.
#'
#' @return A list of two data frames:
#' \itemize{
#'   \item \code{otu_table}: A sample-by-ASV (or sample-by-cluster) abundance
#'     table, with sample identifiers as row names and zeroes for missing
#'     abundances.
#'   \item \code{tax_table}: A taxonomy table containing only ASVs present in
#'     \code{otu_table}, with ASV identifiers as row names.
#' }
#'
#' @details
#' The taxonomy table contains the columns \code{domain}, \code{phylum},
#' \code{class}, \code{order}, \code{genus}, and \code{species}. If a
#' \code{kingdom} column is present in the input taxonomy table, it is also
#' retained.
#'
#' If \code{replicate} is present in the OTU table, the resulting sample
#' identifier is constructed as \code{sample-replicate}. The same
#' transformation is applied to the sample-type table when controls are
#' removed.
#'
#' @examples
#' \dontrun{
#' result <- format_for_vegan(
#'   otu = "otu.csv",
#'   tax = "taxonomy.csv",
#'   sample_type = "sample_type.csv"
#' )
#'
#' otu_table <- result[[1]]
#' tax_table <- result[[2]]
#' }
#'
#' @export

format_for_vegan <- function(otu, tax, rm_coltrol = TRUE, sample_type = NULL, 
                             outfile_motu = NULL, outfile_taxa = NULL, 
                             sep = ",", log_file = NULL){
  
  if (missing(otu)) stop("Argument 'otu' is required")
  if (missing(tax)) stop("Argument 'tax' is required")
  
  # resolve path to log_file
  log_file <- get_log_file(file = log_file)
  # get function name, all arguments and stat time. If log_file == NA, no log
  log <- collect_log(file = log_file)
  
  if(rm_coltrol){
    if(is.null(sample_type)){
      msg <- "sample_type must be specified when rm_coltrol is TRUE"
      stop(msg)
    }else{
      sample_type_df <- read_input(sample_type, sep = sep) 
    }
  }
  
  
  # otu ###################
  
  otu_table <- read_input(otu, sep = sep)
  # make one column with sample-replicate
  if("replicate" %in% colnames(otu_table)){
    otu_table <- otu_table  %>%
      mutate(sample = paste(sample, replicate, sep="-")) %>%
      select(-replicate)
    
    sample_type_df <- sample_type_df %>%
      mutate(sample = paste(sample, replicate, sep = "-"))
  }
  
  if("cluster_id" %in% colnames(otu_table)){ # if cluster_id make output with clusters
    otu_table <- otu_table  %>%
      select(-asv_id, -asv) %>%
      group_by(cluster_id, sample) %>%
      summarise(read_count = sum(read_count), .groups = "drop") %>%
      rename(asv_id = cluster_id) 
  } else {
    otu_table <- otu_table  %>%
      select(-asv) %>%
      group_by(asv_id, sample) %>%
      summarise(read_count = sum(read_count), .groups = "drop")
  }
  
  if(rm_coltrol){
    
    sample_type_df <- sample_type_df %>%
      select(sample, sample_type) %>%
      distinct() %>%
      filter(sample_type == "real")
    
    otu_table <- otu_table %>%
      filter(sample %in%  sample_type_df$sample) 
  }
  
  otu_table <- pivot_wider(otu_table, 
                           names_from = asv_id,
                           values_from = read_count,
                           values_fill = 0)
  otu_table <- as.data.frame(otu_table)
  rownames(otu_table) <- otu_table$sample
  otu_table <- select(otu_table, -sample)
  
  
  # taxa ##########################
  
  tax_table <- read_input(tax, sep = sep)
  tax_table <- as.data.frame(tax_table)
  # remove taxa not in otu
  tax_table <- tax_table %>%
    filter(asv_id %in% colnames(otu_table))
  rownames(tax_table) <- tax_table$asv_id
  
  if("kingdom" %in% colnames(tax_table)){
    tax_table <- tax_table %>%
      select(domain, kingdom, phylum, class, order, genus, species)
  } else {
    tax_table <- tax_table %>%
      select(domain, phylum, class, order, genus, species)
  }
  
  if(!is.null(outfile_motu)){
    check_dir(outfile_motu, is_file = TRUE)
    write.table(otu_table, file = outfile_motu, sep = sep, col.names = NA)
  }
  if(!is.null(outfile_taxa)){
    check_dir(outfile_taxa, is_file = TRUE)
    write.table(tax_table, file = outfile_taxa, sep = sep, col.names = NA)
  }
  
  df_list <- list(otu_table, tax_table)
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(df_list)
}


#' Format OTU, taxonomy, and sample data as a phyloseq object
#'
#' Converts OTU/ASV abundance data and taxonomy data into a
#' \code{\link[phyloseq]{phyloseq}} object. Optionally, sample metadata can
#' be included and control samples can be removed.
#'
#' The OTU table is formatted with taxa as rows and samples as columns.
#' If a \code{cluster_id} column is present in the OTU table, reads are
#' aggregated by cluster and sample; otherwise, reads are aggregated by ASV
#' and sample. Missing abundances are replaced by zero.
#'
#' If a \code{replicate} column is present, the sample identifier is
#' constructed by combining \code{sample} and \code{replicate} with a
#' hyphen. The same transformation is applied to the sample metadata when
#' provided.
#'
#' @param otu Character string giving the path to the OTU/ASV abundance
#'   table, or a data frame. The table
#'   must contain \code{asv_id}, \code{sample}, and \code{read_count}
#'   columns. Optional columns include \code{asv}, \code{replicate}, and
#'   \code{cluster_id}.
#' @param tax Character string giving the path to the taxonomy table,  or a data frame. 
#'   The table must contain an
#'   \code{asv_id} column and the relevant taxonomic ranks.
#' @param samples Optional character string giving the path to the sample
#'   metadata table, or a data frame.
#'   When provided, the table must contain \code{sample} and
#'   \code{sample_type} columns. The metadata are added to the resulting
#'   \code{phyloseq} object as sample data.
#' @param outfile Optional character string specifying the output file
#'   path for the \code{phyloseq} object. If \code{NULL}, the file is not
#'   written to disk.
#' @param rm_control Logical; if \code{TRUE}, samples whose
#'   \code{sample_type} is not \code{"real"} are removed. 
#'   When \code{TRUE}, \code{samples} must be provided.
#' @param sep Character used to separate fields in the input files.
#' @param log_file Character string specifying the path to the CSV log file.
#'   The path is resolved with the following priority:
#'   \enumerate{
#'     \item the `log_file` argument, if explicitly provided by the user;
#'     \item the package-level option/variable storing a default log path
#'       (if set);
#'     \item `"vtamR_log.csv"` in the current working directory, used as a
#'       last resort if neither of the above is set.
#'   }
#'   To disable logging entirely, set `log_file = NA`.
#'
#' @return A \code{\link[phyloseq]{phyloseq}} object containing:
#' \itemize{
#'   \item an OTU table with taxa as rows and samples as columns;
#'   \item a taxonomy table containing taxonomy for the taxa present in the
#'     OTU table; and
#'   \item sample metadata when \code{samples} is provided.
#' }
#'
#' @details
#' The taxonomy table contains the columns \code{domain}, \code{phylum},
#' \code{class}, \code{order}, \code{genus}, and \code{species}. If a
#' \code{kingdom} column is present in the input taxonomy table, it is also
#' retained.
#'
#' Taxa present in the taxonomy table but absent from the OTU table are
#' removed. When sample metadata are provided, only the first metadata row
#' for each sample is retained.
#'
#' The function requires the \pkg{phyloseq} package. If it is not installed,
#' the function stops and provides installation instructions.
#'
#' @examples
#' \dontrun{
#' # Create a phyloseq object without sample metadata
#' physeq <- format_for_phyloseq(
#'   otu = "otu.csv",
#'   tax = "taxonomy.csv"
#' )
#'
#' # Include sample metadata
#' physeq <- format_for_phyloseq(
#'   otu = "otu.csv",
#'   tax = "taxonomy.csv",
#'   samples = "samples.csv"
#' )
#'
#' # Remove control samples
#' physeq <- format_for_phyloseq(
#'   otu = "otu.csv",
#'   tax = "taxonomy.csv",
#'   samples = "samples.csv",
#'   rm_control = TRUE
#' )
#' }
#'
#' @export

format_for_phyloseq <- function(otu, tax, samples = NULL, outfile = NULL,
                                rm_control = FALSE, 
                                log_file = NULL, sep = ",")
{
  
  if (missing(otu)) stop("Argument 'otu' is required")
  if (missing(tax)) stop("Argument 'tax' is required")
  #  if (missing(samples)) stop("Argument 'samples' is required")
  
  ###### test if phyloseq is installed
  if (!requireNamespace("phyloseq", quietly = TRUE) ) {
    stop(
      "Package 'phyloseq' is required for this function.\n",
      "Please install it with:\n",
      "  if (!requireNamespace('BiocManager', quietly = TRUE))\n",
      "    install.packages('BiocManager')\n",
      "  BiocManager::install('phyloseq')",
      call. = FALSE
    )
  }
  
  # get function name, all arguments and stat time
  # resolve path to log_file
  log_file <- get_log_file(file = log_file)
  # get function name, all arguments and stat time. If log_file == NA, no log
  log <- collect_log(file = log_file)
  
  if(rm_control & is.null(samples)){
    msg <- "samples must be specified when rm_control is TRUE and contain a sample and sample_type columns"
    stop(msg)
  }
  if(!is.null(samples)){
    sample_type_df <- read_input(samples, sep = sep) 
  }
  
  ### otu ########################################
  
  otu_mat <- read_input(otu, sep = sep)
  # make one column with sample-replicate
  if("replicate" %in% colnames(otu_mat)){
    otu_mat <- otu_mat  %>%
      mutate(sample = paste(sample, replicate, sep="-")) %>%
      select(-replicate)
    
    if(!is.null(samples)){
      sample_type_df <- sample_type_df %>%
        mutate(sample = paste(sample, replicate, sep = "-"))
    }
  }
  
  if("cluster_id" %in% colnames(otu_mat)){ # if cluster_id make output with clusters
    otu_mat <- otu_mat  %>%
      select(-asv_id, -asv) %>%
      group_by(cluster_id, sample) %>%
      summarise(read_count = sum(read_count), .groups = "drop") %>%
      rename(asv_id = cluster_id) 
  } else {
    otu_mat <- otu_mat  %>%
      select(-asv) %>%
      group_by(asv_id, sample) %>%
      summarise(read_count = sum(read_count), .groups = "drop")
  }
  
  if(rm_control){
    
    sample_type_df <- sample_type_df %>%
      filter(sample_type == "real")
    
    otu_mat <- otu_mat %>%
      filter(sample %in%  sample_type_df$sample) 
  }
  
  otu_mat <- pivot_wider(otu_mat, 
                         names_from = sample,
                         values_from = read_count,
                         values_fill = 0)
  otu_mat <- as.data.frame(otu_mat)
  rownames(otu_mat) <- otu_mat$asv_id
  otu_mat <- select(otu_mat, -asv_id)
  otu_mat <- as.matrix(otu_mat)
  
  ### tax ########################################
  
  tax_mat <- read_input(tax, sep = sep)
  tax_mat <- as.data.frame(tax_mat)
  # remove taxa not in otu
  tax_mat <- tax_mat %>%
    filter(asv_id %in% rownames(otu_mat))
  rownames(tax_mat) <- tax_mat$asv_id
  
  if("kingdom" %in% colnames(tax_mat)){
    tax_mat <- tax_mat %>%
      select(domain, kingdom, phylum, class, order, genus, species)
  } else {
    tax_mat <- tax_mat %>%
      select(domain, phylum, class, order, genus, species)
  }
  tax_mat <- as.matrix(tax_mat)
  
  ### samples ########################################
  if(!is.null(samples)){
    sample_type_df <- sample_type_df %>%
      group_by(sample) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    sample_df <- as.data.frame(sample_type_df)
    rownames(sample_df) <- sample_df$sample
    sample_df <- select(sample_df, -sample)
  }
  
  ### phyloseq ########################################
  if(is.null(samples)){
    OTU = phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
    TAX = phyloseq::tax_table(tax_mat)
    phy_object <- phyloseq::phyloseq(OTU, TAX)
  } else {
    OTU = phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
    TAX = phyloseq::tax_table(tax_mat)
    sample_df = phyloseq::sample_data(sample_df)
    
    phy_object <- phyloseq::phyloseq(OTU, TAX, sample_df)
  }
  
  if(!is.null(outfile)){
    check_dir(outfile, is_file = TRUE)
    saveRDS(phy_object, file = outfile)
  }
  
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(phy_object)
}
