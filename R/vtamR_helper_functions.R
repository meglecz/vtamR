

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
#' @examples
#' \dontrun{
#' my_function <- function(a, b = 10, method = "mean") {
#'     if(!is.null(log_file)){
#'     log <- collect_log()
#'     }
#'   ## rest of the function
#' }
#' }
#'
#' @export
#' 
collect_log <- function() {
  
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
#' @examples
#' \dontrun{
#' my_function <- function(a, b = 10, method = "mean", file=NULL) {
#'
#'   log <- if(!is.null(log_file)){
#'     log <- collect_log()
#'   }
#'   
#'
#'   ## rest of the function
#'
#'   complete_log(log, file=log_file)
#' }
#' }
#'
#' @export
#' 
write_log <- function(log, file = NULL, sep=",") {
  
  if (is.null(file))
    return(invisible(log))
  
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

