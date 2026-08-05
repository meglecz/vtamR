library(dplyr)
library(vtamR)

sep = ","
log = "~/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/vtamR_log_zfzr.csv"
r_versions = "~/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/R_package_versions.csv"
outfile = "~/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/MIEM_bioinformatics.csv"
miem_bioinformatics(log = log, r_versions = r_versions, outfile = outfile)


log = "/home/meglecz/vtamR_benchmark/pipeline_output_11/shark/vtamR/vtamR_log.csv"
r_versions = "/home/meglecz/vtamR_benchmark/pipeline_output_11/shark/vtamR/R_package_versions.csv"
outfile = "/home/meglecz/vtamR_benchmark/pipeline_output_11/shark/vtamR/MIEM_bioinformatics.csv"
miem_bioinformatics(log = log, r_versions = r_versions, outfile = outfile)





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
#' @param log A data.frame containing the vtamR log or the path to a CSV
#'   log file.
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
  log,
  r_versions,
  outfile = "",
  sep = ","){
  
  
  # read input 
  log_df <-read_input(log, sep = sep)
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
  
  ###############################################################################
  ## Parameter table
  ###############################################################################
  
  table_to_text <- function(x){
    
    txt <- paste(
      capture.output(
        print(
          x,
          row.names = FALSE
        )
      ),
      collapse = "\n"
    )
    return(txt)
  }
  
  ###############################################################################
  ## Fill one MIEM field
  ###############################################################################
  
  set_miem_field <- function(miem, field, text) {
    
    i <- match(field, miem$Step)
    
    if (is.na(i))
      stop("Unknown MIEM field: ", field)
    
    miem$Information[i] <- text
    
    return(miem)
    
  }
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
                      "demultiplex_and_trim",
                      "demultiplex_fastq_pairs")
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
#' @noRd

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






