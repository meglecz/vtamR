#' Pool data from different markers
#'
#' Take two or more input files containing filtered results of the same samples 
#' from different but strongly overlapping markers.
#' Files are in long format with asv_id, sample, replicate (optional), 
#' read_count and asv columns.
#'  
#' ASVs identical on their overlapping 
#' regions are pooled into groups, and different ASVs of the same group 
#' are pooled under the centroid (longest ASV of the group). The asv_id are
#' prefixed by the marker, to avoid confounding different ASVs of different 
#' markers, with the same id.
#' Pooling can take the mean of the read counts of the ASV (default), their sum
#' or maximum.
#'  
#' @param ... Data frames with the following variables: 
#' marker, asv_id, sample, replicate (optional), read_count, asv.
#' @param method Character string specifying how read counts should be pooled.
#'   Must be one of "mean", "max" or "sum".
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param asv_with_centroids Character string: csv file name of the output file 
#' containing the same information as the concatenated input files, 
#' completed by centroid_id and centroid columns. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param vsearch_path Character string: path to vsearch executables. 
#' @param num_threads Positive integer: Number of CPUs. If 0, use all available CPUs.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @return Data frame with asv_id, sample, replicate (optional), read_count, asv columns.
#' @examples
#' \dontrun{
#' markers_pooled <- pool_markers(df_mfzr, df_zfzr, method="mean")
#' }
#' @export
#'
pool_markers <- function(..., 
                         method="mean", 
                         outfile="", 
                         asv_with_centroids="", 
                         sep=",", 
                         vsearch_path="vsearch", 
                         num_threads=0,
                         quiet=T
){
  
  #### method
  method <- match.arg(method, c("mean", "max", "sum"))
  fun <- switch(method,
                mean = function(x) mean(x, na.rm = TRUE),
                max  = function(x) max(x, na.rm = TRUE),
                sum  = function(x) sum(x, na.rm = TRUE))
  
  #### num_threads
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  #### make tmp_dir
  tmp_dir <-paste('tmp_pool_markers_', 
                  trunc(as.numeric(Sys.time())), 
                  sample(1:100, 1), 
                  sep='')
  tmp_dir <- file.path(tempdir(), tmp_dir)
  check_dir(tmp_dir)
  
  #### concatenate input df
  # take first, determine repl_bool
  #  df_list <- list(mfzr, zfzr)
  df_list <- list(...)
  df <-  df_list[[1]]
  if("replicate" %in% colnames(df)){
    repl_bool <- TRUE
  }else{
    repl_bool <- FALSE
  }
  # read the other df
  for(i in 2:length(df_list)){
    
    if(repl_bool){
      tmp <- df_list[[i]] %>%
        select(marker, asv_id, sample, replicate, read_count, asv)
    }else{
      tmp <- df_list[[i]] %>%
        select(marker, asv_id, sample, read_count, asv)
    }
    df <- rbind(df, tmp)
  }
  
  
  ###
  # Pool ASVs identical on their overlapping region
  ###
  # add marker to asv_id to avoid incompatibility among asv_id across markers 
  df <- df %>%
    mutate(asv_id = paste(marker, asv_id, sep="_"))
  
  asvs <- df %>%
    group_by(asv_id, asv) %>%
    summarize("rc" = sum(read_count), .groups="drop")
  
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
  ##### run cmd
  args <- c(
    "--cluster_smallmem", fasta,
    "--consout", centroids_file,
    "--blast6out", blastout_file,
    "--id", 1 
  )
  if(num_threads > 0){
    args <- append(args, c("-threads", num_threads))
  }
  if(quiet){
    args <- append(args, c("--quiet"))
  }
  run_system2(vsearch_path, args, quiet=quiet)
  
  ###
  # Make cent data frame with a complete list of ASVs and the centroïd for each of them.
  ###
  # read the ids of centoids, and get the list of centroids
  # >centroid=mfzr_2374;seqs=2
  cent <- read.table(centroids_file)
  colnames(cent) <- c("centroid_id")
  cent <- cent %>%
    filter(grepl(">centroid=", centroid_id)) # keep only fasta definition lines
  cent$centroid_id <- gsub(">centroid=", "", cent$centroid_id)
  cent$nbseq <-   gsub(".+;seqs=", "", cent$centroid_id)
  cent$centroid_id <- gsub(";.+", "", cent$centroid_id)
  cent$nbseq <- as.numeric(cent$nbseq)
  
  # add to centroide the asv_id that are in the same cluster
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
  # add the centroid sequence to each centroid_id in df
  df <- left_join(df, asvs, by=c("centroid_id"="asv_id")) %>%
    select(-length, -rc) %>%
    rename("asv"=asv.x, "centroid"=asv.y) %>%
    arrange(centroid_id, marker)
  # order the columns
  if(repl_bool){
    df <- df %>%
      select(centroid_id,asv_id,marker,sample,replicate,read_count,asv,centroid)
  }else{
    df <- df %>%
      select(centroid_id,asv_id,marker,sample,read_count,asv,centroid)
  }
  
  
  if(repl_bool){
    df_pool <- df %>%
      group_by(centroid_id, sample, replicate) %>%
      summarize("read_count"=round(fun(read_count), digits=0), .groups =  "drop")
  }else{
    df_pool <- df %>%
      group_by(centroid_id, sample) %>%
      summarize("read_count"=round(fun(read_count), digits=0), .groups =  "drop")
  }
  
  
  # add asv column and select columns
  # df_pool is a simple output with the format identical to the read_count_sample dfs,
  # but no info on the asv that has been pooled together
  df_pool <- left_join(df_pool, asvs, by=c("centroid_id" = "asv_id"))
  if(repl_bool){
    df_pool <- df_pool %>%
      select("asv_id"=centroid_id, sample, replicate, read_count, asv) 
  }else{
    df_pool <- df_pool %>%
      select("asv_id"=centroid_id, sample, read_count, asv) 
  }
  
  
  if(asv_with_centroids != ""){
    check_dir(asv_with_centroids, is_file=TRUE)
    write.table(df, file=asv_with_centroids, sep=sep, row.names = F)
  }
  
  
  unlink(tmp_dir, recursive = TRUE)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(df_pool, file=outfile, sep=sep, row.names = F)
  }
  
  return(df_pool)
}

#' Pool Datasets
#' 
#' Deprecated: This function has been replaced bu pool_markers and pool_datasets.
#' 
#' Take several input files, each in long format containing 
#' asv_id, sample, replicate (optional), read_count and asv columns.
#' Pool the different data sets, if all have the same marker.
#'  
#' If more than one markers, ASVs identical on their overlapping 
#' regions are pooled into groups, and different ASVs of the same group 
#' are pooled under the centroid (longest ASV of the group). The asv_id are
#' prefixed by the marker, to avoid confounding different ASVs of different 
#' markers, with the same id.
#' Pooling can take the mean of the read counts of the ASV (default) or their sum.
#'  
#' @param files Data frame with the following variables: file (name of input files), marker.
#' Input files must have asv_id, sample, replicate (optional), read_count and asv columns.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param asv_with_centroids Character string: csv file name of the output file 
#' containing the same information as the concatenated input files, 
#' completed by centroid_id and centroid columns.
#' @param sep Field separator character in input and output csv files.
#' @param mean_over_markers logical: If TRUE, the mean read count is calculated 
#' over different ASVs of each cluster. Sum otherwise.
#' @param vsearch_path Character string: path to vsearch executables. 
#' @param num_threads Positive integer: Number of CPUs. If 0, use all available CPUs.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @return Data frame with asv_id, sample, replicate (optional), read_count, asv columns.
#' @examples
#' \dontrun{
#' files <- data.frame(file=c("vtamR_test/out_mfzr/14_PoolReplicates.csv", 
#'     "vtamR_test/out_zfzr/14_PoolReplicates.csv"),
#'     marker=c("MFZR", "ZFZR"))
#' PoolDatasets(files, vsearch_path=vsearch_path)
#' }
#' @export
PoolDatasets <- function(files, 
                         outfile="", 
                         asv_with_centroids="", 
                         sep=",", 
                         mean_over_markers=T, 
                         vsearch_path="vsearch", 
                         num_threads=0,
                         quiet=T
){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  
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
  
  # Read first file set replicate_col, initialise df
  marker <- files[1, "marker"]
  df <- read.csv(files[1, "file"], sep=sep)
  if("replicate" %in% colnames(df)){
    replicate_col <- TRUE
    df <- df %>%
      select(asv_id, sample, replicate, read_count, asv)
  }else{
    replicate_col <- FALSE
    df <- df %>%
      select(asv_id, sample, read_count, asv)
  }
  df$marker <- rep(marker, nrow(df)) # add maker
  samples <- unique(paste(df$marker, df$sample, sep="."))
  
  
  for(i in 2:nrow(files)){
    file <- files[i, "file"]
    marker <- files[i, "marker"]
    
    tmp <- read.csv(file, sep=sep)
    
    if(replicate_col){
      tmp <- tmp %>%
        select(asv_id, sample, replicate, read_count, asv)
    }else{
      tmp <- tmp %>%
        select(asv_id, sample, read_count, asv)
    }
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
  # complete the asv_id by 
  if(length(marker_list) > 1){ 
    # add marker to asv_id to avoid incompatibility among asv_id across markers
    df <- df %>%
      mutate(asv_id = paste(marker, asv_id, sep="_"))
    # get full list of ASVs
    
    asvs <- df %>%
      group_by(asv_id, asv) %>%
      summarize("rc" = sum(read_count), .groups="drop")
    
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
    ##### run cmd
    args <- c(
      "--cluster_smallmem", fasta,
      "--consout", centroids_file,
      "--blast6out", blastout_file,
      "--id", 1 
    )
    if(num_threads > 0){
      args <- append(args, c("-threads", num_threads))
    }
    if(quiet){
      args <- append(args, c("--quiet"))
    }
    run_system2(vsearch_path, args, quiet=quiet)
    
    ###
    # Make cent data frame with a complete list of ASVs and the centroïd for each of them.
    ###
    # read the ids of centoids, and get the list of centroids
    cent <- read.table(centroids_file)
    colnames(cent) <- c("centroid_id")
    cent <- cent %>%
      filter(grepl(">centroid=", centroid_id)) # keep only fasta definition lines
    cent$centroid_id <- gsub(">centroid=", "", cent$centroid_id)
    cent$nbseq <-   gsub(".+;seqs=", "", cent$centroid_id)
    cent$centroid_id <- gsub(";.+", "", cent$centroid_id)
    #    cent$centroid_id <- as.numeric(cent$centroid_id)
    cent$nbseq <- as.numeric(cent$nbseq)
    
    # add to centroide the asv_id that are in the same cluster
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
    # add the centroid sequence to each centroid_id in df
    df <- left_join(df, asvs, by=c("centroid_id"="asv_id")) %>%
      select(-length, -rc) %>%
      rename("asv"=asv.x, "centroid"=asv.y) %>%
      arrange(centroid_id, marker)
    # order the columns
    if(replicate_col){
      df <- df %>%
        select(centroid_id,asv_id,marker,sample,replicate,read_count,asv,centroid)
    }else{
      df <- df %>%
        select(centroid_id,asv_id,marker,sample,read_count,asv,centroid)
    }
    
    
    
    if(mean_over_markers){
      if(replicate_col){
        df_pool <- df %>%
          group_by(centroid_id, sample, replicate) %>%
          summarize("read_count"=round(mean(read_count), digits=0), .groups =  "drop")
      }else{
        df_pool <- df %>%
          group_by(centroid_id, sample) %>%
          summarize("read_count"=round(mean(read_count), digits=0), .groups =  "drop")
      }
    }else{# sum over markers
      if(replicate_col){
        df_pool <- df %>%
          group_by(centroid_id, sample, replicate) %>%
          summarize("read_count"=sum(read_count), .groups =  "drop" ) 
      }else{
        df_pool <- df %>%
          group_by(centroid_id, sample) %>%
          summarize("read_count"=sum(read_count), .groups =  "drop" ) 
      }
    }
    
    # add asv column and select columns
    # df_pool is a simple output with the format identical to the read_count_sample dfs,
    # but no info on the asv that has been pooled together
    df_pool <- left_join(df_pool, asvs, by=c("centroid_id" = "asv_id"))
    if(replicate_col){
      df_pool <- df_pool %>%
        select("asv_id"=centroid_id, sample, replicate, read_count, asv) 
    }else{
      df_pool <- df_pool %>%
        select("asv_id"=centroid_id, sample, read_count, asv) 
    }
    
    
    if(asv_with_centroids != ""){
      check_dir(asv_with_centroids, is_file=TRUE)
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


#' Filter PCR errors
#' 
#' Filters only ASV and not occurrences
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
#'   occurrences below this value are ignored, to speed up the analyses.
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
#' filtered_read_count_df <- filter_pcr_error_old(
#'   read_count_df,
#'   vsearch_path = vsearch_path,
#'   pcr_error_var_prop = 0.2,
#'   max_mismatch = 2,
#'   by_sample = TRUE,
#'   sample_prop = 0.8
#' )
#' 
#' filtered_read_count_df <- filter_pcr_error_old(
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
filter_pcr_error_old <- function(read_count,
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




#' Filter PCR errors
#' 
#' Min_read_count applied to the total number of reads. Can filter ASv or occurrences.
#'
#' Remove ASVs flagged as potential PCR errors based on sequence similarity
#' (`max_mismatch`) and relative abundance (`pcr_error_var_prop`).
#'
#' An ASV is considered a PCR error when it is highly similar to a more abundant
#' ASV and its abundance ratio is less than or equal to `pcr_error_var_prop`.
#'
#' The analysis can be performed across the full dataset (`by_sample = FALSE`)
#' or independently within each sample (`by_sample = TRUE`).
#'
#' When the analysis is performed across the full dataset (`by_sample = FALSE`),
#' ASVs classified as PCR errors are removed entirely from the dataset.
#'
#' When the analysis is performed sample by sample (`by_sample = TRUE`),
#' `filter_occurrence` determines whether the entire ASV is removed
#' (`filter_occurrence = FALSE`) or only the occurrence is removed
#' (`filter_occurrence = TRUE`).
#'
#' When `filter_occurrence = FALSE`, an ASV is removed only if it is flagged as
#' a PCR error in at least the proportion of samples specified by `sample_prop`.
#'
#' When `filter_occurrence = TRUE`, the ASV is removed only from the samples in
#' which it is flagged as a PCR error and may remain present in other samples.
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
#' @param filter_occurrence Logical. If TRUE, an ASV is removed only from the
#'   samples in which it is flagged as a PCR error and may remain present in other
#'   samples. If FALSE, an ASV is removed entirely if it is flagged as a PCR
#'   error in at least the proportion of samples specified by sample_prop.
#'   Otherwise, none of its occurrences are removed.
#' @param sample_prop Numeric value between 0 and 1 specifying the minimum proportion 
#'   of samples in which an ASV must be flagged as a PCR error (when `by_sample = TRUE`) 
#'   to be removed.
#' @param min_read_count Positive integer specifying the minimum read count threshold; 
#'   ASVs with total read count below this value are ignored, to speed up the analyses.
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
#'   filter_occurrence = FALSE,
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
                             filter_occurrence=TRUE,
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
  if(filter_occurrence & by_sample==FALSE){
    warning("When by_sample==FALSE the filtering eliminates entire ASVs and not occurrences, 
            even is filter_occurrence is TRUE")
  }
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  
  # get unique list of ASV and only ASV with the total read count >=  min_read_count
  unique_asv_df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count = sum(read_count)) %>%
    filter(read_count >= min_read_count) %>%
    ungroup()
  
  if(by_sample){
    
    sample_list <- unique(read_count_df$sample)
    unique_asv_sample <- read_count_df %>%
      filter(asv %in% unique_asv_df$asv) %>% # use only asv with >= read_count
      group_by(asv, sample) %>%
      summarize(read_count = sum(read_count), .groups="drop")%>%
      ungroup()
    
    if(filter_occurrence){ # by_sample=TRUE, filter_occurrence = TRUE
      
      pcr_flags <- data.frame(
        asv = character(),
        sample = character(),
        PCRerror = numeric()
      )
      
      # loop over samples
      for(sample_loc in sample_list){
        # get unique list of ASVs with their total read_count in the sample
        sample_df <- unique_asv_sample %>%
          filter(sample == sample_loc) 
        
        # flag PCR errors; 
        # add one column to sample_df for each sample with 1 if ASV is flagged in the sample, 
        # 0 otherwise
        sample_df <- flag_pcr_error(sample_df, 
                                    vsearch_path=vsearch_path, 
                                    num_threads=num_threads,
                                    pcr_error_var_prop=pcr_error_var_prop, 
                                    max_mismatch=max_mismatch,
                                    quiet=quiet
        )
        sample_df <- sample_df %>%
          select(asv, sample, PCRerror)
        pcr_flags <- rbind(pcr_flags, sample_df)
      }
      # delete occurrences flagged as PCRerror
      read_count_df <- left_join(read_count_df, pcr_flags, by=c("sample", "asv")) %>%
        filter(PCRerror == 0 | is.na(PCRerror)) %>%
        select(-PCRerror)
      
    }else { # by_sample=TRUE, filter_occurrence = FALSE
      
      # loop over samples
      for(sample_loc in sample_list){
        # get unique list of ASVs with their total read_count in the sample
        sample_df <- unique_asv_sample %>%
          filter(sample == sample_loc) 
        
        # flag PCR errors; 
        # add one column to sample_df with 1 if ASV is flagged as a PCR error in the sample, 
        # 0 otherwise
        sample_df <- flag_pcr_error(sample_df, 
                                    vsearch_path=vsearch_path, 
                                    num_threads=num_threads,
                                    pcr_error_var_prop=pcr_error_var_prop, 
                                    max_mismatch=max_mismatch,
                                    quiet=quiet
        )
        sample_df <- sample_df %>%
          select(asv, PCRerror)
        
        # add a column for for each sample to unique_asv_df, 
        # with 1 if ASV is flagged in the sample, 0 otherwise
        unique_asv_df <- left_join(unique_asv_df, sample_df, by = "asv")
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
    } # end by_sample=TRUE, filter_occurrence ==FALSE
  } else { # end by_sample, Filter the all at once
    
    unique_asv_df <- flag_pcr_error(unique_asv_df, 
                                    vsearch_path=vsearch_path, 
                                    num_threads=num_threads,
                                    pcr_error_var_prop=pcr_error_var_prop, 
                                    max_mismatch=max_mismatch)
    
    
    unique_asv_df <- unique_asv_df %>%
      filter(PCRerror == 1)
    # eliminate potential PCRerrors from read_count_df
    read_count_df <- read_count_df %>%
      filter(!asv %in% unique_asv_df$asv)
  } # end by_sample = FALSE
  
  ###### Print output
  if(!is.null(outfile)){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}



#' Filter chimeric sequences
#' 
#' Filters ASV not occurrences
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


#' Summarize results for a MIEM checklist
#' 
#' Depending on the functions used and their order during the data analysis with vtamR,
#' the information that can or should be reported in the MIEM file may vary.
#' This function provides a comprehensive summary of the analysis by using
#' the information recorded in the log file.
#'
#' In particular, the function identifies the preprocessing and filtering
#' functions used during the analysis and extracts the number of reads and
#' ASVs present in the output of each step. It also summarizes read counts
#' by sample type, taxonomic assignments, and control classification
#' results when the corresponding information is available.
#'
#' The resulting files provide users with the information needed to complete
#' the MIEM checklist. Users can select and report the values that are most
#' appropriate for their MIEM checklist based on the functions and analyses
#' performed.
#'
#' @param log Character string giving the path to the MIEM log file.
#' @param outdir Character string giving the directory in which result
#'   summary files will be written. The directory is created if necessary.
#' @param sampleinfo Character string giving the path to the sample
#'   information file containing `sample` and `sample_type` columns. 
#'   If `NULL`, the function attempts to identify the
#'   most recently recorded `sampleinfo`, `fastqinfo`, or `fastainfo` file
#'   from the log, that usually contains this infromation.
#' @param taxa Character string giving the path to the taxonomic assignment
#'   file. If `NULL`, the function attempts to identify the taxonomic
#'   assignment output from the log.
#' @param asv Character string giving the path to an ASV read-count file.
#'   If `NULL`, the function attempts to identify the output of the last
#'   filtering step (excluding `cluster_asv`).
#' @param motu Character string giving the path to the mOTU file produced
#'   by `cluster_asv`. This argument is only required when clustering has
#'   been performed and the mOTU output cannot be identified from the log.
#' @param sep Character string used to separate fields in the input and
#'   output files. Defaults to `","`.
#'
#' @return The function writes summary files to `outdir`. It does not return a data frame.
#'
#' @details
#' The function performs the following analyses:
#'
#' \enumerate{
#'   \item Counts reads in the raw FASTQ input files when the corresponding
#'   information is recorded in the log.
#'   \item Counts reads after preprocessing steps.
#'   \item Counts reads and ASVs after each filtering step.
#'   \item Summarizes read counts by sample type for the input and final
#'     ASV datasets.
#'   \item Counts ASVs or mOTUs assigned to each major taxonomic rank.
#'   \item Extracts false-positive and false-negative occurrences from
#'   control classification results.
#' }
#'
#' When possible, input files are automatically identified from the log.
#' Explicitly providing `sampleinfo`, `taxa`, `asv`, or `motu` allows the
#' user to override this automatic detection.
#'
#' The following output files are written to `outdir`:
#'
#' \itemize{
#'   \item `Read_ASV_count_by_step.csv`: read and ASV counts after each
#'   processing step.
#'   \item `Read_count_by_sample.csv`: read-count summary statistics by
#'   sample type and dataset (before the first filtering and after the last filtering 
#'   step)
#'   \item `ASV_count_by_taxonomomic_rank.csv`: number of ASVs assigned to
#'   each taxonomic rank, when applicable.
#'   \item `mOTU_count_by_taxonomomic_rank.csv`: number of mOTUs assigned
#'   to each taxonomic rank when clustering has been performed.
#'   \item `false_positives_and_nevatives.csv`: false-positive and
#'   false-negative occurrences identified from control classification
#'   results.
#' }
#'
#' @examples
#' \dontrun{
#' miem_results(
#'   log = "path/to/miem_log.csv",
#'   outdir = "path/to/results"
#' )
#'
#' # Provide specific input files when they cannot be identified
#' # automatically from the log.
#' miem_results(
#'   log = "path/to/miem_log.csv",
#'   outdir = "path/to/results",
#'   sampleinfo = "path/to/sampleinfo.csv",
#'   taxa = "path/to/taxonomy.csv",
#'   asv = "path/to/asv.csv",
#'   motu = "path/to/motu.csv"
#' )
#' }
#'
#' @export
miem_results <- function(log, outdir, sampleinfo = NULL, taxa = NULL, asv = NULL, motu = NULL, sep = ","){
  
  read_count_steps <- data.frame(
    function_name = as.character(),
    read_count = as.numeric(),
    asv_count = as.numeric()
  )
  outdir <- check_dir(outdir)
  
  # read input 
  log_df <- read_input(log, sep = sep)
  
  ## Number of raw input reads #########################
  read_count_steps <- count_input_read_count(log_df, read_count_steps)
  
  ## Preprocessing steps #########################
  read_count_steps <- count_preprocessing_read_count(log_df, read_count_steps)
  
  ## Filtering steps #########################
  results <- count_filtering_read_count(log_df, read_count_steps)
  read_count_steps <- results[[1]]
  filter_df <- results[[2]]
  
  outfile <- file.path(outdir, "Read_ASV_count_by_step.csv")
  write.table(read_count_steps, file = outfile, sep = sep, row.names = FALSE)
  
  ######## Read count by sample and by sample type
  ## if sampleinfo = NULL get the last sampleinfo or fastainfo or fastqinfo used which is a file and not a df
  # get the sample type from sample info
  # get the input and output of filter_dfing, separate samples by type
  # get the minimum, maximum, mean, median, for each sample_type and dataset
  # columns:
  # input_real, input_mock, input_negative, filtered_real, filtered_mock, filtered_negative
  # rows: min, max, mean, median, number of samples
  
  # get sampleinfo if not provided by the user
  if(is.null(sampleinfo)){
    arguments <- c(
      "sampleinfo",
      "fastqinfo",
      "fastainfo"
    )
    sampleinfo <- log_df %>%
      filter(argument_name %in% arguments) %>%
      select(value) %>%
      filter(!grepl("<data.frame>", value)) %>% # get line, with filename and not data frame
      last() %>%
      pull()
    if(is.na(sampleinfo)){
      msg <- paste0("Sampleinfo file is not recorded in the log file. Please provide the filename using sampleinfo argument")
      stop(msg)
    }
  } 
  
  
  # select the data file which is the input of filtering and output of filtering (sample, ...)
  input_data_file <- filter_df$value[1]
  # get read count ... by sample for the input of filtering steps
  read_count_sample = extract_sample_read_count(read_count = input_data_file, sampleinfo = sampleinfo, dataset = "input_filter", sep = ",")
  
  # if user defines the output of any of the filtering steps (typically the last), get the same counts for that data set
  if(!is.null(asv)){ 
    tmp = extract_sample_read_count(read_count = asv, sampleinfo = sampleinfo, dataset = "user_defined_asv", sep = ",")
    read_count_sample <- rbind(read_count_sample, tmp)
  }else{ # guess the last asv filter
    tmp <- filter_df %>%
      filter(function_name != "cluster_asv")
    asv <- filter_df$value[nrow(tmp)]
    
    last_filter_name <-  filter_df$function_name[nrow(tmp)]
    last_filter_name <- paste0("output_", last_filter_name)
    tmp = extract_sample_read_count(read_count = asv, sampleinfo = sampleinfo, dataset = last_filter_name, sep = ",")
    read_count_sample <- rbind(read_count_sample, tmp)
    
  }
  
  outfile <- file.path(outdir, "Read_count_by_sample.csv")
  write.table(read_count_sample, file = outfile, sep = sep, row.names = FALSE)
  
  ######## get the results of taxassign
  # If taxassign before clustering, count assigned ASV by ASV, else only for mOTUs
  # Count the number of mOUTs/ ASV mOTU assigned to taxa at each major taxonomic level
  
  if(is.null(taxa)){ # user did not provide taxa
    functions <- c("assign_taxonomy_ltg",
                   "assign_taxonomy_rdp")
    arguments <- c("outfile")
    
    taxa <- get_log_entries(
      log_df,
      functions = functions,
      arguments = arguments,
      remove_duplicates = TRUE) %>%
      filter(!grepl("<data.frame>", value)) %>%
      first() %>%
      pull()
    
    if(is.na(taxa)){
      msg <- paste0("Taxonomic assignment file is not recorded in the log file. Please provide the filename using the taxa argument\n")
      stop(msg)
    }
  }
  
  # read taxa and make columns asv_id, rank_index, taxonomic_level
  taxa_df <- format_taxa(taxa, sep)
  
  # get the order of the functions
  tmp <- log_df %>%
    select(function_name, start_time) %>%
    distinct() %>%
    select(function_name) %>%
    mutate(order = rownames(.))
  
  taxassign_index <- tmp %>%
    filter(function_name %in% c("assign_taxonomy_ltg", "assign_taxonomy_rdp")) %>%
    first() %>%
    pull()
  
  cluster_index <- tmp %>%
    filter(function_name == "cluster_asv") %>%
    last() %>%
    pull()
  
  if(is.na(cluster_index) || cluster_index == "NA"){ # no clustering, work only with ASV
    asv_by_rank <- count_taxassing_by_rank(asv, taxa_df, sep = ",")
    outfile <- file.path(outdir, "ASV_count_by_taxonomomic_rank.csv")
    write.table(asv_by_rank, file = outfile, sep = sep, row.names = FALSE)
  } else if (!is.null(motu)){ # motu file was defined by the user, even if cluster_asv is not in log
    motu_by_rank <- count_taxassing_by_rank(motu, taxa_df, sep = ",")
    outfile <- file.path(outdir, "mOTU_count_by_taxonomomic_rank.csv")
    write.table(motu_by_rank, file = outfile, sep = sep, row.names = FALSE)
  } else { # clustering has been done and motu not defined by the user
    ##### get motu filename
    motu <- get_log_entries(
      log_df,
      functions = "cluster_asv",
      arguments = "outfile",
      remove_duplicates = TRUE) %>%
      select(value) %>%
      last() %>%
      pull()
    
    if(is.na(motu)){ # data frame instead of file
      msg <- paste0("The mOTU file is not recorded in the log file. Please provide the filename using the motu argument\n")
      warning(msg)
    }else{ 
      motu_by_rank <- count_taxassing_by_rank(motu, taxa_df, sep = ",")
      outfile <- file.path(outdir, "mOTU_count_by_taxonomomic_rank.csv")
      write.table(motu_by_rank, file = outfile, sep = sep, row.names = FALSE)
      
      if(cluster_index > taxassign_index){ # taxassing before clustering, use ASVs can be classified
        asv_by_rank <- count_taxassing_by_rank(asv, taxa_df, sep = ",")
        outfile <- file.path(outdir, "ASV_count_by_taxonomomic_rank.csv")
        write.table(asv_by_rank, file = outfile, sep = sep, row.names = FALSE)
      }
    }
  }
  
  outfile <- file.path(outdir, "ASV_count_by_taxonomomic_rank.csv")
  write.table(asv_by_rank, file = outfile, sep = sep, row.names = FALSE)
  
  ####### get the number of FP, FN TP 
  
  known_occurrences <- get_log_entries(
    log_df,
    functions = "classify_control_occurrences",
    arguments = "known_occurrences",
    remove_duplicates = TRUE) %>%
    select(value) %>% last() %>%  pull()
  
  false_negatives <- get_log_entries(
    log_df,
    functions = "classify_control_occurrences",
    arguments = "false_negatives",
    remove_duplicates = TRUE) %>%
    select(value) %>% last() %>%  pull()
  
  performance_metrics <- get_log_entries(
    log_df,
    functions = "classify_control_occurrences",
    arguments = "performance_metrics",
    remove_duplicates = TRUE) %>%
    select(value) %>% last() %>%  pull()
  
  fp <- read_input(known_occurrences, sep = sep) %>%
    filter(action == "delete") %>%
    mutate(occurrence_type = "FP") %>%
    select(occurrence_type, sample, asv_id, asv)
  
  fn <- read_input(false_negatives, sep = sep) 
  if(!"asv_id" %in% colnames(fn)){
    mutate(asv_id = NA)
  }
  fn <- fn %>%
    mutate(occurrence_type = "FN") %>%
    select(occurrence_type, sample, asv_id, asv)
  
  tmp <- rbind(fp, fn)
  
  performance_metrics_df <- read_input(performance_metrics, sep = sep)
  outfile <- file.path(outdir, "false_positives_and_nevatives.csv")
  write.table(tmp, file = outfile, sep = sep, row.names = FALSE)
  
}







#' Count reads in raw FASTQ input files
#'
#' Identifies the raw FASTQ input files recorded in the log and calculates
#' the total number of reads across all input files. The function searches
#' the log for preprocessing functions that can start directly from raw
#' FASTQ files and uses the corresponding `fastqinfo` and `fastq_dir`
#' arguments to locate the input files.
#'
#' @param log_df A data frame containing the parsed log information. It must
#'   contain the function names, argument names, and values required by
#'   [get_log_entries()].
#' @param read_count_steps A data frame containing read-count results from
#'   previous processing steps. It must contain `function_name`,
#'   `read_count`, and `asv_count` columns.
#' @param sep Character string used to separate fields in the `fastqinfo`
#'   input file. Defaults to `","`.
#'
#' @return A data frame containing the read-count results with an additional
#'   row named `"RAW INPUT FASTQ"` when raw FASTQ input information is found
#'   in the log. The row contains the total number of reads across all raw
#'   FASTQ files and `NA` for the ASV count.
#'
#' @details
#' The function searches for the `fastqinfo` and `fastq_dir` arguments
#' associated with the `demultiplex_and_trim_fastq` and
#' `merge_fastq_pairs` functions.
#'
#' The `fastqinfo` file is expected to contain a `fastq_fw` column listing
#' the forward FASTQ files. Each file is located relative to `fastq_dir`,
#' and the number of reads is calculated using [count_reads()] with
#' `file_type = "fastq"`.
#'
#' If the log does not contain information about a function starting from
#' raw FASTQ files, a warning is issued and `read_count_steps` is returned
#' unchanged.
#'
#' @keywords internal
count_input_read_count <- function(log_df, read_count_steps, sep = ","){
  
  ## Number of reads in the input #########################
  # These are the functions that can start for the raw fastq files
  preprocess_functs <- c("demultiplex_and_trim_fastq", "merge_fastq_pairs") 
  args <- c("fastqinfo", "fastq_dir")
  
  tmp <- get_log_entries(
    log_df,
    functions = preprocess_functs,
    arguments = args,
    remove_duplicates = TRUE) %>% 
    group_by(argument_name) %>% # get first fastq_dir and fastqinfo
    summarize(value = first(value), .groups = "drop")
  
  if(nrow(tmp) > 0){
    
    fastq_dir <- tmp %>%
      filter(argument_name == "fastq_dir") %>%
      pull(value)
    if(is.na(fastq_dir) || is.null(fastq_dir) || fastq_dir == "NULL"){
      fastq_dir = "."
    }
    
    fastqinfo <- tmp %>%
      filter(argument_name == "fastqinfo") %>%
      pull(value)
    
    if(is.na(fastqinfo) || is.null(fastqinfo) || fastqinfo == "NULL" || grepl("<data.frame>", fastqinfo)){
      f <- paste0(preprocess_functs, collapse = ", ")
      msg <- paste0( "The log file does not contain information about functions starting from ", 
                     "raw FASTQ files (", f, ") or their info file. ", 
                     "You can use the count_reads_in_dir() function to count the number of reads ", 
                     "in the input file." )
      warning(msg)
    } else{
      
      fastq_files <- read.table(fastqinfo, header = TRUE, sep = sep) %>%
        select(fastq_fw) %>%
        distinct()
      
      total_read_count_input = 0
      for(i in 1:length(fastq_files$fastq_fw) ){
        file <- file.path(fastq_dir, fastq_files$fastq_fw[i])
        n <- count_reads(file, file_type="fastq")
        total_read_count_input = total_read_count_input + n
      }
      
      read_count_steps <- read_count_steps %>%
        add_row(function_name = "RAW INPUT FASTQ", read_count = total_read_count_input, asv_count = NA)
    }
    
  } else {
    f <- paste0(preprocess_functs, collapse = ", ")
    msg <- paste0( "The log file does not contain information about functions starting from ", 
                   "raw FASTQ files (", f, "). ", 
                   "You can use the count_reads_in_dir() function to count the number of reads ", 
                   "in the input file." )
    warning(msg)
  }
  return(read_count_steps)
}


#' Count reads after preprocessing steps
#'
#' Extracts preprocessing steps recorded in the log and calculates the total
#' number of reads remaining after each preprocessing step. The function
#' identifies preprocessing output directories from the `outdir` arguments
#' recorded in the log and retrieves read counts from the corresponding
#' `info.csv` files.
#'
#' @param log_df A data frame containing the parsed log information. It must
#'   contain the function names, argument names, and values required by
#'   [get_log_entries()].
#' @param read_count_steps A data frame containing read-count results from
#'   previous processing steps. It must contain `function_name`,
#'   `read_count`, and `asv_count` columns.
#' @param sep Character string used to separate fields in the preprocessing
#'   output files. Defaults to `","`.
#'
#' @return A data frame containing the read counts for the preprocessing
#'   steps, appended to the input `read_count_steps` data frame. The returned
#'   data frame contains the columns `function_name`, `read_count`, and
#'   `asv_count`.
#'
#' @details
#' The function searches the log for the following preprocessing functions:
#' `demultiplex_and_trim_fastq`, `demultiplex_and_trim_fasta`,
#' `merge_fastq_pairs`, `subsample_fasta`, and `trim_primers`.
#'
#' For each preprocessing step, the output directory is obtained from the
#' `outdir` argument recorded in the log. The function then identifies the
#' corresponding `info.csv` file and sums its `read_count` column.
#'
#' If no preprocessing functions are found in the log, a warning is issued
#' and the input `read_count_steps` data frame is returned unchanged.
#'
#' @keywords internal
count_preprocessing_read_count<- function(log_df, read_count_steps, sep = ","){
  
  # all preprocess functions
  functions <- c("demultiplex_and_trim_fastq", # outdir, fastqinfo.csv
                 "demultiplex_and_trim_fasta", # outdir, sampleinfo.csv
                 "merge_fastq_pairs", # outdir, fastainfo.csv
                 "subsample_fasta", # outdir, fastainfo.csv
                 "trim_primers") # outdir, sampleinfo.csv
  arguments <- c("outdir") # fastqinfo, fatsainfo, sorterinfo are in th outdir
  
  preprocess <- get_log_entries(
    log_df,
    functions = functions,
    arguments = arguments,
    remove_duplicates = TRUE) %>%
    mutate("read_count" = NA,
           "asv_count" = NA)
  
  if(nrow(preprocess) == 0){
    f <- paste0(functions, collapse = ", ")
    msg <- paste0( "The log file does not contain information about preprocessing functions ", 
                   "raw FASTQ files (", f, "). ", 
                   "You can use the count_reads_in_dir() function to count the number of reads." )
    warning(msg)
    return(read_count_steps)
  } else {
    
    for(i in 1: nrow(preprocess)){
      
      dir <- preprocess$value[i]
      files <- list.files(dir, pattern = "info\\.csv$") # get the name of the info file
      if(!is.null(files[1])){
        
        info_file <- file.path(dir, files[1])
        info_df <- read.table(file = info_file, sep = sep, header = TRUE) 
        
        info_df <- info_df %>%
          select((ncol(.) - 1):ncol(.)) %>% # keep the last filename column and the read_count
          distinct() # get unique list
        
        read_count <- sum(info_df$read_count)
        preprocess[i,"read_count"] <- read_count
      }
    }
    
    preprocess <- preprocess %>%
      select(function_name, read_count, asv_count)
    
    read_count_steps <- rbind(read_count_steps, preprocess)
    return(read_count_steps)
  }
}


#' Count reads and ASVs after filtering steps
#'
#' Extracts filtering steps recorded in the log and calculates the total
#' number of reads and ASVs remaining after each filtering step. The
#' filtering output files are identified from the `outfile` arguments
#' recorded in the log.
#'
#' @param log_df A data frame containing the parsed log information. It must
#'   contain the function names, argument names, and corresponding values
#'   required by [get_log_entries()].
#' @param read_count_steps A data frame containing read-count results from
#'   previous processing steps. It must contain `function_name`,
#'   `read_count`, and `asv_count` columns.
#' @param sep Character string used to separate fields in the filtering
#'   output files. Defaults to `","`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{read_count_steps}{A data frame containing the number of reads and
#'   ASVs remaining at each processing step.}
#'   \item{filter}{A data frame containing the filtering steps extracted
#'   from the log, together with their corresponding read and ASV counts.}
#' }
#'
#' @details
#' The function searches the log for known filtering functions and their
#' `outfile` arguments. For each filtering output file, reads are summed
#' across ASVs and the number of distinct ASVs is counted.
#'
#' If no filtering functions are found in the log, a warning is issued and
#' the input `read_count_steps` is returned without being completed .
#'
#' @keywords internal

count_filtering_read_count <- function(log_df, read_count_steps, sep = ","){
  
  # all filtering functions
  functions <- c("dereplicate",
                 "denoise_by_swarm",
                 "denoise_by_swarm",
                 "filter_asv_global",
                 "filter_contaminant",
                 "filter_chimera",
                 "filter_stop_codon",
                 "filter_contaminant",
                 "filter_indel",
                 "filter_pcr_error",
                 "filter_min_replicate",
                 "filter_occurrence_read_count",
                 "filter_occurrence_sample",
                 "filter_occurrence_variant",
                 "filter_pcr_error",
                 "pool_filters",
                 "filter_replicate",
                 "cluster_asv"
  )
  
  arguments <- c("outfile")
  
  filter_df <- get_log_entries(
    log_df,
    functions = functions,
    arguments = arguments,
    remove_duplicates = TRUE) %>%
    mutate("read_count" = NA,
           "asv_count" = NA)
  
  
  if(nrow(filter_df) == 0){
    f <- paste0(functions, collapse = ", ")
    msg <- paste0( "The log file does not contain information about filtering functions (", 
                   f, "). ")
    warning(msg)
    l <- list(read_count_steps, filter_df)
    return(l)
  }
  
  for(i in 1: nrow(filter_df)){
    file <- filter_df$value[i]
    print(file)
    if(!is.null(file) && file != "NULL"){
      info_df <- read.table(file, sep = sep, header = TRUE) %>%
        group_by(asv) %>%
        summarize(read_count = sum(read_count), .groups = "drop")
      
      filter_df[i,"read_count"] <- sum(info_df$read_count)
      filter_df[i,"asv_count"] <- nrow(info_df)
    }
  }
  
  tmp <- filter_df %>%
    select(function_name, read_count, asv_count)
  
  read_count_steps <- rbind(read_count_steps, tmp)
  
  l <- list(read_count_steps, filter_df)
  
  return(l)
  
}


#' Summarize read counts by sample type
#'
#' Calculates summary statistics for read counts across samples, grouped by
#' sample type. Read counts are first summed for each sample and then
#' summarized within each sample type (read, negative, mock).
#'
#' @param read_count Character string giving the path to the input read-count
#'   file. The file must contain at least `sample` and `read_count` columns.
#' @param sampleinfo Character string giving the path to the sample information
#'   file. The file must contain `sample` and `sample_type` columns.
#' @param dataset Character string identifying the dataset represented by the
#'   read counts. Defaults to `"input_to_filtering"`.
#' @param sep Character string used to separate fields in the input files.
#'   Defaults to `","`.
#'
#' @return A data frame containing read-count summary statistics for each
#'   sample type. The returned data frame includes the dataset name,
#'   minimum and maximum read counts, mean and median read counts, and
#'   number of samples.
#'
#' @keywords internal

extract_sample_read_count <- function(read_count, sampleinfo, dataset = "input_to_filtering", sep = ","){
  
  
  if(!is.null(sampleinfo) && !is.na(sampleinfo) && sampleinfo != "NULL" && !is.null(read_count) && !is.na(read_count) && read_count != "NULL"){
    sampleinfo_df <- read_input(sampleinfo, sep = sep) %>%
      select(sample, sample_type) %>%
      distinct()
    
    read_count_sample <- read_input(read_count, sep = sep) %>%
      group_by(sample) %>%
      summarise(read_count = sum(read_count), .groups = "drop") %>%
      left_join(sampleinfo_df, by = "sample") %>%
      group_by(sample_type) %>%
      summarize(
        minimum_read_count = min(read_count),
        maximum_read_count = max(read_count),
        mean_read_count = round(mean(read_count), digits = 0),
        median_read_count = round(median(read_count), digits = 0),
        number_of_samples = n()
      ) %>%
      mutate(dataset = dataset) %>%
      select(dataset, everything())
  }else{
    read_count_sample <- data.frame(
      minimum_read_count = numeric(),
      maximum_read_count = numeric(),
      mean_read_count = numeric(),
      median_read_count = numeric(),
      number_of_samples = numeric(),
      dataset = character()
    )
  }
  
  return(read_count_sample)
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
  
  if( !is.na(read_count) && !is.null(read_count) && read_count != "NULL" && read_count != "" && !grepl("<data.frame>", read_count)){
    
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
  }
  return(asv_by_rank)
}


table_to_text <- function(x){
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  m <- format(x, justify = "left", na.encode = FALSE)
  
  header <- names(m)
  widths <- pmax(nchar(header), vapply(m, function(col) max(nchar(col)), integer(1)))
  
  pad <- function(s, w) formatC(s, width = -w)  # left-justify, pad right
  
  header_line <- paste(mapply(pad, header, widths), collapse = "  ")
  body_lines  <- do.call(paste, c(
    lapply(seq_along(m), function(i) pad(m[[i]], widths[i])),
    sep = "  "
  ))
  
  paste(c(header_line, body_lines), collapse = "\n")
}

#' Create a phyloseq object from OTU, taxonomy, and sample data
#'
#' Constructs a \code{phyloseq} object from three input tables: an OTU
#' abundance table, a taxonomy table, and a sample metadata table. Input
#' objects are read using \code{\link{read_input}} and may therefore be
#' provided either as CSV files or data frames.
#'
#' The OTU table is transformed into a taxa-by-sample matrix. If a
#' \code{replicate} column is present, it is combined with \code{sample}
#' to create a unique sample identifier. If a \code{cluster_id} column is
#' present, reads are aggregated at the cluster level; otherwise, reads are
#' aggregated at the ASV level.
#'
#' Missing abundance values are replaced by zero. The taxonomy table is
#' converted to a matrix containing the standard taxonomic ranks from domain
#' to species. Sample metadata are reduced to one row per sample, using the
#' first occurrence when multiple rows are present.
#'
#' @param otu OTU/ASV abundance data. Can be a data frame or CSV file path.
#' The table must contain \code{asv_id},
#'   \code{sample}, and \code{read_count} columns. An optional \code{asv},
#'   \code{replicate}, or \code{cluster_id} column may also be present.
#' @param tax Taxonomy data. Can be a data frame or CSV file path.
#'  The table must contain an \code{asv_id} column
#'   and the taxonomic ranks \code{domain}, \code{kingdom} (optional), \code{phylum},
#'   \code{class}, \code{order}, \code{genus}, and \code{species}.
#' @param samples Sample metadata (optional).  Can be a data frame or CSV file path.
#'  The table must contain a \code{sample} column,
#'   which is used as the row names of the sample metadata.
#' @param sep Character used as the field separator when reading input files.
#'
#' @return A \code{phyloseq} object containing:
#' \itemize{
#'   \item an OTU table with taxa as rows and samples as columns;
#'   \item a taxonomy table containing domain through species classifications;
#'   \item sample metadata with one row per sample.
#' }
#'
#' @details
#' When \code{replicate} is present in the OTU table, the sample identifier is
#' constructed as \code{sample-replicate}. When \code{cluster_id} is present,
#' \code{asv_id} and \code{asv} are removed and read counts are summed by
#' cluster and sample. Otherwise, read counts are summed by ASV and sample.
#'
#' The resulting OTU matrix is converted to a matrix with taxa as rows, as
#' required by \code{phyloseq::otu_table()} with \code{taxa_are_rows = TRUE}.
#'
#' @examples
#' \dontrun{
#' phy_object <- make_phyloseq(
#'   otu = "otu.csv",
#'   tax = "taxonomy.csv",
#'   samples = "samples.csv"
#' )
#' 
#' phy_object
#' }
#'
#' @export
make_phyloseq <- function(otu,  tax, samples = NULL, sep = ","){
  
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
  
  ### otu ########################################
  
  otu_mat <- read_input(otu, sep = sep)
  # make one column with sample-replicate
  if("replicate" %in% colnames(otu_mat)){
    otu_mat <- otu_mat  %>%
      mutate(sample = paste(sample, replicate, sep="-")) %>%
      select(-replicate)
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
  
  otu_mat <- pivot_wider(otu_mat, 
                         names_from = sample,
                         values_from = read_count)
  otu_mat[is.na(otu_mat)] <- 0
  otu_mat <- as.data.frame(otu_mat)
  rownames(otu_mat) <- otu_mat$asv_id
  otu_mat <- select(otu_mat, -asv_id)
  otu_mat <- as.matrix(otu_mat)
  
  ### tax ########################################
  
  tax_mat <- read_input(tax, sep = sep)
  tax_mat <- as.data.frame(tax_mat)
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
    samples <- read_input(samples, sep = sep) %>%
      group_by(sample) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    sample_df <- as.data.frame(samples)
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
  return(phy_object)
}


######################################################

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


count_fastq_records <- function(fastq, chunk_lines = 4e6) {
  
  con <- open_any(fastq, "rt")
  on.exit(close(con))
  total_lines <- 0L
  while (length(lines <- readLines(con, n = chunk_lines)) > 0) {
    total_lines <- total_lines + length(lines)
  }
  if (total_lines %% 4 != 0) {
    warning(sprintf("%s: line count (%d) not a multiple of 4 - malformed FASTQ?",
                    fastq, total_lines), call. = FALSE)
  }
  total_lines %/% 4L
}

#### will be replaced by count_reads
count_fasta_headers <- function(fasta, chunk_lines = 1e6) {
  con <- open_any(fasta, "rt")
  on.exit(close(con))
  total <- 0L
  while (length(lines <- readLines(con, n = chunk_lines)) > 0) {
    total <- total + sum(startsWith(lines, ">"))
  }
  total
}


###########################################
## copy-whole-file helper (used when n >= total)
#.copy_fastq <- function(fastq, outfile) fast_copy(fastq, outfile)

#' Random-subsample a single FASTQ file
#' ### replaced by random_sample_fastq_pair
random_sample_fastq <- function(fastq, outfile, n = 1e6, randseed = NULL,
                                quiet = TRUE, chunk_records = 2.5e5,
                                pigz_path="pigz", compress_method = "R",
                                num_threads = 0) {
  s <- sample_fastq_indices(fastq, n = n, randseed = randseed, quiet = quiet)
  if (is.null(s$keep_idx)) {
    warning(sprintf("WARNING: %s contains %d records.\nThe input file is copied to output.",
                    fastq, s$total), call. = FALSE)
    fast_copy(fastq, outfile,  
              pigz_path=pigz_path, 
              compress_method = compress_method,
              num_threads = num_threads)
    return(invisible(s$total))
  }
  keep_lookup <- logical(s$total); keep_lookup[s$keep_idx] <- TRUE
  if (!quiet) cat("Extracting sampled records.\n")
  extract_fastq_records(fastq, outfile, keep_lookup, s$total, chunk_records)
  invisible(n)
}


########################################################
count_reads_original <- function(file, file_type="fastq", chunk_lines=1e5, fast_count = TRUE){
  
  if (endsWith(file, ".zip")) {
    stop("File compression type is not supported.")
  }
  
  if(is_linux() & fast_count ){
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
    print("WARNING: This command is quicker on linux-like systems using fast_count == TRUE.")
    
    if(file_type == "fasta"){ # can deal with compressed and uncompressed files
      
      con <- open_any(file, "rt")
      on.exit(close(con))
      seq_count <- 0L
      while (length(lines <- readLines(con, n = chunk_lines)) > 0) {
        seq_count <- seq_count + sum(startsWith(lines, ">"))
      }
    }else { # fastq and others
      
      con <- open_any(file, "rt")
      on.exit(close(con))
      total_lines <- 0L
      while (length(lines <- readLines(con, n = chunk_lines)) > 0) {
        total_lines <- total_lines + length(lines)
      }
      if (total_lines %% 4 != 0) {
        warning(sprintf("%s: line count (%d) not a multiple of 4 - malformed FASTQ?",
                        file, total_lines), call. = FALSE)
      }
      
      if(file_type == "fastq"){
        seq_count <- total_lines %/% 4L
      }else{
        seq_count <- total_lines
        msg <- paste(file_type, "is neither fasta nor fastq. 
                     The number of lines in file will be returned for", file)
        print(msg)
      }
      
    }
    return(seq_count)
  } # end non-linux-like
}
