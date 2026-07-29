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