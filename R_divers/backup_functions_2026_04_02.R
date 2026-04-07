

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
