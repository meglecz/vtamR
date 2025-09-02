#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first if_else
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym := !!
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL

#' Pairwise identity
#' 
#' Align all pairs of asv with at least min_id similarity and
#' make a data frame with pairs of asv and their percentage of identity. 
#' asv pairs with identity bellow min_id are not listed.
#' 
#' @param asv Data frame or csv file with asv and asv_id columns.
#' @param min_id Real: Bellow this percentage of identity do not align asv
#' @param vsearch_path Character string: path to vsearch executables.
#' @param num_threads Positive integer: Number of CPUs. If 0, use all available CPUs.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with asv pairs and their percentage of identity. 
#' asv pairs with identity bellow min_id are not listed.
#' @examples 
#' \dontrun{
#' identity_df <- PairwiseIdentity(asv, 
#'                                 min_id = 0.8, 
#'                                 vsearch_path=vsearch, 
#'                                 num_threads=8)
#' }
#' @export
#' 
PairwiseIdentity <- function(asv, 
                             min_id = 0.8, 
                             vsearch_path=vsearch, 
                             num_threads=0,
                             outfile="",
                             sep=",",
                             quiet=TRUE
                             ){
  
  # can accept df or file as an input
  if(is.character(asv)){
    # read known occurrences
    asv_df <- read.csv(asv, header=T, sep=sep)
  }else{
    asv_df <- asv
  }
  
  # change to percent
  min_id_perc <- min_id * 100 
  
  #### get unique asv list and make fasta file
  asv_df <- asv_df %>%
    select(asv, asv_id) %>%
    distinct()
  t <- check_one_to_one_relationship(asv_df) # stop execution, if FALSE
  
  fasta <- file.path(tempdir(), "asv.fasta")
  write_fasta_df(asv_df, outfile=fasta, read_count=FALSE)
  
  ### run vsearch allpairs_global
  vsearch_out <- file.path(tempdir(), "vsearch_out.tsv")
  cmd <- paste(vsearch_path, "--allpairs_global", fasta, "--userout", vsearch_out, 
               '--userfields "query+target+id+ids+alnlen+aln"', 
               "--id", min_id, "--threads", num_threads, sep=" ")
  if(!quiet){
    print(cmd)
  }
  system(cmd, show.output.on.console = FALSE)
  
  #### read vsearch result
  if(file.exists(vsearch_out) && file.size(vsearch_out) > 0){
    # read vsearch results
    results_vsearch<- read.csv(vsearch_out, header = FALSE, sep="\t")
    colnames(results_vsearch) <- c("query","target","identity","nb_ids","aln_length","aln")
    # none of the values easily outputted by vsearch take into the external gaps as a diff 
    # => correct this, based on the alnlen and the number of identities
    # in a correctly filtered data set this correction should not make a big difference, but it does 
    # in unfiltered data sets, with strongly variable asv length.
    results_vsearch$id_bis <- round(results_vsearch$nb_ids / nchar(results_vsearch$aln) * 100, digits=1)
    
    results_vsearch <- results_vsearch %>%
      select(query, target, identity=id_bis) %>%
      filter(identity >= min_id_perc)
    results_vsearch$query <- as.numeric(results_vsearch$query)
    results_vsearch$target <- as.numeric(results_vsearch$target)
  } else {
    results_vsearch <- data.frame(query=numeric(),
                                  target= numeric(),
                                  identity= numeric())
  }
  
  unlink(fasta)
  unlink(vsearch_out)
  
  #### write csv
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(out_df, file = outfile,  row.names = F, sep=sep)
  }
  return(results_vsearch)
}


#' Pairwise identity density plot for different Swarm d values
#' 
#' Cluster by swarm all asv with a range of d values.
#' For each d, make a density plot of pairwise percentage of identities between 
#' asv, using different colors for identities between asv of the same or different
#' clusters.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv, sample, replicate (optional), read_count.
#' @param swarm_d_min Positive integer: Minimum value of d for Swarm.
#' @param swarm_d_max Positive integer: Maximum value of d for Swarm.
#' @param swarm_d_increment Positive integer: increase d by swarm_d_increment between
#' swarm_d_min and swarm_d_max.
#' @param min_id Real: Bellow this percentage of identity asv pairs are not aligned
#' and their identity is not plotted.
#' @param vsearch_path Character string: path to vsearch executable.
#' @param swarm_path Character string: path to swarm executable. 
#' @param num_threads Positive integer: Number of CPUs.
#' @param outfile Character string: png file name to print the output plot. 
#' If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns A density plot pairwise percentage of identities.
#' @examples 
#' \dontrun{
#' plot <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
#'                                       swarm_d_min=2, 
#'                                       swarm_d_max=12,
#'                                       swarm_d_increment=2,
#'                                       min_id = 0.8, 
#'                                       vsearch_path="vsearch", 
#'                                       swarm_path="swarm",
#'                                       num_threads=8,
#'                                       outfile="density_plot.png")
#' }
#' @export
PairwiseIdentityPlotPerSwarmD <- function(read_count, 
                                          swarm_d_min=1, 
                                          swarm_d_max=15,
                                          swarm_d_increment=1,
                                          min_id = 0.8, 
                                          vsearch_path="vsearch", 
                                          swarm_path="swarm",
                                          num_threads=0,
                                          outfile="", 
                                          sep=",", 
                                          quiet=TRUE
){
  
  
  ##### make df if read_count is file
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # system2 cannot use ~ as a home
  vsearch_path <- path.expand(vsearch_path)
  swarm_path <- path.expand(swarm_path)
  
  #####
  # make pairwise_id df with query, target, identity
  if(!quiet){
    print("Calculating pairwise identities")
  }
  pairwise_id <- PairwiseIdentity(read_count_df, min_id = 0.8, vsearch_path=vsearch_path, quiet=TRUE, num_threads=0)
  
  #####
  # make a df with unique asv, asv_id and readcount (sum)
  asv_df <- read_count_df %>%
    group_by(asv, asv_id) %>%
    summarize(rc = sum(read_count), .groups="drop")
  
  #####
  # make a fasta file with unique asv format adapted to swarm
  input_swarm <- file.path(tempdir(), "swarm_input.fasta")
  writeLines(paste(">", asv_df$asv_id, "_", 
                   asv_df$rc, "\n", 
                   asv_df$asv, 
                   sep=""), 
             input_swarm)
  
  #####
  # initialize data frame (cluster: same/different)
  pairwise_id_d <- data.frame(identity = numeric(),
                              cluster = character(),
                              swarm_d= factor())
  
  #####
  # for each d
  for(d in seq(swarm_d_min, swarm_d_max, by=swarm_d_increment)){
    if(!quiet){
    cat("Running swarm with d=", d)
    }
    
    #####
    # run swarm
    
    # make tmp dir separately for each d
    tmp_dir <-paste('tmp_swarm_', d, '_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    tmp_dir <- file.path(tempdir(), tmp_dir)
    check_dir(tmp_dir)
    
    # clusters.txt each line is a cluster, with asv_ids separated  by space
    clusters <- file.path(tmp_dir, "clusters.txt")
    
    # Build argument vector
    args <- c(
      "-d", d,
      "-o", clusters,
      input_swarm
    )
    if(num_threads > 0){
      args <- append(args, c("-t", num_threads), after = 2)
    }
    
    if (!quiet) {
      # Show the full command that will be run
      cat("Running command:\n")
      cat(swarm_path, paste(shQuote(args), collapse = " "), "\n")
      
      system2(
        command = swarm_path,
        args = args,
        stdout = "",
        stderr = ""
      )
    }else{
      output <- suppressWarnings(system2(
        command = swarm_path,
        args = args,
        stdout = TRUE,
        stderr = TRUE
      ))
      
      # Extract only error/warning lines
      errors_only <- grep("error|warning|fail", output, ignore.case = TRUE, value = TRUE)
      cat(errors_only, sep = "\n")
    }
    

    #####
    # make a data frame with cluster_id and asv_id columns, 
    # where asv_id has all swarm input asv_id, 
    # and  cluster_id is the name of the cluster they belong to
    
    # read.table is unpredictable, when the number of element is variable among lines. Use a more complicated, but more sure solution.
    # Read the file line by line
    lines <- readLines(clusters)
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
    
    #  repeat each cluster_id as many times as columns
    # transpose and flatten to make asv_id
    # This will produce some lines with NA for asv_id. Filter them out afterwards.
    cluster_df <- data.frame(cluster_id = rep(cluster_df$V1, 
                                              each = ncol(cluster_df)),
                             asv_id = as.vector(t(cluster_df[,]))) 
    # delete lines with NA values in asv_id
    cluster_df <- cluster_df %>%
      filter(!is.na(asv_id))
    # delete read_counts from id
    cluster_df$cluster_id <- as.numeric(
      sub("_[0-9]+", "", cluster_df$cluster_id ))
    cluster_df$asv_id <- as.numeric(
      sub("_[0-9]+", "", cluster_df$asv_id ))
    
    #####
    # add to pairwise_id the clusters of each query and target and define 
    # if they are in the same or different clusters
    pairwise_id_local <- left_join(pairwise_id, cluster_df, by=c("query"="asv_id")) %>%
      rename(cluster_id_query = cluster_id)
    pairwise_id_local <- left_join(pairwise_id_local, cluster_df, by=c("target"="asv_id")) %>%
      rename(cluster_id_target = cluster_id) 
    
    # add within/between column, and the d for each line
    pairwise_id_local <- pairwise_id_local %>%
      mutate(cluster = ifelse(cluster_id_query == cluster_id_target, "same", "different")) %>%
      mutate(swarm_d=paste0("swarm's d: ",d)) %>%
      select(identity, cluster, swarm_d)
    # add lines to the overall df
    pairwise_id_d <- rbind(pairwise_id_d, pairwise_id_local)
    
    unlink(tmp_dir, recursive = TRUE)
  } # end for
  unlink(input_swarm)
  
  ####
  # Make the plot
  # fix order of d
  pairwise_id_d$swarm_d <- factor(pairwise_id_d$swarm_d, levels = unique(pairwise_id_d$swarm_d))
  
  p <-ggplot(pairwise_id_d, aes(x = identity, fill = cluster, color = cluster)) +
    geom_density(adjust = 1.5, alpha = 0.4) +
    scale_x_continuous(limits = c(80, 100)) +
    facet_wrap(~swarm_d, scales = "free_y") +  # one subplot per d
    ggtitle("Pairwise percent identity between ASV of the \n same or different cluster according Swarm's d") +
    theme(plot.title = element_text(size=12, hjust=0.5)) 
  
  
  ### 
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    png(filename=outfile) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  return(p)
  
} # end function



#' Pairwise identity density plot for different % identity thresholds of clustering
#' 
#' Cluster all asv using cluster_size algorithm of vsearch with a range of identity values.
#' 
#' For each identity, make a density plot of pairwise percentage of identities between 
#' asv, using different colors for identities between asv of the same or different
#' clusters.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv, sample, replicate (optional), read_count.
#' @param identity_min Real; Minimum identity threshold between asv and centroid.
#' @param identity_max Real; Maximum identity threshold between asv and centroid.
#' @param identity_increment Real; Identity thresholds vary from identity_min 
#' to identity_max by identity_increment.
#' @param min_id Real: Bellow this percentage of identity asv pairs are not aligned
#' and their identity is not plotted.
#' @param vsearch_path Character string: path to vsearch executable.
#' @param num_threads Positive integer: Number of CPUs.
#' @param outfile Character string: png file name to print the output plot. 
#' If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns A density plot of pairwise percentage of identities.
#' @examples 
#' \dontrun{
#' plot <- PairwiseIdentityPlotPerClusterIdThreshold(read_count_df, 
#'                                       identity_min=0.9, 
#'                                       identity_max=0.99,
#'                                       identity_increment=0.01,
#'                                       min_id = 0.8, 
#'                                       vsearch_path="vsearch", 
#'                                       num_threads=8,
#'                                       outfile="density_plot.png")
#' }
#' @export
PairwiseIdentityPlotPerClusterIdThreshold <- function(read_count, 
                                          identity_min=0.9, 
                                          identity_max=0.99,
                                          identity_increment=0.01,
                                          min_id = 0.8, 
                                          vsearch_path="vsearch", 
                                          num_threads=0,
                                          outfile="", 
                                          sep=",", 
                                          quiet=TRUE){
  
  
  ##### make df if read_count is file
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # system2 cannot use ~ as a hone
  vsearch_path <- path.expand(vsearch_path)
  
  
  #####
  # make pairwise_id df with query, target, identity
  if(!quiet){
    print("Calculating pairwise identities")
  }
  pairwise_id <- PairwiseIdentity(read_count_df, min_id = 0.8, vsearch_path=vsearch_path, quiet=quiet, num_threads=0)
  
  #####
  ### get unique list of asv with read_count and asv_id
  asv_df <- read_count_df %>%
    group_by(asv, asv_id) %>%
    summarize(read_count = sum(read_count), .groups="drop")%>%
    arrange(desc(read_count))
  
  #####
  # make a fasta file with unique asv format adapted to cluster_size of vsearch
  input_cluster_size <- file.path(tempdir(), "input_cluster_size.fasta")
  # make fasta file with abundances
  write_fasta_rc(asv_df, input_cluster_size)
  
  #####
  # initialize data frame (cluster: same/different)
  pairwise_id_d <- data.frame(identity = numeric(),
                              cluster = character(),
                              identity_threshold= factor())
  
  #####
  # for each  identity threshold
  for(d in seq(identity_min, identity_max, by=identity_increment)){
    
    if (!quiet) {
      cat("Running cluster_size of vsearch with identity threshold", d)
    }
    
    #####
    # run vsearch
    
    # make tmp dir separately for each d
    tmp_dir <-paste('tmp_clustersize_', d, '_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    tmp_dir <- file.path(tempdir(), tmp_dir)
    check_dir(tmp_dir)
    
    # output of clustering column & and 2: "merged_id","centroid_id"
    blast6_file <- file.path(tmp_dir, "clusters.txt")
    
    # Build argument vector
    args <- c(
      "--cluster_size", input_cluster_size,
      "--blast6out", blast6_file,
      "--id", d
    )
    if(num_threads > 0){
      args <- append(args, c("--threads", num_threads))
    }
    
    # run cluster_size
    if (!quiet) {
      # Show the full command that will be run
      cat("Running command:\n")
      cat(vsearch_path, paste(shQuote(args), collapse = " "), "\n")
      
      system2(
        command = vsearch_path,
        args = args,
        stdout = "",
        stderr = ""
      )
    }else{
      output <- suppressWarnings(system2(
        command = vsearch_path,
        args = args,
        stdout = TRUE,
        stderr = TRUE
      ))
      # Extract only error/warning lines
      errors_only <- grep("error|warning|fail", output, ignore.case = TRUE, value = TRUE)
      cat(errors_only, sep = "\n")
    } # end run cluster_size
    
    file_info <- file.info(blast6_file)
    if(file_info$size == 0){ # No output of clustering
      # Delete the temp directory
      unlink(outdir_tmp, recursive = TRUE)
      df_centroids <- data.frame("merged_id"=as.numeric(),
                                 "cluster_id"=as.numeric())
    } else{
      # read clustering results to df_centroids
      # merged_id (include to a cluster), centroid_id (most abundant asv_id of the cluster)
      # if sequence is not in a cluster or if it is a centroid, asv not in blast6_file
      df_centroids <- read_blast6out(blast6_file)
      df_centroids <- df_centroids %>%
        rename(cluster_id = centroid_id)
    }
  
    #####
    # add to pairwise_id the clusters of each query and target and define 
    # if they are in the same or different clusters

    # Add cluster of query
    pairwise_id_local <- left_join(pairwise_id, df_centroids, by=c("query"="merged_id")) %>%
      rename(cluster_id_query = cluster_id)
    # Add cluster of target
    pairwise_id_local <- left_join(pairwise_id_local, df_centroids, by=c("target"="merged_id")) %>%
      rename(cluster_id_target = cluster_id) 
    # if cluster_id_query or cluster_id_target is NA, add query or target
    pairwise_id_local <- pairwise_id_local %>%
      mutate(cluster_id_query = if_else(is.na(cluster_id_query), query, cluster_id_query)) %>%
      mutate(cluster_id_target = if_else(is.na(cluster_id_target), target, cluster_id_target))
    
    # add within/between column, and the d for each line
    pairwise_id_local <- pairwise_id_local %>%
      mutate(cluster = ifelse(cluster_id_query == cluster_id_target, "same", "different")) %>%
      mutate(identity_threshold= d) %>%
      select(identity, cluster, identity_threshold)
    # add lines to the overall df
    pairwise_id_d <- rbind(pairwise_id_d, pairwise_id_local)
    
    unlink(tmp_dir, recursive = TRUE)
  } # end for
  unlink(input_cluster_size)
  
  ####
  # Make the plot
  # fix order of d
  pairwise_id_d$identity_threshold <- factor(pairwise_id_d$identity_threshold, levels = unique(pairwise_id_d$identity_threshold))

  
  p <-ggplot(pairwise_id_d, aes(x = identity, fill = cluster, color = cluster)) +
     geom_density(adjust = 1.5, alpha = 0.4) +
     scale_x_continuous(limits = c(80, 100)) +
     facet_wrap(~identity_threshold, scales = "free_y") +  # one subplot per d
     ggtitle("Pairwise percent identity between ASV of the \n same or different cluster according identity threshold of clustering") +
    theme(plot.title = element_text(size=12, hjust=0.5)) 
  
  ### 
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    png(filename=outfile) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  return(p)
  
} # end function




