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
#' @param min_id Real: [0-1] Bellow this identity do not align asv
#' @param vsearch_path Character string: path to vsearch executable.
#' @param num_threads Positive integer: Number of CPUs. If 0, use all available CPUs.
#' @param outfile Character string: csv file name to print the output data frame.
#' If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with asv pairs and their percentage of identity. 
#' asv pairs with identity bellow min_id are not listed. Colums: query, target, identity 
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
  # Build argument vector
  args <- c(
    "--allpairs_global", fasta,
    "--userout", vsearch_out,
#    "--userfields", '"query+target+id+ids+alnlen+aln"',
    "--userfields", "query+target+id+ids+alnlen+aln",
    "--id", min_id,
    "--threads", num_threads
  )
  
  run_system2(vsearch_path, args, quiet=quiet)
  
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
#' @param num_threads Positive integer: Number of CPUs.
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



#' Cluster all input ASV by cluster_size function of Vsearch
#' 
#' Cluster all input ASV and return a data frame with asv_id, cluster_id
#' columns
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, asv, read_count.
#' @param identity Real;[0-1] Identity threshold for clustering.
#' @param vsearch_path Character string: path to vsearch executable. 
#' @param num_threads Positive integer: Number of CPUs.
#' @param outfile Character string: name of the output csv file with asv_id, 
#' and cluster_id columns
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the following columns: asv_id, cluster_id
#' @examples 
#' \dontrun{
#' cluster_df <- GetClusterIdVsearch(read_count_df, 
#'                                       identity=0.97,
#'                                       vsearch_path="vsearch",
#'                                       num_threads=8)
#' }
#' @export
GetClusterIdVsearch <- function(read_count, 
                          identity=0.97, 
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
  
  # if called from a ClusterASV, without specifying the path, it has a "" value
  if(vsearch_path==""){
    vsearch_path<- "vsearch"
  }
  
  #####
  # make a df with unique asv, asv_id and readcount (sum)
  asv_df <- read_count_df %>%
    group_by(asv, asv_id) %>%
    summarize(read_count = sum(read_count), .groups="drop")%>%
    ungroup() %>%
    arrange(desc(read_count))
  
  ##### make tmp dir and input files
  tmp_dir <-paste('tmp_cluster_size_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  tmp_dir <- file.path(tempdir(), tmp_dir)
  check_dir(tmp_dir)
  
  # make a fasta file with unique asv format adapted to swarm
  input_fas <- file.path(tmp_dir, "input.fasta")
  # make fasta file with abundances
  write_fasta_rc(asv_df, input_fas)
  # Outfile name. Each line is a cluster, with asv_ids separated  by space
  outfile <- file.path(tmp_dir, "out_cluster_size.txt")
  
  # Build argument vector
  args <- c(
    "--cluster_size", input_fas,
    "--blast6out", outfile,
    "--id", identity
  )
  if(num_threads > 0){
    args <- append(args, c("--threads", num_threads))
  }
  
  run_system2(vsearch_path, args, quiet=quiet)
  
  file_info <- file.info(outfile)
  if(file_info$size == 0){ # No output of clustering => cluster_id is the same as asv_id
    cluster_df <- asv_df %>%
      select(asv_id) %>%
      mutate(cluster_id = asv_id)
  } else{
    # read clustering results to df_centroids
    # merged_id (include to a cluster), centroid_id (most abundant asv_id of the cluster)
    # if sequence is not in a cluster or if it is a centroid, asv not in blast6_file
    cluster_df1 <- read_blast6out(outfile)
    cluster_df1 <- cluster_df1 %>%
      rename(cluster_id = centroid_id) %>%
      rename(asv_id = merged_id)
    cluster_df <- left_join(asv_df, cluster_df1, by= "asv_id") %>%
      select(asv_id, cluster_id) %>%
      mutate(cluster_id = ifelse(is.na(cluster_id), asv_id, cluster_id))
  }
  
  unlink(tmp_dir, recursive = TRUE)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(cluster_df, file = outfile,  row.names = F, sep=sep)
  }
  return(cluster_df)
}


#' Make a density plot: pairwise percent identities between ASVs within and across 
#' clusters using different Swarm's d for clustering.
#' 
#' Cluster by swarm all ASV with a range of d values.
#' For each clustering (d), make a density plot of pairwise percentage of 
#' identities between ASVs, using different colors for identities between ASV
#' of the same or different clusters.
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
#' @param outfile Character string: output csv file name with the following columns:
#' pairwise_asv_identity, cluster, clustering_parameter; If empty, no file is written.
#' @param plotfile Character string: png file name for the output plot; 
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
#'                                       plotfile="density_plot.png")
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
                                          plotfile="", 
                                          sep=",", 
                                          quiet=TRUE
){
  
  
  ##### make df if read_count is file
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  #####
  # make pairwise_id df with query, target, identity
  if(!quiet){
    print("Calculating pairwise identities")
  }
  pairwise_id <- PairwiseIdentity(read_count_df, min_id = 0.8, vsearch_path=vsearch_path, quiet=TRUE, num_threads=0)
  
  
  #####
  # initialize data frame (cluster: same/different)
  pairwise_id_final <- data.frame(pairwise_asv_identity = numeric(),
                              cluster = character(),
                              clustering_parameter= factor())
  
  #####
  # for each d
  for(d in seq(swarm_d_min, swarm_d_max, by=swarm_d_increment)){
    
    if(!quiet){
      cat("Running swarm with d=", d, "\n")
    }
    #### Run swarm 
    # cluster_df : avs_id, clsuter_id
    cluster_df <- GetClusterIdSwarm(read_count_df, swarm_d=d,fastidious=FALSE,
                  swarm_path=swarm_path, num_threads=num_threads, quiet=quiet)
    #####
    # add to pairwise_id the clusters of each query and target and define 
    # if they are in the same or different clusters
    pairwise_id_local <- left_join(pairwise_id, cluster_df, by=c("query"="asv_id")) %>%
      rename(cluster_id_query = cluster_id)
    pairwise_id_local <- left_join(pairwise_id_local, cluster_df, by=c("target"="asv_id")) %>%
      rename(cluster_id_target = cluster_id) 
    
    # add same/different column, and the d for each line
    pairwise_id_local <- pairwise_id_local %>%
      mutate(cluster = ifelse(cluster_id_query == cluster_id_target, "same", "different")) %>%
#      mutate(clustering_parameter=paste0("swarm's d: ",d)) %>%
      mutate(clustering_parameter=d) %>%
      select(pairwise_asv_identity=identity, cluster, clustering_parameter)
    # add lines to the overall df
    pairwise_id_final <- rbind(pairwise_id_final, pairwise_id_local)
  } # end for
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(pairwise_id_final, file = outfile,  row.names = F, sep=sep)
  }
  
  pairwise_id_final <- pairwise_id_final %>%
    mutate(clustering_parameter = paste0("Swarm d = ",clustering_parameter))
  ####
  # Make the plot
  # fix order of d
  pairwise_id_final$clustering_parameter <- factor(pairwise_id_final$clustering_parameter, levels = unique(pairwise_id_final$clustering_parameter))
  
  p <-ggplot(pairwise_id_final, aes(x = pairwise_asv_identity, fill = cluster, color = cluster)) +
    geom_density(adjust = 1.5, alpha = 0.4) +
    scale_x_continuous(limits = c(80, 100)) +
    facet_wrap(~clustering_parameter, scales = "free_y") +  # one subplot per d
    ggtitle("Pairwise Percent Identity Between ASVs \n Within and Across Clusters (by Swarm's d)") +
    xlab("Pairwise Percent Identity Between ASVs") +
    ylab("Density") +
    theme(plot.title = element_text(size=12, hjust=0.5))
  
  ### 

  
  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  return(p)
  
} # end function

#' Make a density plot: pairwise percent identities between ASVs within and across 
#' clusters using different clustering identity thresholds.
#' 
#' Cluster all ASVs using cluster_size algorithm of vsearch with a range of identity values.
#' 
#' For each clustering identity threshold, make a density plot of pairwise 
#' percentage of identities between ASVs, using different colors for identities 
#' between ASVs of the same or different clusters.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv, sample, replicate (optional), read_count.
#' @param identity_min Real; Minimum identity threshold between asv and centroid.
#' @param identity_max Real; Maximum identity threshold between asv and centroid.
#' @param identity_increment Real; Identity thresholds vary from identity_min 
#' to identity_max by identity_increment.
#' @param min_id Real;[0-1] Bellow this pairwise identity asv pairs are not aligned
#' and their identity is not plotted.
#' @param vsearch_path Character string: path to vsearch executable.
#' @param num_threads Positive integer: Number of CPUs.
#' @param outfile Character string: csv file name with the following columns: 
#' pairwise_asv_identity, cluster, clustering_parameter; If empty, no file is written.
#' @param plotfile Character string: png file name for the output plot; 
#' If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns A density plot of pairwise percentage of identities.
#' @examples 
#' \dontrun{
#' plot <- PairwiseIdentityPlotPerClusterIdentityThreshold(read_count_df, 
#'                                       identity_min=0.9, 
#'                                       identity_max=0.99,
#'                                       identity_increment=0.01,
#'                                       min_id = 0.8, 
#'                                       vsearch_path="vsearch", 
#'                                       num_threads=8,
#'                                       plotfile="density_plot.png")
#' }
#' @export
PairwiseIdentityPlotPerClusterIdentityThreshold <- function(read_count, 
                                          identity_min=0.9, 
                                          identity_max=0.99,
                                          identity_increment=0.01,
                                          min_id = 0.8, 
                                          vsearch_path="vsearch", 
                                          num_threads=0,
                                          outfile="", 
                                          plotfile="",
                                          sep=",", 
                                          quiet=TRUE){
  
  
  ##### make df if read_count is file
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  #####
  # make pairwise_id df with query, target, identity
  if(!quiet){
    print("Calculating pairwise identities")
  }
  pairwise_id <- PairwiseIdentity(read_count_df, min_id = 0.8, vsearch_path=vsearch_path, quiet=quiet, num_threads=0)
  
  #####
  # initialize data frame (cluster: same/different)
  pairwise_id_final <- data.frame(pairwise_asv_identity = numeric(),
                              cluster = character(),
                              clustering_parameter= factor())
  
  #####
  # for each  identity threshold
  for(d in seq(identity_min, identity_max, by=identity_increment)){
    
    if (!quiet) {
      cat("Running cluster_size of vsearch with identity threshold", d, "\n")
    }
    
    #####
    # run vsearch
    cluster_size <- GetClusterIdVsearch(read_count_df, 
                                                 identity=d, 
                                                 vsearch_path=vsearch_path, 
                                                 num_threads=num_threads, 
                                                 quiet=TRUE)
  
    #####
    # add to pairwise_id the clusters of each query and target and define 
    # if they are in the same or different clusters

    # Add cluster of query
    pairwise_id_local <- left_join(pairwise_id, cluster_size, by=c("query"="asv_id")) %>%
      rename(cluster_id_query = cluster_id)
    # Add cluster of target
    pairwise_id_local <- left_join(pairwise_id_local, cluster_size, by=c("target"="asv_id")) %>%
      rename(cluster_id_target = cluster_id) 

    # add within/between column, and the d for each line
    pairwise_id_local <- pairwise_id_local %>%
      mutate(cluster = ifelse(cluster_id_query == cluster_id_target, "same", "different")) %>%
#      mutate(clustering_parameter= paste0("Clustering with ", d, "identity")) %>%
      mutate(clustering_parameter= d) %>%
      select(pairwise_asv_identity=identity, cluster, clustering_parameter)
    # add lines to the overall df
    pairwise_id_final <- rbind(pairwise_id_final, pairwise_id_local)
  } # end for
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(pairwise_id_final, file = outfile,  row.names = F, sep=sep)
  }
  
  pairwise_id_final <- pairwise_id_final %>%
    mutate(clustering_parameter = paste0("Cluster = ",clustering_parameter))
  ####
  # Make the plot
  # fix order of d
  pairwise_id_final$clustering_parameter <- factor(pairwise_id_final$clustering_parameter, levels = unique(pairwise_id_final$clustering_parameter))

  
  p <-ggplot(pairwise_id_final, aes(x = pairwise_asv_identity, fill = cluster, color = cluster)) +
     geom_density(adjust = 1.5, alpha = 0.4) +
     scale_x_continuous(limits = c(80, 100)) +
     facet_wrap(~clustering_parameter, scales = "free_y") +  # one subplot per d
     xlab("Pairwise Percent Identity Between ASVs") +
     ylab("Density") +
     ggtitle("Pairwise Percent Identity Between ASVs \n Within and Across Clusters \n using different identity thresholds for clustering") +
    theme(plot.title = element_text(size=12, hjust=0.5)) 
  
  ### 

  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  return(p)
  
} # end function

#' Classify clusters based on taxonomic agreement among their ASVs
#'
#' Classifies each cluster according to the taxonomic assignment at a given level:
#' 
#' - closed: All ASVs in the cluster are assigned to the same taxon, 
#' and all ASVs of that taxon belong exclusively to this cluster.  
#' - open: All ASVs in the cluster are assigned to the same taxon, 
#' but some ASVs of that taxon are found in other clusters.  
#' - hybrid: the cluster contains ASVs assigned to more than one taxon.  
#'   
#' 
#' @param cluster Data frame or csv file with the following variables: 
#' asv_id, cluster_id
#' @param taxa Data frame or CSV file with the following variables:  
#' asv_id, one columns per taxonomic level. (e.g. domain, phylum, class, order, 
#' family, genus, species). 
#' Missing values within the lineage are not allowed. If the taxonomic assignment 
#' is incomplete due to low resolution, higher-level taxa should be recorded as `NA`.
#' @param taxlevels Character vector with names of the taxonomic levels 
#' @param outfile Character string: name of the output csv file
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the following columns: cluster_id, classification
#' for each taxonomic level
#' @examples 
#' \dontrun{
#' plot <- ClassifyClusters(read_count_df, 
#'                                       swarm_d=7, 
#'                                       vsearch_path="vsearch",
#'                                       num_threads=8)
#' }
#' @export
ClassifyClusters <- function(cluster, taxa, outfile="", sep=",", quiet=TRUE, 
                             taxlevels=c("domain", "phylum", "class", "order","family", "genus", "species")
                             ){
  
  ##### make df if read_count is file
  if(is.character(cluster)){
    cluster_df <- read.csv(cluster, header=T, sep=sep)
  }else{
    cluster_df <- cluster
  }
  cluster_df <- cluster_df %>%
    select(asv_id, cluster_id)
  
  ##### make df if read_count is file
  if(is.character(taxa)){
    taxa_df <- read.csv(taxa, header=T, sep=sep)
  }else{
    taxa_df <- taxa
  }
  # taxa_df <- asv_tax
  taxa_df <- taxa_df %>%
    select(asv_id, all_of(taxlevels))
  
  taxa_df <- left_join(taxa_df, cluster_df, by="asv_id")
  
  # Data frame to complete at each iteration, for each taxlevel
  class <- data.frame("cluster_id" = unique(taxa_df$cluster_id))
  
  for(t in 1:length(taxlevels)){
    tl = taxlevels[t]
    
    # select the taxlevel, delete NA
    taxa_df_tl <- taxa_df %>%
      rename(taxon = all_of(tl)) %>% # rename the actual taxlevel to taxon
      select(cluster_id, taxon) %>%
      filter(!is.na(taxon))

    # group_by cluster
    taxa_df_tl <- taxa_df_tl %>%
      group_by(cluster_id) %>% # one line per cluster
      summarise(unique_taxa = list(unique(taxon)), # list of unique values of the taxon
                .groups = "drop") %>%  
      rowwise() %>% # for each cluster
      mutate(
        n_unique_taxa = length(unique_taxa),
      )
    
    # of all taxa in a vector instead of lists
    taxa_simple <- unlist(taxa_df_tl$unique_taxa)

    # classify
    taxa_df_tl <- taxa_df_tl %>%
      mutate(
        classification = case_when(
          # --- Case 1: only one unique taxon
          n_unique_taxa == 1 ~ {
            taxon <- unique_taxa[[1]]
            n = sum(taxa_simple == taxon, na.rm=TRUE) # number of clusters where 
                                                      # the taxon is present
            if(n == 1) "closed" else "open"
            },
          TRUE ~ "hybrid"
          )
      ) %>%
      select(cluster_id, classification)
    
    new_name <- paste("classification", tl, sep = "_")
    class <- left_join(class, taxa_df_tl, by="cluster_id") %>%
    rename(!!new_name := classification)
    
  }
  return(class)
}



#' Plot Cluster Classification according to taxa
#'
#' Cluster ASV with using different clustering parameters, than classify each 
#' cluster at each cluster setting :
#'  
#' - closed: ASVs in the cluster are assigned to the same taxon, 
#' and all ASVs of that taxon belong exclusively to this cluster.  
#' - open: ASVs in the cluster are assigned to the same taxon, 
#' but some ASVs of that taxon are found in other clusters.  
#' - hybrid: the cluster contains ASVs assigned to more than one taxon.  
#' 
#' Plot the number of clusters in each class for each setting.
#'   
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, asv, read_count
#' @param taxa Data frame or CSV file with the following variables:  
#' asv_id, one columns per taxonomic level. (e.g. domain, phylum, class, order, 
#' family, genus, species). 
#' Missing values within the lineage are not allowed. If the taxonomic assignment 
#' is incomplete due to low resolution, higher-level taxa should be recorded as `NA`.
#' @param clustering_method character; swarm or vsearch
#' @param cluster_params numerical vector of either swarm's d to be used, or 
#' percentage of percentage of identity threshold for the cluster_size algorithm of vsearch.
#' @param vsearch_path Character string: path to vsearch executable.
#' @param swarm_path Character string: path to swarm executable.
#' @param taxlevels Character vector with names of the taxonomic levels to be classed and plotted
#' @param outfile Character string: name of the output csv file with the following 
#' columns: classification, number_of_clusters, taxlevel, clustering_parameter;
#' If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns A connected scatterplot of the number of clusters in different classes 
#' (open, closed, hybrid)
#' @examples 
#' \dontrun{
#' plot <- PlotClusterClasstification(read_count, 
#'                                    taxa,
#'                                    clustering_method = "swarm"
#'                                    cluster_params = c(2, 4, 6, 8, 10)
#'                                    swarm_path= swarm_path,
#'                                    taxlevels = c(species, genus),
#'                                    num_threads=8)
#' }
#' @export

PlotClusterClasstification <- function(read_count, taxa, 
                           clustering_method="swarm", 
                           cluster_params=c(2,4,6,8,10), 
                           vsearch_path="vsearch", 
                           swarm_path="swarm", 
                           taxlevels= c("species", "genus"),
                           outfile="",
                           sep= ",",
                           quiet = TRUE){
  
  ##### make df if read_count is file
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  read_count_df <- read_count_df %>%
    group_by(asv_id, asv) %>%
    summarize(read_count=sum(read_count), .groups="drop") 
  
  ##### make df if taxa is a file
  if(is.character(taxa)){
    taxa_df <- read.csv(taxa, header=T, sep=sep)
  }else{
    taxa_df <- taxa
  }
  # taxa_df <- asv_tax
  taxa_df <- taxa_df %>%
    select(asv_id, all_of(taxlevels))
  
  clusters_count <- data.frame(
    "classification" = character(),
    "number_of_clusters" = numeric(),
    "taxlevel" = character(),
    "clustering_parameter" = numeric())

  #### for each d or % id
  for(i in cluster_params){
    ### cluster
      if(clustering_method =="swarm"){
        cluster_df <- GetClusterIdSwarm(read_count_df, 
                                       swarm_d=i, 
                                       swarm_path=swarm_path, 
                                       num_threads=num_threads, 
                                       quiet=quiet)
      }else{
        cluster_df <- GetClusterIdVsearch(read_count_df, 
                                                   identity=i/100, 
                                                   vsearch_path=vsearch_path, 
                                                   num_threads=num_threads, 
                                                   quiet=quiet)
      }
      
      ### classify
      classification <- ClassifyClusters(cluster_df, taxa_df, quiet=quiet, 
                                          taxlevels=taxlevels)
      #### count the number of clusters in each class (open, closed, hybrid, NA)
      for(tl in taxlevels){
        col_name <- paste("classification", tl, sep="_")
  
          tmp <- classification %>%
            count(!!sym(col_name)) %>%
            rename("classification" = !!sym(col_name), "number_of_clusters"=n) %>%
            mutate("taxlevel" = tl, "clustering_parameter" = i) 
          
          clusters_count <- rbind(clusters_count, tmp)
      }
  } # end for each cluster_params
  
  clusters_count$taxlevel <- factor(clusters_count$taxlevel, levels = unique(clusters_count$taxlevel))
  
  if(clustering_method == "swarm"){
    plot_title = "Number of different cluster types according the the taxonomic assignment of ASVs"
    x_label = "Swarm's d used for clustering"
  }else{
    plot_title = "Number of different cluster types according the the taxonomic assignment of ASVs"
    x_label = "% identity threshold for cluster_size algorithm of vsearch"
  }
  p <- ggplot(clusters_count %>% filter(!is.na(classification)), aes(x=clustering_parameter, y = number_of_clusters, group=classification, color = classification)) +
    geom_line() +
    geom_point() +
    facet_wrap(~taxlevel, scales = "free_y") + # one subplot per tax level
    ggtitle("Number of different cluster types \n according the the taxonomic assignment of ASVs")  +
    ylab("Number of clusters") +
    xlab(x_label) +
    theme(plot.title = element_text(size=12, hjust=0.5)) 
  
  ### 
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(clusters_count, file = outfile,  row.names = F, sep=sep)
  }
  return(p)

}


#' Pool ASVs of the same cluster
#' 
#' Pools variants of the same cluster and sums read counts of the 
#' underlying ASVs.
#' 
#' @param read_count_df Data frame with the following variables: 
#' asv, sample, replicate (optional), read_count.
#' @param cluster_df Data frame with the following variables: 
#' asv_id, cluster_id.
#' @returns Data frame: same structure as the input (read_count_df), but ASVs of 
#' the same cluster pooled to one row.
#' @examples
#' \dontrun{
#' pool_by_cluster(read_count=read_count_df, cluster_df=cluster_df)
#' }
#' @export
#' 
pool_by_cluster <- function(read_count_df, 
                          cluster_df){
  
  ### check if one to one relationship between asv and asv_id
  t <- check_one_to_one_relationship(read_count_df)
  
  ### get unique list of asv and asv_id
  asv_unique <- read_count_df%>%
    select(asv_id, asv) %>%
    distinct()
  
  # replace asv_id by cluster_id in read_count_df
  read_count_df <- left_join(read_count_df, cluster_df,  by= "asv_id") %>%
    select(-asv_id, -asv)
  
  
  if("replicate" %in% colnames(read_count_df)){ # if replicates in input, keep replicates
    
    read_count_df <- read_count_df %>%
      group_by(cluster_id, sample, replicate) %>%
      summarize(read_count=sum(read_count), .groups="drop") %>%
      rename(asv_id=cluster_id)
    
  }else{ # if no replicate column
    read_count_df <- read_count_df %>%
      group_by(cluster_id, sample) %>%
      summarize(read_count=sum(read_count), .groups="drop") %>%
      rename(asv_id=cluster_id)
  }
  
  read_count_df <- left_join(read_count_df, asv_unique, by="asv_id")
  return(read_count_df)
}

#' Cluster ASVs by Swarm or Vsearch
#' 
#' This function runs Swarm or the cluster_size command of Vserach 
#' on the ASVs in the input data frame.
#' Each ASV is assigned to a cluster, and the function provides 
#' two possible output formats, controlled by the argument group.
#' 
#' If group = TRUE, the function aggregates the ASVs belonging to the same 
#' cluster into a single row. In this case, the asv_id and asv columns contain 
#' the identifier and the sequence of the cluster's centroid. 
#' read_count is summed over ASVs, samples and replicates are unchanged.
#' 
#' If group = FALSE, the function returns the original input data frame 
#' with an additional column: cluster_id. Each row still corresponds to one ASV.
#' 
#' The clustering can be run on the whole data set at once, 
#' or sample by sample (by_sample)
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv, sample, replicate (optional), read_count.
#' @param group Logical; If TRUE, ASVs of the same cluster are pooled to one row.
#' if FALSE, input Df is simply completed by a cluser_id column. 
#' @param by_sample Logical: run clustering separately for each sample.
#' @param method Character string: program to use for clustering: swarm or vsearch 
#' @param path Character string: path to swarm or vsearch executables. 
#' @param swarm_d Positive integer: d parameter for swarm (if applicable).
#' Maximum number of differences allowed between two ASVs, 
#' meaning that two ASVs will be grouped if they have d (or less) differences.
#' @param fastidious Logical: When clustering with swarm and working with d = 1, 
#' perform a second clustering pass to reduce the number of small clusters.
#' @param identity real: (0-1), the identity threshold used for
#' clustering when the cluster_size algorithm of vsearch is used.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param num_threads Positive integer: Number of CPUs.
#' @param sep Field separator character in input and output csv files.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns read_count_df: same structure as the input, but ASVs of 
#' the same cluster pooled to one row if group==TRUE, or rows kept as in the input 
#' and an additional cluster_id column is added if group==FALSE.
#' @examples
#' \dontrun{
#' read_count_df <- ClusterASV(read_count=read_count_df, group=TRUE,method="vsearch",
#' by_sample=TRUE, path=swarm_path, 
#' num_threads=4)
#' }
#' @export
#' 
ClusterASV <- function(read_count, 
                      group = TRUE,
                      by_sample=FALSE,
                      method = "swarm",
                      path="", 
                      num_threads=0, 
                      swarm_d=1, 
                      fastidious=T, 
                      identity = 0.97,
                      outfile="", 
                      sep=",", 
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
    if(group==FALSE){
      out_df$cluster_id <- numeric()
    }
    
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
      
      # get cluster_id for all ASV
      if(method == "swarm"){
        cluster_df <- GetClusterIdSwarm(df_sample, 
                                        swarm_d=swarm_d, 
                                        swarm_path=path, 
                                        num_threads=num_threads, 
                                        quiet=quiet)
      }else{
        cluster_df <- GetClusterIdVsearch(df_sample, 
                                          identity=identity, 
                                          vsearch_path=path, 
                                          num_threads=num_threads, 
                                          quiet=quiet)
      }
      
      # modify input df according to group
      if(group){ # pool ASV by cluster
        df_sample <- pool_by_cluster(df_sample, cluster_df)
      }
      else{ # add cluster_id column
        df_sample <- left_join(df_sample, cluster_df, by="asv_id")
      }
      # add output of the sample to the total data frame    
      out_df <- rbind(out_df, df_sample)
    }
  }else{ # run swarm for all samples together
    # get cluster_id for all ASV
    if(method == "swarm"){
      cluster_df <- GetClusterIdSwarm(read_count_df, 
                                      swarm_d=swarm_d, 
                                      swarm_path=path, 
                                      num_threads=num_threads, 
                                      quiet=quiet)
    }else{
      cluster_df <- GetClusterIdVsearch(read_count_df, 
                                        identity=identity, 
                                        vsearch_path=path, 
                                        num_threads=num_threads, 
                                        quiet=quiet)
    }
    
    if(group){# pool ASV by cluster
      out_df <- pool_by_cluster(read_count_df, cluster_df)
    }
    else{  # add cluster_id column
      out_df <- left_join(read_count_df, cluster_df, by="asv_id")
    }
    
  }
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(out_df, file = outfile,  row.names = F, sep=sep)
  }
  return(out_df)
}


