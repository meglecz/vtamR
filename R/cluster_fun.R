#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first 
#' @importFrom dplyr if_else case_when 
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat
#' @importFrom ggplot2 scale_x_continuous facet_wrap ggtitle xlab ylab geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom grDevices dev.off png
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym := !!
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL

#' Pairwise Identity
#' 
#' Computes pairwise sequence identity between ASVs.
#' 
#' All ASV pairs with sequence similarity greater than or equal to `min_id`
#' are aligned, and their percentage identity is reported in a data frame.
#' Pairs with identity below `min_id` are not included in the output.
#' 
#' @param asv Data frame or CSV file containing `asv` and `asv_id` columns.
#' @param min_id Numeric value between 0 and 1. Minimum identity threshold;
#' ASV pairs with identity below this value are ignored.
#' @param vsearch_path Character string: path to the VSEARCH executable.
#' @param num_threads Positive integer: number of CPUs to use. If 0, all available CPUs are used.
#' @param outfile Character string: name of the output CSV file. If empty, no file is written.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: If TRUE, suppress informational messages and display only warnings or errors.
#' 
#' @return A data frame listing ASV pairs and their sequence identity, with columns:
#' `query`, `target`, `identity`.
#' 
#' @examples
#' \dontrun{
#' identity_df <- pairwise_identity(
#'   asv,
#'   min_id = 0.8,
#'   vsearch_path = vsearch,
#'   num_threads = 8
#' )
#' }
#' @export
#' 
pairwise_identity <- function(asv, 
                              min_id = 0.8, 
                              vsearch_path=vsearch, 
                              num_threads=0,
                              outfile="",
                              sep=",",
                              quiet=TRUE
){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  
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
  t <- check_one_to_one(asv_df) # stop execution, if FALSE
  
  fasta <- file.path(tempdir(), "asv.fasta")
  write_fasta_with_counts(asv_df, outfile=fasta, read_count=FALSE)
  
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


#' Cluster ASVs using Swarm
#' 
#' Clusters all input ASVs using the Swarm algorithm and returns a data frame
#' mapping each ASV to a cluster.
#' 
#' @param read_count Data frame or CSV file containing the following columns:
#' `asv_id`, `asv`, `read_count`.
#' @param swarm_d Positive integer: clustering distance parameter (`d`) for Swarm.
#' @param fastidious Logical: if TRUE and `swarm_d = 1`, performs a second clustering
#' pass to reduce the number of small clusters.
#' @param swarm_path Character string: path to the Swarm executable.
#' @param num_threads Positive integer: number of CPUs to use. If 0, all available CPUs are used.
#' @param outfile Character string: name of the output CSV file containing
#' `asv_id` and `cluster_id`. If empty, no file is written.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and display only warnings or errors.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#' 
#' @return A data frame with columns: `asv_id`, `cluster_id`.
#' 
#' @examples
#' \dontrun{
#' cluster_df <- cluster_swarm(
#'   read_count_df,
#'   swarm_d = 7,
#'   swarm_path = "swarm",
#'   num_threads = 8
#' )
#' }
#' @export
#' 
cluster_swarm <- function(read_count, 
                          swarm_d=1, 
                          fastidious=FALSE,
                          swarm_path="swarm", 
                          num_threads=0, 
                          outfile="", 
                          sep=",", 
                          quiet=TRUE,
                          log_file=NULL){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # get function name, all arguments and stat time
  if(!is.null(log_file)){
    log <- collect_log()
  }
  ##### make df if read_count is file
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  # if called from a cluster_asv, without specifying the path, it has a "" value
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
  out_uclust <- file.path(tmp_dir, "out_swarm_uclust.tsv")
  out_swarm <- file.path(tmp_dir, "out_swarm.txt")
  
  ##### run swarm
  # Build argument vector
  args <- c(
    "-d", swarm_d,
    "-u", out_uclust,
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
  #1. Record type: S, H, or C.
  #2. Cluster number (ID) (zero-based).
  #3. Centroid length (S), query length (H), or cluster size (C).
  #4. Percentage of similarity with the centroid sequence (H), or set to ’*’ (S, C).
  #5. Match orientation + or - (H), or set to ’*’ (S, C).
  #6. Not used, always set to ’*’ (S, C) or to zero (H).
  #7. Not used, always set to ’*’ (S, C) or to zero (H).
  #8. set to ’*’ (S, C) or, for H, compact representation of the pairwise alignment using
  #the CIGAR format (Compact Idiosyncratic Gapped Alignment Report): M
  #(match), D (deletion) and I (insertion). The equal sign ’=’ indicates that the query
  #is identical to the centroid sequence.
  #9. Header of the query sequence (H), or of the centroid sequence (S, C).
  #10. Header of the centroid sequence (H), or set to ’*’ (S, C)
  
  cluster_df <- read.table(out_uclust, header=FALSE, sep="\t")
  cluster_df <- cluster_df[,c(1,9,10)]
  colnames(cluster_df) <- c("record_type", "asv_id", "cluster_id")
  
  hits <- cluster_df %>%
    filter(record_type == "H") %>%
    select(-record_type) %>%
    mutate(asv_id = as.numeric(sub("_[0-9]+$", "", asv_id))) %>%
    mutate(cluster_id = as.numeric(sub("_[0-9]+$", "", cluster_id)))
  
  centroids <- cluster_df %>%
    filter(record_type == "C") %>%
    select(-record_type) %>%
    mutate(asv_id = as.numeric(sub("_[0-9]+$", "", asv_id))) %>%
    mutate(cluster_id = asv_id)
  
  
  cluster_df <- rbind(hits, centroids) %>%
    arrange(cluster_id, asv_id)
  
  unlink(tmp_dir, recursive = TRUE)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(cluster_df, file = outfile,  row.names = F, sep=sep)
  }
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(cluster_df)
}


#' Cluster ASVs using VSEARCH
#' 
#' Clusters all input ASVs using the VSEARCH `--cluster_size` method and
#' returns a data frame mapping each ASV to a cluster identifier.
#' 
#' @param read_count Data frame or CSV file containing the following columns:
#' `asv_id`, `asv`, `read_count`.
#' @param identity Numeric value between 0 and 1. Sequence identity threshold used for clustering.
#' @param vsearch_path Character string: path to the VSEARCH executable.
#' @param num_threads Positive integer: number of CPUs to use. If 0, all available CPUs are used.
#' @param outfile Character string: name of the output CSV file containing
#' `asv_id` and `cluster_id`. If empty, no file is written.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and display only warnings or errors.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#' 
#' @return A data frame with columns: `asv_id`, `cluster_id`.
#' 
#' @examples
#' \dontrun{
#' cluster_df <- cluster_vsearch(
#'   read_count_df,
#'   identity = 0.97,
#'   vsearch_path = "vsearch",
#'   num_threads = 8
#' )
#' }
#' @export
#' 
cluster_vsearch <- function(read_count, 
                            identity=0.97, 
                            vsearch_path="vsearch", 
                            num_threads=0, 
                            outfile="", 
                            sep=",", 
                            quiet=TRUE,
                            log_file=NULL){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # get function name, all arguments and stat time
  if(!is.null(log_file)){
    log <- collect_log()
  }
  ##### make df if read_count is file
  if(is.character(read_count)){
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # if called from a cluster_asv, without specifying the path, it has a "" value
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
  write_fasta_with_counts(asv_df, input_fas, read_count = TRUE)
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
    cluster_df1 <- read_vsearch_cluster_size_outmft6(outfile)
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
  
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(cluster_df)
}

#' Density plot of pairwise ASV identities across Swarm clusters
#' 
#' Computes pairwise sequence identities between ASVs clustered using Swarm
#' across a range of `d` values, and generates density plots comparing identity
#' distributions within clusters and between clusters.
#' 
#' For each Swarm clustering result, pairwise ASV identities are computed.
#' Identities below `min_id` are not included in the analysis or plot.
#' Density distributions are then generated for within-cluster and between-cluster
#' comparisons, using different colors.
#' 
#' @param read_count Data frame or CSV file containing the following columns:
#' `asv`, `sample`, `replicate` (optional), `read_count`.
#' @param swarm_d_min Positive integer: minimum Swarm `d` value to test.
#' @param swarm_d_max Positive integer: maximum Swarm `d` value to test.
#' @param swarm_d_increment Positive integer: increment step between `d` values.
#' @param min_id Numeric value between 0 and 1. Minimum identity threshold;
#' ASV pairs below this value are excluded from the analysis.
#' @param vsearch_path Character string: path to the VSEARCH executable.
#' @param swarm_path Character string: path to the Swarm executable.
#' @param num_threads Positive integer: number of CPUs to use. If 0, all available CPUs are used.
#' @param outfile Character string: name of the output CSV file containing
#' pairwise ASV identities and clustering metadata. If empty, no file is written.
#' @param plotfile Character string: name of the output PNG file for the density plot.
#' If empty, no file is written.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and display only warnings or errors.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#' 
#' @return A density plot object showing distributions of pairwise ASV identities
#' within and between Swarm clusters across tested `d` values.
#' 
#' @examples
#' \dontrun{
#' plot <- plot_pairwise_identity_swarm(
#'   read_count_df,
#'   swarm_d_min = 2,
#'   swarm_d_max = 12,
#'   swarm_d_increment = 2,
#'   min_id = 0.8,
#'   vsearch_path = "vsearch",
#'   swarm_path = "swarm",
#'   num_threads = 8,
#'   plotfile = "density_plot.png"
#' )
#' }
#' @export
#' 
plot_pairwise_identity_swarm <- function(read_count, 
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
                                         quiet=TRUE,
                                         log_file=NULL
){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # get function name, all arguments and stat time
  if(!is.null(log_file)){
    log <- collect_log()
  }
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
  pairwise_id <- pairwise_identity(read_count_df, min_id = 0.8, vsearch_path=vsearch_path, quiet=TRUE, num_threads=0)
  
  
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
    cluster_df <- cluster_swarm(read_count= read_count_df, 
                                swarm_d=d,
                                fastidious=FALSE,
                                swarm_path=swarm_path, 
                                num_threads=num_threads, 
                                quiet=quiet)
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
    theme(plot.title = element_text(size=12, hjust=0.5),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 8)
    )
  
  ### 
  
  
  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile, width = 2000, height = 1500, res = 300)
    print(p) # print plot to file
    dev.off()
  }
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(p)
  
} # end function

#' Density plot of pairwise ASV identities across VSEARCH clustering thresholds
#' 
#' Generates density plots of pairwise ASV sequence identities within and
#' between clusters obtained using VSEARCH clustering across a range of
#' identity thresholds.
#' 
#' ASVs are clustered using the VSEARCH `--cluster_size` algorithm for each
#' identity threshold between `identity_min` and `identity_max`.
#' For each clustering result, pairwise ASV identities are computed and
#' visualized as density distributions, separating within-cluster and
#' between-cluster comparisons.
#' 
#' Identity values below `min_id` are excluded from the analysis and are not plotted.
#' 
#' @param read_count Data frame or CSV file containing the following columns:
#' `asv`, `sample`, `replicate` (optional), `read_count`.
#' @param identity_min Numeric value between 0 and 1. Minimum clustering identity threshold.
#' @param identity_max Numeric value between 0 and 1. Maximum clustering identity threshold.
#' @param identity_increment Numeric value between 0 and 1. Step size used to vary
#' clustering identity thresholds from `identity_min` to `identity_max`.
#' @param min_id Numeric value between 0 and 1. Minimum pairwise identity;
#' ASV pairs below this threshold are excluded from alignment and plotting.
#' @param vsearch_path Character string: path to the VSEARCH executable.
#' @param num_threads Positive integer: number of CPUs to use. If 0, all available CPUs are used.
#' @param outfile Character string: name of the output CSV file containing
#' `pairwise_asv_identity`, cluster assignment, and clustering parameter.
#' If empty, no file is written.
#' @param plotfile Character string: name of the output PNG file for the plot.
#' If empty, no file is written.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and display only warnings or errors.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#' 
#' @return A density plot of pairwise ASV identity distributions across clustering thresholds.
#' 
#' @examples
#' \dontrun{
#' plot <- plot_pairwise_identity_vsearch(
#'   read_count_df,
#'   identity_min = 0.9,
#'   identity_max = 0.99,
#'   identity_increment = 0.01,
#'   min_id = 0.8,
#'   vsearch_path = "vsearch",
#'   num_threads = 8,
#'   plotfile = "density_plot.png"
#' )
#' }
#' @export
#' 
plot_pairwise_identity_vsearch <- function(read_count, 
                                           identity_min=0.9, 
                                           identity_max=0.99,
                                           identity_increment=0.01,
                                           min_id = 0.8, 
                                           vsearch_path="vsearch", 
                                           num_threads=0,
                                           outfile="", 
                                           plotfile="",
                                           sep=",", 
                                           quiet=TRUE,
                                           log_file=NULL){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  # get function name, all arguments and stat time
  if(!is.null(log_file)){
    log <- collect_log()
  }
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
  pairwise_id <- pairwise_identity(asv = read_count_df, 
                                   min_id = 0.8, 
                                   vsearch_path=vsearch_path, 
                                   quiet=quiet, 
                                   num_threads=0
                                   )
  
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
    cluster_size <- cluster_vsearch(read_count = read_count_df, 
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
    theme(plot.title = element_text(size=12, hjust=0.5),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          strip.text = element_text(size = 10)
    )
  
  ### 
  
  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile, width = 2000, height = 1500, res = 300)
    print(p) # print plot to file
    dev.off()
  }
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(p)
  
} # end function


#' Classify clusters based on taxonomic agreement among ASVs
#'
#' Assigns each cluster a classification based on the taxonomic consistency
#' of the ASVs it contains at one or more taxonomic levels.
#'
#' Clusters are classified as:
#' - **closed**: all ASVs in the cluster belong to the same taxon, and that taxon
#'   is not present in any other cluster.
#' - **open**: all ASVs in the cluster belong to the same taxon, but that taxon
#'   is found in multiple clusters.
#' - **hybrid**: the cluster contains ASVs assigned to more than one taxon.
#'
#' @param cluster Data frame or CSV file with the following columns:
#' `asv_id`, `cluster_id`.
#' @param taxa Data frame or CSV file containing taxonomic assignments with one
#' column per taxonomic level (e.g. domain, phylum, class, order, family, genus, species).
#' Missing values are allowed and should be represented as `NA` when higher-level
#' resolution is unavailable.
#' @param taxlevels Character vector specifying the taxonomic levels to evaluate
#' (e.g. `"species"`, `"genus"`, `"family"`).
#' @param outfile Character string: name of the output CSV file. If empty, no file is written.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and display only warnings or errors.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#'
#' @return A data frame with columns:
#' `cluster_id` and one column per taxonomic level containing the classification
#' (`closed`, `open`, or `hybrid`).
#'
#' @examples
#' \dontrun{
#' df <- classify_clusters(
#'   cluster = cluster_df,
#'   taxa = taxa_df,
#'   taxlevels = c("species", "genus")
#' )
#' }
#' @export
#' 
classify_clusters <- function(cluster, taxa, outfile="", 
                              sep=",", quiet=TRUE, 
                              taxlevels=c("domain", "phylum", "class", "order","family", "genus", "species"),
                              log_file=NULL
){
  # get function name, all arguments and stat time
  if(!is.null(log_file)){
    log <- collect_log()
  }
  
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
        classification = dplyr::case_when(
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
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(class)
}


#' Plot cluster classification according to taxonomy
#'
#' Clusters ASVs using different clustering parameters and classifies each
#' cluster based on taxonomic consistency at selected taxonomic levels.
#'
#' Clusters are classified as:
#' - **closed**: all ASVs in the cluster belong to the same taxon, and that taxon
#'   is not present in any other cluster.
#' - **open**: all ASVs in the cluster belong to the same taxon, but that taxon
#'   is found in multiple clusters.
#' - **hybrid**: the cluster contains ASVs assigned to more than one taxon.
#'
#' The function then plots the number of clusters in each class for each clustering
#' parameter and taxonomic level.
#'
#' @param read_count Data frame or CSV file containing the following columns:
#' `asv_id`, `asv`, `read_count`.
#' @param taxa Data frame or CSV file containing taxonomic assignments with one
#' column per taxonomic level (e.g. domain, phylum, class, order, family, genus, species).
#' Missing values should be represented as `NA` when no assignment is available.
#' @param clustering_method Character string; clustering method to use:
#' `"swarm"` or `"vsearch"`.
#' @param cluster_params Numeric vector of clustering parameters:
#' Swarm `d` values or VSEARCH identity thresholds (values between 0 and 1).
#' @param vsearch_path Character string: path to the VSEARCH executable.
#' @param swarm_path Character string: path to the Swarm executable.
#' @param taxlevels Character vector specifying taxonomic levels to classify and plot
#' (e.g. `"species"`, `"genus"`).
#' @param num_threads Positive integer: number of CPUs to use. If 0, all available CPUs are used.
#' @param outfile Character string: name of the output CSV file containing columns:
#' `classification`, `number_of_clusters`, `taxlevel`, `clustering_parameter`.
#' If empty, no file is written.
#' @param plotfile Character string: name of the output PNG file for the plot.
#' If empty, no file is written.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and display only warnings or errors.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#'
#' @return A connected scatterplot showing the number of clusters in each class
#' (closed, open, hybrid) across clustering parameters and taxonomic levels.
#'
#' @examples
#' \dontrun{
#' plot <- plot_cluster_classification(
#'   read_count,
#'   taxa,
#'   clustering_method = "swarm",
#'   cluster_params = c(2, 4, 6, 8, 10),
#'   swarm_path = swarm_path,
#'   taxlevels = c("species", "genus"),
#'   num_threads = 8
#' )
#' }
#' @export
#' 
plot_cluster_classification <- function(read_count, taxa, 
                                        clustering_method="swarm", 
                                        cluster_params=c(2,4,6,8,10), 
                                        vsearch_path="vsearch", 
                                        swarm_path="swarm", 
                                        taxlevels= c("species", "genus"),
                                        outfile="",
                                        plotfile="",
                                        sep= ",",
                                        num_threads=0,
                                        quiet = TRUE,
                                        log_file=NULL){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
#  collect_log()
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
      cluster_df <- cluster_swarm(read_count = read_count_df, 
                                  swarm_d=i, 
                                  swarm_path=swarm_path, 
                                  num_threads=num_threads, 
                                  quiet=quiet)
    }else{
      cluster_df <- cluster_vsearch(read_count = read_count_df, 
                                    identity=i, 
                                    vsearch_path=vsearch_path, 
                                    num_threads=num_threads, 
                                    quiet=quiet)
    }
    
    ### classify
    classification <- classify_clusters(cluster = cluster_df, taxa = taxa_df, quiet=quiet, 
                                        taxlevels=taxlevels)
    #### count the number of clusters in each class (open, closed, hybrid, NA)
    for(tl in taxlevels){
      col_name <- paste("classification", tl, sep="_")
      
      tmp <- classification %>%
        dplyr::count(!!sym(col_name)) %>%
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
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  
  
  ### 
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(clusters_count, file = outfile,  row.names = F, sep=sep)
  }
  
  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile, width = 2000, height = 1500, res = 300) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  # add end_time and runtime, print
#  write_log(log, file=log_file)
  return(p)
  
}

#' Pool ASVs within clusters
#'
#' Aggregates ASVs belonging to the same cluster and sums their read counts.
#'
#' All ASVs assigned to the same `cluster_id` are merged, and their read counts
#' are summed for each sample (and replicate, if present).
#'
#' @param read_count_df Data frame containing the following columns:
#' `asv`, `sample`, `replicate` (optional), `read_count`.
#' @param cluster_df Data frame containing the following columns:
#' `asv_id`, `cluster_id`.
#'
#' @return A data frame with the same structure as `read_count_df`, where ASVs
#' belonging to the same cluster are pooled into a single entry per sample
#' (and replicate, if applicable).
#'
#' @examples
#' \dontrun{
#' pooled_df <- pool_by_cluster(
#'   read_count = read_count_df,
#'   cluster_df = cluster_df
#' )
#' }
#' @export
#' 
pool_by_cluster <- function(read_count_df, 
                            cluster_df){
  
  ### check if one to one relationship between asv and asv_id
  t <- check_one_to_one(read_count_df)
  
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

#' Cluster ASVs using Swarm or VSEARCH
#'
#' This function clusters ASVs into molecular operational taxonomic units (mOTUs)
#' or performs denoising depending on the chosen method and parameters.
#' When using Swarm with `d = 1`, the result is equivalent to denoising.
#'
#' The function can use either Swarm or VSEARCH (`cluster_size`) to cluster ASVs
#' in the input data set. Each ASV is assigned to a cluster, and the output format
#' depends on the `group` argument.
#'
#' Two output modes are available:
#'
#' - **group = TRUE**: ASVs belonging to the same cluster are merged into a single row.
#'   The `asv_id` and `asv` columns correspond to the cluster centroid, and read counts
#'   are summed across all ASVs in the cluster for each sample (and replicate, if present).
#'
#' - **group = FALSE**: The original data frame is returned with an additional
#'   `cluster_id` column. Each row corresponds to an individual ASV.
#'
#' Clustering can be performed globally or independently for each sample using
#' the `by_sample` argument.
#'
#' Note: when `by_sample = TRUE` and `group = FALSE`, the same `asv_id` may be
#' assigned different `cluster_id` values across samples.
#'
#' @param read_count Data frame or CSV file containing the following columns:
#' `asv`, `sample`, `replicate` (optional), `read_count`.
#' @param group Logical; if TRUE, ASVs in the same cluster are aggregated into a single row.
#' If FALSE, the input data frame is returned with an added `cluster_id` column.
#' @param by_sample Logical: if TRUE, clustering is performed separately for each sample.
#' @param method Character string specifying the clustering method: `"swarm"` or `"vsearch"`.
#' @param path Character string: path to the Swarm or VSEARCH executable.
#' @param swarm_d Positive integer: Swarm `d` parameter (maximum number of differences
#' allowed between ASVs to cluster them together).
#' @param fastidious Logical: if TRUE and `swarm_d = 1`, performs an additional
#' clustering pass to reduce small clusters.
#' @param identity Numeric value between 0 and 1. Identity threshold used for
#' VSEARCH `cluster_size` clustering.
#' @param outfile Character string: name of the output CSV file. If empty, no file is written.
#' @param num_threads Positive integer: number of CPUs to use. If 0, all available CPUs are used.
#' @param sep Field separator character used in input and output CSV files.
#' @param quiet Logical: if TRUE, suppress informational messages and display only warnings or errors.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#'
#' @return A data frame. If `group = TRUE`, ASVs are merged by cluster;
#' otherwise, the original structure is preserved with an additional `cluster_id` column.
#'
#' @examples
#' \dontrun{
#' read_count_df <- cluster_asv(
#'   read_count = read_count_df,
#'   group = TRUE,
#'   method = "vsearch",
#'   by_sample = TRUE,
#'   path = swarm_path,
#'   num_threads = 4
#' )
#' }
#' @export
#' 
cluster_asv <- function(read_count, 
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
                        quiet=T,
                        log_file=NULL
){
  
  if(num_threads == 0){
    num_threads <- parallel::detectCores()
  }
  
  # get function name, all arguments and stat time
  if(!is.null(log_file)){
    log <- collect_log()
  }
  
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
        cluster_df <- cluster_swarm(read_count = df_sample, 
                                    swarm_d=swarm_d, 
                                    swarm_path=path, 
                                    num_threads=num_threads, 
                                    quiet=quiet)
      }else{
        cluster_df <- cluster_vsearch(read_count = df_sample, 
                                      identity=identity, 
                                      vsearch_path=path, 
                                      num_threads=num_threads, 
                                      quiet=quiet)
      }
      
      # modify input df according to group
      if(group){ # pool ASV by cluster
        df_sample <- pool_by_cluster(read_count_df = df_sample, cluster_df = cluster_df)
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
      cluster_df <- cluster_swarm(read_count = read_count_df, 
                                  swarm_d=swarm_d, 
                                  swarm_path=path, 
                                  num_threads=num_threads, 
                                  quiet=quiet)
    }else{
      cluster_df <- cluster_vsearch(read_count = read_count_df, 
                                    identity=identity, 
                                    vsearch_path=path, 
                                    num_threads=num_threads, 
                                    quiet=quiet)
    }
    
    if(group){# pool ASV by cluster
      out_df <- pool_by_cluster(read_count_df = read_count_df, cluster_df = cluster_df)
    }
    else{  # add cluster_id column
      out_df <- left_join(read_count_df, cluster_df, by="asv_id")
    }
    
  }
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(out_df, file = outfile,  row.names = F, sep=sep)
  }
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(out_df)
}


#' Denoise data using SWARM
#'
#' This function performs **denoising** using SWARM clustering.
#'
#' Clustering can be applied to the entire dataset or independently for each
#' sample using the `by_sample` argument.
#'
#' By default, read counts of ASVs belonging to the same cluster are summed
#' within each sample–replicate combination. Optionally, an additional
#' post-processing step (`split_clusters = TRUE`) can be applied to further
#' refine clusters. In this step, clusters are split based on the relative
#' abundance of ASVs compared to the cluster centroid (see `?split_swarm_clusters`).
#' This option is only available when clustering is performed per sample
#' (`by_sample = TRUE`).
#'
#' @param read_count Data frame or CSV file with the following variables:
#'   `asv`, `sample`, `replicate` (optional), `read_count`.
#' @param by_sample Logical; if `TRUE`, clustering is performed separately for
#'   each sample.
#' @param split_clusters Logical; if `TRUE`, apply post-processing to split
#'   clusters based on highly abundant ASVs relative to the centroid.
#' @param min_abundance_ratio Numeric; used when `split_clusters = TRUE`.
#'   Minimum ratio of an ASV's read count relative to the centroid read count
#'   required for the ASV to define a new cluster.
#' @param min_read_count Numeric; used when `split_clusters = TRUE`.
#'   Minimum absolute read count required for an ASV to be considered for splitting.
#' @param swarm_path Character string: path to SWARM executable.
#' @param swarm_d Positive integer; SWARM `d` parameter. Maximum number of
#'   differences allowed between two ASVs for clustering (ASVs are grouped if
#'   they differ by ≤ d positions).
#' @param fastidious Logical; when `swarm_d = 1`, perform an additional clustering
#'   pass to reduce the number of small clusters.
#' @param num_threads Positive integer: number of CPUs to use. If 0, use all available CPUs.
#' @param sep Field separator character for input and output CSV files.
#' @param quiet Logical; if `TRUE`, suppress informational messages and show only warnings/errors.
#' @param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
#'
#' @return A data frame with the same structure as the input, where ASVs
#' belonging to the same cluster are pooled (read counts summed) within each
#' sample–replicate combination. If `split_clusters = TRUE`, clusters may be
#' further subdivided according to the splitting criteria.
#'
#' @examples
#' \dontrun{
#' read_count_df <- denoise_by_swarm(
#'   read_count = read_count_df,
#'   by_sample = TRUE,
#'   split_clusters = TRUE
#' )
#' }
#' @export
#' 
denoise_by_swarm <- function(read_count, 
                             by_sample=FALSE,
                             split_clusters=FALSE,
                             min_abundance_ratio = 0.2,
                             min_read_count = 10,
                             swarm_path="swarm", 
                             num_threads=0, 
                             swarm_d=1, 
                             fastidious=TRUE, 
                             outfile="", 
                             sep=",", 
                             quiet=TRUE,
                             log_file=NULL
){
  
  # get function name, all arguments and stat time
  if(!is.null(log_file)){
    log <- collect_log()
  }
  
  if (split_clusters & (!by_sample | swarm_d != 1 | !fastidious)) {
    stop(
      "ERROR: When split_clusters is TRUE:\n",
      "- by_sample must be TRUE\n",
      "- swarm_d must be 1\n",
      "- fastidious must be TRUE"
    )
  }
  
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
  
  if(by_sample){
    # make an empty output df, with the same columns and variable types as read_count_df
    # can dela with df with out without replicates
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
      
      # get cluster_id for all ASV
      cluster_df <- cluster_swarm(read_count = df_sample, 
                                  swarm_d=swarm_d, 
                                  swarm_path=swarm_path, 
                                  num_threads=num_threads, 
                                  quiet=quiet)
      
      # pool ASV by cluster
      if(split_clusters){
        df_sample <- split_swarm_clusters(df_sample,
                                          cluster_df,
                                          min_abundance_ratio = min_abundance_ratio,
                                          min_read_count = min_read_count
        )
      } else{
        df_sample <- pool_by_cluster(df_sample, cluster_df)
      }
      # add output of the sample to the total data frame    
      out_df <- rbind(out_df, df_sample)
    }# end for each sample
  }else{ # run swarm for all samples together
    # get cluster_id for all ASV
    cluster_df <- cluster_swarm(read_count = read_count_df, 
                                swarm_d=swarm_d, 
                                swarm_path=swarm_path, 
                                num_threads=num_threads, 
                                quiet=quiet)
    out_df <- pool_by_cluster(read_count_df, cluster_df)
  }
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(out_df, file = outfile,  row.names = F, sep=sep)
  }
  # add end_time and runtime, print
  write_log(log, file=log_file)
  return(out_df)
}

#' Split SWARM Clusters
#'
#' During SWARM denoising (including when `fastidious = TRUE` and `d = 1`),
#' ASVs differing by a single nucleotide can be grouped into the same cluster,
#' even when they have comparable abundances. This function refines such
#' clusters by identifying ASVs whose read counts exceed both a minimum
#' relative abundance threshold (`min_abundance_ratio`) and a minimum absolute
#' read count (`min_read_count`).
#'
#' For each cluster, ASVs meeting both criteria are retained as independent
#' units and each defines a new cluster. The total number of reads in the
#' original cluster is then redistributed among these selected ASVs while
#' preserving their relative abundances.
#'
#' If no ASVs meet the criteria, reads are summed across ASVs within each
#' cluster–sample–replicate combination, and the cluster is retained as a
#' single unit.
#'
#' @param read_count Data frame or CSV file with the following variables:
#'   `asv`, `sample`, `replicate` (optional), `read_count`.
#' @param cluster Data frame or CSV file with the following variables:
#'   `asv_id`, `cluster_id`.
#' @param min_abundance_ratio Numeric; minimum ratio of an ASV's read count
#'   relative to the cluster centroid required for the ASV to define a new
#'   cluster (e.g., 0.2 corresponds to 20% of centroid abundance).
#' @param min_read_count Numeric; minimum absolute read count required for an
#'   ASV to be considered for splitting.
#' @param outfile Character string: name of the output CSV file. If empty, no
#'   file is written.
#' @param sep Field separator character for input and output CSV files.
#' @param quiet Logical; if `TRUE`, suppress informational messages and show
#'   only warnings or errors.
#'
#' @return A data frame with the same structure as the input. For each
#' cluster–sample–replicate combination, ASVs are either:
#' \itemize{
#'   \item split into multiple clusters defined by ASVs that meet the
#'   `min_abundance_ratio` and `min_read_count` criteria, with read counts
#'   redistributed while preserving relative abundances, or
#'   \item pooled into a single row (reads summed across ASVs) if no ASV meets
#'   the splitting criteria.
#' }
#'
#' @examples
#' \dontrun{
#' read_count_df <- split_swarm_clusters(
#'   read_count = read_count_df,
#'   min_abundance_ratio = 0.1,
#'   min_read_count = 10
#' )
#' }
#' @export
#' 
split_swarm_clusters <-function(read_count,
                                cluster,
                                min_abundance_ratio = 0.1,
                                min_read_count = 10,
                                outfile = "",
                                sep = ",",
                                quiet = TRUE
){
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  if(is.character(cluster)){
    # read known occurrences
    cluster_df <- read.csv(cluster, header=T, sep=sep)
  }else{
    cluster_df <- cluster
  }
  
  read_count_df <- left_join(read_count_df, cluster_df, by="asv_id")
  
  ###############################
  # prepare list of clusters to split, and for each of them the list of ASV to keep apart
  clusters_to_split <- read_count_df %>%
    group_by(asv_id, cluster_id) %>%
    summarise(read_count = sum(read_count), .groups = "drop") %>%
    # add centroid read count
    group_by(cluster_id) %>%
    mutate(centroid_read_count = read_count[asv_id == cluster_id]) %>%
    ungroup() %>%
    # keep only asv, that have at least 20% of the centroids reads
    filter(read_count / centroid_read_count > min_abundance_ratio & 
             read_count > min_read_count &
             asv_id != cluster_id) %>%
    select(asv_id, cluster_id) 
  
  # get a list of centroids of the clusters_to_split
  centroids <- data.frame(
    asv_id = unique(clusters_to_split$cluster_id),
    cluster_id = unique(clusters_to_split$cluster_id)
  )
  # add them to clusters_to_split
  clusters_to_split <- rbind(centroids, clusters_to_split) %>%
    arrange(cluster_id)
  
  ###############################
  
  # for all clusters that should be split up among selected ASVs,
  # modify the original read count of the selected ASV to
  #  - keep their proportions among different replicates
  #  - Their sum should equal the total number of read of the cluster.
  # In other words, split of the total number of reads of a cluster among selected
  # ASV, and keep the proportions of the original read counts.
  #
  # Replace the cluster_id by the asv_id
  
  cluster_selected <- read_count_df %>% # add the total number of reads of the cluster
    group_by(cluster_id) %>%
    mutate(cluster_rc_all = sum(read_count)) %>%
    ungroup() %>%
    # keep only asv, that should be kept apart
    filter(asv_id %in% clusters_to_split$asv_id) %>%
    # get the total number of reads of the selected ASVs of the cluster
    group_by(cluster_id) %>%
    mutate(cluster_rc_selected = sum(read_count)) %>%
    ungroup() %>%
    # multiply all read count by the proportion of original and selected read_count of the cluster
    mutate(
      read_count_cis = round(
        as.numeric(read_count) * cluster_rc_all / cluster_rc_selected,
        digits = 0
      )
    )
  #    mutate(read_count_cis = round(read_count * cluster_rc_all / cluster_rc_selected, digits=0))
  # Check if the sum of the modified read count equals of the original total read of the cluster
  #  group_by(cluster_id) %>%
  #  mutate(read_count_cis_sum = sum(read_count_cis)) %>%
  #  ungroup() %>%
  
  if( "replicate" %in% colnames(read_count_df)){
    cluster_selected <- cluster_selected %>%
      select(asv_id, sample, replicate, read_count = read_count_cis, asv)
  }else{
    cluster_selected <- cluster_selected %>%
      select(asv_id, sample, read_count = read_count_cis, asv)
  }
  
  ###############################
  # make a list of asv_id asv for the centroids 
  cluster_seq <- read_count_df %>%
    filter(asv_id == cluster_id) %>%
    select(asv_id, asv) %>%
    distinct()
  
  ###############################
  unmodified_clusters <- read_count_df %>%
    filter(!cluster_id %in% clusters_to_split$cluster_id)
  
  if( "replicate" %in% colnames(read_count_df)){
    unmodified_clusters <- unmodified_clusters %>%
      group_by(cluster_id, sample, replicate) %>%
      summarize(read_count = sum(read_count), .groups = "drop") %>%
      select(asv_id = cluster_id, sample, replicate, read_count) %>%
      left_join(cluster_seq, by="asv_id")
    
  }else{
    unmodified_clusters <- unmodified_clusters %>%
      group_by(cluster_id, sample) %>%
      summarize(read_count = sum(read_count), .groups = "drop") %>%
      select(asv_id = cluster_id, sample, read_count) %>%
      left_join(cluster_seq, by="asv_id")
  }
  
  ###############################
  # pool modified and unmodified clusters
  grouped_swarm <- rbind(cluster_selected, unmodified_clusters)
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(grouped_swarm, file = outfile,  row.names = F, sep=sep)
  }
  
  return(grouped_swarm)
}

