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
  
  # system2 cannot use ~ as a hone
  vsearch_path <- path.expand(vsearch_path)
  
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
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the following columns: asv_id, cluster_id
#' @examples 
#' \dontrun{
#' plot <- cluster_swarm(read_count_df, 
#'                      swarm_d=7, 
#'                      swarm_path="swarm",
#'                      num_threads=8)
#' }
#' @export
cluster_swarm <- function(read_count, 
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
  
  # system2 cannot use ~ as a home
  swarm_path <- path.expand(swarm_path)
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



#' Cluster all input ASV by cluster_size function of vsearch
#' 
#' Cluster all input ASV and return a data frame with asv_id, cluster_id
#' columns
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, asv, read_count.
#' @param identity Real; Identity threshold for clustering.
#' @param vsearch_path Character string: path to vsearch executable. 
#' @param num_threads Positive integer: Number of CPUs.
#' @param outfile Character string: name of the output csv file
#' @param sep Field separator character in input and output csv files.
#' @param quiet Logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the following columns: asv_id, cluster_id
#' @examples 
#' \dontrun{
#' plot <- cluster_vsearch_cluster_size(read_count_df, 
#'                                       identity=0.97,
#'                                       vsearch_path="vsearch",
#'                                       num_threads=8)
#' }
#' @export
cluster_vsearch_cluster_size <- function(read_count, 
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
  
  # system2 cannot use ~ as a home
  vsearch_path <- path.expand(vsearch_path)

  #####
  # make a df with unique asv, asv_id and readcount (sum)
  asv_df <- read_count_df %>%
    group_by(asv, asv_id) %>%
    summarize(read_count = sum(read_count), .groups="drop")%>%
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
#' @param outfile Character string: output csv file name with the following columns:
#' identity, cluster, cluster_criterium; If empty, no file is written.
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
  
  #####
  # make pairwise_id df with query, target, identity
  if(!quiet){
    print("Calculating pairwise identities")
  }
  pairwise_id <- PairwiseIdentity(read_count_df, min_id = 0.8, vsearch_path=vsearch_path, quiet=TRUE, num_threads=0)
  
  
  #####
  # initialize data frame (cluster: same/different)
  pairwise_id_final <- data.frame(identity = numeric(),
                              cluster = character(),
                              cluster_criterium= factor())
  
  #####
  # for each d
  for(d in seq(swarm_d_min, swarm_d_max, by=swarm_d_increment)){
    
    if(!quiet){
      cat("Running swarm with d=", d, "\n")
    }
    #### Run swarm 
    # cluster_df : avs_id, clsuter_id
    cluster_df <- cluster_swarm(read_count_df, swarm_d=d,fastidious=FALSE,
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
      mutate(cluster_criterium=paste0("swarm's d: ",d)) %>%
      select(identity, cluster, cluster_criterium)
    # add lines to the overall df
    pairwise_id_final <- rbind(pairwise_id_final, pairwise_id_local)
  } # end for
  
  ####
  # Make the plot
  # fix order of d
  pairwise_id_final$cluster_criterium <- factor(pairwise_id_final$cluster_criterium, levels = unique(pairwise_id_final$cluster_criterium))
  
  p <-ggplot(pairwise_id_final, aes(x = identity, fill = cluster, color = cluster)) +
    geom_density(adjust = 1.5, alpha = 0.4) +
    scale_x_continuous(limits = c(80, 100)) +
    facet_wrap(~cluster_criterium, scales = "free_y") +  # one subplot per d
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
#' @param outfile Character string: csv file name with the following columns: 
#' identity, cluster, cluster_criterium; If empty, no file is written.
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
  
  #####
  # make pairwise_id df with query, target, identity
  if(!quiet){
    print("Calculating pairwise identities")
  }
  pairwise_id <- PairwiseIdentity(read_count_df, min_id = 0.8, vsearch_path=vsearch_path, quiet=quiet, num_threads=0)
  
  #####
  # initialize data frame (cluster: same/different)
  pairwise_id_final <- data.frame(identity = numeric(),
                              cluster = character(),
                              cluster_criterium= factor())
  
  #####
  # for each  identity threshold
  for(d in seq(identity_min, identity_max, by=identity_increment)){
    
    if (!quiet) {
      cat("Running cluster_size of vsearch with identity threshold", d, "\n")
    }
    
    #####
    # run vsearch
    cluster_size <- cluster_vsearch_cluster_size(read_count_df, 
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
      mutate(cluster_criterium= paste0("Clustering with ", d, "indentity")) %>%
      select(identity, cluster, cluster_criterium)
    # add lines to the overall df
    pairwise_id_final <- rbind(pairwise_id_final, pairwise_id_local)
  } # end for
  
  ####
  # Make the plot
  # fix order of d
  pairwise_id_final$cluster_criterium <- factor(pairwise_id_final$cluster_criterium, levels = unique(pairwise_id_final$cluster_criterium))

  
  p <-ggplot(pairwise_id_final, aes(x = identity, fill = cluster, color = cluster)) +
     geom_density(adjust = 1.5, alpha = 0.4) +
     scale_x_continuous(limits = c(80, 100)) +
     facet_wrap(~cluster_criterium, scales = "free_y") +  # one subplot per d
     ggtitle("Pairwise percent identity between ASV of the \n same or different cluster according identity threshold of clustering") +
    theme(plot.title = element_text(size=12, hjust=0.5)) 
  
  ### 
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(pairwise_id_final, file = outfile,  row.names = F, sep=sep)
  }
  return(p)
  
} # end function

#' Classify clusters based on taxonomic agreement among their ASVs
#'
#' Classifies each cluster according to the taxonomic assignment at a given level:
#' 
#' - closed: ASVs in the cluster are assigned to the same taxon, 
#' and all ASVs of that taxon belong exclusively to this cluster.  
#' - open: ASVs in the cluster are assigned to the same taxon, 
#' but some ASVs of that taxon are found in other clusters.  
#' - hybrid: the cluster contains ASVs assigned to more than one taxon.  
#'   
#' 
#' @param cluster Data frame or csv file with the following variables: 
#' asv_id, cluter_id
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
#' @returns Data frame with the following columns: cluster_id, class for each taxonomic level
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
    
    # seletc the taxlevel, delete NA
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
#' @param clustering_method character; [swarm/cluster_size]
#' @param cluster_params numerical vector of either swarm's d to be used, or 
#' percentage of identity threshold for the cluster_siz algorithm of vsearch.
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
        cluster_df <- cluster_swarm(read_count_df, 
                                       swarm_d=i, 
                                       swarm_path=swarm_path, 
                                       num_threads=num_threads, 
                                       quiet=quiet)
      }else{
        cluster_df <- cluster_vsearch_cluster_size(read_count_df, 
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
  








