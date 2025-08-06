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
  
  #### get unique asv list and make fasta file
  asv_df <- asv_df %>%
    select(asv, asv_id) %>%
    distinct()
  t <- check_one_to_one_relationship(asv_df) # stop execution, if FALSE
  
  fasta <- file.path(tempdir(), "asv.fasta")
  print(fasta)
  write_fasta_df(asv_df, outfile=fasta, read_count=FALSE)
  
  ### run vsearch allpairs_global
  vsearch_out <- file.path(tempdir(), "vsearch_out.tsv")
  cmd <- paste(vsearch_path, "--allpairs_global", fasta, "--userout", vsearch_out, 
               '--userfields "query+target+id+ids+alnlen+aln"', 
               "--id", min_id, "--threads", num_threads, sep=" ")
  if(!quiet){
    print(cmd)
  }
  system(cmd)
  
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
      mutate(identity = as.numeric(identity)) %>%
      filter(identity >= min_id)
    results_vsearch$query <- as.numeric(results_vsearch$query)
    results_vsearch$target <- as.numeric(results_vsearch$target)
  } else {
    results_vsearch <- data.frame(query=numeric(),
                                  target= numeric(),
                                  identity= numeric())
  }
  
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
  
  #####
  # make pairwise_id df with query, target, identity
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
    print(d)
    
    #####
    # run swarm
    
    # make tmp dir separatelit for each d
    tmp_dir <-paste('tmp_swarm_', d, '_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    tmp_dir <- file.path(tempdir(), tmp_dir)
    check_dir(tmp_dir)
    
    # clusters.txt each line is a cluster, with asv_ids separated  by space
    clusters <- file.path(tmp_dir, "clusters.txt")
    swarm <- paste(swarm_path, " -d ", d, " -o ", clusters, sep="")
    if(num_threads > 0){ # if num_threads have been specified
      swarm <- paste(swarm, " -t ", num_threads, sep="")
    }
    swarm <- paste(swarm, input_swarm, sep=" ")
    if(!quiet){
      print(swarm)
    }
    system(swarm)
    
    #####
    # make a data frame with cluster_id and asv_id columns, 
    # where asv_id has all swarm input asv_id, 
    # and  cluster_id is the name of the cluster they belong to
    
    # read.table is unpredictable. Use a more complicated, but more sure solution.
    # Read the file line by line
    lines <- readLines(clusters)
    # Split each line by whitespace
    split_lines <- strsplit(lines, "[[:space:]]+")
    # Determine the maximum number of fields
    max_cols <- max(sapply(split_lines, length))
    # Convert to a data frame and fill empty cells by NA
#    cluster_df <- as.data.frame(do.call(rbind, lapply(split_lines, function(x) c(x, rep(NA, max_cols - length(x))))),
#                                stringsAsFactors = FALSE)
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
    # add lines to the oveall df
    pairwise_id_d <- rbind(pairwise_id_d, pairwise_id_local)
  } # end for
  
  bad_rows <- pairwise_id_d %>%
    filter(!is.finite(identity) | identity < 80 | identity > 100)
  
  print(bad_rows)
  
  
  ####
  # Make the plot
  # fix order of d
  pairwise_id_d$swarm_d <- factor(pairwise_id_d$swarm_d, levels = unique(pairwise_id_d$swarm_d))

  p <-ggplot(pairwise_id_d, aes(x = identity, fill = cluster)) +
    geom_density(adjust = 1.5, alpha = 0.4) +
    scale_x_continuous(limits = c(80, 100)) +
    facet_wrap(~swarm_d, scales = "free_y") +  # one subplot per d
    ggtitle("Pairwise percent idnetity between ASV ofthe same or different cluster according Swarm's d") +
    theme_minimal()
  
  ### 
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    png(filename=outfile) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  return(p)

} # end function


