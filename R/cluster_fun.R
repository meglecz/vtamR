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

#' Pairwise identity
#' 
#' Align all pairs of asv with at least min_id similarity and
#' make a data frame with pairs of asv and their percentage of identity. 
#' asv pairs with identity bellow min_id are not listed.
#' 
#' @param asv data frame or csv file with asv and asv_id columns.
#' @param min_id bellow this percentage of identity do not aligne asv
#' @param vsearch_path Character string: path to vsearch executables.
#' @param num_threads Positive integer: Number of CPUs. If 0, use all available CPUs.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns data frame with asv pairs and their percentage of identity. 
#' asv pairs with identity bellow min_id are not listed.
#' @examples 
#' \dontrun{
#' check_dir(path="data")
#' }
#' @export
#' 
PairwiseIdentity <- function(asv, min_id = 0.8, vsearch_path=vsearch, quiet=TRUE, num_threads=0){
  
  # can accept df or file as an input
  if(is.character(asv)){
    # read known occurrences
    asv_df <- read.csv(asv, header=T, sep=sep)
  }else{
    asv_df <- asv
  }
  
  asv_df <- asv_df %>%
    select(asv, asv_id) %>%
    distinct()
  t <- check_one_to_one_relationship(asv_df) # stop execution, if FALSE
  
  fasta <- file.path(tempdir(), "asv.fasta")
  print(fasta)
  write_fasta_df(asv_df, outfile=fasta, read_count=FALSE)
  
  vsearch_out <- file.path(tempdir(), "vsearch_out.tsv")
  
  cmd <- paste(vsearch_path, "--allpairs_global", fasta, "--userout", vsearch_out, 
               '--userfields "query+target+id+ids+alnlen+aln"', "--id",  min_id, "--threads", num_threads, sep=" ")
  if(!quiet){
    print(cmd)
  }
  system(cmd)
  
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
      filter(identity >= min_id)
  } else {
    results_vsearch <- data.frame(query=numeric(),
                                  target= numeric(),
                                  identity= numeric())
  }
  return(results_vsearch)
}


#' Pairwise identity Density plot for different Swarm d values
#' 
#' Cluster by swarm all asv with a range of d values.
#' For each d, make a density plot of pairwise percentage of identities between 
#' asv of the same cluster (within) or different clusters (between).
#' 
#' @param asv data frame or csv file with asv and asv_id and read_count columns.
#' @param min_id bellow this percentage of identity do not aligne asv
#' @param vsearch_path Character string: path to vsearch executables.
#' @param num_threads Positive integer: Number of CPUs. If 0, use all available CPUs.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv, sample, replicate (optional), read_count.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param swarm_path Character string: path to swarm executables. 
#' @param num_threads Positive integer: Number of CPUs.
#' @param swarm_d_min Positive integer: Minimum value of d for Swarm.
#' @param swarm_d_max Positive integer: Maximum value of d for Swarm.
#' @param swarm_d_increment Positive integer: increase d by swarm_d_increment between
#' swarm_d_min and swarm_d_max.
#' @param sep Field separator character in input and output csv files.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' 
#' show warnings or errors.
#' @returns data frame with asv pairs and their percentage of identity. 
#' asv pairs with identity bellow min_id are not listed.
#' @examples 
#' \dontrun{
#' check_dir(path="data")
#' }
#' @export
#' 
#PairwiseIdentityPlotPerSwarmD <- function(asv, min_id = 0.8, vsearch_path=vsearch, quiet=TRUE, num_threads=0){
PairwiseIdentityPlotPerSwarmD <- function(read_count, 
                                          outfile="", 
                                          swarm_path="swarm", 
                                          num_threads=0, 
                                          swarm_d_min=3, 
                                          swarm_d_max=15,
                                          swarm_d_increment=3,
                                          sep=",", 
                                          quiet=T ){
  
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  # make df with query, target, identity
  pairwise_id <- PairwiseIdentity(read_count_df, min_id = 0.8, vsearch_path=vsearch_path, quiet=TRUE, num_threads=0)
  
  # make a df with unique asv, asv_id and readcount (sum)
  asv_df <- read_count_df %>%
    group_by(asv, asv_id) %>%
    summarize(rc = sum(read_count), .groups="drop")
  # make a fasta file with these
  input_swarm <- file.path(tempdir(), "swarm_input.fasta")
  writeLines(paste(">", asv_df$asv_id, "_", 
                   asv_df$rc, "\n", 
                   asv_df$asv, 
                   sep=""), 
             input_swarm)
  
  for(d in seq(swarm_d_min, swarm_d_max, by=swarm_d_increment)){
    print(d)
    
    tmp_dir <-paste('tmp_swarm_', d, '_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    tmp_dir <- file.path(tempdir(), tmp_dir)
    check_dir(tmp_dir)
    
    # swarm on the same fasta with different D
    # clusters.txt each line is a cluster, with asv_ids separated  by space
    clusters <- file.path(tmp_dir, "clusters.txt")
    swarm <- paste(swarm_path, 
                   " -d ", d,
                   " -o ", clusters, 
                   sep=""
    )
    if(num_threads > 0){ # if num_threads have been specified
      swarm <- paste(swarm, 
                     " -t ", num_threads, 
                     sep=""
      )
    }
    swarm <- paste(swarm, input_swarm, sep=" ")
    if(!quiet){
      print(swarm)
    }
    system(swarm)
    
    
    
    # get df with asv_id cluster_id
    # add to pairwise_id the clusters of eaxh quary and target and define if it is within or between cluster
    # prot %id with betwee, and within coloration
    
  }
  
}

# PairwiseIdentityPlotPerSwarmD(read_count_df, num_threads=0, swarm_d_min=3, swarm_d_max=15, swarm_d_increment=3)

#results_vsearch <- PairwiseIdentity(read_count_df, vsearch_path=vsearch_path, num_threads=0, quiet=FALSE, min_id=0.8)


