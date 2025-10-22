

## load
library(vtamR)
library(dplyr)
library(ggplot2)

setwd("/home/meglecz/vtamR/")
library("devtools")
library("roxygen2")
load_all(".")
roxygenise()
usethis::use_roxygen_md()

### set up
cutadapt_path <- "~/miniconda3/envs/vtam/bin/cutadapt"
vsearch_path <- "~/miniconda3/envs/vtam/bin/vsearch"
blast_path <- "~/miniconda3/envs/vtam/bin/blastn"
swarm_path <- "swarm"
sep <- ","
outdir <- "~/vtamR_demo_out"

fastq_dir <- system.file("extdata/demo/fastq", package = "vtamR")
fastqinfo <-  system.file("extdata/demo/fastqinfo.csv", package = "vtamR")
mock_composition <-  system.file("extdata/demo/mock_composition.csv", package = "vtamR")
asv_list <-  system.file("extdata/demo/asv_list.csv", package = "vtamR")
taxonomy <- system.file("extdata/db_test/taxonomy_reduced.tsv", package = "vtamR")
blast_db <- system.file("extdata/db_test", package = "vtamR")
blast_db <- file.path(blast_db, "COInr_reduced")

### Merge
merged_dir <- file.path(outdir, "merged")
fastainfo_df <- Merge(fastqinfo, 
                      fastq_dir=fastq_dir, 
                      vsearch_path=vsearch_path, 
                      outdir=merged_dir,
                      fastq_maxee=1,
                      fastq_maxns=0,
                      fastq_allowmergestagger=F)

### demultiplex
sorted_dir <- file.path(outdir, "sorted")
sortedinfo_df <- SortReads(fastainfo_df, 
                           fasta_dir=merged_dir, 
                           outdir=sorted_dir, 
                           check_reverse=TRUE, 
                           cutadapt_path=cutadapt_path, 
                           vsearch_path=vsearch_path,
                           tag_to_end = T,
                           primer_to_end=T)

###############
### dereplicate
###############
updated_asv_list <- file.path(outdir, "updated_asv_list.tsv")
sortedinfo <- file.path(sorted_dir, "sortedinfo.csv")
read_count_df <- Dereplicate(sortedinfo, 
                             dir=sorted_dir, 
                             asv_list=asv_list,
                             updated_asv_list = updated_asv_list
)

### stat
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())

stat_df <- GetStat(read_count_df, stat_df, stage="input_sample_replicate", params=NA)

#### swarm
by_sample <- TRUE
d=1
fastidious= TRUE
quiet=TRUE
outfile <- file.path(outdir, "filter", "2_Swarm_by_sample.csv")

read_count_df <- ClusterASV(read_count_df,
                            method = "swarm",
                            swarm_d=d,
                            fastidious=fastidious,
                            by_sample=by_sample,
                            group = TRUE,
                            path=swarm_path
                            )


stat_df <- GetStat(read_count_df, stat_df, stage="Swarm", params=by_sample)

#### LFNglobalReadCount
global_read_count_cutoff = 2

read_count_df <- LFNglobalReadCount(read_count_df, 
                                    cutoff=global_read_count_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNglobalReadCount", params=NA)

#### FilterIndel
read_count_df <- FilterIndel(read_count_df)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterIndel", params=NA)

### FilterCodonStop
genetic_code = 5
read_count_df <- FilterCodonStop(read_count_df, 
                                 genetic_code=genetic_code)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterCodonStop", params=NA)

### FilterExternalContaminant
conta_file <- file.path(outdir, "tmp", "external_contamination.csv")
read_count_df <- FilterExternalContaminant(read_count_df, 
                          sample_types=sortedinfo, 
                          conta_file=conta_file)

stat_df <- GetStat(read_count_df, stat_df, stage="FilterExternalContaminant", params=NA)

### FilterChimera
abskew=2
by_sample = T
sample_prop = 0.8
read_count_df <- FilterChimera(read_count_df, 
                               vsearch_path=vsearch_path, 
                               by_sample=by_sample, 
                               sample_prop=sample_prop, 
                               abskew=abskew)

stat_df <- GetStat(read_count_df, stat_df, stage="FilterChimera", params=NA)

#### FilterRenkonen
cutoff <- 0.4
read_count_df <- FilterRenkonen(read_count_df, 
                                cutoff=cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterRenkonen", params=cutoff)

### FilterPCRerror
pcr_error_var_prop <- 0.05
max_mismatch <- 2
read_count_df <- FilterPCRerror(read_count_df, 
                                vsearch_path=vsearch_path, 
                                pcr_error_var_prop=pcr_error_var_prop, 
                                max_mismatch=max_mismatch)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterPCRerror", params=NA)

### LFNsampleReplicate
lfn_sample_replicate_cutoff <- 0.004
read_count_df <- LFNsampleReplicate(read_count_df, 
                                    cutoff=lfn_sample_replicate_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNsampleReplicate", params=NA)

### FilterMinReplicate
min_replicate_number <- 2
read_count_df <- FilterMinReplicate(read_count_df, 
                                    cutoff=min_replicate_number)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterMinReplicate", params=NA)


### LFNvariant
lnf_variant_cutoff = 0.001
read_count_df_lnf_variant <- LFNvariant(read_count_df, 
                                        cutoff=lnf_variant_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNvariant", params=NA)

### LFNreadCount
lfn_read_count_cutoff <- 10
read_count_df_lfn_read_count <- LFNreadCount(read_count_df, 
                                             cutoff=lfn_read_count_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNreadCount", params=NA)




### Combine results
read_count_df <- PoolFilters(read_count_df_lfn_read_count, 
                             read_count_df_lnf_variant)
stat_df <- GetStat(read_count_df, stat_df, stage="Combine results", params=NA)

# delete temporary data frames
rm(read_count_df_lfn_read_count)
rm(read_count_df_lnf_variant)

### FilterMinReplicate
min_replicate_number <- 2
read_count_df <- FilterMinReplicate(read_count_df, 
                                    cutoff=min_replicate_number)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterMinReplicate", params=NA)


plot_swarm <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      outfile="13_FilterMinReplicate.png")




#' Generate ASV-Specific Cutoff Values
#'
#' In some datasets, certain variants may occur in many samples with high read counts.  
#' As a result, reads due to tag-jump contamination may remain unfiltered when using  
#' a fixed cutoff value in the *LFNvariant* function.  
#'
#' This function computes variant-specific cutoff values for use in *LFNvariant*,  
#' targeting ASVs known to have false-positive occurrences 
#' (detected by MakeKnownOccurrences). 
#'
#' For each ASV, the function:
#' - Identifies all its false-positive occurrences.
#' - Takes the maximum read count among these false positives.
#' - Divides this value by the total number of reads of the ASV across the dataset,  
#'   or within replicates if `by_replicate = TRUE`.  
#'
#' Because some false positives may not result from tag-jump contamination, these  
#' cutoff values can be excessively high. Therefore:
#' - Filter the dataset as much as possible **before** using ASV-specific cutoffs.
#' - Set a reasonable upper limit for ASV-specific cutoffs using the `max_cutoff` parameter.  
#'
#' @param read_count Data frame or csv file with the following variables: 
#' asv, sample, replicate, read_count.
#' @param max_cutoff Numeric; the maximum allowed cutoff value.
#' @param mock_composition Data frame or csv file with columns: 
#' sample, action (keep/tolerate), asv.
#' @param habitat_proportion Numeric. Value between 0 and 1: for each asv, if the proportion 
#' of reads in a habitat is below this cutoff,
#' it is considered as a false positive in all samples of the habitat.
#' @param by_replicate Logical; if `TRUE`, compute cutoffs separately by replicates.
#' @param outfile Character string: output file. If empty, no file is written. 
#' @param sep Field separator character in input and output csv files.
#' @returns Data frame with the following columns: 
#' asv_id, replicate (if by_replicate), cutoff
#' @seealso [LFNvariant()]
#' @examples
#' \dontrun{
#' ASVspecificCutoff(read_count=read_count_samples_df, 
#'     mock_composition="data/mock_composition.csv"
#'     )
#' }
#' @export
#'
ASVspecificCutoff <- function(read_count, 
                              max_cutoff=0.05,
                              mock_composition="",
                              habitat_proportion=0.5,
                              by_replicate=FALSE, 
                              outfile="", 
                              sep=",")  {
  
  # can accept df or file as an input
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  ### MakeKnownOccurrences to make known_occurrences_df
  results <- MakeKnownOccurrences(read_count_df, 
                                  sortedinfo=sortedinfo, 
                                  mock_composition=mock_composition,
                                  habitat_proportion=habitat_proportion)

  known_occurrences_df <- results[[1]]
  
  # total number of read by asv, ou asv.replicate
  if(by_replicate){
    asv_total_rc <- read_count_df %>%
      group_by(asv_id, replicate) %>%
      summarize(total_rc = sum(read_count), .groups="drop")
  }else{
    asv_total_rc <- read_count_df %>%
      group_by(asv_id) %>%
      summarize(total_rc = sum(read_count))%>%
      ungroup()
  }
  
  # list of asv-sample FP
  delete_occurrences_df <- known_occurrences_df %>%
    filter(action=="delete") %>%
    select(sample,asv_id) %>%
    distinct() %>%
    mutate(asv_sample = paste(asv_id, sample, sep="."))
  
  # 
  if(by_replicate){
    asv_spec_cutoff_df <- read_count_df %>%
      mutate(asv_sample = paste(asv_id, sample, sep=".")) %>%
      filter(asv_sample %in% delete_occurrences_df$asv_sample) %>%
      group_by(asv_id, replicate) %>%
      filter(read_count==max(read_count))%>%
      ungroup() %>%
      left_join(asv_total_rc, by=c("asv_id", "replicate")) %>%
      mutate(cutoff_asv_spec = read_count/total_rc) %>%
      select(asv_id, replicate, cutoff_asv_spec)
      
  }else{
    asv_spec_cutoff_df <- read_count_df %>%
      mutate(asv_sample = paste(asv_id, sample, sep=".")) %>%
      filter(asv_sample %in% delete_occurrences_df$asv_sample) %>%
      group_by(asv_id) %>%
      filter(read_count==max(read_count))%>%
      ungroup() %>%
      left_join(asv_total_rc, by=c("asv_id")) %>%
      mutate(cutoff_asv_spec= read_count/total_rc)%>%
      select(asv_id, cutoff_asv_spec)
  }
  
  # adjust too high values to max_cutoff
  asv_spec_cutoff_df <- asv_spec_cutoff_df %>%
    mutate(cutoff_asv_spec = if_else(
      cutoff_asv_spec > max_cutoff, max_cutoff, cutoff_asv_spec))

  # write to outfile
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(asv_spec_cutoff_df, file=outfile, row.names = F, sep=sep)
  }
  
  return(asv_spec_cutoff_df)
}



ASVspecificCutoff_df <- ASVspecificCutoff(read_count_df,  mock_composition=mock_composition,
                              by_replicate=TRUE, 
                              outfile="tmp/ASVspecificCutoff_by_replicate_false.csv")



#' LFNvariant2
#' 
#' This function filters out false positives present dut to tag-jump or light
#' intersame contamination.
#' 
#' If by_replicate is FALSE: Eliminate occurrences where the 
#' (read_count/read_count of the asv in the data set) is less than cutoff.
#' If by_replicate is TRUE: Eliminate occurrences where the 
#' (read_count/read_count of the asv in its replicate) is less than cutoff.
#' 
#' Issues a warning if the total read count of an ASV has been reduced 
#' bellow min_read_count_prop, since it can indicate a to high cutoff value.
#' 
#' By default, the same cutoff value is applied for all variants. However, it is
#' also possible to use variant specific cutoffs, present in asv_specific_cutoffs
#' data frame or csv file.
#' 
#' @param read_count Data frame or csv file with the following variables: 
#' asv_id, sample, replicate, read_count, asv.
#' @param cutoff Numeric. Value between 0 and 1: minimum proportion of the read count of
#'  an occurrence within all reads of the asv or asv-replicate. Bellow this cutoff
#'  the occurrence is deleted.
#' @param  asv_specific_cutoffs a data fral or csv file with the following columns.
#' asv_id, replicate (if by_replicate), asv_specific_cutoff
#' @param by_replicate logical: Compare read count of the occurrence to the 
#' read counts of the ASV-replicate.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param sep Field separator character in input and output csv files.
#' @param min_read_count_prop Numeric. Value between 0 and 1: If the proportion of the read count 
#' of a variant in the output compared to the input is less then 
#' min_read_count_prop, prints out a warning, since it suggest a 
#' to high cutoff value
#' @returns Filtered read_count_df data frame.
#' @examples
#' \dontrun{
#' filtered_read_count_df <- LFNvariant(read_count_df, cutoff=0.005, min_read_count_prop=0.8)
#' }
#' @export
#' 
LFNvariant2 <- function(read_count, 
                       cutoff=NULL, 
                       asv_specific_cutoff = NULL,
                       by_replicate=TRUE, 
                       outfile="", 
                       sep=",", 
                       min_read_count_prop=0.7){
  
  #### get read_count_df
  if(is.character(read_count)){
    # read known occurrences
    read_count_df <- read.csv(read_count, header=T, sep=sep)
  }else{
    read_count_df <- read_count
  }
  
  ##### check coherence of parameters
  if(is.null(asv_specific_cutoff)){
    if(is.null(cutoff)){
      stop("ERROR: cutoff and asv_specific_cutoff are both NULL, Please, specify at least one of them.")
    }
  }else{
    # make asv_specific_cutoff_df
    if (is.character(asv_specific_cutoff)){
      # read known occurrences
      asv_specific_cutoff_df <- read.csv(asv_specific_cutoff, header=T, sep=sep)
    }else{
      asv_specific_cutoff_df <- asv_specific_cutoff
    }
    
    if(!("replicate" %in% colnames(asv_specific_cutoff_df)) & by_replicate==TRUE){
      stop("ERROR: When by_replicate id TRUE, asv_specific_cutoff should have a replicate column.")
    }
    if("replicate" %in% colnames(asv_specific_cutoff_df) & by_replicate==FALSE){
      stop("ERROR: When by_replicate id FALSE, asv_specific_cutoff should not have a replicate column.")
    }
  }
  

  #### input read count and sample count for leater comparaison
  asvs <- read_count_df %>%
    group_by(asv_id) %>%
    summarize("sample_count_input" = length(sample), "read_count_input"=sum(read_count)) %>%
    filter(read_count_input > 10) %>%
    ungroup()
  
  #### make df with asv total read count 
  if(by_replicate){
    sum_by_asv <- read_count_df %>%
      group_by(asv_id,replicate) %>%
      summarize(asv_sum = sum(read_count), .groups="drop")
  } else{
    sum_by_asv <- read_count_df %>%
      group_by(asv_id) %>%
      summarize(asv_sum = sum(read_count)) %>%
      ungroup()
  }
  
  #### Simple case of fixed cutoff
  if(is.null(asv_specific_cutoff)){
    if(by_replicate){
      read_count_df <- left_join(read_count_df, sum_by_asv, by=c("asv_id", "replicate")) %>%
        filter(read_count/asv_sum >= cutoff)
    } else{
      read_count_df <- left_join(read_count_df, sum_by_asv, by=c("asv_id")) %>%
        filter(read_count/asv_sum >= cutoff)
    }
  }

  #### ASV specific cutoff
  # get sum_by_asv with sum read_cout of asv and threshold
  if(!is.null(asv_specific_cutoff)){

    # add cutoff_asv_spec from input asv_specific_cutoff
    if(by_replicate){
      sum_by_asv <- left_join(sum_by_asv, asv_specific_cutoff_df, by=c("asv_id", "replicate"))
    }else{
      sum_by_asv <- left_join(sum_by_asv, asv_specific_cutoff_df, by=c("asv_id"))
    }
    
    if(is.null(cutoff)){ # no fix cutoff
      sum_by_asv <- sum_by_asv %>%
        mutate(cutoff_asv_spec = if_else(is.na(cutoff_asv_spec), 0, cutoff_asv_spec))
    }
    else{
      ### add fixed cutoff, when not specified in the input asv_specific_cutoff_df
      sum_by_asv <- sum_by_asv %>%
        mutate(cutoff_asv_spec = if_else(is.na(cutoff_asv_spec), cutoff, cutoff_asv_spec))
    }

    ### filter
    if(by_replicate){
      read_count_df <- left_join(read_count_df, sum_by_asv, by = c("asv_id", "replicate"))
    }else{
      read_count_df <- left_join(read_count_df, sum_by_asv, by = c("asv_id"))
    }
    read_count_df <- read_count_df %>%
      filter(read_count/asv_sum >= cutoff_asv_spec) %>%
      select(-cutoff_asv_spec, -asv_sum )
  }

  ###
  # Check if filter do not eliminate occurrences with relatively high readcount 
  ###
  asvs_output <- read_count_df %>%
    group_by(asv_id) %>%
    summarize("sample_count_output" = length(sample), 
              "read_count_output"=sum(read_count)) %>%
    ungroup()
  # join sample and read counts before and after filtering
    asvs <- left_join(asvs, asvs_output, by="asv_id")
  asvs$sample_prop <- asvs$sample_count_output / asvs$sample_count_input
  asvs$read_count_prop <- asvs$read_count_output / asvs$read_count_input
  asvs <- asvs %>%
    filter(read_count_prop<min_read_count_prop) %>%
    arrange(read_count_prop, sample_prop) %>%
    select("asv_id", "read_count_input", "read_count_output", "read_count_prop",
           "sample_count_input", "sample_count_output", "sample_prop")
  
  if(nrow(asvs > 0)){
    cat("WARNING: The following ASVs have lost a high proportion of their 
          reads during this filtering step. 
          The cutoff value of LFNvariant function might need to be reduced.")
    print(asvs)
  }
  
  
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}


tmp <- LFNvariant2(read_count_df, 
                        cutoff=0.8, 
                        asv_specific_cutoff = ASVspecificCutoff_df,
                        by_replicate=TRUE, 
                        outfile="", 
                        sep=",", 
                        min_read_count_prop=0.7)











### MakeKnownOccurrences performance_metrics
results <- MakeKnownOccurrences(read_count_samples_df, 
                                sortedinfo=sortedinfo_df, 
                                mock_composition=mock_composition)


updated_asv_list <- file.path(outdir, "updated_asv_list_end.tsv")
UpdateASVlist(asv_list1 = read_count_samples_df,
              asv_list2 =asv_list, 
              outfile=updated_asv_list
)

### TaxAssign
asv_tax <- TaxAssign(asv=read_count_samples_df, 
                     taxonomy=taxonomy, 
                     blast_db=blast_db, 
                     blast_path=blast_path, 
                     num_threads=num_threads)


plot <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      num_threads=num_threads,
                                      outfile="2_swarm.png")


### WriteASVtable
outfile=file.path(outdir, "Final_asvtable_with_TaxAssign.csv")
asv_table_df <- WriteASVtable(read_count_samples_df, 
                              outfile=outfile, 
                              asv_tax=asv_tax, 
                              sortedinfo=sortedinfo_df, 
                              add_empty_samples=T, 
                              add_sums_by_sample=T, 
                              add_sums_by_asv=T, 
                              add_expected_asv=T, 
                              mock_composition=mock_composition)


#####################
#####################
#####################
# make mOTU

#####################
### mOTU with swarm

stat_df <- GetStat(read_count_samples_df, stat_df, stage="PoolReplicates", params=NA)

by_sample <- FALSE
d = 7
read_count_df_swarm_motu <- Swarm(read_count_samples_df, 
                       swarm_path=swarm_path, 
                       swarm_d=d,
                       fastidious=FALSE,
                       num_threads=num_threads, 
                       by_sample=by_sample)

stat_df <- GetStat(read_count_df_swarm_motu, stat_df, stage="swarm_motu_7", params=d)

### mOTU with ClusterSize
identity <- 0.97 
read_count_df_clustersize_motu <- ClusterSize(read_count_samples_df,
                                     id=identity, 
                                     vsearch_path=vsearch_path,
                                     num_threads=num_threads,
                                     by_sample=FALSE)

stat_df <- GetStat(read_count_df_clustersize_motu, stat_df, stage="clustersize_motu_7", params=identity)



### ClusterSize
identity <- 0.97 
read_count_samples_df <- ClusterSize(read_count_samples_df,
                                     id=identity, 
                                     vsearch_path=vsearch_path,
                                     by_sample=FALSE)

### TaxAssign
asv_tax <- TaxAssign(asv=read_count_samples_df, 
                     taxonomy=taxonomy, 
                     blast_db=blast_db, 
                     blast_path=blast_path, 
                     num_threads=num_threads)

### WriteASVtable
outfile=file.path(outdir, "Final_asvtable_with_TaxAssign.csv")
asv_table_df <- WriteASVtable(read_count_samples_df, 
                              outfile=outfile, 
                              asv_tax=asv_tax, 
                              sortedinfo=sortedinfo_df, 
                              add_empty_samples=T, 
                              add_sums_by_sample=T, 
                              add_sums_by_asv=T, 
                              add_expected_asv=T, 
                              mock_composition=mock_composition)




#################
#################
################
# test ClusterSize after dereplicate
################

#### ClusterSize by sample
## sample_replicate
by_sample <- TRUE
read_count_df_ClusterSize_sample_replicate <- ClusterSize(read_count_df, 
                                                          id= 0.97,
                                                          vsearch_path=vsearch_path, 
                                                          num_threads=num_threads, 
                                                          by_sample=by_sample)


stat_df <- GetStat(read_count_df_ClusterSize_sample_replicate, stat_df, stage="ClusterSize_sample_replicate", params=by_sample)

## sample
read_count_df_ClusterSize_sample <- ClusterSize(read_count_sample_df, 
                                          id= 0.97,
                                          vsearch_path=vsearch_path, 
                                          num_threads=num_threads, 
                                          by_sample=by_sample)


stat_df <- GetStat(read_count_df_ClusterSize_sample, stat_df, stage="ClusterSize_sample", params=by_sample)

#### ClusterSize by sample=FALSE

## sample_replicate
by_sample <- FALSE
read_count_df_ClusterSize_all_sample_replicate <- ClusterSize(read_count_df, 
                                                              id= 0.97,
                                                              vsearch_path=vsearch_path, 
                                                              num_threads=num_threads, 
                                                              by_sample=by_sample)

stat_df <- GetStat(read_count_df_ClusterSize_all_sample_replicate, stat_df, stage="ClusterSize_all_sample_replicate", params=by_sample)

## sample
read_count_df_ClusterSize_all_sample <- ClusterSize(read_count_sample_df, 
                                             id= 0.97,
                                             vsearch_path=vsearch_path, 
                                       num_threads=num_threads, 
                                       by_sample=by_sample)


stat_df <- GetStat(read_count_df_ClusterSize_all_sample, stat_df, stage="ClusterSize_all_sample", params=by_sample)





################
# test swarm after dereplicate
################
#### Swarm by sample
## sample_replicate
by_sample <- TRUE
read_count_df_swam_sample_replicate <- Swarm(read_count_df, 
                       swarm_path=swarm_path, 
                       num_threads=num_threads, 
                       by_sample=by_sample
)
stat_df <- GetStat(read_count_df_swam_sample_replicate, stat_df, stage="Swarm_sample_replicate", params=by_sample)

## sample
read_count_df_swam_sample <- Swarm(read_count_sample_df, 
                                   swarm_path=swarm_path, 
                                   num_threads=num_threads, 
                                   by_sample=by_sample)


stat_df <- GetStat(read_count_df_swam_sample, stat_df, stage="Swarm_sample", params=by_sample)

#### Swarm by sample=FALSE

## sample_replicate
by_sample <- FALSE
read_count_df_swam_all_sample_replicate <- Swarm(read_count_df, 
                                             swarm_path=swarm_path, 
                                             num_threads=num_threads, 
                                             by_sample=by_sample
)
stat_df <- GetStat(read_count_df_swam_all_sample_replicate, stat_df, stage="Swarm_all_sample_replicate", params=by_sample)

## sample
read_count_df_swam_all_sample <- Swarm(read_count_sample_df, 
                                   swarm_path=swarm_path, 
                                   num_threads=num_threads, 
                                   by_sample=by_sample)


stat_df <- GetStat(read_count_df_swam_all_sample, stat_df, stage="Swarm_all_sample", params=by_sample)



plot <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      num_threads=num_threads,
                                      outfile="density_plot_1_15_3.png")

tmp <- PairwiseIdentity(read_count_df, 
                             min_id = 0.8, 
                             num_threads=num_threads,
                             vsearch_path=vsearch_path)
