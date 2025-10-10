

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




#' Make ASV specific cutoff values
#' 
#' For all ASVs that have a know false positive occurrences, 
#' calculate an ASV specific cutoff value for the LFNvariant filter.
#' 
#' For each ASV take the maximum read count among its false positive occurrences
#' and divide it by the total number of reads of the ASV. #' 
#'  
#' @param read_count Data frame or csv file with the following variables: 
#' asv, sample, replicate, read_count.
#' @param mock_composition Data frame or csv file with columns: 
#' sample, action (keep/tolerate), asv.
#' @param sep Field separator character in input and output csv files.
#' @param out Character string: output file. If empty, no file is written. 
#' @returns Data frame with the following columns: 
#' asv_id, cutoff
#' @examples
#' \dontrun{
#' make_missing_occurrences(read_count_samples=read_count_samples_df, 
#'     mock_composition="data/mock_composition.csv"
#'     )
#' }
#' @export
#'
ASVspecificCutoff <- function(read_count, 
                              by_replicate=FALSE, 
                              min_read_count_prop=0.7,
                              mock_composition="",
                              sortedinfo="",
                              habitat_proportion,
                              outfile="", 
                              sep=",", 
                              outfile="")  {
  
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
    asv_spec_cutoff <- read_count_df %>%
      mutate(asv_sample = paste(asv_id, sample, sep=".")) %>%
      filter(asv_sample %in% delete_occurrences_df$asv_sample) %>%
      group_by(asv_id, replicate) %>%
      filter(read_count==max(read_count))%>%
      ungroup() %>%
      left_join(asv_total_rc, by=c("asv_id", "replicate")) %>%
      mutate(asv_specific_cutoff= read_count/total_rc) %>%
      select(asv_id, replicate, asv_specific_cutoff)
      
  }else{
    asv_spec_cutoff <- read_count_df %>%
      mutate(asv_sample = paste(asv_id, sample, sep=".")) %>%
      filter(asv_sample %in% delete_occurrences_df$asv_sample) %>%
      group_by(asv_id) %>%
      filter(read_count==max(read_count))%>%
      ungroup() %>%
      left_join(asv_total_rc, by=c("asv_id")) %>%
      mutate(asv_specific_cutoff= read_count/total_rc)%>%
      select(asv_id, asv_specific_cutoff)
  }

  # write to outfile
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(asv_spec_cutoff, file=outfile, row.names = F, sep=sep)
  }
  
  return(asv_spec_cutoff)
}




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
