

## load
library(vtamR)
library(dplyr)
library(gplot2)

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
num_threads <- 0
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
                      fastq_allowmergestagger=F,
                      num_threads=num_threads
                      
)

### demultiplex

sorted_dir <- file.path(outdir, "sorted")
sortedinfo_df <- SortReads(fastainfo_df, 
                           fasta_dir=merged_dir, 
                           outdir=sorted_dir, 
                           check_reverse=TRUE, 
                           cutadapt_path=cutadapt_path, 
                           vsearch_path=vsearch_path,
                           num_threads=num_threads
)

###############
### dereplicate
###############
updated_asv_list <- file.path(outdir, "updated_asv_list.tsv")
read_count_df <- Dereplicate(sortedinfo_df, 
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
by_sample <- FALSE

read_count_df <- Swarm(read_count_df, 
                       swarm_path=swarm_path, 
                       num_threads=num_threads, 
                       by_sample=by_sample)
stat_df <- GetStat(read_count_df, stat_df, stage="Swarm", params=NA)

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

### FilterChimera
abskew=2
by_sample = T
sample_prop = 0.8
read_count_df <- FilterChimera(read_count_df, 
                               vsearch_path=vsearch_path, 
                               num_threads=num_threads,
                               by_sample=by_sample, 
                               sample_prop=sample_prop, 
                               abskew=abskew)

#### FilterRenkonen
cutoff <- 0.4
read_count_df <- FilterRenkonen(read_count_df, 
                                cutoff=cutoff)

### FilterPCRerror
pcr_error_var_prop <- 0.05
max_mismatch <- 2
read_count_df <- FilterPCRerror(read_count_df, 
                                vsearch_path=vsearch_path, 
                                num_threads=num_threads,
                                pcr_error_var_prop=pcr_error_var_prop, 
                                max_mismatch=max_mismatch)

### LFNsampleReplicate
lfn_sample_replicate_cutoff <- 0.004
read_count_df <- LFNsampleReplicate(read_count_df, 
                                    cutoff=lfn_sample_replicate_cutoff)

### FilterMinReplicate
min_replicate_number <- 2
read_count_df <- FilterMinReplicate(read_count_df, 
                                    cutoff=min_replicate_number)

### LFNvariant
lnf_variant_cutoff = 0.001
read_count_df_lnf_variant <- LFNvariant(read_count_df, 
                                        cutoff=lnf_variant_cutoff)
### LFNreadCount
lfn_read_count_cutoff <- 10
read_count_df_lfn_read_count <- LFNreadCount(read_count_df, 
                                             cutoff=lfn_read_count_cutoff)
### Combine results
read_count_df <- PoolFilters(read_count_df_lfn_read_count, 
                             read_count_df_lnf_variant)

# delete temporary data frames
rm(read_count_df_lfn_read_count)
rm(read_count_df_lnf_variant)

### FilterMinReplicate
min_replicate_number <- 2
read_count_df <- FilterMinReplicate(read_count_df, 
                                    cutoff=min_replicate_number)

#### PoolReplicates
read_count_samples_df <- PoolReplicates(read_count_df)

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
                                      swarm_d_increment=1,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      num_threads=0,
                                      outfile="density_plot_1_15.png")

tmp <- PairwiseIdentity(read_count_df, 
                             min_id = 0.8, 
                             vsearch_path=vsearch_path)
