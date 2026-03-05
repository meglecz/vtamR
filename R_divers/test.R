

## load
library(vtamR)
library(dplyr)
library(ggplot2)
library(rRDP)
library(rRDPData)

setwd("/home/meglecz/vtamR/")
setwd("C:/Users/emese/vtamR")
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
pigz_path <- "pigz"
sep <- ","
outdir <- "~/vtamR_demo_out"

fastq_dir <- system.file("extdata/demo/fastq", package = "vtamR")
fastqinfo <-  system.file("extdata/demo/fastqinfo.csv", package = "vtamR")
mock_composition <-  system.file("extdata/demo/mock_composition.csv", package = "vtamR")
asv_list <-  system.file("extdata/demo/asv_list.csv", package = "vtamR")
taxonomy <- system.file("extdata/db_test/taxonomy_reduced.tsv", package = "vtamR")
blast_db <- system.file("extdata/db_test", package = "vtamR")
blast_db <- file.path(blast_db, "COInr_reduced")

taxonomy_COInr <- "/home/meglecz/mkCOInr/COInr/COInr_for_vtam_2025_05_23_dbV5/COInr_for_vtam_taxonomy.tsv"




### Merge
merged_dir_uncompress <- file.path(outdir, "merged_uncompress")
fastainfo_df_uncompress <- Merge(fastqinfo, 
                      fastq_dir=fastq_dir, 
                      compress_method="pigz",
                      pigz_path = pigz_path,
                      vsearch_path=vsearch_path, 
                      outdir=merged_dir_uncompress,
                      fastq_maxee=1,
                      fastq_maxns=0,
                      fastq_allowmergestagger=F,
                      compress=FALSE)


outdir_RandomSeq <- file.path(outdir, "RandomSeqR_gz")
fastainfo_randomSeq <-RandomSeq(fastainfo_df_uncompress, 
           n=40000,
           fasta_dir=merged_dir_uncompress,
           outdir=outdir_RandomSeq, 
           randseed=0, 
           quiet=TRUE)


### demultiplex
demultiplexed_dir <- file.path(outdir, "demultiplexed")
sampleinfo_df <- SortReads(fastainfo_randomSeq, 
                           fasta_dir=outdir_RandomSeq, 
                           outdir=demultiplexed_dir, 
                           check_reverse=TRUE, 
                           cutadapt_path=cutadapt_path, 
                           vsearch_path=vsearch_path,
                           tag_to_end = T,
                           primer_to_end=T)

###############
### dereplicate
###############
updated_asv_list <- file.path(outdir, "updated_asv_list.tsv")
read_count_df <- Dereplicate(sampleinfo_df, 
                             dir=demultiplexed_dir, 
                             input_asv_list=asv_list,
                             output_asv_list = updated_asv_list)

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


read_count_df_backup <- read_count_df

### FilterExternalContaminant
conta_file <- file.path(outdir, "tmp", "external_contamination.csv")
read_count_df <- FilterExternalContaminant(read_count_df, 
                                           sampleinfo=sampleinfo_df,
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



#####################################
PoolReplicates_csv <- file.path(outdir, "PoolReplicates", "PoolReplicates_mean.csv")
read_count_sample_mean <- PoolReplicates(read_count_df, digits=0, outfile = PoolReplicates_csv)
read_count_sample_fun_mean <- PoolReplicates(read_count_df, digits=0, outfile = PoolReplicates_csv)
read_count_sample_fun_max <- PoolReplicates(read_count_df, method="max", digits=0, outfile = PoolReplicates_csv)
read_count_sample_fun_sum <- PoolReplicates(read_count_df, method="sum", digits=0, outfile = PoolReplicates_csv)

identical(read_count_sample_fun_max, read_count_sample_fun_mean)


read_count_cluster_df <- ClusterASV(read_count_df,
                            method = "swarm",
                            swarm_d=7,
                            fastidious=fastidious,
                            by_sample=FALSE,
                            group = FALSE,
                            path=swarm_path
                            )

read_count_cluster_df_max <- PoolReplicates(read_count_cluster_df, method="max")
#####################################

tmp_replicate <- WriteASVtable(read_count_df, 
                          sampleinfo=sampleinfo_df, 
                          pool_replicates=FALSE,
                          add_sums_by_asv=TRUE
                     )
tmp_sample <- WriteASVtable(read_count_df, 
                               sampleinfo=sampleinfo_df, 
                               pool_replicates=TRUE,
                               add_sums_by_asv=TRUE
                            )

tmp_sample_new <- WriteASVtable(read_count_df, 
                            sampleinfo=sampleinfo_df, 
                            pool_replicates=TRUE,
                            method="max",
                            add_sums_by_asv=TRUE,
                            add_sums_by_sample=FALSE
)
#####################################

opt <- OptimizeLFNreadCountLFNvariant(read_count_df, 
                               outdir=file.path(outdir, "optimize"),
                               known_occurrences = NULL, 
                               mock_composition = mock_composition_df,
                               sampleinfo = sampleinfo_df,
                               habitat_proportion = 0.5,
                               min_replicate_number=2, 
                               quiet=T
                               )

mock_composition_df <- read.csv(mock_composition)

mock_composition_df[7,] <- c("tpos1", "keep", "TCTATACCTTATTTTCGGCGCAATTTCAGGTATTGCAGGTACCGCTTTATCTCTTTACATTCGAATTACTTTATCTCAACCTAATGGTAATTTTTTAGAATACAATCACCACTTTTATAATGTGATTATAACGGGTCACGCTCTTCTTATGATTTTTTTCATGGTAATGCCAATCTTGATT", NA, NA)

read_count_samples <- PoolReplicates(read_count_df)
missing <- make_missing_occurrences(read_count_samples, mock_composition_df, quiet = FALSE)


#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' 
#'     missing_asv_text <- paste(capture.output(print(df)), collapse = "\n")
warning(
  paste0(
    "\n  Some expected ASVs are missing from the mock samples.\n",
    "----------------------------------------------------------\n",
    missing_asv_text,
    "\n----------------------------------------------------------\n"
  ),
  call. = FALSE
)
#####################################
  



#### FilterRenkonen
cutoff <- 0.4
read_count_df <- FilterRenkonen(read_count_df, 
                                cutoff=cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterRenkonen", params=cutoff)


### FilterPCRerror

mock_composition_df <- read.csv(mock_composition)

mock_composition_df[5,"asv"] <- "TTTTTTTTTTTT"
mock_composition_df[6,"asv"] <- "AAAAAAAAAa"

dir_opt = file.path(outdir, "OptimizeLFNreadCountLFNvariant") 
opt_rc_var <- OptimizeLFNreadCountLFNvariant(read_count_df, 
                                           outdir = dir_opt,
                                           known_occurrences = NULL, 
                                           mock_composition = mock_composition_df,
                                           sampleinfo = sampleinfo_df,
                                           habitat_proportion = 0.5,
                                           sep=",",
                                           min_lfn_read_count_cutoff=10, 
                                           max_lfn_read_count_cutoff=100, 
                                           increment_lfn_read_count_cutoff=5, 
                                           min_lnf_variant_cutoff=0.001, 
                                           max_lnf_variant_cutoff=0.01, 
                                           increment_lnf_variant_cutoff=0.001, 
                                           by_replicate=FALSE, 
                                           min_replicate_number=2, 
                                           quiet=T)

opt_rc_var <- OptimizeLFNreadCountLFNvariant(read_count_df, 
                                             outdir = dir_opt,
                                             known_occurrences = "/home/meglecz/vtamR_demo_out/OptimizeLFNreadCountLFNvariant/known_occurrences.csv", 
                                             sep=",",
                                             min_lfn_read_count_cutoff=10, 
                                             max_lfn_read_count_cutoff=100, 
                                             increment_lfn_read_count_cutoff=5, 
                                             min_lnf_variant_cutoff=0.001, 
                                             max_lnf_variant_cutoff=0.01, 
                                             increment_lnf_variant_cutoff=0.001, 
                                             by_replicate=FALSE, 
                                             min_replicate_number=2, 
                                             quiet=T)




optPCR <- OptimizePCRerror(read_count=read_count_df, 
                 mock_composition=mock_composition_df, 
                 vsearch_path=vsearch_path, 
                 max_mismatch=2, 
                 min_read_count=5)

optLFN_sample <- OptimizeLFNsampleReplicate(read_count_df, mock_composition=mock_composition_df)

optLFN_var <- OptimizeLFNreadCountLFNvariant(read_count_df, 
                                           known_occurrences, 
                                           sep=",",
                                           outfile="", 
                                           min_lfn_read_count_cutoff=10, 
                                           max_lfn_read_count_cutoff=100, 
                                           increment_lfn_read_count_cutoff=5, 
                                           min_lnf_variant_cutoff=0.001, 
                                           max_lnf_variant_cutoff=0.01, 
                                           increment_lnf_variant_cutoff=0.001, 
                                           by_replicate=FALSE, 
                                           min_replicate_number=2, 
                                           quiet=T
)


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

read_count_df_backup <- read_count_df

plot_swarm <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      outfile="13_FilterMinReplicate.png")





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


###################################
###################################
# test TaxAssignRDP
file <- "/home/meglecz/vtamR/tmp/test_rdp/9_PoolReplicates_TAS2_16S.csv"
read_count_16S <- read.csv(file)

taxa <- TaxAsssignRDP(asv=read_count_16S, confidence=0.7, max_memory=8, rm_chloroplast=FALSE)

plot_vsearch <- PairwiseIdentityPlotPerClusterIdentityThreshold(read_count_16S,
                                                                identity_min=0.9,
                                                                identity_max=0.99,
                                                                identity_increment=0.01,
                                                                min_id = 0.8, 
                                                                vsearch_path=vsearch_path)

scatterplot_vsearch <- PlotClusterClasstification(
  read_count=read_count_16S,
  taxa=taxa,
  clustering_method="vsearch",
  cluster_params=c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99),
  vsearch_path=vsearch_path,
  taxlevels= c("species", "genus"))



read_count_samples_df_ClusterSize <- ClusterASV(read_count_16S,
                                                method = "vsearch",
                                                identity = 0.97,
                                                group = TRUE,
                                                by_sample=FALSE,
                                                path=vsearch_path)

dim(read_count_16S)
dim(read_count_samples_df_ClusterSize)

meta_info <- "/home/meglecz/vtamR/tmp/test_rdp/fastqinfo.csv"
mock <- "/home/meglecz/vtamR/tmp/test_rdp/mock_composition.csv"
out <- "/home/meglecz/vtamR/tmp/test_rdp/asv_table.csv"
asv_table_df <- WriteASVtable(read_count_16S, 
                              asv_tax=taxa, 
                              sortedinfo=meta_info, 
                              pool_replicates=TRUE,
                              add_empty_samples=T, 
                              add_sums_by_sample=T, 
                              add_sums_by_asv=T, 
                              add_expected_asv=T, 
                              outfile=out,
                              mock_composition=mock
)
