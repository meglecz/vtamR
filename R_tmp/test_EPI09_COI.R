

## load
library(vtamR)
#library(dplyr)
#library(ggplot2)

setwd("/home/meglecz/vtamR/")
#library("devtools")
#library("roxygen2")
#load_all(".")
#roxygenise()
#usethis::use_roxygen_md()

### set up
cutadapt_path <- "~/miniconda3/envs/vtam/bin/cutadapt"
vsearch_path <- "~/miniconda3/envs/vtam/bin/vsearch"
blast_path <- "~/miniconda3/envs/vtam/bin/blastn"
swarm_path <- "swarm"
sep <- ","
outdir <- "~/vtamR_test_EPI09_COI"

fastq_dir <- "/home/meglecz/vtamR_large_files/EPI09"
fastqinfo <-  "/home/meglecz/vtamR_large_files/EPI09/metainfo/fastqinfo_Epi09_COI.csv"
mock_composition <-  "/home/meglecz/vtamR_large_files/EPI09/metainfo/mock_composition_EPI09_COI.csv"
taxonomy <- "/home/meglecz/mkCOInr/COInr/COInr_for_vtam_2025_05_23_dbV5/COInr_for_vtam_taxonomy.tsv"
blast_db <- "/home/meglecz/mkCOInr/COInr/COInr_for_vtam_2025_05_23_dbV5/COInr_for_vtam"
asv_list <- "/home/meglecz/vtamR_large_files/EPI09/metainfo/asv_list_COI.csv"

### Merge
t1 <- proc.time()
merged_dir <- file.path(outdir, "merged")
fastainfo_df <- Merge(fastqinfo, 
                      fastq_dir=fastq_dir, 
                      vsearch_path=vsearch_path, 
                      outdir=merged_dir,
                      fastq_maxee=1,
                      fastq_maxns=0,
                      fastq_allowmergestagger=F, quiet=FALSE)

t2 <- proc.time()
t <- t2 - t1
time_df <- data.frame(
  Step = "Merge",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE)

### demultiplex
t1 <- proc.time()
demultiplexed_dir <- file.path(outdir, "demultiplexed")
sampleinfo_df <- SortReads(fastainfo_df, 
                           fasta_dir=merged_dir, 
                           outdir=demultiplexed_dir, 
                           check_reverse=TRUE, 
                           cutadapt_path=cutadapt_path, 
                           vsearch_path=vsearch_path,
                           tag_to_end = T,
                           primer_to_end=T)
t2 <- proc.time()
t <- t2 - t1

time_df <- rbind(time_df, data.frame(
  Step = "SortReads",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

###############
### dereplicate
###############
t1 <- proc.time()
updated_asv_list <- file.path(outdir, "updated_asv_list.tsv")
demultiplexed_dir <- file.path(outdir, "demultiplexed")
sampleinfo <- file.path(demultiplexed_dir, "sampleinfo.csv")
outfile <- file.path(outdir, "filter", "1_input_dereplicated.csv")
read_count_df <- Dereplicate(sampleinfo, 
                             dir=demultiplexed_dir, 
                             input_asv_list=asv_list,
                             output_asv_list = updated_asv_list,
                             outfile=outfile)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "Dereplicate",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

### stat
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())

stat_df <- GetStat(read_count_df, stat_df, stage="input_sample_replicate", params=NA)

#### swarm
t1 <- proc.time()
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
                            path=swarm_path,
                            outfile=outfile
                            )


stat_df <- GetStat(read_count_df, stat_df, stage="Swarm", params=by_sample)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "ClusterASV_SWARM",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

#### LFNglobalReadCount
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "3_LFNglobalReadCount.csv")
global_read_count_cutoff = 2

read_count_df <- LFNglobalReadCount(read_count_df, 
                                    cutoff=global_read_count_cutoff,
                                    outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNglobalReadCount", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "LFNglobalReadCount",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

#### FilterIndel
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "4_FilterIndel.csv")
read_count_df <- FilterIndel(read_count_df, 
                             outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterIndel", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "FilterIndel",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

### FilterCodonStop
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "5_FilterCodonStop.csv")
genetic_code = 5
read_count_df <- FilterCodonStop(read_count_df, 
                                 genetic_code=genetic_code,
                                 outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterCodonStop", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "FilterCodonStop",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

### FilterExternalContaminant
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "6_FilterExternalContaminant.csv")
conta_file <- file.path(outdir, "filter", "external_contamination.csv")
read_count_df <- FilterExternalContaminant(read_count_df, 
                          sample_types=sampleinfo, 
                          conta_file=conta_file,
                          outfile=outfile)

stat_df <- GetStat(read_count_df, stat_df, stage="FilterExternalContaminant", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "FilterExternalContaminant",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

### FilterChimera
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "7_FilterChimera.csv")
abskew=2
by_sample = T
sample_prop = 0.8
read_count_df <- FilterChimera(read_count_df, 
                               vsearch_path=vsearch_path, 
                               by_sample=by_sample, 
                               sample_prop=sample_prop, 
                               abskew=abskew,
                               outfile=outfile)

stat_df <- GetStat(read_count_df, stat_df, stage="FilterChimera", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "FilterChimera",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

#### FilterRenkonen
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "8_FilterRenkonen.csv")
cutoff <- 0.4
read_count_df <- FilterRenkonen(read_count_df, 
                                cutoff=cutoff,
                                outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterRenkonen", params=cutoff)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "FilterRenkonen",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

### OptimizePCRerror
t1 <- proc.time()
outfile <- file.path(outdir, "optimize", "OptPCRerror.csv")
OptPCR <- OptimizePCRerror(read_count_df, 
                             mock_composition=mock_composition, 
                             vsearch_path= vsearch_path, 
                             outfile=outfile, 
                             max_mismatch=1, 
                             min_read_count=10)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "OptimizePCRerror",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

### FilterPCRerror
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "9_FilterPCRerror.csv")
pcr_error_var_prop <- 0.1
max_mismatch <- 1
read_count_df <- FilterPCRerror(read_count_df, 
                                vsearch_path=vsearch_path, 
                                pcr_error_var_prop=pcr_error_var_prop, 
                                max_mismatch=max_mismatch)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterPCRerror", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "FilterPCRerror",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

#### OptimizeLFNsampleReplicate
t1 <- proc.time()
outfile <- file.path(outdir, "optimize", "OptimizeLFNsampleReplicate.csv")
optLFN_sample <- OptimizeLFNsampleReplicate(read_count_df, 
                                            mock_composition=mock_composition,
                                            outfile=outfile)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "OptimizeLFNsampleReplicate",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))


### LFNsampleReplicate
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "10_LFNsampleReplicate.csv")
lfn_sample_replicate_cutoff <- 0.004
read_count_df <- LFNsampleReplicate(read_count_df, 
                                    cutoff=lfn_sample_replicate_cutoff,
                                    outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNsampleReplicate", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "LFNsampleReplicate",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

#### OptimizeLFNreadCountLFNvariant
t1 <- proc.time()
dir_opt <- file.path(outdir, "optimize")
OptLFNreadCountLFNvariant <- OptimizeLFNreadCountLFNvariant(read_count_df, 
                                           outdir = dir_opt,
                                           mock_composition = mock_composition,
                                           sampleinfo = sampleinfo_df,
                                           habitat_proportion = 0.5,
                                           by_replicate=TRUE, 
                                           min_replicate_number=2)

t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "OptimizeLFNreadCountLFNvariant",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

### LFNvariant
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "11_LFNvariant.csv")
lnf_variant_cutoff = 0.001
read_count_df_lnf_variant <- LFNvariant(read_count_df, 
                                        cutoff=lnf_variant_cutoff,
                                        outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNvariant", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "LFNvariant",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

### LFNreadCount
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "12_LFNreadCount.csv")
lfn_read_count_cutoff <- 10
read_count_df_lfn_read_count <- LFNreadCount(read_count_df, 
                                             cutoff=lfn_read_count_cutoff,
                                             outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNreadCount", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "LFNreadCount",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))




### Combine results
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "13_PoolFilters.csv")
read_count_df <- PoolFilters(read_count_df_lfn_read_count, 
                             read_count_df_lnf_variant,
                             outfile=outfile)
# delete temporary data frames
rm(read_count_df_lfn_read_count)
rm(read_count_df_lnf_variant)

stat_df <- GetStat(read_count_df, stat_df, stage="Combine results", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "PoolFilters",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))


### FilterMinReplicate
t1 <- proc.time()
outfile <- file.path(outdir, "filter", "14_FilterMinReplicate.csv")
min_replicate_number <- 2
read_count_df <- FilterMinReplicate(read_count_df, 
                                    cutoff=min_replicate_number,
                                    outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterMinReplicate", params=NA)
t2 <- proc.time()
t <- t2 - t1
time_df <- rbind(time_df, data.frame(
  Step = "PoolFilters",
  user = t["user.self"],
  system = t["sys.self"],
  elapsed = t["elapsed"],
  stringsAsFactors = FALSE
))

outfile <- file.path(outdir, "stat_time.csv")
write.csv(time_df, file=outfile, row.names = FALSE)

outfile <- file.path(outdir, "stat_count.csv")
write.csv(stat_df, file=outfile, row.names = FALSE)





plot_swarm <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      outfile="13_FilterMinReplicate.png")






ASVspecificCutoff_df_round <- ASVspecificCutoff(read_count_df,  mock_composition=mock_composition,
                              by_replicate=FALSE, 
                              outfile="tmp/ASVspecificCutoff_by_replicate_false.csv")






tmp2 <- LFNvariant2(read_count_df, 
                        asv_specific_cutoff = NULL,
                        cutoff=0.01,
                        by_replicate=FALSE, 
                        outfile="", 
                        sep=",", 
                        min_read_count_prop=0.7)

tmp <- LFNvariant(read_count_df, 
                    cutoff=0.01,
                    by_replicate=FALSE, 
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
