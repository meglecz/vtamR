if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install("meglecz/vtamR@develop")

environment_vtamR_yml_path <- system.file("environment_vtamR.yml", package = "vtamR")
cat(environment_vtamR_yml_path)



################################################"

### Win
setwd("C:/Users/emese/vtamR")
cutadapt_path <- "C:/Users/Public/cutadapt"
vsearch_path <- "C:/Users/Public/vsearch-2.23.0-win-x86_64/bin/vsearch"
blast_path <- "C:/Users/Public/blast-2.16.0+/bin/blastn"
swarm_path <- "C:/Users/Public/swarm-3.1.5-win-x86_64/bin/swarm"
num_threads <- 4
sep <- ","
outdir <- "C:/Users/emese/vtamR_demo"


### set up Endoume
cutadapt_path <- "/home/emese/anaconda3/envs/vtamRenv/bin/cutadapt"
vsearch_path <- "/home/emese/anaconda3/envs/vtamRenv/bin/vsearch"
blast_path <- "/home/emese/anaconda3/envs/vtamRenv/bin/blastn"
swarm_path <- "/home/emese/anaconda3/envs/vtamRenv/bin/swarm"
num_threads <- 8
sep <- ","
outdir <- "~/vtamR_demo_out"
setwd("/home/emese/vtamR/")

library(dplyr)
library(ggplot2)
library("devtools")
library("roxygen2")
load_all(".")
roxygenise()
usethis::use_roxygen_md()


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
sortedinfo <- file.path(outdir, "sorted", "sortedinfo.csv")
sorted_dir <- file.path(outdir, "sorted")
updated_asv_list <- file.path(outdir, "updated_asv_list.tsv")
read_count_df <- Dereplicate(sortedinfo, 
                             dir=sorted_dir, 
                             asv_list=asv_list,
                             updated_asv_list = updated_asv_list
)

is_grouped_df(read_count_df)


### stat
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())

stat_df <- GetStat(read_count_df, stat_df, stage="input_sample_replicate", params=NA)

#### swarm denoising
by_sample <- TRUE
group <- TRUE
method <- "swarm"
swarm_d <- 1
fastidious <- T
read_count_df1 <- ClusterASV(read_count_df, 
                        group=group,
                        method=method,
                        swarm_d=swarm_d,
                        fastidious=fastidious,
                        path=swarm_path, 
                       num_threads=num_threads, 
                       by_sample=by_sample)
par = paste(by_sample, group, swarm_d, fastidious, method, sep=";")
stat_df <- GetStat(read_count_df1, stat_df, stage="ClusterASV", params=par)

is_grouped_df(read_count_df1)

#### LFNglobalReadCount
global_read_count_cutoff = 2

read_count_df <- LFNglobalReadCount(read_count_df, 
                                    cutoff=global_read_count_cutoff)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNglobalReadCount", params=global_read_count_cutoff)

#### FilterIndel
read_count_df <- FilterIndel(read_count_df)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterIndel", params=NA)

### FilterCodonStop
genetic_code = 5
read_count_df <- FilterCodonStop(read_count_df, 
                                 genetic_code=genetic_code)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterCodonStop", params=NA)

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
par= paste(abskew, by_sample, sample_prop, sep=";")
stat_df <- GetStat(read_count_df, stat_df, stage="FilterChimera", params=par)

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


### Make mOTU
plot <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                      swarm_d_min=1, 
                                      swarm_d_max=15,
                                      swarm_d_increment=3,
                                      min_id = 0.8, 
                                      vsearch_path=vsearch_path, 
                                      swarm_path=swarm_path,
                                      num_threads=num_threads)

plot_vsearch <- PairwiseIdentityPlotPerClusterIdentityThreshold(read_count_df, 
                                                                identity_min=0.9, 
                                                                identity_max=0.99,
                                                                identity_increment=0.01,
                                                                min_id = 0.8, 
                                                                vsearch_path=vsearch_path, 
                                                                num_threads=num_threads)

group=TRUE
by_sample=FALSE
method="vsearch"
identity= 0.95
read_count_df_grouped <- ClusterASV(read_count_df, 
                            group=group,  
                            by_sample=by_sample, 
                            identity=identity,
                            method="vsearch",
                            path=vsearch_path, 
                            num_threads=num_threads)

par <- paste(group, by_sample, method, identity)
stat_df <- GetStat(read_count_df_grouped, stat_df, stage="ClusterASVs", params=par)

group=FALSE
by_sample=FALSE
method="vsearch"
identity= 0.95
read_count_df_asv <- ClusterASV(read_count_df, 
                                    group=group,  
                                    by_sample=by_sample, 
                                    identity=identity,
                                    method="vsearch",
                                    path=vsearch_path, 
                                    num_threads=num_threads)

par <- paste(group, by_sample, method, identity)
stat_df <- GetStat(read_count_df_asv, stat_df, stage="ClusterASVs", params=par)


### MakeKnownOccurrences performance_metrics
results <- MakeKnownOccurrences(read_count_df_asv, 
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






cluster_df <- GetClusterIdSwarm(read_count_df, 
                                swarm_d=5, 
                                fastidious=FALSE,
                                swarm_path=swarm_path, 
                                num_threads=num_threads, 
                                quiet=FALSE)

asv_tax <- TaxAssign(asv=read_count_df, 
                     taxonomy=taxonomy, 
                     blast_db=blast_db, 
                     blast_path=blast_path, 
                     num_threads=num_threads,
                     quiet=FALSE)

cluster_class <- ClassifyClusters(cluster_df, asv_tax, outfile="", sep=",", quiet=TRUE, 
                                  taxlevels=c("domain", "phylum", "class", "order","family", "genus", "species"))

classplot <- PlotClusterClasstification(read_count_df, asv_tax, 
                                        clustering_method="vsearch", 
                                        cluster_params=c(90, 95, 97), 
                                        vsearch_path=vsearch_path, 
                                        swarm_path=swarm_path, 
                                        taxlevels= c("species", "genus"),
                                        quiet = FALSE)



pairID <- PairwiseIdentity(read_count_df, 
                           min_id = 0.8, 
                           vsearch_path=vsearch_path, 
                           num_threads=num_threads,
                           quiet=TRUE)

p <- PairwiseIdentityPlotPerSwarmD(read_count_df, 
                                   swarm_d_min=2, 
                                   swarm_d_max=16,
                                   swarm_d_increment=2,
                                   min_id = 0.8, 
                                   vsearch_path=vsearch_path, 
                                   swarm_path=swarm_path,
                                   num_threads=num_threads,
                                   quiet=TRUE)

