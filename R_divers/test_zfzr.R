
computer = "Endoume" # Bombyx/Endoume/Win
if(computer == "Bombyx"){
  setwd("/home/meglecz/vtamR/")
  cutadapt_path <- "~/miniconda3/envs/vtam/bin/cutadapt"
  vsearch_path <- "~/miniconda3/envs/vtam/bin/vsearch"
  blast_path <- "~/miniconda3/envs/vtam/bin/blastn"
  swarm_path <- "swarm"
  pigz_path <- "pigz"
  sep <- ","
  outdir <- "~/vtamR_demo_out_mfzr"
}else if(computer == "Endoume"){
  setwd("/home/emese/vtamR")
  cutadapt_path <- "cutadapt"
  vsearch_path <- "vsearch"
  blast_path <- "blastn"
  swarm_path <- "swarm"
  pigz_path <- "pigz"
  sep <- ","
  outdir <- "~/vtamR_demo_out_mfzr"
}else{
  setwd("C:/Users/emese/vtamR")
  cutadapt_path <- "cutadapt"
  vsearch_path <- "vsearch"
  blast_path <- "blastn"
  swarm_path <- "swarm"
  pigz_path <- "pigz"
  sep <- ","
  outdir <- "C:/Users/emese/vtamR_demo_out_mfzr"
}

## load
library(vtamR)
library(dplyr)
library(ggplot2)
library(rRDP)
library(rRDPData)

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
mock_ncbi_fasta <- system.file("extdata/demo/mock_ncbi.fasta", package = "vtamR")

#zfzr
#fastqinfo <-  system.file("extdata/demo/fastqinfo_marker2.csv", package = "vtamR")


#taxonomy_COInr <- "/home/meglecz/mkCOInr/COInr/COInr_for_vtam_2025_05_23_dbV5/COInr_for_vtam_taxonomy.tsv"


### Merge
merged_dir <- file.path(outdir, "1_merged")
fastainfo_df <- Merge(
  fastqinfo, 
  fastq_dir=fastq_dir, 
  vsearch_path=vsearch_path, 
  outdir=merged_dir,
  fastq_maxee=1,
  fastq_maxns=0,
  fastq_allowmergestagger=F,
  quiet=TRUE)

random_dir <- file.path(outdir, "2_RandomSeq")
fastainfo_random_seq <-RandomSeq(
  fastainfo_df, 
  n=40000,
  fasta_dir=merged_dir,
  outdir=random_dir, 
  randseed=0, 
  quiet=TRUE
  )


### demultiplex
demultiplexed_dir <- file.path(outdir, "3_demultiplexed")
sampleinfo_df <- SortReads(
  fastainfo_random_seq, 
  fasta_dir=random_dir, 
  outdir=demultiplexed_dir, 
  check_reverse=TRUE, 
  cutadapt_path=cutadapt_path, 
  vsearch_path=vsearch_path,
  quiet=TRUE
  )

###############
### dereplicate
###############
updated_asv_list <- file.path(outdir, "updated_asv_list.csv")
outfile <- file.path(outdir, "4_filter", "1_Input.csv")
read_count_df <- Dereplicate(
  sampleinfo_df, 
  outfile=outfile,
  dir=demultiplexed_dir
  )


### stat
stat_df <- data.frame(
  parameters=character(),
  asv_count=integer(),
  read_count=integer(),
  sample_count=integer(),
  sample_replicate_count=integer()
  )

stat_df <- GetStat(read_count_df, stat_df, stage="Input", params=NA)

#### swarm
by_sample <- TRUE
split_clusters <- TRUE
outfile <- file.path(outdir, "4_filter", "2_Swarm.csv")
read_count_df <- denoise_by_swarm(
  read_count_df,
  outfile=outfile,
  by_sample=by_sample,
  split_clusters=split_clusters,
  min_abundance_ratio = 0.2,
  min_read_count = 10,
  swarm_path=swarm_path,
  quiet=TRUE)

param <- paste("by_sample:", by_sample, "split_clusters:", split_clusters, sep=" ")
stat_df <- GetStat(read_count_df, stat_df, stage="Swarm", params=param)

#### LFNglobalReadCount
global_read_count_cutoff = 2
outfile <- file.path(outdir, "4_filter", "3_LFNglobalReadCount.csv")
read_count_df <- LFNglobalReadCount(
  read_count_df, 
  outfile=outfile,
  cutoff=global_read_count_cutoff)

param <- paste("global_read_count_cutoff:", global_read_count_cutoff, sep=" ")
stat_df <- GetStat(read_count_df, stat_df, stage="LFNglobalReadCount", params=NA)

#### FilterIndel
outfile <- file.path(outdir, "4_filter", "4_FilterIndel.csv")
read_count_df <- FilterIndel(
  read_count_df,
  outfile=outfile
  )

stat_df <- GetStat(read_count_df, stat_df, stage="FilterIndel", params=NA)

### FilterCodonStop
outfile <- file.path(outdir, "4_filter", "5_FilterCodonStop.csv")
genetic_code = 5
read_count_df <- FilterCodonStop(
  read_count_df, 
  genetic_code=genetic_code,
  outfile=outfile
  )

param <- paste("genetic_code:", genetic_code, sep=" ")
stat_df <- GetStat(read_count_df, stat_df, stage="FilterCodonStop", params=param)

### FilterExternalContaminant
outfile <- file.path(outdir, "4_filter", "6_FilterExternalContaminant.csv")
conta_file <- file.path(outdir, "4_filter", "Potential_contamination.csv")
read_count_df <- FilterExternalContaminant(
  read_count_df, 
  sampleinfo=sampleinfo_df,
  conta_file=conta_file,
  outfile=outfile
  )

stat_df <- GetStat(read_count_df, stat_df, stage="FilterExternalContaminant", params=NA)

### FilterChimera
outfile <- file.path(outdir, "4_filter", "7_FilterChimera.csv")
read_count_df <- FilterChimera(
  read_count_df, 
  vsearch_path=vsearch_path, 
  outfile=outfile
  )

stat_df <- GetStat(read_count_df, stat_df, stage="FilterChimera", params=NA)

### FilterRenkonen
outfile <- file.path(outdir, "4_filter", "8_FilterRenkonen.csv")
renkonen_distance_quantile <- 0.9
read_count_df <- FilterRenkonen(
  read_count_df, 
  renkonen_distance_quantile = renkonen_distance_quantile,
  outfile = outfile
  )

param <- paste("renkonen_distance_quantile:", renkonen_distance_quantile, sep=" ")
stat_df <- GetStat(read_count_df, stat_df, stage="FilterRenkonen", params=param)

### Make Mock_composition

outdir_mock <- file.path(outdir, "mock_composition")

mock_template <- MakeMockCompositionLTG(
  read_count=read_count_df,
  fas=mock_ncbi_fasta,
  taxonomy=taxonomy,
  sampleinfo = sampleinfo_df,
  outdir= outdir_mock,
  blast_path=blast_path
  )

mock_composition <- file.path(outdir_mock, "mock_composition_template_to_check.csv")

### LFNsampleReplicate

#### Suggest 
outfile = file.path(outdir, "4_filter", "optimize", "OptimizeLFNsampleReplicate.csv")
  
OptimizeLFNsampleReplicate_df <- OptimizeLFNsampleReplicate(
  read_count=read_count_df,
  mock_composition=mock_composition,
  outfile=outfile
  )

head(OptimizeLFNsampleReplicate_df)

#### LFNsampleReplicate
    
lfn_sample_replicate_cutoff <- 0.001
outfile <- file.path(outdir, "4_filter", "9_LFNsampleReplicate.csv")
  
read_count_df <- LFNsampleReplicate(
  read_count_df, 
  cutoff=lfn_sample_replicate_cutoff, 
  outfile=outfile
  )

param <- paste("cutoff:", lfn_sample_replicate_cutoff, sep=" ")
stat_df <- GetStat(read_count_df, stat_df, stage="LFNsampleReplicate", params=param)
  

  
### FilterMinReplicate 1

min_replicate_number <- 2
outfile <- file.path(outdir, "4_filter", "11_FilterMinReplicate.csv")
  
read_count_df <- FilterMinReplicate(
  read_count_df, 
  cutoff=min_replicate_number, 
  outfile=outfile
  )

param <- paste("cutoff:", min_replicate_number, sep=" ")
stat_df <- GetStat(read_count_df, stat_df, stage="FilterMinReplicate", params=param)
  

### LFNvariant and LFNreadCount

#### Suggest
outdir_optimize = file.path(outdir, "4_filter", "optimize")
              
OptimizeLFNreadCountLFNvariant_df <- OptimizeLFNreadCountLFNvariant(
  read_count_df,
  mock_composition = mock_composition,
  sampleinfo = sampleinfo_df,
  habitat_proportion = 0.5,
  outdir= outdir_optimize, 
  min_replicate_number=2
)

head(OptimizeLFNreadCountLFNvariant_df)

### LFNvariant

outfile <- file.path(outdir, "4_filter", "12_LFNvariant.csv")
lnf_variant_cutoff = 0.001
              
read_count_df_lnf_variant <- LFNvariant(
  read_count_df, 
  cutoff=lnf_variant_cutoff, 
  outfile=outfile
  )

param <- paste("cutoff:", lnf_variant_cutoff, sep=" ")         
stat_df <- GetStat(read_count_df_lnf_variant, stat_df, stage="LFNvariant", params=param)

#### LFNreadCount

lfn_read_count_cutoff <- 25
outfile <- file.path(outdir, "4_filter", "13_LFNreadCount.csv")
              
read_count_df_lfn_read_count <- LFNreadCount(
  read_count_df, 
  cutoff=lfn_read_count_cutoff, 
  outfile=outfile
  )

param <- paste("cutoff:", lfn_read_count_cutoff, sep=" ")        
stat_df <- GetStat(read_count_df_lfn_read_count, stat_df, stage="LFNreadCount", params=param)

####  Combine results

outfile <- file.path(outdir, "4_filter", "14_poolLFN.csv")
read_count_df <- PoolFilters(
  read_count_df_lfn_read_count, read_count_df_lnf_variant,    
  outfile=outfile
  )

stat_df <- GetStat(read_count_df, stat_df, stage="PoolFilters", params=NA)
              
# delete temporary data frames
rm(read_count_df_lfn_read_count)
rm(read_count_df_lnf_variant)
              

### FilterMinReplicate 2

min_replicate_number <- 2
outfile <- file.path(outdir, "4_filter", "15_FilterMinReplicate.csv")

read_count_df <- FilterMinReplicate(
  read_count_df, 
  cutoff=min_replicate_number, 
  outfile=outfile
  )

param <- paste("cutoff:", min_replicate_number, sep=" ")  
stat_df <- GetStat(read_count_df, stat_df, stage="FilterMinReplicate", params=param)

### Pool replicates

outfile <- file.path(outdir, "4_filter", "16_pool_replicates.csv")

read_count_sample_df <- PoolReplicates(
  read_count_df,
  method = "mean",   # Aggregate replicates using the mean read count
  digits = 0,        # Round aggregated read counts to integers
  outfile = outfile
)

### Get performance metrics

missing_occurrences <- file.path(outdir, "4_filter", "Missing_occurrences.csv")
performance_metrics <- file.path(outdir, "4_filter", "Performance_metrics.csv")
known_occurrences <- file.path(outdir, "4_filter", "Known_occurrences.csv")

results <- MakeKnownOccurrences(
  read_count_sample_df, 
  sampleinfo=sampleinfo_df, 
  mock_composition=mock_composition, 
  known_occurrences=known_occurrences, 
  missing_occurrences=missing_occurrences,
  performance_metrics=performance_metrics
  )

# give explicit names to the 3 output data frames
known_occurrences_df <- results[[1]]
missing_occurrences_df <- results[[2]]
performance_metrics_df <- results[[3]]

print(performance_metrics_df)

## Taxonomic Assignment

outfile <- file.path(outdir, "4_filter", "TaxAssignLTG.csv")
asv_tax <- TaxAssignLTG(
  asv=read_count_sample_df, 
  taxonomy=taxonomy, 
  blast_db=blast_db, 
  blast_path=blast_path, # can be omitted if BLAST is in the PATH
  outfile=outfile
  )

