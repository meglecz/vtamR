

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

### set up windows
cutadapt_path <- "cutadapt"
vsearch_path <- "vsearch"
blast_path <- "blastn"
swarm_path <- "swarm"
pigz_path <- "pigz"
sep <- ","
outdir <- "C:/Users/emese/vtamR_demo_out_marker2"

fastq_dir <- system.file("extdata/demo/fastq", package = "vtamR")
fastqinfo <-  system.file("extdata/demo/fastqinfo_marker2.csv", package = "vtamR")
#mock_composition <-  system.file("extdata/demo/mock_composition.csv", package = "vtamR")
#asv_list <-  system.file("extdata/demo/asv_list.csv", package = "vtamR")
asv_list <- "C:/Users/emese/vtamR_demo_out/updated_asv_list.csv"
taxonomy <- system.file("extdata/db_test/taxonomy_reduced.tsv", package = "vtamR")
blast_db <- system.file("extdata/db_test", package = "vtamR")
blast_db <- file.path(blast_db, "COInr_reduced")


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
updated_asv_list <- file.path(outdir, "updated_asv_list_2markers.tsv")
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

write.csv(read_count_df, file=file.path(outdir, "zfzr.csv"), row.names = FALSE)

zfzr <- read_count_df %>%
  mutate(marker = "zfzr")

mfzr <- read.csv2("C:/Users/emese/vtamR_demo_out/mfzr.csv") %>%
  mutate(marker="mfzr")

df <- pool_markers(mfzr, zfzr, 
             method ="mean", 
             outfile = file.path(outdir, "pooled", "2markers.csv"),
             asv_with_centroids  = file.path(outdir, "pooled", "2markers_asv_with_centroids .csv")
             )



