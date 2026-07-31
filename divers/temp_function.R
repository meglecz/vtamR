setwd("~/vtamR")

library(vtamR)
library(dplyr)

outdir <- "vtamR_demo_out/zfzr_plate1"
#set_log_file(file.path(outdir, "vtamR_log.csv"))
#set_log_file(NULL)



fastq_dir <- system.file("extdata/demo/fastq", package = "vtamR")
fastqinfo <- system.file("extdata/demo/fastqinfo_zfzr_plate1.csv", package = "vtamR")
mock_ncbi_fasta <- system.file("extdata/demo/mock_ncbi.fasta", package = "vtamR")
asv_list <- system.file("extdata/demo/ASV_list_with_IDs.csv", package = "vtamR") # all ASV of MFZR
taxonomy <- system.file("extdata/db_test/taxonomy_reduced.tsv", package = "vtamR")
blast_db <- system.file("extdata/db_test", package = "vtamR")
blast_db <- file.path(blast_db, "COInr_reduced")


merged_dir <- file.path(outdir, "1_merged")

fastainfo_df <- merge_fastq_pairs(
  fastqinfo=fastqinfo,
  fastq_dir=fastq_dir,
  outdir=merged_dir,
  fastq_maxee=1,
  fastq_maxns=0,
  fastq_allowmergestagger=T
)


demultiplexed_dir <- file.path(outdir, "3_demultiplexed")

sampleinfo_df <- demultiplex_and_trim(
  fastainfo=fastainfo_df,
  fasta_dir=merged_dir,
  outdir=demultiplexed_dir,
  check_reverse=TRUE,
  cutadapt_minimum_length = 150,
  cutadapt_maximum_length = 165
)



# get function name, all arguments and stat time
log <- collect_log()

# add end_time and runtime, print
write_log(log)

@param log_file Character string specifying the path to the CSV log file.
#'   If `NULL`, no log file is written.
,
log_file = NULL

log <- collect_log()

write_log(log, file=log_file)


compute_renkonen_distances(
  read_count_df ==> correcte to read_count
  
  filter_replicate
  Error in data.frame(function_name = rep(fun_name, length(args)), argument_name = names(args), : arguments imply differing number of rows: 4, 5
                      
suggest_sample_cutoff
rror in UseMethod("mutate") : pas de méthode pour 'mutate' applicable pour un objet de classe "function"

pool_filters
plot_cluster_classification
