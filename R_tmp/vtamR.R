install.packages("devtools")
install.packages("roxygen2")
install.packages("seqinr")
install.packages("tidyr")
#install.packages("Rtools")


library("devtools")
library("roxygen2")
library("seqinr") # splitseq for FilterCodonStop
library("dplyr")
library("tidyr") # gather for read_asv_table; pivot_wider in WriteAsVtable and stat_sample !!sym
#library("utils") # to handle zipped files
library("ggplot2") 
#library("Rtools")
#install_github("meglecz/vtamR")
#devtools::install_github("CalCOFI/rCRUX", build_vignettes = TRUE)
#sudo apt-get install r-cran-devtools
#devtools::install_github("meglecz/vtamR", build_vignettes = TRUE)

vignette("rd")
vignette("rd-other")
vignette("rd-formatting")
vignette("reuse")
vignette("namespace")

#library("Biostrings")

computer <- "Windows" # Bombyx/Endoume/Windows
if(computer == "Bombyx"){
  vtam_dir <- "~/vtamR"
  cutadapt_path="/home/meglecz/miniconda3/envs/vtam_2/bin/"
  vsearch_path = ""
  blast_path="~/ncbi-blast-2.11.0+/bin/" # bombyx
  swarm_path <- ""
  db_path="~/mkLTG/COInr_for_vtam_2022_05_06_dbV5/"
      fastq_dir <- "vtamR_test/data/"
      fastqinfo <- "vtamR_test/data/fastqinfo_zfzr.csv"
      outdir <- "vtamR_test/out_zfzr/"
      mock_composition <- "vtamR_test/data/mock_composition_zfzr.csv"
      asv_list <- "vtamR_test/data/asv_list.csv"
  #      fastq_dir <- "/home/meglecz/vtamR_large_files/fastq/"
      #      fastqinfo <- "/home/meglecz/vtamR_large_files/user_input/fastqinfo_mfzr.csv"
      #     outdir <- "/home/meglecz/vtamR_large_files/out/"
     #     mock_composition <- "/home/meglecz/vtamR_large_files/user_input/mock_composition_mfzr.csv"
     #     asv_list <- "/home/meglecz/vtamR_large_files/user_input/asv_list.csv"

  num_threads=8
  compress = T
} else if (computer == "Endoume"){
  vtam_dir <- "~/vtamR"
  cutadapt_path="/home/emese/miniconda3/bin/"
  vsearch_path = "/home/emese/miniconda3/bin/"
  blast_path= "" # deactivate conda
  swarm_path <- ""
  db_path= "~/mkCOInr/COInr/COInr_for_vtam_2023_05_03_dbV5/"
  #    fastq_dir <- "vtamR_test/data/"
  #     fastqinfo <- "vtamR_test/data/fastqinfo_mfzr.csv"
  #     outdir <- "vtamR_test/out_mfzr/"
  #     mock_composition <- "vtamR_test/data/mock_composition_mfzr.csv"
  #     asv_list <- "vtamR_test/data/asv_list_zfzr.csv"
      fastq_dir <- "~/vtamR_large_data"
      fastqinfo <- "~/vtamR_large_data/metadata/fastqinfo_Sea18_IIICBR_vtamR.csv"
      outdir <- "/home/emese/vtamR_large_data/out/"
     mock_composition <- "~/vtamR_large_data/metadata/mock_composition_Sea18_IIICBR_vtamR.csv"
      asv_list <- "~/vtamR_large_data/metadata/asv_list.csv"
  num_threads=8
  compress = T
}else if (computer == "Windows"){
  vtam_dir <- "C:/Users/emese/vtamR/"
  cutadapt_path="C:/Users/Public/"
  vsearch_path = "C:/Users/Public/vsearch-2.23.0-win-x86_64/bin/"
  blast_path="C:/Users/Public/blast-2.14.1+/bin/"
  swarm_path <- "C:/Users/Public/swarm-3.1.4-win-x86_64/bin/"
  db_path="C:/Users/Public/COInr_for_vtam_2023_05_03_dbV5/"
#  fastq_dir <- "C:/Users/emese/vtamR_private/fastq/"
  fastq_dir <- "vtamR_test/data/"
  fastqinfo <- "vtamR_test/data/fastqinfo_zfzr.csv"
  outdir <- "vtamR_test/out_zfzr/"
  mock_composition <- "vtamR_test/data/mock_composition_zfzr.csv"
#  asv_list <- "vtamR_test/data/asv_list.csv"
  asv_list <- "vtamR_test/out_mfzr/ASV_list_with_IDs.csv"
  num_threads=4
  compress = F
}
sep=","
setwd(vtam_dir)

taxonomy=paste(db_path, "COInr_for_vtam_taxonomy.tsv", sep="")
blast_db=paste(db_path, "COInr_for_vtam", sep="")

ltg_params = "vtamR_test/data/ltg_params.csv"
ltg_params_df = data.frame( pid=c(100,97,95,90,85,80),
                            pcov=c(70,70,70,70,70,70),
                            phit=c(70,70,70,70,70,70),
                            taxn=c(1,1,2,3,4,4),
                            seqn=c(1,1,2,3,4,4),
                            refres=c("species","species","species","genus","family","family"),
                            ltgres=c("species","species","species","species", "genus","genus")
)

ltg_params_df = data.frame( pid=c(100,97,95,90,85,80),
                            pcov=c(70,70,70,70,70,70),
                            phit=c(70,70,70,70,70,70),
                            taxn=c(1,1,2,3,4,4),
                            seqn=c(1,1,2,3,4,4),
                            refres=c(8,8,8,7,6,6),
                            ltgres=c(8,8,8,8,7,7)
)

# load local packages
load_all(".")
roxygenise() 
usethis::use_roxygen_md()

####
# Make TrimPrimer
###
fastainfo <- "/home/meglecz/vtamR/vtamR_test/data/fastainfo_notag.csv"
fasta_dir <- "vtamR_test/out/test_sortreads/no_tag/"
outdir <- "vtamR_test/out/test_sortreads/primer_trimmed/"
outdir <- check_dir(dir=outdir)
check_reverse <- T
primer_to_end <-F
cutadapt_error_rate <- 0.1 # -e in cutadapt
cutadapt_minimum_length <- 50 # -m in cutadapt
cutadapt_maximum_length <- 500 # -M in cutadapt
compress <- T

fastainfo_df <- TrimPrimer(fastainfo, fasta_dir=fasta_dir, outdir=outdir, compress=compress, cutadapt_path=cutadapt_path, vsearch_path=vsearch_path, check_reverse=check_reverse, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length)


###
# Check coherence (often) user-made input files
###
# file_type [fastqinfo/fastainfo/sortedinfo/mock_composition/known_occurrences/read_count/read_count_sample/asv_list]
CheckFileinfo(file=fastqinfo, dir=fastq_dir, file_type="fastqinfo", sep=sep)
CheckFileinfo(file=mock_composition, file_type="mock_composition", sep=sep)
CheckFileinfo(file=asv_list, file_type="asv_list", sep=sep)

###
# Test major functions
###
Test_MergeSortReads(test_dir="vtamR_test/", vsearch_path=vsearch_path, cutadapt_path=cutadapt_path, delete_tmp=T)
Test_Filters(test_dir="vtamR_test/", swarm_path=swarm_path, vsearch_path=vsearch_path, sep=sep, delete_tmp=T)
Test_MakeKnownOccurrences(test_dir="vtamR_test/", sep=sep, delete_tmp=T)
Test_Optimize(test_dir="vtamR_test/", vsearch_path=vsearch_path, delete_tmp=T)
Test_TaxAssign(test_dir="vtamR_test/", sep=sep, blast_path=blast_path, num_threads=num_threads)


####

# define input filenames
#fasta_dir <- "local/mfzr/sorted/"
#sortedinfo <- "local/user_input/fileinfo_mfzr_eu.csv"


#fasta_dir <- "/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/"
#sortedinfo <-"/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/fileinfo_vtamr.csv"
#mock_composition <- "/home/meglecz/vtamR/local/user_input/mock_composition_mfzr_prerun.csv"
#sep="\t"

#fasta_dir <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/"
#sortedinfo <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/fileinfo_vtamr.csv"


# create the output directory and check the the slash at the end
outdir <- check_dir(dir=outdir)
fastq_dir <- check_dir(dir=fastq_dir)
# Measure runtime using system.time()
start_time <- Sys.time()  # Record the start time
# define stat data frame that will be completed with counts after each step
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())


###
### Merge
###
fastq_ascii <- 33
fastq_maxdiffs <- 10
fastq_maxee <- 1
fastq_minlen <- 50
fastq_maxlen <- 500
fastq_minmergelen <- 50
fastq_maxmergelen <-500
fastq_maxns <- 0
fastq_truncqual <- 10
fastq_minovlen <- 50
fastq_allowmergestagger <- F
merged_dir <- paste(outdir, "merged/", sep="")
compress = F
# read fastqinfo
fastainfo_df <- Merge(fastqinfo, fastq_dir=fastq_dir, vsearch_path=vsearch_path, outdir=merged_dir, fastq_ascii=fastq_ascii, fastq_maxdiffs=fastq_maxdiffs, fastq_maxee=fastq_maxee, fastq_minlen=fastq_minlen, fastq_maxlen=fastq_maxlen, fastq_minmergelen=fastq_minmergelen, fastq_maxmergelen=fastq_maxmergelen, fastq_maxns=fastq_maxns, fastq_truncqual=fastq_truncqual, fastq_minovlen=fastq_minovlen, fastq_allowmergestagger=fastq_allowmergestagger, sep=sep, compress=compress)

#test_file <- '/home/emese/vtamR_large_data/out/merged/Sea18_COI_R1_S8_R1_001.fasta' 
#cmd <- paste("gzip -k", test_file, sep=" ")
#system(cmd)
###
### RandomSeq
# RandomSeq is about 5 times quicker than RandomSeq_windows, and it is more memory efficient.
# However, RandomSeq does not work on Windows
###
randomseq_dir = paste(outdir, "random_seq/", sep="")
fastainfo <- paste(merged_dir, "fastainfo.csv", sep="")
#fastainfo_df <- read.csv(file=fastainfo, header=T, sep=sep)
compress = T
fastainfo_df <- RandomSeq(fastainfo, fasta_dir=merged_dir, outdir=randomseq_dir, vsearch_path=vsearch_path, n=10000, randseed=0, compress=compress, sep=sep)
#fastainfo_df <- RandomSeqWindows(fastainfo_df, fasta_dir=merged_dir, outdir=randomseq_dir, n=1000, randseed=0, compress=compress, sep=sep)


###
### SortReads
###
sorted_dir <- paste(outdir, "sorted/", sep="")
check_reverse <- T
tag_to_end <- F
primer_to_end <-F
cutadapt_error_rate <- 0.1 # -e in cutadapt
cutadapt_minimum_length <- 50 # -m in cutadapt
cutadapt_maximum_length <- 500 # -M in cutadapt
compress <- T
fastainfo <- paste(randomseq_dir, "fastainfo.csv", sep="")
sortedinfo_df <- SortReads(fastainfo_df, fasta_dir=randomseq_dir, outdir=sorted_dir, cutadapt_path=cutadapt_path, vsearch_path=vsearch_path, check_reverse=check_reverse, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)


###############
randomseq_dir = paste(outdir, "random_seq_after_sortreads/", sep="")
fastainfo_df <- RandomSeq(sortedinfo_df, fasta_dir=sorted_dir, outdir=randomseq_dir, vsearch_path=vsearch_path, n=100, randseed=0, compress=compress, sep=sep)
#fastainfo_df <- RandomSeqWindows(fastainfo_df, fasta_dir=merged_dir, outdir=randomseq_dir, n=1000, randseed=0, compress=compress, sep=sep)
###############


###
### Read input fasta files, dereplicate reads to ASV, and count the number of reads of each ASV in each sample-replicate, add a unique id for ASVs, can take into account ASVs from earlier analyses
###
outfile <- paste(outdir, "1_before_filter.csv", sep="")
sortedinfo_df <- read.csv(paste(sorted_dir, "sortedinfo.csv", sep =""), sep=sep)
updated_asv_list <- sub("\\.", "_updated_2024_02_29.", asv_list) # add date to the name of the input asv_list to get a file name for the updated_file
read_count_df <- Dereplicate(sortedinfo_df, dir=sorted_dir, outfile=outfile, sep=sep, asv_list=asv_list, updated_asv_list=updated_asv_list)
# make stat counts
stat_df <- GetStat(read_count_df, stat_df, stage="Input", params=NA)

###
# Run swarm
###
read_count_df <- read.csv(paste(outdir, "1_before_filter.csv", sep=""), sep=sep)

swarm_d <- 1
fastidious <- TRUE
by_sample <- FALSE
quiet <- F
outfile <- paste(outdir, "2_Swarm_by_sample.csv", sep="")
read_count_df <- Swarm(read_count_df, outfile=outfile, swarm_path=swarm_path, num_threads=num_threads, swarm_d=swarm_d, fastidious=fastidious, sep=sep, by_sample=by_sample, quiet=quiet)
params <- paste(swarm_d, fastidious, by_sample, sep=";")
stat_df <- GetStat(read_count_df, stat_df, stage="swarm_by_sample", params=params)

###
### LFNglobalReadCount
###
# Eliminate variants with less than global_read_count_cutoff reads in the dataset
global_read_count_cutoff = 2
outfile <- paste(outdir, "3_LFNglobalReadCount.csv", sep="")
read_count_df <- LFNglobalReadCount(read_count_df, cutoff=global_read_count_cutoff, outfile=outfile)
stat_df <- GetStat(read_count_df, stat_df, stage="LFNglobalReadCount", params=global_read_count_cutoff)

###
### LFN_filters
###
# LFNreadCount
lfn_read_count_cutoff <- 10
outfile <- paste(outdir, "4_LFNreadCount.csv", sep="")
read_count_df_lfn_read_count <- LFNreadCount(read_count_df, cutoff=lfn_read_count_cutoff, outfile=outfile, sep=sep)
stat_df <- GetStat(read_count_df_lfn_read_count, stat_df, stage="LFNreadCount", params=lfn_read_count_cutoff)


# LFNsampleReplicate (by column)
lfn_sample_replicate_cutoff <- 0.001
outfile <- paste(outdir, "5_LFNsampleReplicate.csv", sep="")
read_count_df_lnf_sample_replicate <- LFNsampleReplicate(read_count_df, cutoff=lfn_sample_replicate_cutoff, outfile=outfile, sep=sep)
stat_df <- GetStat(read_count_df_lnf_sample_replicate, stat_df, stage="LFNsampleReplicate", params=lfn_sample_replicate_cutoff)


# LFNvariant (by line)
lnf_variant_cutoff = 0.001
by_replicate = TRUE
outfile <- paste(outdir, "6_LFNvariant_replicate.csv", sep="")
min_read_count_prop = 0.7
read_count_df_lnf_variant <- LFNvariant(read_count_df, cutoff=lnf_variant_cutoff, by_replicate, outfile=outfile, sep=sep, min_read_count_prop=min_read_count_prop)
param_values <- paste(lnf_variant_cutoff, by_replicate, sep=";")
stat_df <- GetStat(read_count_df_lnf_variant, stat_df, stage="LFNvariant", params=param_values)


# pool the results of the different filterLFN to one data frame; keep only occurrences that passed all filters
outfile <- paste(outdir, "7_poolLFN.csv", sep="")
read_count_df <- PoolFilters(read_count_df_lfn_read_count, read_count_df_lnf_variant, read_count_df_lnf_sample_replicate, outfile=outfile, sep=sep)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterLFN")
# delete temporary data frames
read_count_df_lfn_read_count <- NULL
read_count_df_lnf_variant <- NULL
read_count_df_lnf_sample_replicate <- NULL


###
### keep repeatable occurrences
###
min_replicate_number <- 2
outfile <- paste(outdir, "8_FilterMinReplicate.csv", sep="")
read_count_df <- FilterMinReplicate(read_count_df, min_replicate_number, outfile=outfile, sep=sep)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterMinReplicate", params=min_replicate_number)

###
### FilerIndel
###
outfile <- paste(outdir, "9_FilterIndel.csv", sep="")
read_count_df <- FilterIndel(read_count_df, outfile=outfile, sep=sep)
stat_df <- GetStat(read_count_df, stat_df, stage="FilterIndel")

###
### FilterCodonStop
###
outfile <- paste(outdir, "10_FilterCodonStop.csv", sep="")
genetic_code = 5
read_count_df <- FilterCodonStop(read_count_df, outfile=outfile, genetic_code=genetic_code, sep=sep)
stat_df <- GetStat(read_count_df, stat_df, stage="FilerCodonStop", params=genetic_code)


###
### FilerPCRerror
###
pcr_error_var_prop <- 0.1
max_mismatch <- 2
by_sample <- T
sample_prop <- 0.8
outfile <- paste(outdir, "11_FilterPCRerror.csv", sep="")
read_count_df <- FilterPCRerror(read_count_df, outfile=outfile, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, sep=sep)
params <- paste(pcr_error_var_prop, max_mismatch, by_sample, sample_prop, sep=";")
stat_df <- GetStat(read_count_df, stat_df, stage="FilerPCRerror", params=params)

###
### FilterChimera
###
abskew=2
by_sample = T
sample_prop = 0.8
outfile <- paste(outdir, "12_FilterChimera.csv", sep="")
read_count_df <- FilterChimera(read_count_df, outfile=outfile, vsearch_path=vsearch_path, by_sample=by_sample, sample_prop=sample_prop, abskew=abskew, sep=sep)
params <- paste(abskew, by_sample, sample_prop, sep=";")
stat_df <- GetStat(read_count_df, stat_df, stage="FilterChimera", params=params)



###
# Avoid rerunning the longest step
###
read_count_df <- read.csv(paste(outdir, "12_FilterChimera.csv", sep =""), sep=sep)
sorted_dir <- paste(outdir, "sorted/", sep="")
sortedinfo <- paste(sorted_dir, "sortedinfo.csv", sep="")
stat_df <- GetStat(read_count_df, stat_df, stage="FilterChimera", params=params)

###
### FilterRenkonen
###
# Renkonen index:
# PS = summ(min(p1i, p2i))
# p1i = number of reads for variant i in replicate 1 / number of reads in replicate 1
renkonen_distance_quantile = 0.9
outfile <- paste(outdir, "13_FilterRenkonen.csv", sep="")
read_count_df <- FilterRenkonen(read_count_df, outfile=outfile, renkonen_distance_quantile=renkonen_distance_quantile, sep=sep)
stat_df <- GetStat(read_count_df, stat_df, stage="FilerRenkonen", params=renkonen_distance_quantile)

###
### PoolReplicates
###
digits = 0
outfile <- paste(outdir, "14_PoolReplicates.csv", sep="")
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, outfile=outfile, sep=sep)
stat_df <- GetStat(read_count_samples_df, stat_df, stage="PoolReplicates")

###
### TaxAssign
###
outfile <- paste(outdir, "TaxAssign.csv", sep="")
asv_tax <- TaxAssign(asv=read_count_samples_df, ltg_params=ltg_params, taxonomy=taxonomy, blast_db=blast_db, blast_path=blast_path, outfile=outfile, num_threads=num_threads)


###
### print output files
###
# write sequence and variant counts after each step
write.csv(stat_df, file = paste(outdir, "stat_steps.csv", sep=""))
# long format, each line corresponds to an occurrence (); if csv files are written at each step, this is not usefull
#write.csv(read_count_samples_df, file = paste(outdir, "16_Final_asvtable_long.csv", sep=""), row.names=F)
# wide format (ASV table), samples are in columns, ASVs in lines
outfile=paste(outdir, "Final_asvtable.csv", sep="")
sortedinfo <- paste(sorted_dir, "sortedinfo.csv", sep ="")
WriteASVtable(read_count_samples_df, outfile=outfile, sortedinfo=sortedinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)
# write ASV table completed by taxonomic assignments
outfile=paste(outdir, "Final_asvtable_with_taxassign.csv", sep="")
WriteASVtable(read_count_samples_df, outfile=outfile, asv_tax=asv_tax, sortedinfo=sortedinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)



# start optimize from almost unfiltered data (after eliminating ASV with low global reads count)
file <- paste(outdir, "2_Swarm_by_sample.csv", sep="")
read_count_df <- read.csv(file, sep=sep)
dim(read_count_df)
###
### OptimizePCRerror
###
outfile <- paste(outdir, "OptimizePCRerror.csv", sep="")
OptimizePCRError_df <- OptimizePCRerror(read_count=read_count_df, mock_composition=mock_composition, sep=sep, outfile=outfile, max_mismatch=1, min_read_count=10)

###
### OptimizeLFNsampleReplicate
###
outfile = paste(outdir, "OptimizeLFNsampleReplicate.csv", sep="")
OptimizeLFNsampleReplicate_df <- OptimizeLFNsampleReplicate(read_count=read_count_df, mock_composition=mock_composition, sep=sep, outfile=outfile)

###
### Make known occurrences
###
file <- paste(outdir, "14_PoolReplicates.csv", sep="")
read_count_samples_df <- read.csv(file, sep=sep)
sortedinfo <- paste(sorted_dir, "sortedinfo.csv", sep ="")
known_occurrences <- paste(outdir, "known_occurrences.csv", sep= "")
missing_occurrences <- paste(outdir, "missing_occurrences.csv", sep= "")
performance_metrics <- paste(outdir, "performance_metrics.csv", sep= "")
habitat_proportion= 0.5 # for each asv, if the proportion of reads in a habitat is below this cutoff, is is considered as an artifact in all samples of the habitat
results <- MakeKnownOccurrences(read_count_samples = read_count_samples_df, sortedinfo=sortedinfo, mock_composition=mock_composition, sep=sep, known_occurrences=known_occurrences, missing_occurrences=missing_occurrences, performance_metrics=performance_metrics, habitat_proportion=habitat_proportion)
known_occurrences_df <- results[[1]]
missing_occurrences <- results[[2]]
performance_metrics <- results[[3]]

###
### OptimizeLFNReaCountAndLFNvariant
###

min_lfn_read_count_cutoff=10
max_lfn_read_count_cutoff=100
increment_lfn_read_count_cutoff=10
min_lnf_variant_cutoff=0.001
max_lnf_variant_cutoff=0.01
increment_lnf_variant_cutoff=0.01
by_replicate=FALSE
max_mismatch=1
by_sample=T
sample_prop=0.8
min_replicate_number=2
outfile = paste(outdir, "OptimizeLFNreadCountLFNvariant.csv", sep="")
OptimizeLFNreadCountLFNvariant_df <- OptimizeLFNreadCountLFNvariant(read_count=read_count_df, known_occurrences=known_occurrences_df, sep=sep, outfile= outfile, min_lfn_read_count_cutoff=lfn_read_count_cutoff, max_lfn_read_count_cutoff=max_lfn_read_count_cutoff, increment_lfn_read_count_cutoff=increment_lfn_read_count_cutoff, min_lnf_variant_cutoff=min_lnf_variant_cutoff, max_lnf_variant_cutoff=max_lnf_variant_cutoff, increment_lnf_variant_cutoff=increment_lnf_variant_cutoff, by_replicate=by_replicate, min_replicate_number=min_replicate_number)


##################
# Helper functions
##################

###
# check the a particular value of a feature (asv, asv_id, sample, replicate) in all intermediate output files
###
tmp <- HistoryBy(dir=outdir, feature="asv_id", value="915", sep=sep)

###
# group lines by grouped_by (asv, asv_id, sample, replicate) variable and summarize a feature (asv, asv_id, sample, replicate, read_count) for each group
###
outfile <- paste(outdir, "read_count_by_sample.csv", sep="")
read_count_by_sample <- SummarizeBy(dir=outdir, sep=sep, outfile=outfile, feature="read_count", grouped_by="sample")



###
# The safest way to homogenize ASVs and asv_ids among different data sets is updating the asv_list systematically at the Dereplicate step.
# This will keep all ASVs ever seen in any of the data sets you have analyzed.
# However, if you have many large data sets, the file will grow quickly that can cause memory issues. 
# A quite safe solution is to update the asv_list after swarm, since swarm has eliminated the majority of the rare ASVs that will be unlikely to be seen among validated ASVs in further runs.
###

updated_asv_list <- sub("\\.", "_updated_2024_02_19_after_swarm.", asv_list) # add date to the name of the input asv_list to get a file name for the updated_file
UpdateASVlist(read_count_df, asv_list=asv_list, outfile=updated_asv_list)

###
# Pool different data sets
###
files <- data.frame(file=c("vtamR_test/out_mfzr/14_PoolReplicates.csv", "vtamR_test/out_zfzr/14_PoolReplicates.csv"),
                    marker=c("MFZR", "ZFZR"))
outfile <- paste(outdir, "Pooled_datasets.csv", sep="") # pooled data sets; ASVs are grouped in identical in their overlapping region
asv_with_centroids <- paste(outdir, "Pooled_datasets_asv_with_centroids.csv", sep="") # original ASVs + centroids for each of them
read_count_pool <- PoolDatasets(files, outfile=outfile, asv_with_centroids=asv_with_centroids, sep=sep, mean_over_markers=T, vsearch_path=vsearch_path)

###
# count reads in fasta of fastq
###

dir <- "vtamR_test/data"
df <- CountReadsDir(dir, pattern="_fw.fastq.gz", file_type="fastq", outfile="tmp.csv", sep=",")


###
# Graphs
###

#https://r-graph-gallery.com
#https://www.data-to-viz.com/graph/barplot.html


file <- paste(outdir, "12_FilterChimera.csv", sep="")
read_count_df <- read.csv(file, sep=sep)
###
# Barplot_ReadCountBySample
###
sortedinfo <- paste(sorted_dir, "sortedinfo.csv", sep ="")
p <- Barplot_ReadCountBySample(read_count_df=read_count_df, sample_replicate=F, sample_types=sortedinfo, sep=sep, x_axis_label_size=8 )
print(p)
p <- Barplot_ReadCountBySample(read_count_df=read_count_df, sample_replicate=T, sep=sep, x_axis_label_size=4 )
print(p)

###
# Histogram_ReadCountByVariant
###
p <- Histogram_ReadCountByVariant(read_count_df, min_read_count=10, binwidth=1000)
print(p)

###
# graph_renkonen_distances
###
renkonen_within_df <- MakeRenkonenDistances(read_count_df, compare="within")
renkonen_all_df <- MakeRenkonenDistances(read_count_df, compare="all")
p <- Barplot_RenkonenDistance(renkonen_within_df, sample_types=sortedinfo, sep=sep, x_axis_label_size=6)
print(p)
p <- DensityPlot_RenkonenDistance(renkonen_within_df)
print(p)
  

  
start_time <- Sys.time() 
end_time <- Sys.time()  # Record the end time
runtime <- end_time - start_time  # Calculate the run time
print(runtime)



**Set general parameters Windows**

  
```{r set_path, eval=TRUE}
cutadapt_path <- "C:/Users/Public"
vsearch_path <- "C:/Users/Public/vsearch-2.23.0-win-x86_64/bin"
blast_path <- "C:/Users/Public/blast-2.14.1+/bin/"
swarm_path <- "C:/Users/Public/swarm-3.1.4-win-x86_64/bin"
vtam_dir <- "C:/Users/emese/vtamR"
setwd(vtam_dir)
num_threads <- 4
sep <- ","
```

**Set general parameters Bombyx**
  ```{r, eval=TRUE}
cutadapt_path="/home/meglecz/miniconda3/envs/vtam_2/bin/"
vsearch_path = ""
blast_path="~/ncbi-blast-2.11.0+/bin/"
swarm_path <- ""
num_threads <- 8
sep <- ","
```


**Set general parameters Endoume**
  ```{r, eval=TRUE}
cutadapt_path="/home/emese/miniconda3/bin/"
vsearch_path = ""
blast_path=""
swarm_path <- ""
num_threads <- 8
sep <- ","
```

* Adapt the path to third party programs according to your installation. If the program is in your path, it can be empty (e.g. `cutadapt_path=""`)
* `num_threads` is the number of CPUs for multithreaded programs
* `sep` is the separator used in csv files


**TaxAssign reference data base**
  
  You can download a ready to use, comprehensive COI database from [OSF](https://osf.io/vrfwz/) ([Meglécz, 2023](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756)). Set the path to the database and the accompanying taxonomy file.

```{r, eval=TRUE}
taxonomy <- "C:/Users/Public/COInr_for_vtam_2023_05_03_dbV5/COInr_for_vtam_taxonomy.tsv"
blast_db <- "C:/Users/Public/COInr_for_vtam_2023_05_03_dbV5/COInr_for_vtam"
```

**TaxAssign reference data base Bombyx**
  ```{r, eval=TRUE}
taxonomy <- "~/mkLTG/COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam_taxonomy.tsv"
blast_db <- "~/mkLTG/COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam"
```


**TaxAssign reference data base Endoume**
  ```{r, eval=TRUE}
taxonomy <- "~/mkCOInr/COInr/COInr_for_vtam_2023_05_03_dbV5/COInr_for_vtam_taxonomy.tsv"
blast_db <- "~/mkCOInr/COInr/COInr_for_vtam_2023_05_03_dbV5/COInr_for_vtam"
```

* `taxonomy` CSV file file with taxonomic information
* `blast_db` BLAST database
* Details are in [Reference database for taxonomic assignments section](#reference-database-for-taxonomic-assignments)
  
