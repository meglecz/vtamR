install.packages("devtools")
install.packages("roxygen2")
install.packages("seqinr")
install.packages("tidyr")

library("devtools")
library("roxygen2")
library("seqinr") # splitseq for FilterCodonStop
library("dplyr")
library("tidyr") # gather for read_asv_table; pivot_wider in write_asvtable
library("utils") # to handle zipped files
#library("Biostrings")


computer <- "Bombyx" # Bombyx/Endoume/Windows
if(computer == "Bombyx"){
  vtam_dir <- "~/vtamR"
  cutadapt_path="/home/meglecz/miniconda3/envs/vtam_2/bin/"
  vsearch_path = ""
  blast_path="~/ncbi-blast-2.11.0+/bin/" # bombyx
  swarm_path <- ""
  db_path="~/mkLTG/COInr_for_vtam_2022_05_06_dbV5/"
    fastqdir <- "vtamR_test/data/"
    fastqinfo <- "vtamR_test/data/fastqinfo_zfzr.csv"
    outdir <- "vtamR_test/out_zfzr/"
    mock_composition <- "vtamR_test/data/mock_composition_zfzr.csv"
    asv_list <- "vtamR_test/data/asv_list_zfzr.csv"
  #fastqdir <- "/home/meglecz/vtamR_large_files/fastq/"
  #fastqinfo <- "/home/meglecz/vtamR_large_files/user_input/fastqinfo_mfzr.csv"
  #outdir <- "/home/meglecz/vtamR_large_files/out/"
  #mock_composition <- "local/user_input/mock_composition_mfzr_prerun.csv"

  num_threads=8
  compress = T
} else if (computer == "Endoume"){
  vtam_dir <- "~/vtamR"
  cutadapt_path="/home/emese/miniconda3/bin/"
  vsearch_path = "/home/emese/miniconda3/bin/"
  blast_path= "" # deactivate conda
  swarm_path <- ""
  db_path= "/home/emese/mkCOInr/COInr/COInr_for_vtam_2023_05_03_dbV5/"
#  fastqdir <- "local/fastq/"
  fastqdir <- "vtamR_test/data/"
  fastqinfo <- "vtamR_test/data/fastqinfo_mfzr_gz.csv"
  outdir <- "vtamR_test/out_zfzr/"
  num_threads=8
  compress = T
}else if (computer == "Windows"){
  vtam_dir <- "C:/Users/emese/vtamR/"
  cutadapt_path="C:/Users/Public/"
  vsearch_path = "C:/Users/Public/vsearch-2.23.0-win-x86_64/bin/"
  blast_path="C:/Users/Public/blast-2.14.1+/bin/"
  swarm_path <- "C:/swarm-3.1.4-win-x86_64/bin/"
  db_path="C:/Users/Public/COInr_for_vtam_2023_05_03_dbV5/"
#  fastqdir <- "C:/Users/emese/vtamR_private/fastq/"
  fastqdir <- "vtamR_test/data/"
  fastqinfo <- "vtamR_test/data/fastqinfo_zfzr_gz.csv"
  outdir <- "vtamR_test/out_zfzr/"
  mock_composition <- "vtamR_test/data/mock_composition_zfzr_eu.csv"
  num_threads=4
  compress = F
}
sep=","
setwd(vtam_dir)


taxonomy=paste(db_path, "COInr_for_vtam_taxonomy.tsv", sep="")
blast_db=paste(db_path, "COInr_for_vtam", sep="")


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



#setwd("D:/vtamR")
# load local packages
load_all(".")
roxygenise() # Builds the help files
usethis::use_roxygen_md() # rebuild the help files


###
# Test major functions
###
#test_merge_and_sortreads(test_dir="vtamR_test/", vsearch_path=vsearch_path, cutadapt_path=cutadapt_path)
#test_filters(test_dir="vtamR_test/", vsearch_path=vsearch_path, sep=sep)
#test_make_known_occurrences(test_dir="vtamR_test/", sep=sep)
#test_taxassign(test_dir="vtamR_test/", sep=sep, blast_path=blast_path, blast_db=blast_db, taxonomy=taxonomy, num_threads=num_threads)
#test_optimize(test_dir="vtamR_test/", vsearch_path=vsearch_path, sep=sep)


####
# define input filenames
#fastadir <- "local/mfzr/sorted/"
#fileinfo <- "local/user_input/fileinfo_mfzr_eu.csv"


#fastadir <- "/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/"
#fileinfo <-"/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/fileinfo_vtamr.csv"
#mock_composition <- "/home/meglecz/vtamR/local/user_input/mock_composition_mfzr_prerun.csv"
#sep="\t"

#fastadir <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/"
#fileinfo <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/fileinfo_vtamr.csv"


# create the output directory and check the the slash at the end
outdir <- check_dir(dir=outdir)
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
compress = T
# read fastqinfo
fastqinfo_df <- read.csv(fastqinfo, header=T, sep=sep)
fastainfo_df <- Merge(fastqinfo_df=fastqinfo_df, fastqdir=fastqdir, vsearch_path=vsearch_path, outdir=merged_dir, fastq_ascii=fastq_ascii, fastq_maxdiffs=fastq_maxdiffs, fastq_maxee=fastq_maxee, fastq_minlen=fastq_minlen, fastq_maxlen=fastq_maxlen, fastq_minmergelen=fastq_minmergelen, fastq_maxmergelen=fastq_maxmergelen, fastq_maxns=fastq_maxns, fastq_truncqual=fastq_truncqual, fastq_minovlen=fastq_minovlen, fastq_allowmergestagger=fastq_allowmergestagger, sep=sep, compress=compress)


###
### RandomSeq
###
randomseq_dir = paste(outdir, "random_seq/", sep="")
#fastainfo <- paste(merged_dir, "fastainfo_gz.csv", sep="")
#fastainfo_df <- read.csv(file=fastainfo, header=T, sep=sep)
compress = T
RandomSeq(fastainfo_df, fasta_dir=merged_dir, outdir=randomseq_dir, vsearch_path=vsearch_path, n=10000, randseed=0, compress=compress)

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
compress <- F
sortedinfo_df <- SortReads(fastainfo_df=fastainfo_df, fastadir=randomseq_dir, outdir=sorted_dir, cutadapt_path=cutadapt_path, vsearch_path=vsearch_path, check_reverse=check_reverse, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)


###
### Read input fasta files, dereplicate reads to ASV, and count the number of reads of each ASV in each sample-replicate, add a unique id for ASVs, can take into account ASVs from earlier analyses
###
outfile <- paste(outdir, "1_before_filter.csv", sep="")
sortedinfo_df <- read.csv(paste(sorted_dir, "sortedinfo.csv", sep =""), sep=sep)
updated_asv_list <- sub("\\.", "_updated_2024_02_19.", asv_list) # add date to the name of the input asv_list to get a file name for the updated_file
read_count_df <- read_fastas_from_fileinfo(sortedinfo_df, dir=sorted_dir, outfile=outfile, sep=sep, asv_list=asv_list, updated_asv_list=updated_asv_list)
read_count_df_backup <- read_count_df
read_count_df <- read_count_df_backup
# make stat counts
stat_df <- get_stat(read_count_df, stat_df, stage="Input", params=NA)


###
# Run swarm
###
swarm_d <- 1
fastidious <- TRUE
by_sample <- TRUE
outfile <- paste(outdir, "2_Swarm.csv", sep="")
read_count_df <- swarm(read_count_df, outfile=outfile, swarm_path=swarm_path, num_threads=num_threads, swarm_d=swarm_d, fastidious=fastidious, write_csv=T, sep=sep, by_sample=by_sample)
params <- paste(swarm_d, fastidious, by_sample, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="swarm", params=params)


###
# The safest way to homogenize ASVs and asv_ids among different datasets is updating the asv_list sytematically at the read_fastas_from_fileinfo step.
# This will keep all ASVs ever seen in any of the data sets you have analysed.
# However, if you have many large data sets, the file will grow quickly that can cause memory issues. 
# A quite safe solution is to updates the asv_list after swarm, since swarm has eliminated the majority of the rare ASVs that will be unlikely to be seen among validates ASVs in further runs.
###
#updated_asv_list <- sub("\\.", "_updated_2024_02_19_after_swarm.", asv_list) # add date to the name of the input asv_list to get a file name for the updated_file
#update_asv_list(read_count_df, asv_list=asv_list, outfile=updated_asv_list)


###
### LFN_global_read_count
###
# Eliminate variants with less than global_read_count_cutoff reads in the dataset
global_read_count_cutoff = 2
outfile <- paste(outdir, "3_LFN_global_read_count.csv", sep="")
read_count_df <- LFN_global_read_count(read_count_df, global_read_count_cutoff, outfile=outfile, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="LFN_global_read_count", params=global_read_count_cutoff)

###
### LFN_filters
###
# LFN_read_count
lfn_read_count_cutoff <- 10
outfile <- paste(outdir, "4_LFN_read_count.csv", sep="")
read_count_df_lfn_read_count <- LFN_read_count(read_count_df, cutoff=lfn_read_count_cutoff, outfile=outfile, sep=sep)
stat_df <- get_stat(read_count_df_lfn_read_count, stat_df, stage="LFN_read_count", params=lfn_read_count_cutoff)


# LFN_sample_replicate (by column)
lfn_sample_replicate_cutoff <- 0.001
outfile <- paste(outdir, "5_LFN_sample_replicate.csv", sep="")
read_count_df_lnf_sample_replicate <- LFN_sample_replicate(read_count_df, cutoff=lfn_sample_replicate_cutoff, outfile=outfile, sep=sep)
stat_df <- get_stat(read_count_df_lnf_sample_replicate, stat_df, stage="LFN_sample_replicate", params=lfn_sample_replicate_cutoff)


# LFN_sample_variant (by line)
lnf_variant_cutoff = 0.001
by_replicate = TRUE
outfile <- paste(outdir, "6_LFN_variant_replicate.csv", sep="")
read_count_df_lnf_variant <- LFN_variant(read_count_df, cutoff=lnf_variant_cutoff, by_replicate, outfile=outfile, sep=sep)
param_values <- paste(lnf_variant_cutoff, by_replicate, sep=";")
stat_df <- get_stat(read_count_df_lnf_variant, stat_df, stage="LFN_variant", params=param_values)


# pool the results of the different filterLFN to one data frame; keep only occurrences that passed all filters
outfile <- paste(outdir, "7_pool_LFN.csv", sep="")
read_count_df <- pool_LFN(read_count_df_lfn_read_count, read_count_df_lnf_variant, read_count_df_lnf_sample_replicate, outfile=outfile, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterLFN")
# delete temporary data frames
read_count_df_lfn_read_count <- NULL
read_count_df_lnf_variant <- NULL
read_count_df_lnf_sample_replicate <- NULL


###
### keep repeatable occurrences
###
min_replicate_number <- 2
outfile <- paste(outdir, "8_FilterMinReplicateNumber.csv", sep="")
read_count_df <- FilterMinReplicateNumber(read_count_df, min_replicate_number, outfile=outfile, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterMinReplicateNumber", params=min_replicate_number)

###
### FilerIndel
###
outfile <- paste(outdir, "9_FilterIndel.csv", sep="")
read_count_df <- FilterIndel(read_count_df, outfile=outfile, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterIndel")

###
### FilerCodonStop
###
outfile <- paste(outdir, "10_FilterCodonStop.csv", sep="")
genetic_code = 5
read_count_df <- FilterCodonStop(read_count_df, outfile=outfile, genetic_code=genetic_code, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilerCodonStop", params=genetic_code)


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
stat_df <- get_stat(read_count_df, stat_df, stage="FilerPCRerror", params=params)

###
### FilterChimera
###
abskew=2
by_sample = T
sample_prop = 0.8
outfile <- paste(outdir, "12_FilterChimera.csv", sep="")
read_count_df <- FilterChimera(read_count_df, outfile=outfile, vsearch_path=vsearch_path, by_sample=by_sample, sample_prop=sample_prop, abskew=abskew, sep=sep)
params <- paste(abskew, by_sample, sample_prop, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="FilterChimera", params=params)

###
### FilerRenkonen
###
# Renkonen index:
# PS = summ(min(p1i, p2i))
# p1i = number of reads for variant i in replicate 1 / number of reads in replicate 1
renkonen_distance_quantile = 0.9
outfile <- paste(outdir, "13_FilterRenkonen.csv", sep="")
read_count_df <- FilterRenkonen(read_count_df, outfile=outfile, renkonen_distance_quantile=renkonen_distance_quantile, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilerRenkonen", params=renkonen_distance_quantile)

###
### PoolReplicates
###
digits = 0
outfile <- paste(outdir, "14_PoolReplicates.csv", sep="")
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, outfile=outfile, sep=sep)


###
### TaxAssign
###
outfile <- paste(outdir, "15_TaxAssign.csv", sep="")
asv_tax <- TaxAssign(df=read_count_samples_df, ltg_params_df=ltg_params_df, taxonomy=taxonomy, blast_db=blast_db, blast_path=blast_path, outfile=outfile, num_threads=num_threads)


###
### print output files
###
# write sequence and variant counts after each step
write.csv(stat_df, file = paste(outdir, "stat_steps.csv", sep=""))
# long format, each line corresponds to an occurrence (); if csv files are written at each step, this is not usefull
#write.csv(read_count_samples_df, file = paste(outdir, "16_Final_asvtable_long.csv", sep=""), row.names=F)
# wide format (ASV table), samples are in columns, ASVs in lines
outfile=paste(outdir, "Final_asvtable.csv", sep="")
fileinfo <- paste(sorted_dir, "sortedinfo.csv", sep ="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)
# write ASV table completed by taxonomic assignments
outfile=paste(outdir, "Final_asvtable_with_taxassign.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, asv_tax=asv_tax, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)



# start optimize from almost unfiltered data (after eliminating ASV with low global reads count)
file <- paste(outdir, "2_Swarm.csv", sep="")
read_count_df <- read.csv(file, sep=sep)
dim(read_count_df)
###
### OptimizePCRError
###
outfile <- paste(outdir, "OptimizePCRError.csv", sep="")
OptimizePCRError_df <- OptimizePCRError(read_count_df, mock_composition=mock_composition, sep=sep, outfile=outfile, max_mismatch=1, min_read_count=10)

###
### OptimizeLFNsampleReplicate
###
outfile = paste(outdir, "OptimizeLFNsampleReplicate.csv", sep="")
OptimizeLFNsampleReplicate_df <- OptimizeLFNsampleReplicate(read_count_df, mock_composition=mock_composition, sep=sep, outfile=outfile)

###
### Make known occurrences
###
file <- paste(outdir, "14_PoolReplicates.csv", sep="")
read_count_samples_df <- read.csv(file, sep=sep)
fileinfo <- paste(sorted_dir, "sortedinfo.csv", sep ="")
known_occurrences <- paste(outdir, "known_occurrences.csv", sep= "")
missing_occurrences <- paste(outdir, "missing_occurrences.csv", sep= "")
habitat_proportion= 0.5 # for each asv, if the proportion of reads in a habitat is below this cutoff, is is considered as an artifact in all samples of the habitat
TP_df <- make_known_occurrences(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, known_occurrences=known_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=habitat_proportion)


###
### OptimizeLFNReaCountAndLFNvariant
###
min_replicate_number=2
lfn_sample_replicate_cutoff=0.004
pcr_error_var_prop=0.1

min_lfn_read_count_cutoff=10
max_lfn_read_count_cutoff=100
increment_lfn_read_count_cutoff=5
min_lnf_variant_cutoff=0.001
max_lnf_variant_cutoff=0.05
increment_lnf_variant_cutoff=0.001
by_replicate=FALSE
vsearch_path=""
max_mismatch=1
by_sample=T
sample_prop=0.8
min_replicate_number=2
OptimizeLFNReadCountAndLFNvariant(optimize_read_count_df, known_occurrences=known_occurrences, sep=sep, outdir=optimize_dir, min_lfn_read_count_cutoff=lfn_read_count_cutoff, max_lfn_read_count_cutoff=max_lfn_read_count_cutoff, increment_lfn_read_count_cutoff=increment_lfn_read_count_cutoff, min_lnf_variant_cutoff=min_lnf_variant_cutoff, max_lnf_variant_cutoff=max_lnf_variant_cutoff, increment_lnf_variant_cutoff=increment_lnf_variant_cutoff, by_replicate=by_replicate, lfn_sample_replicate_cutoff=lfn_sample_replicate_cutoff, pcr_error_var_prop=pcr_error_var_prop, vsearch_path=vsearch_path, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, min_replicate_number=min_replicate_number)

end_time <- Sys.time()  # Record the end time
runtime <- end_time - start_time  # Calculate the run time
print(runtime)



###
# PoolReplicates
###
digits = 0
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)
###
# TaxAssign
###
asv_tax <- TaxAssign(df=read_count_samples_df, ltg_params_df=ltg_params_df, taxonomy=taxonomy, blast_db=blast_db, blast_path=blast_path, outdir=outdir, num_threads=num_threads)
# write the list of ASV and their taxonomic assignment
write.csv(asv_tax, file = paste(outdir, "taxa.csv", sep=""), row.names = F)
fileinfo <- "/home/meglecz/vtamR/vtamR_test/out_zfzr/sorted/sortedinfo.csv"
write_asvtable(read_count_samples_df, outfile="vtamR_test/out_zfzr/asvtable_swarm_zfzr.csv", fileinfo=fileinfo, asv_tax=asv_tax, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)


###
# Pool different datasets
###
files <- data.frame(file=c("vtamR_test/out_mfzr/PoolReplicates.csv", "vtamR_test/out_zfzr/PoolReplicates.csv"),
                    marker=c("MFZR", "ZFZR"))

files <- data.frame(file=c("vtamR_test/out_mfzr/PoolReplicates.csv"),
                    marker=c("MFZR"))

read_count_pool <- pool_datasets(files, outdir=outdir, sep=sep, mean_over_markers=T, write_csv=T)



