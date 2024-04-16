install.packages("devtools")
install.packages("roxygen2")
install.packages("seqinr")
install.packages("tidyr")

library("devtools")
library("roxygen2")
library("seqinr") # splitseq for FilterCodonStop
library("dplyr")
library("tidyr") # gather for read_asv_table; pivot_wider in WriteAsVtable and stat_sample !!sym
#library("utils") # to handle zipped files
library("ggplot2") 

#library("Biostrings")

computer <- "Bombyx" # Bombyx/Endoume/Windows
if(computer == "Bombyx"){
  vtam_dir <- "~/vtamR"
  cutadapt_path="/home/meglecz/miniconda3/envs/vtam_2/bin/"
  vsearch_path = ""
  blast_path="~/ncbi-blast-2.11.0+/bin/" # bombyx
  swarm_path <- ""
  db_path="~/mkLTG/COInr_for_vtam_2022_05_06_dbV5/"
#      fastq_dir <- "vtamR_test/data/"
#      fastqinfo <- "vtamR_test/data/fastqinfo_zfzr.csv"
#      outdir <- "vtamR_test/out_zfzr/"
      mock_composition <- "vtamR_test/data/mock_composition_zfzr.csv"
      asv_list <- "vtamR_test/data/asv_list_zfzr.csv"
  #      fastq_dir <- "/home/meglecz/vtamR_large_files/fastq/"
      #      fastqinfo <- "/home/meglecz/vtamR_large_files/user_input/fastqinfo_mfzr.csv"
      #     outdir <- "/home/meglecz/vtamR_large_files/out/"
     #     mock_composition <- "/home/meglecz/vtamR_large_files/user_input/mock_composition_mfzr.csv"
     #     asv_list <- "/home/meglecz/vtamR_large_files/user_input/asv_list.csv"
      fastq_dir <- "/home/meglecz/Bureau/carp/data/"
      fastqinfo <- "/home/meglecz/Bureau/carp/data/fastainfo_carp.csv"
      outdir <- "/home/meglecz/Bureau/carp/out/"

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
  swarm_path <- "C:/swarm-3.1.4-win-x86_64/bin/"
  db_path="C:/Users/Public/COInr_for_vtam_2023_05_03_dbV5/"
#  fastq_dir <- "C:/Users/emese/vtamR_private/fastq/"
  fastq_dir <- "vtamR_test/data/"
  fastqinfo <- "vtamR_test/data/fastqinfo_mfzr.csv"
  outdir <- "vtamR_test/out_mfzr/"
  mock_composition <- "vtamR_test/data/mock_composition_mfzr.csv"
  asv_list <- "vtamR_test/data/asv_list.csv"
  num_threads=4
  compress = F
}
sep=","
setwd(vtam_dir)

taxonomy=paste(db_path, "COInr_for_vtam_taxonomy.tsv", sep="")
blast_db=paste(db_path, "COInr_for_vtam", sep="")


ltg_params_df = data.frame( pid=c(100,98,95,90,85,80),
                            pcov=c(70,70,70,70,70,70),
                            phit=c(70,70,70,70,70,70),
                            taxn=c(1,1,2,3,4,4),
                            seqn=c(1,1,2,3,4,4),
                            refres=c("species","species","species","genus","family","family"),
                            ltgres=c("species","species","species","species", "genus","genus")
)

ltg_params_df = data.frame( pid=c(100,98,95,90,85,80),
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


# create the output directory and check the the slash at the end
outdir <- check_dir(dir=outdir)
fastq_dir <- check_dir(dir=fastq_dir)
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
fastq_allowmergestagger <- T
quiet <- F
merged_dir <- paste(outdir, "merged/", sep="")
compress = TRUE
# read fastqinfo
fastainfo_df <- Merge(fastqinfo, fastq_dir=fastq_dir, vsearch_path=vsearch_path, outdir=merged_dir, fastq_ascii=fastq_ascii, fastq_maxdiffs=fastq_maxdiffs, fastq_maxee=fastq_maxee, fastq_minlen=fastq_minlen, fastq_maxlen=fastq_maxlen, fastq_minmergelen=fastq_minmergelen, fastq_maxmergelen=fastq_maxmergelen, fastq_maxns=fastq_maxns, fastq_truncqual=fastq_truncqual, fastq_minovlen=fastq_minovlen, fastq_allowmergestagger=fastq_allowmergestagger, sep=sep, compress=compress, quiet=quiet)


####
# Make TrimPrimer
###
fastainfo <- "/home/meglecz/Bureau/carp/out/merged/fastainfo.csv"
fasta_dir <- "/home/meglecz/Bureau/carp/out/merged/"
sorted_dir<- "/home/meglecz/Bureau/carp/out/sorted/"
sortedinfo_df <- check_dir(dir=outdir)
check_reverse <- F
primer_to_end <-F
cutadapt_error_rate <- 0.1 # -e in cutadapt
cutadapt_minimum_length <- 50 # -m in cutadapt
cutadapt_maximum_length <- 500 # -M in cutadapt
compress <- F
quiet <- F

fastainfo_df <- TrimPrimer(fastainfo, fasta_dir=fasta_dir, outdir=sorted_dir, compress=compress, cutadapt_path=cutadapt_path, vsearch_path=vsearch_path, check_reverse=check_reverse, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, quiet=quiet)



###
### Read input fasta files, dereplicate reads to ASV, and count the number of reads of each ASV in each sample-replicate, add a unique id for ASVs, can take into account ASVs from earlier analyses
###
outfile <- paste(outdir, "1_before_filter.csv", sep="")
sortedinfo_df <- read.csv(paste(sorted_dir, "sortedinfo.csv", sep =""), sep=sep)
updated_asv_list <- paste(outdir, "asv_list.csv", sep="") # add date to the name of the input asv_list to get a file name for the updated_file
read_count_df <- read_fastas_from_sortedinfo(sortedinfo_df, dir=sorted_dir, outfile=outfile, sep=sep, updated_asv_list=updated_asv_list, quiet=quiet)
# make stat counts
stat_df <- get_stat(read_count_df, stat_df, stage="Input", params=NA)

###
# Run swarm
###
quiet <- F
swarm_d <- 1
fastidious <- TRUE
by_sample <- FALSE
outfile <- paste(outdir, "2.csv", sep="")
read_count_df <- Swarm(read_count_df, outfile=outfile, swarm_path=swarm_path, num_threads=num_threads, swarm_d=swarm_d, fastidious=fastidious, sep=sep, by_sample=by_sample, quiet=quiet)
params <- paste(swarm_d, fastidious, by_sample, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="swarm_by_sample", params=params)

###
### LFN_global_read_count
###
# Eliminate variants with less than global_read_count_cutoff reads in the dataset
global_read_count_cutoff = 2
outfile <- paste(outdir, "3_LFN_global_read_count.csv", sep="")
read_count_df <- LFN_global_read_count(read_count_df, cutoff=global_read_count_cutoff, outfile=outfile)
stat_df <- get_stat(read_count_df, stat_df, stage="LFN_global_read_count", params=global_read_count_cutoff)

###
### FilterChimera
###
quiet <- F
abskew=4
by_sample = F
sample_prop = 0.8
quiet = F
outfile <- paste(outdir, "4_FilterChimera.csv", sep="")
read_count_df <- FilterChimera(read_count_df, outfile=outfile, vsearch_path=vsearch_path, by_sample=by_sample, sample_prop=sample_prop, abskew=abskew, sep=sep, quiet=quiet)
params <- paste(abskew, by_sample, sample_prop, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="FilterChimera", params=params)


###
### PoolReplicates
###
digits = 0
outfile <- paste(outdir, "5_PoolReplicates.csv", sep="")
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, outfile=outfile, sep=sep)
stat_df <- get_stat(read_count_samples_df, stat_df, stage="PoolReplicates")

###
### TaxAssign
###
outfile <- paste(outdir, "TaxAssign.csv", sep="")
asv_tax <- TaxAssign(asv=read_count_samples_df, ltg_params=ltg_params_df, taxonomy=taxonomy, blast_db=blast_db, blast_path=blast_path, outfile=outfile, num_threads=num_threads)

###
### print output files
###
outfile <- paste(outdir, "ASV_table_with_TaxAssign.csv", sep="")
asv_table_df <- WriteASVtable(read_count_samples_df, outfile=outfile, asv_tax=asv_tax, sortedinfo=sortedinfo_df, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, sep=sep)

