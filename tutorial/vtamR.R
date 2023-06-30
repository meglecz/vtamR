install.packages("devtools")
install.packages("roxygen2")
install.packages("seqinr")
install.packages("tidyr")

library("devtools")
library("roxygen2")
library("seqinr")
library("dplyr")
library("tidyr")
#library("Biostrings")

setwd("~/vtamR")
#setwd("D:/vtamR")
# load local packages
load_all(".")
roxygenise() # Builds the help files
usethis::use_roxygen_md() # rebuild the help files ?

####
# define input filenames
fastqdir <- "local/small_test"
fileinfo <- "local/user_input/fileinfo_small.csv"


fastqdir <- "local/mfzr/sorted/"
fileinfo <- "local/user_input/fileinfo_mfzr.csv"

fastqdir <- "/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/"
fileinfo <-"/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/fileinfo_vtamr.csv"

fastqdir <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/"
fileinfo <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/fileinfo_vtamr.csv"

# create the output directory and check the the slash at the end
outdir <- check_dir(dir="local/out")

# Measure runtime using system.time()
start_time <- Sys.time()  # Record the start time

# define stat data frame that will be completed with counts after each step
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())

# read input fasta files in fileinfo, demultiplex and count the number of reads in each plate-sample-replicate
read_count_df <- read_fastas_from_fileinfo(file=fileinfo, dir=fastqdir, write_csv=F, outdir=outdir)
# make stat counts
stat_df <- get_stat(read_count_df, stat_df, stage="Input", params=NA)

temp_df <- read_count_df
read_count_df <- temp_df

###
### FilterChimera
###
vsearch_path = ""
abskew=2
by_sample = T
sample_prop = 0.8
read_count_df <- FilterChimera(read_count_df, write_csv=F, outdir=outdir, vsearch_path=vsearch_path, by_sample=by_sample, sample_prop=sample_prop, abskew=abskew)
params <- paste(abskew, by_sample, sample_prop, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="FilterChimera", params=params)


###
### FilerPCRerror
###
pcr_error_var_prop = 0.1
max_mismatch = 1
by_sample = T
sample_prop = 0.8
vsearch_path = ""
read_count_df <- FilterPCRerror(read_count_df, write_csv=F, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop)
params <- paste(pcr_error_var_prop, max_mismatch, by_sample, sample_prop, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="FilerPCRerror", params=params)


###
### LFN_global_read_count
###
# Eliminate variants with less than global_read_count_cutoff reads in the dataset
global_read_count_cutoff = 2
read_count_df <- LFN_global_read_count(read_count_df, global_read_count_cutoff, write_csv=F, outdir=outdir)
stat_df <- get_stat(read_count_df, stat_df, stage="LFN_global_read_count", params=global_read_count_cutoff)


###
### LFN_filters
###
# LFN_read_count
lfn_read_count_cutoff = 10
read_count_df_lfn_read_count <- LFN_read_count(read_count_df, lfn_read_count_cutoff, write_csv=F, outdir = outdir)
stat_df <- get_stat(read_count_df_lfn_read_count, stat_df, stage="LFN_read_count", params=lfn_read_count_cutoff)


# LFN_sample_replicate (by column)
lfn_sample_replicate_cutoff = 0.001
read_count_df_lnf_sample_replicate <- LFN_sample_replicate(read_count_df, lfn_sample_replicate_cutoff, write_csv=F, outdir = outdir)
stat_df <- get_stat(read_count_df_lnf_sample_replicate, stat_df, stage="LFN_sample_replicate", params=lfn_sample_replicate_cutoff)


# LFN_sample_variant (by line)
lnf_variant_cutoff = 0.001
by_replicate = TRUE
read_count_df_lnf_variant <- LFN_variant(read_count_df, lnf_variant_cutoff, by_replicate, write_csv=F, outdir = outdir)
param_values <- paste(lnf_variant_cutoff, by_replicate, sep=";")
stat_df <- get_stat(read_count_df_lnf_variant, stat_df, stage="LFN_variant", params=param_values)


# pool the results of the different filterLFN to one data frame; keep only occurrences that passed all filters
read_count_df <- pool_LFN(read_count_df_lfn_read_count, read_count_df_lnf_variant, read_count_df_lnf_sample_replicate, write_csv=T, outdir = outdir)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterLFN")
# delete temporary data frames
read_count_df_lfn_read_count <- NULL
read_count_df_lnf_variant <- NULL
read_count_df_lnf_sample_replicate <- NULL

###
### keep repeatable occurrences
###
min_relicate_number = 2
read_count_df <- FilterMinReplicateNumber(read_count_df, min_relicate_number, write_csv=F, outdir = outdir)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterMinReplicateNumber", params=min_relicate_number)

###
### FilerIndel
###
read_count_df <- FilterIndel(read_count_df, write_csv=F, outdir=outdir)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterIndel")


###
### FilerCodonStop
###
genetic_code = 5
read_count_df <- FilterCodonStop(read_count_df, write_csv=F, outdir=outdir, genetic_code=genetic_code)
stat_df <- get_stat(read_count_df, stat_df, stage="FilerCodonStop", params=genetic_code)


###
### print output files
###
write.csv(stat_df, file = paste(outdir, "count_stat.csv", sep=""))
write.csv(read_count_df, file = paste(outdir, "read_count_final.csv", sep=""))

end_time <- Sys.time()  # Record the end time
runtime <- end_time - start_time  # Calculate the runtime
print(runtime)

####################################################
####################################################
####################################################
####################################################

# info divers

### run blast with seqinr
library("seqinr")
# Load the query sequence from a file or define it directly
query_sequence <- read.fasta(file = "/home/meglecz/vtamR/local/small_test/small_test1_1.fas")
# Perform local sequence search using blast()
results <- blast(query = query_sequence, database = "/home/meglecz/vtamR/local/small_test/small_test1_1.fas", type = "DNA")
# Print the results
print(results)


if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")
library("DECIPHER")

# Load the query sequence from a file or define it directly
query_sequence <- readDNAStringSet("/home/meglecz/vtamR/local/small_test/small_test1_1.fas")

# Load the reference sequences from a file or define them directly
reference_sequences <- readDNAStringSet("database.fasta")

# Perform local sequence search using blastn()
results <- blastn(query_sequence, reference_sequences)

# Print the results
print(results)

### get the list of function in a package
functions_seqinr <- ls("package:seqinr", all.names = TRUE)
print(functions_seqinr)


### Define the command to run the third-party program
command <- "path/to/program"
args <- c("arg1", "arg2")
# Run the command using system2()
output <- system2(command, args, stdout = TRUE, stderr = TRUE)
# Print the output
print(output)


command <- "path/to/program arg1 arg2"
# Run the command using system()
system(command)



# start blast from R
myPipe <- pipe( "blastall -p blastp -i text.fasta -d data.fasta" )
results <- read.table( myPipe )
colnames( blastResults ) <- c( "QueryID",  "SubjectID", "Perc.Ident",
                               "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                               "S.start", "S.end", "E", "Bits" )


# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html



# Create a connection to an external command (e.g., "ls" command on Unix/Linux)
cmd <- pipe("ls", "r")  # Open the "ls" command for reading

# Read the output from the command
output <- readLines(cmd)

# Close the connection
close(cmd)

# Print the output
print(output)



makeblastdb -in local/out/small/tmp_1687945245/test.fas -dbtype nucl
blastn –db local/out/small/tmp_1687945245/test.fas –query local/out/small/tmp_1687945245/test.fas –outfmt 6 –out local/out/small/tmp_1687945245/test_blastout.out

pipe_vsearch <- pipe(vsearch)
results <- read.table( pipe_vsearch )

blast <- "blastn -db local/out/small/tmp_1687945245/test.fas -query local/out/small/tmp_1687945245/test.fas -outfmt 6"
myPipe <- pipe(blast)
results <- read.table( myPipe )
colnames( results ) <- c( "QueryID",  "SubjectID", "Perc.Ident",
                          "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                          "S.start", "S.end", "E", "Bits" )

