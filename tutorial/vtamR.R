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
vsearch_path = ""
blast_path="/home/meglecz/ncbi-blast-2.11.0+/bin/"
taxonomy="/home/meglecz/mkLTG/COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam_taxonomy.tsv"
blast_db="/home/meglecz/mkLTG/COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam"
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
usethis::use_roxygen_md() # rebuild the help files ?

####
# define input filenames
#fastadir <- "local/small_test"
#fileinfo <- "local/user_input/fileinfo_small.csv"


fastadir <- "local/mfzr/sorted/"
fileinfo <- "local/user_input/fileinfo_mfzr_eu.csv"
mock_composition <- "local/user_input/mock_composition_mfzr_eu.csv"
sep=";"


#fastadir <- "/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/"
#fileinfo <-"/home/meglecz/vtam_benchmark_local/vtam_fish/sorted_mfzr/fileinfo_vtamr.csv"
#mock_composition <- "/home/meglecz/vtamR/local/user_input/mock_composition_mfzr_prerun.csv"
#sep="\t"

#fastadir <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/"
#fileinfo <- "/home/meglecz/vtam_benchmark_local/vtam_bat/fasta/fileinfo_vtamr.csv"

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
read_count_df <- read_fastas_from_fileinfo(file=fileinfo, dir=fastadir, write_csv=F, outdir=outdir, sep=sep)
# make stat counts
stat_df <- get_stat(read_count_df, stat_df, stage="Input", params=NA)

#temp_df <- read_count_df
#read_count_df <- temp_df



###
### LFN_global_read_count
###
# Eliminate variants with less than global_read_count_cutoff reads in the dataset
global_read_count_cutoff = 60
read_count_df <- LFN_global_read_count(read_count_df, global_read_count_cutoff, write_csv=T, outdir=outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="LFN_global_read_count", params=global_read_count_cutoff)

###
### PoolReplicates
###
digits = 0
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)

###
### TaxAssign
###

TaxAssign <- function(read_count_samples_df, ltg_params_df, taxonomy=taxonomy, blast_db=blast_db, blast_path=blast_path){
  
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste(outdir, 'tmp_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  outdir_tmp <- check_dir(outdir_tmp)
  
  # get unique list of ASVs 
  seqs <- unique(read_count_samples_df$asv)
  blast_out <- paste(outdir_tmp, 'blast.out', sep="")
  # run blast
  run_blast(seqs, blast_db=blast_db, blast_path=blast_path, out=blast_out, qcov_hsp_perc=min(ltg_params_df$pcov), perc_identity=min(ltg_params_df$pid), num_threads=8)

  # read taxonomy file; quote="" is important, snce some of the taxon names have quites and this should be ignored
  tax_df <- read.delim(taxonomy, header=T, sep="\t", fill=T, quote="")
  # make data frame with old taxids as line numbers and taxids in a colums
  old_taxid <- tax_df %>%
    filter(!is.na(old_tax_id)) %>%
    select(tax_id, old_tax_id)

  # delete old_tax_ids from tax_df and make taxids unique
  tax_df <- tax_df %>%
    select(-old_tax_id)
  tax_df <- unique(tax_df)

  blast_res <- read.delim(blast_out, header=F, sep="\t", fill=T, quote="")
  colnames(blast_res) <- c("qseqid","sseqid","pident","length","qcovhsp","staxids","evalue") 
  blast_res <- blast_res %>%
    select(qseqid, pident, qcovhsp, staxids)
  
  # replace old taxids (if any) by up to date ones 
  blast_res <- left_join(blast_res, old_taxid, by=c("staxids" = "old_tax_id"))
  blast_res$staxids[which(!is.na(blast_res$tax_id))] <- blast_res$tax_id[which(!is.na(blast_res$tax_id))]
  blast_res <- blast_res %>%
    select(-tax_id)
  
  # add taxonomy info
  blast_res <- left_join(blast_res, tax_df, by=c("staxids" = "tax_id"))  
  
  for(i in 1:length(seqs)){ # go through all sequences 
    df <- blast_res %>%
      filter(qseqid==i)
    
    for(p in 1:nrow(ltg_params_df)){ # for each pid
      pidl <- ltg_params_df[p,"pid"]
      pcovl <- ltg_params_df[p,"pcov"]
      phitl <- ltg_params_df[p,"phit"]
      taxnl <- ltg_params_df[p,"taxn"]
      seqnl <- ltg_params_df[p,"seqn"]
      refresl <- ltg_params_df[p,"refres"]
      ltgresl <- ltg_params_df[p,"ltgres"]
      
      # filter the blastres according to  pid, pcov, refres
      df_intern <- df %>%
        filter(pident>=pidl) %>%
        filter(qcovhsp>=pcovl) %>%
        filter(taxlevel>=refresl)
      
      # check if enough tata and seq among validated hits
      tn <- length(unique(df_intern$staxids))
      if(tn >= taxnl & nrow(df_intern) > seqnl ){
        #        make_ltg(df_intern$staxids)
        break
      }
    }
  }
  
}

make_ltg <- function(taxids){

}


###
### LFN_filters
###
# LFN_read_count
read_count_df_lfn_read_count <- LFN_read_count(read_count_df, lfn_read_count_cutoff, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df_lfn_read_count, stat_df, stage="LFN_read_count", params=lfn_read_count_cutoff)


# LFN_sample_replicate (by column)
read_count_df_lnf_sample_replicate <- LFN_sample_replicate(read_count_df, lfn_sample_replicate_cutoff, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df_lnf_sample_replicate, stat_df, stage="LFN_sample_replicate", params=lfn_sample_replicate_cutoff)


# LFN_sample_variant (by line)
lnf_variant_cutoff = 0.001
by_replicate = TRUE
read_count_df_lnf_variant <- LFN_variant(read_count_df, lnf_variant_cutoff, by_replicate, write_csv=T, outdir = outdir, sep=sep)
param_values <- paste(lnf_variant_cutoff, by_replicate, sep=";")
stat_df <- get_stat(read_count_df_lnf_variant, stat_df, stage="LFN_variant", params=param_values)


# pool the results of the different filterLFN to one data frame; keep only occurrences that passed all filters
read_count_df <- pool_LFN(read_count_df_lfn_read_count, read_count_df_lnf_variant, read_count_df_lnf_sample_replicate, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterLFN")
# delete temporary data frames
read_count_df_lfn_read_count <- NULL
read_count_df_lnf_variant <- NULL
read_count_df_lnf_sample_replicate <- NULL

###
### keep repeatable occurrences
###
read_count_df <- FilterMinReplicateNumber(read_count_df, min_replicate_number, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterMinReplicateNumber", params=min_replicate_number)

###
### FilerPCRerror
###
read_count_df <- FilterPCRerror(read_count_df, write_csv=T, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, sep=sep)
params <- paste(pcr_error_var_prop, max_mismatch, by_sample, sample_prop, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="FilerPCRerror", params=params)

###
### FilterChimera
###
vsearch_path = ""
abskew=2
by_sample = T
sample_prop = 0.8
read_count_df <- FilterChimera(read_count_df, write_csv=T, outdir=outdir, vsearch_path=vsearch_path, by_sample=by_sample, sample_prop=sample_prop, abskew=abskew, sep=sep)
params <- paste(abskew, by_sample, sample_prop, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="FilterChimera", params=params)

###
### FilerRenkonen
###
# Renkonen index:
# PS = summ(min(p1i, p2i))
# p1i = number of reads for variant i in replicate 1 / number of reads in replicate 1
renkonen_distance_quantile = 0.9
read_count_df <- FilerRenkonen(read_count_df, write_csv=T, outdir=outdir, renkonen_distance_quantile=renkonen_distance_quantile, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilerRenkonen", params=renkonen_distance_quantile)

###
### FilerIndel
###
read_count_df <- FilterIndel(read_count_df, write_csv=T, outdir=outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterIndel")

###
### FilerCodonStop
###
genetic_code = 5
read_count_df <- FilterCodonStop(read_count_df, write_csv=T, outdir=outdir, genetic_code=genetic_code, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilerCodonStop", params=genetic_code)

###
### PoolReplicates
###
digits = 0
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)

###
### print output files
###

write.csv(stat_df, file = paste(outdir, "count_stat.csv", sep=""))
write.csv(read_count_samples_df, file = paste(outdir, "Final_asvtable_long.csv", sep=""))
outfile=paste(outdir, "Final_asvtable.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)







###
### PoolReplicates
###
digits = 0
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)


###
### Make known occurrences
###
known_occurrences <- paste(outdir, "known_occurrences.csv", sep= "")
missing_occurrences <- paste(outdir, "missing_occurrences.csv", sep= "")
habitat_proportion= 0.5 # for each asv, if the proportion of reads in a habitat is below this cutoff, is is considered as an artifact in all samples of the habitat
make_known_occurrences(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=known_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=habitat_proportion)


###
### OptimizePCRError
###
optimize_dir = paste(outdir, "optimize", sep="")
OptimizePCRError(read_count_df, mock_composition=mock_composition, sep=sep, outdir=optimize_dir, min_read_count=10)


###
### OptimizeLFNsampleReplicate
###
optimize_dir = paste(outdir, "optimize", sep="")
OptimizeLFNsampleReplicate(read_count_df, mock_composition=mock_composition, sep=sep, outdir=optimize_dir)


###
### OptimizeLFNReaCountAndLFNvariant
###
lfn_read_count_cutoff=10
lnf_variant_cutoff=0.001
by_replicate=T

lfn_sample_replicate_cutoff=0.003

pcr_error_var_prop=0.1
max_mismatch=1
sample_prop=0.8
by_sample=T

min_replicate_number=2

optimize_dir = paste(outdir, "optimize", sep="")

OptimizeLFNReaCountAndLFNvariant(read_count_df, known_occurrences=known_occurrences, sep=sep, outdir=optimize_dir, min_lfn_read_count_cutoff=lfn_read_count_cutoff, min_lnf_variant_cutoff=lnf_variant_cutoff, by_replicate=by_replicate, lfn_sample_replicate_cutoff=lfn_sample_replicate_cutoff, pcr_error_var_prop=pcr_error_var_prop, vsearch_path=vsearch_path, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, min_replicate_number=min_replicate_number)





end_time <- Sys.time()  # Record the end time
runtime <- end_time - start_time  # Calculate the run time
print(runtime)

####################################################
####################################################
####################################################
####################################################

# info divers

#' write_fasta_seq_as_id
#' 
#' Write a fasta file using the sequences as ID
#'  
#' @param sequences list of sequences
#' @param filename output file name, including path
#' @export
#' 
write_fasta_seq_as_id <- function(sequences, filename) {
  # Open the file for writing
  file <- file(filename, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(sequences)) {
    header <- paste0(">", sequences[[i]])
    writeLines(c(header, sequences[[i]], ""), file)
  }
  # Close the file
  close(file)
}



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

