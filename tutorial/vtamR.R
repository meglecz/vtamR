install.packages("devtools")
install.packages("roxygen2")
install.packages("seqinr")
install.packages("tidyr")

library("devtools")
library("roxygen2")
library("seqinr")
library("dplyr")
library("tidyr")
library("utils") # to handle zipped files
#library("Biostrings")

setwd("~/vtamR")
cutadapt_path="/home/meglecz/miniconda3/envs/vtam_2/bin/"
vsearch_path = ""
blast_path="~/ncbi-blast-2.11.0+/bin/" # bombyx
#blast_path="" # endoume deactivate conda
#db_path="~/mkCOInr/COInr/COInr_for_vtam_2023_05_03_dbV5/" # Endoume
db_path="~/mkLTG/COInr_for_vtam_2022_05_06_dbV5/" # Bombyx
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
usethis::use_roxygen_md() # rebuild the help files ?

####
# define input filenames
#fastadir <- "local/small_test"
#fileinfo <- "local/user_input/fileinfo_small.csv"



fastqdir <- "/home/meglecz/vtam_test/example/fastq/"
fastqinfo <- "local/user_input/fastqinfo_mfzr_eu.csv"
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
compress="gz" # "gz" or "zip" for compressing output files; no comprssion by default
merged_dir <- paste(outdir, "merged", sep="")
# read fastqinfo
fastqinfo_df <- read.csv(fastqinfo, header=T, sep=sep)
fastainfo_df <- Merge(fastqinfo_df=fastqinfo_df, fastqdir=fastqdir, vsearch_path=vsearch_path, outdir=merged_dir, fastq_ascii=fastq_ascii, fastq_maxdiffs=fastq_maxdiffs, fastq_maxee=fastq_maxee, fastq_minlen=fastq_minlen, fastq_maxlen=fastq_maxlen, fastq_minmergelen=fastq_minmergelen, fastq_maxmergelen=fastq_maxmergelen, fastq_maxns=fastq_maxns, fastq_truncqual=fastq_truncqual, fastq_minovlen=fastq_minovlen, fastq_allowmergestagger=fastq_allowmergestagger, sep=sep, compress=compress)


sorted_dir <- paste(outdir, "sorted", sep="")
check_reverse <- T
tag_to_end <- F
primer_to_end <-F
cutadapt_error_rate <- 0.1 # -e in cutadapt
cutadapt_minimum_length <- 50 # -m in cutadapt
cutadapt_maximum_length <- 500 # -M in cutadapt

fileinfo_df <- SortReads(fastainfo_df=fastainfo_df, fastadir=merged_dir, outdir=sorted_dir, cutadapt_path=cutadapt_path, check_reverse=check_reverse, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)

SortReads <- function(fastainfo_df, fastadir, outdir="", cutadapt_path="" ,check_reverse=F, tag_to_end=F, primer_to_end=F, cutadapt_error_rate=0.1,cutadapt_minimum_length=50,cutadapt_maximum_length=500, sep=",",  compress=0){
  # run on strand +
  if(check_reverse){
    # use +strand, output to sorted_dir, uncompressed
#    outdir = sorted_dir
    fileinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df, fastadir=merged_dir, outdir=outdir, cutadapt_path=cutadapt_path, check_reverse=F, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=0)

    #### use - strand
    # swap fw and rv tags and primers
    fastainfo_df_tmp <- fastainfo_df %>%
      select(tag_fw_tmp = tag_rv, tag_rv_tmp = tag_fw, primer_fw_tmp = primer_rv, primer_rv_tmp = primer_fw, plate, marker, sample, sample_type,habitat, replicate, fasta) %>%
      select(tag_fw = tag_fw_tmp, tag_rv = tag_rv_tmp, primer_fw = primer_fw_tmp, primer_rv = primer_rv_tmp, plate, marker, sample, sample_type,habitat, replicate, fasta)
    
    # make temp dir 
    outdir <- check_dir(outdir)
    rc_dir <-paste(outdir, 'rc_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    rc_dir <- check_dir(rc_dir)
    # run sortreads on for reverse strand
    fileinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df_tmp, fastadir=merged_dir, outdir=rc_dir, cutadapt_path=cutadapt_path, check_reverse=F, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=0)
    
    
    rc_dir = "/home/meglecz/vtamR/local/out/sorted/rc_169642856935/"
    outdir = "/home/meglecz/vtamR/local/out/sorted/"
    # get list of files demultiplexed on - strand
    files <- list.files(path = rc_dir)
    # Filter the files based on the motif using regular expressions
    files <- grep(pattern = "\\.fasta", x = files, value = TRUE)
    for(i in 1:length(files)){
      i = 1
      plus <- paste(outdir, files[i], sep="")
      minus <- paste(rc_dir, files[i], sep="")
      # open plus file for append
      plus_connection <- file(plus, "a")
      # open minus file for read
      minus_connection <- file(minus, "r")

      # Read the file line by line
#      l <- 0
      while (T) {
        def_line <- readLines(minus_connection, n = 1)
        seq <- readLines(minus_connection, n = 1)
 #       l <- l+2
        rc <- reverse_complement(seq)
        writeLines(c(def_line, rc), con = plus_connection)
      }
      close(minus_connection)
      close(plus_connection)
    }
    
    #!!!!! reverse demultiplexed sequneces
    #!!!!! pool the strand + and - results
  }
  else{
    fileinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df, fastadir=merged_dir, outdir=outdir, cutadapt_path=cutadapt_path, check_reverse=F, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)
  }
  
}
  

SortReads_no_reverse <- function(fastainfo_df, fastadir, outdir="", cutadapt_path="" ,check_reverse=F, tag_to_end=F, primer_to_end=F, cutadapt_error_rate=0.1,cutadapt_minimum_length=50,cutadapt_maximum_length=500, sep=",",  compress=0){
  # do the complete job of demultiplexing and trimming of input file without checking the reverse sequences
  # Makes gz files, but adapt this to windows, later
  # !!!!!!!!!!!!! OK, but deal with no_reverse


  # upper case for all primers and tags
  fastainfo_df$tag_fw <- toupper(fastainfo_df$tag_fw)
  fastainfo_df$tag_rv <- toupper(fastainfo_df$tag_rv)
  fastainfo_df$primer_fw <- toupper(fastainfo_df$primer_fw)
  fastainfo_df$primer_rv <- toupper(fastainfo_df$primer_rv)
  # make a column for output filenames
  fastainfo_df$filename <- NA
  
  # check dirs and make temp dir
  outdir <- check_dir(outdir)
  fastadir<- check_dir(fastadir)

  # get unique list of input fasta files
  fastas <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(fastas)){ # for each input fasta (each of them can be multi run or multi marker)
    # make temp dir for tag file and tag-trimmed files
    tmp_dir <-paste(outdir, 'tmp_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    tmp_dir <- check_dir(tmp_dir)
    
    # select lines in fastainfo_df that corresponds to a given input fasta file
    fasta_file <- fastas[i]
    df <- fastainfo_df %>%
      filter(fasta==fasta_file)
    
    # make a tags.fasta file with all tag combinations of the fasta to be demultiplexed
    tag_file <- make_adapater_fasta(fastainfo_df, fasta_file=fasta_file, tag_to_end=tag_to_end, outdir=tmp_dir)
    # add path
    fasta_file <- paste(fastadir, fasta_file, sep="")
    # demultiplex fasta, write output to tmp file
    demultiplex_cmd = paste(cutadapt_path, "cutadapt --cores=0 -e 0 --no-indels --trimmed-only -g file:", tag_file," -o ", tmp_dir, "tagtrimmed-{name}.fasta ", fasta_file, sep="")
    print(demultiplex_cmd)
    system(demultiplex_cmd)
    
    # get marker list; different markers can have the same tag combinations, so the tagtrimmed files can contain ore than one marker
    markers <- unique(df$marker)
    for(m in 1:length(markers)){ # for each marker trim the sequences from primers
#      m=1
      # select lines for a given input fasta and marker
      df_marker <- df %>%
        filter(marker==markers[m])
      
      # for a given marker, ter is only one primer combination
      primer_fwl <- df_marker[1,"primer_fw"]
      primer_rvl <- df_marker[1,"primer_rv"]
      primer_rvl_rc <- reverse_complement(primer_rvl)
      
      for(f in 1:nrow(df_marker)){# go through each demultiplexed, tagtrimmed file and trim primers
#        f=1
        outfilename <- paste(df_marker[f,"plate"], df_marker[f,"marker"], df_marker[f,"sample"], df_marker[f,"replicate"], sep="-")
        outfilename <- paste(outfilename, ".fasta", sep="")
        if(compress=="gz"){
          outfilename <- paste(outfilename, ".gz", sep="")
        }
        if(compress=="zip"){ # !!!! this might not be operational => work on it whan using windows
          outfilename <- paste(outfilename, ".zip", sep="")
        }
        # complete fastainfo_df with output fasta name
        fastainfo_df$filename[which(fastainfo_df$plate==df_marker[f,"plate"] & fastainfo_df$marker==df_marker[f,"marker"] & fastainfo_df$sample==df_marker[f,"sample"] & fastainfo_df$replicate==df_marker[f,"replicate"])] <- outfilename
        # add path to output file
        primer_trimmed_file <- paste(outdir, outfilename, sep="")
        tag_trimmed_file <- paste(tmp_dir, "tagtrimmed-", df_marker[f,"tag_fw"], "-", df_marker[f,"tag_rv"], ".fasta", sep="")
        if(primer_to_end){
          primer_trim_cmd <- paste(cutadapt_path, "cutadapt --cores=0 -e ",cutadapt_error_rate ," --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length ," --maximum-length ", cutadapt_maximum_length, " -g ^", primer_fwl, "...", primer_rvl_rc, "$ --output ", primer_trimmed_file, " ", tag_trimmed_file, sep="")
        }
        else{
          primer_trim_cmd <- paste(cutadapt_path, "cutadapt --cores=0 -e ",cutadapt_error_rate ," --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length ," --maximum-length ", cutadapt_maximum_length, " -g '", primer_fwl, ";min_overlap=",nchar(primer_fwl),"...", primer_rvl_rc,  ";min_overlap=",nchar(primer_rvl_rc),"' --output ", primer_trimmed_file, " ", tag_trimmed_file, sep="")
        }
        print(primer_trim_cmd)
        system(primer_trim_cmd)
      } # end tagtrimmed within marker
    }# end marker
    # delete the tpp dir wit the tagtrimmed filese
    unlink(tmp_dir, recursive = TRUE)
  }# end fasta
   
  # make sortedinfo file
  fastainfo_df <- fastainfo_df %>%
    select(-tag_fw, -primer_fw, -tag_rv, -primer_rv, -fasta)
  write.table(fastainfo_df, file = paste(outdir, "fileinfo.csv", sep=""),  row.names = F, sep=sep) 
  return(fastainfo_df)
}






# get list of tagtrimmed files from the outdir
tagtrimmed_files <- list.files(path = outdir)
# Filter the files based on the motif using regular expressions
tagtrimmed_files <- grep(pattern = "^tagtrimmed\\.", x = tagtrimmed_files, value = TRUE)




# read input fasta files in fileinfo, demultiplex and count the number of reads in each plate-sample-replicate
read_count_df <- read_fastas_from_fileinfo(file=fileinfo, dir=fastadir, write_csv=F, outdir=outdir, sep=sep)
# make stat counts
stat_df <- get_stat(read_count_df, stat_df, stage="Input", params=NA)



###
### LFN_global_read_count
###
# Eliminate variants with less than global_read_count_cutoff reads in the dataset
global_read_count_cutoff = 2
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
asv_tax <- TaxAssign(df=read_count_samples_df, ltg_params_df=ltg_params_df, taxonomy=taxonomy, blast_db=blast_db, blast_path=blast_path, outdir=outdir)
# write the list of ASV and their taxonomic assignment
write.csv(asv_tax, file = paste(outdir, "taxa.csv", sep=""), row.names = F)
# write ASV table completed by taxonomic assignments
outfile=paste(outdir, "Final_asvtable_with_taxassign.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, asv_tax=asv_tax, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)



###
### LFN_filters
###
# LFN_read_count
lfn_read_count_cutoff <- 10
read_count_df_lfn_read_count <- LFN_read_count(read_count_df, cutoff=lfn_read_count_cutoff, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df_lfn_read_count, stat_df, stage="LFN_read_count", params=lfn_read_count_cutoff)


# LFN_sample_replicate (by column)
lfn_sample_replicate_cutoff <- 0.001
read_count_df_lnf_sample_replicate <- LFN_sample_replicate(read_count_df, cutoff=lfn_sample_replicate_cutoff, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df_lnf_sample_replicate, stat_df, stage="LFN_sample_replicate", params=lfn_sample_replicate_cutoff)


# LFN_sample_variant (by line)
lnf_variant_cutoff = 0.001
by_replicate = TRUE
read_count_df_lnf_variant <- LFN_variant(read_count_df, cutoff=lnf_variant_cutoff, by_replicate, write_csv=T, outdir = outdir, sep=sep)
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
min_replicate_number <- 2
read_count_df <- FilterMinReplicateNumber(read_count_df, min_replicate_number, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterMinReplicateNumber", params=min_replicate_number)

###
### FilerPCRerror
###
pcr_error_var_prop <- 0.1
max_mismatch <- 1
by_sample <- T
sample_prop <- 0.8
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
write_asvtable(read_count_samples_df, asv_tax=asv_tax, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)



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


compressed_file <- compress_file("/home/meglecz/vtamR/local/out/sorted/plate1-MFZR-14ben01-1.fasta", compress="gz", remove_uncompressed=T)
compress_file <- function(file, compress="gz", remove_uncompressed=T){
  #### !!!!! zip is not working correctly in minux, since it produces and archive with the embedded directory structure in the path => correct it when using windows
  outfile <- file
  if(compress == "gz"){
    # Specify the path for the gzipped output file
    file_gz <- paste(file, ".gz", sep="")
    # Open the existing uncompressed file for reading
    file_content <- readBin(file, "raw", file.info(file)$size)
    # Create a gzipped copy of the file
    gz <- gzfile(file_gz, "wb")
    writeBin(file_content, gz)
    close(gz)
    if(remove_uncompressed){
      file.remove(file)
    }
    outfile <- file_gz
  }
  if(compress == "zip"){
    file_zip <- paste(file, ".zip", sep="")
    zip(file_zip, file, recursive = FALSE)
    if(remove_uncompressed){
      file.remove(file)
    }
    outfile <- file_zip
  }
  
  return(outfile)
}


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

