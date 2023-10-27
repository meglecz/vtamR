#' test_merge_and_sortreads
#' 
#' Compare the Merge and SortReads (using default values of vtam) of vtamR to precoputed file obtained by vtam
#'  
#' @param vtam_dir directory of vtamR
#' @param vsearch_path path to vsearch executables
#' @param cutadapt_path path to cutadapt executables
#' @export
#'

test_merge_and_sortreads <- function(vtam_dir="~/vtamR", vsearch_path="", cutadapt_path=""){
  
  merge_pass <- F
  sortreads_pass <- F
  
  backup_wd <- getwd()
  setwd(vtam_dir)

  fastqinfo_df <- read.csv("vtamR_test/data/fastqinfo_mfzr_uncompressed.csv", header=T, sep=sep)
  fastqdir <- "vtamR_test/data/"
  outdir <- "vtamR_test/out/"
  outdir <- check_dir(outdir)
  
  ###
  # run Merge using the same parameters as vtam
  fastq_ascii <- 33 #
  fastq_maxdiffs <- 10 #
  fastq_maxee <- 1 #
  fastq_minlen <- 50 #
  fastq_maxlen <- 90000 #
  fastq_minmergelen <- 100 #
  fastq_maxmergelen <-500 #
  fastq_maxns <- 0 #
  fastq_truncqual <- 10 #
  fastq_minovlen <- 50 #
  fastq_allowmergestagger <- F #
  sep <- ";"
  compress <- 0
  merged_dir <- paste(outdir, "merged/", sep="")
  fastainfo_df <- Merge(fastqinfo_df=fastqinfo_df, fastqdir=fastqdir, vsearch_path=vsearch_path, outdir=merged_dir, fastq_ascii=fastq_ascii, fastq_maxdiffs=fastq_maxdiffs, fastq_maxee=fastq_maxee, fastq_minlen=fastq_minlen, fastq_maxlen=fastq_maxlen, fastq_minmergelen=fastq_minmergelen, fastq_maxmergelen=fastq_maxmergelen, fastq_maxns=fastq_maxns, fastq_truncqual=fastq_truncqual, fastq_minovlen=fastq_minovlen, fastq_allowmergestagger=fastq_allowmergestagger, sep=sep, compress=compress)
  
  ### compare results to precomputed files by vtam
  vtam_out <- "vtamR_test/vtam/merged/"
  vtamfiles <- list.files(path = vtam_out, pattern = "\\.fasta$")
  vtamRfiles <- list.files(path = merged_dir, pattern = "\\.fasta$")
  
  
  if(length(vtamfiles) == length(vtamRfiles)){# same number of files
    
    for(vtamf in vtamfiles){ # go through all files in vtamf
      vtamRf <- paste(merged_dir, vtamf, sep="")
      vtamf <- paste(vtam_out, vtamf, sep="")
      if(file.exists(vtamRf)){ # corresponding file exists for vtamR
      
        vtamRseq <- read.fasta(vtamRf, seqonly = T)
        vtamseq <- read.fasta(vtamf, seqonly = T)
        
        if(!(identical(vtamRseq, vtamseq))){# the sequences differ between the two files
          
          setwd(backup_wd)
          stop(paste(vtamRf, "and", vtamf, "are not identical"))
        }
      }else{ # no corresponding to to vtam output
        setwd(backup_wd)
        stop(paste(vtamRf, "does not exists"))
      }
    }
    merge_pass <- T
  }else{
    setwd(backup_wd)
    stop("Different number of output files for vtam and vtamR")
  }
  
  ###
  # run Sortreads using the same parameters as vtam
  fastainfo <- paste(merged_dir, "fastainfo.csv", sep="")
  fastainfo_df <- read.csv(fastainfo, header=T, sep=sep)
  sorted_dir <- paste(outdir, "sorted/", sep="")
  check_reverse <- T
  tag_to_end <- F
  primer_to_end <-F
  cutadapt_error_rate <- 0.1 # -e in cutadapt
  cutadapt_minimum_length <- 50 # -m in cutadapt
  cutadapt_maximum_length <- 500 # -M in cutadapt
  compress <- "0"
  
  fileinfo_df <- SortReads(fastainfo_df=fastainfo_df, fastadir=merged_dir, outdir=sorted_dir, cutadapt_path=cutadapt_path, vsearch_path=vsearch_path, check_reverse=check_reverse, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)
  vtamR_csv <-  paste(sorted_dir, "fastainfo.csv", sep="")
  ### compare output
  vtam_out <-  "vtamR_test/vtam/sorted/"
  vtam_csv <-  "vtamR_test/vtam/sorted/sortedinfo.tsv"
  fastainfo_vtam_df <- read.csv(vtam_csv, header=T, sep="\t")
  fastainfo_vtam_df <- fastainfo_vtam_df %>% rename(plate = run)
  
  df <- full_join(fastainfo_vtam_df, fileinfo_df, by=c("plate", "marker", "sample", "replicate"))
  
  # check if all output files are present
  if(any(is.na(df$sortedfasta)) | any(is.na(df$filename))){
    setwd(backup_wd)
    vtam_missing <- df %>%
      filter(is.na(sortedfasta))
    print(vtam_missing)
    
    vtamR_missing <- df %>%
      filter(is.na(filename))
    print(vtamR_missing)
    stop("Some output files are missing")
  }
  
  # 
  for(i in 1:nrow(df)){
    vtamRf <- paste(sorted_dir, df$filename[i], sep="")
    vtamf <- paste(vtam_out, df$sortedfasta[i], sep="")
    print(vtamRf)
    
    vtamRseq <- read_fasta_seq(filename=vtamRf, dereplicate=T)
    vtamRseq <- vtamRseq %>% arrange(asv)
    vtamseq <- read_fasta_seq(filename=vtamf, dereplicate=T)
    vtamseq <- vtamseq %>% arrange(asv)

    if(!(identical(vtamRseq, vtamseq))){# the sequences differ between the two files
      setwd(backup_wd)
      stop(paste(vtamRf, "and", vtamf, "are not identical"))
    }
  }
  sortreads_pass <- T
  setwd(backup_wd)
  
  if(merge_pass){
    print("PASS: The results of merge are identical for vtam and vtamR")
  }
  if(sortreads_pass){
    print("PASS: The results of sortreads are identical for vtam and vtamR")
  }
  
  setwd(backup_wd)
}

#' test_filters
#' 
#' Compare run different fintering staps on a test input and compare results to expected (precomputed) output
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param vsearch_path path to vsearch executables
#' @param cutadapt_path path to cutadapt executables
#' @param sep separator of the csv files
#' @export
#'


test_filters <- function(test_dir="~/vtamR/vtamR_test/", vsearch_path="", cutadapt_path="", sep=","){
  
  test_dir <- check_dir(test_dir)
  outdir <- paste(test_dir, "out", sep="")
  outdir <- check_dir(outdir)
  test_input_file <- paste(test_dir, "test/test_file.csv", sep="")
  
  ### make input df
  input_df <- read_asv_table(filename=test_input_file, sep=sep)
  
  
  #LFN_global_read_count
  global_read_count_cutoff = 50
  global_read_count_cutoff_df <- LFN_global_read_count(input_df, global_read_count_cutoff, write_csv=F, outdir=outdir, sep=sep)
  global_read_count_cutoff_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_global_read_count50_out.csv", sep=""), sep=sep)
  comp_LFN_global_read_count <- compare_df(global_read_count_cutoff_df, global_read_count_cutoff_exp_df, step="LFN_global_read_count")
  
  
  ### LFN_filters
  # LFN_read_count
  lfn_read_count_cutoff <- 10
  lfn_read_count_df <- LFN_read_count(input_df, cutoff=lfn_read_count_cutoff, write_csv=F, outdir = outdir, sep=sep)
  lfn_read_count_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_min_read_count_out.csv", sep=""), sep=sep)
  comp_LFN_read_count <- compare_df(lfn_read_count_df, lfn_read_count_exp_df, step="LFN_read_count")
  
  # LFN_sample_replicate (by column)
  input_df_tmp <- input_df %>%
    select(-seq_id)
  lfn_sample_replicate_cutoff <- 0.001
  lnf_sample_replicate_df <- LFN_sample_replicate(input_df_tmp, cutoff=lfn_sample_replicate_cutoff, write_csv=F, outdir = outdir, sep=sep)
  lnf_sample_replicate_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_sample_replicate_out.csv", sep=""), sep=sep)
  comp_LFN_sample_replicate <- compare_df(lnf_sample_replicate_df, lnf_sample_replicate_exp_df, step="LFN_sample_replicate")
  
  # LFN_variant_replicate (by line)
  lnf_variant_cutoff = 0.002
  by_replicate = TRUE
  lnf_variant_replicate_df <- LFN_variant(input_df_tmp, cutoff=lnf_variant_cutoff, by_replicate, write_csv=F, outdir = outdir, sep=sep)
  lnf_variant_replicate_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_variant_replicate002_out.csv", sep=""), sep=sep)
  comp_LFN_variant_replicate <- compare_df(lnf_variant_replicate_df, lnf_variant_replicate_exp_df, step="LFN_variant_replicate")
  
  # LFN_variant (by line)
  lnf_variant_cutoff = 0.002
  by_replicate = FALSE
  lnf_variant_df <- LFN_variant(input_df_tmp, cutoff=lnf_variant_cutoff, by_replicate, write_csv=F, outdir = outdir, sep=sep)
  lnf_variant_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_variant002_out.csv", sep=""), sep=sep)
  comp_LFN_variant <- compare_df(lnf_variant_df, lnf_variant_exp_df, step="LFN_variant")
  
  
  # pool the results of the different filterLFN to one data frame; keep only occurrences that passed all filters
  lfn_pool_df <- pool_LFN(lfn_read_count_df, lnf_sample_replicate_df, lnf_variant_df, write_csv=F, outdir = outdir, sep=sep)
  lnf_pool_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_pool_LFN_out.csv", sep=""), sep=sep)
  comp_LFN_variant <- compare_df(lfn_pool_df, lnf_pool_exp_df, step="pool_LFN")
  
  
  ### keep repeatable occurrences
  min_replicate_number <- 2
  FilterMinReplicateNumber_df <- FilterMinReplicateNumber(input_df, min_replicate_number, write_csv=F, outdir = outdir, sep=sep)
  FilterMinReplicateNumber_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_repeat_out.csv", sep=""), sep=sep)
  comp_FilterMinReplicateNumber <- compare_df(FilterMinReplicateNumber_df, FilterMinReplicateNumber_exp_df, step="FilterMinReplicateNumber")
  
  
  ### FilterPCRerror
  pcr_error_var_prop <- 0.1
  max_mismatch <- 1
  by_sample <- T
  sample_prop <- 0.3
  FilterPCRerror_df1 <- FilterPCRerror(input_df, write_csv=F, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, sep=sep)
  FilterPCRerror_exp_df1 = read_asv_table(filename=paste(test_dir, "test/test_file_pcr1_03_out.csv", sep=""), sep=sep)
  comp_FilterPCRerror1 <- compare_df(FilterPCRerror_df1, FilterPCRerror_exp_df1, step="FilterPCRerror")
  
  pcr_error_var_prop <- 0.1
  max_mismatch <- 1
  by_sample <- T
  sample_prop <- 0.6
  FilterPCRerror_df1_6 <- FilterPCRerror(input_df, write_csv=F, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, sep=sep)
  FilterPCRerror_exp_df1_6 = read_asv_table(filename=paste(test_dir, "test/test_file_pcr1_06_out.csv", sep=""), sep=sep)
  comp_FilterPCRerror1_6 <- compare_df(FilterPCRerror_df1_6, FilterPCRerror_exp_df1_6, step="FilterPCRerror")
  
  
  pcr_error_var_prop <- 0.1
  max_mismatch <- 2
  by_sample <- T
  sample_prop <- 0.3
  FilterPCRerror_df2 <- FilterPCRerror(input_df, write_csv=F, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, sep=sep)
  FilterPCRerror_exp_df2 = read_asv_table(filename=paste(test_dir, "test/test_file_pcr2_03_out.csv", sep=""), sep=sep)
  comp_FilterPCRerror_2 <- compare_df(FilterPCRerror_df2, FilterPCRerror_exp_df2, step="FilterPCRerror")
  
  ### FilterChimera
  test_input_file_cim <- paste(test_dir, "test/test_file2.csv", sep="")
  input_df_chim <- read_asv_table(filename=test_input_file_cim, sep=sep)
  abskew=10
  by_sample = T
  sample_prop = 0.3
  FilterChimera_10_03_df <- FilterChimera(input_df_chim, write_csv=F, outdir=outdir, vsearch_path=vsearch_path, by_sample=by_sample, sample_prop=sample_prop, abskew=abskew, sep=sep)
  FilterChimera_10_03_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_chimera_out_03_10.csv", sep=""), sep=sep)
  comp_FilterChimera_10_03 <- compare_df(FilterChimera_10_03_df, FilterChimera_10_03_exp_df, step="FilterChimera")
  
  abskew=10
  by_sample = T
  sample_prop = 0.6
  FilterChimera_10_06_df <- FilterChimera(input_df_chim, write_csv=F, outdir=outdir, vsearch_path=vsearch_path, by_sample=by_sample, sample_prop=sample_prop, abskew=abskew, sep=sep)
  FilterChimera_10_06_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_chimera_out_06_10.csv", sep=""), sep=sep)
  comp_FilterChimera_10_06 <- compare_df(FilterChimera_10_06_df, FilterChimera_10_06_exp_df, step="FilterChimera")
  
  ### FilterRenkonen
  renkonen_distance_quantile = 0.9
  FilterRenkonen_df <- FilterRenkonen(input_df, write_csv=F, outdir=outdir, renkonen_distance_quantile=renkonen_distance_quantile, sep=sep)
  FilterRenkonen_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_FilerRenkonen_09_out.csv", sep=""), sep=sep)
  comp_FilterRenkonen <- compare_df(FilterRenkonen_df, FilterRenkonen_exp_df, step="FilerRenkonen")
  
  renkonen_distance_quantile = 0.8
  FilterRenkonen_df <- FilterRenkonen(input_df, write_csv=F, outdir=outdir, renkonen_distance_quantile=renkonen_distance_quantile, sep=sep)
  FilterRenkonen_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_FilerRenkonen_08_out.csv", sep=""), sep=sep)
  comp_FilterRenkonen <- compare_df(FilterRenkonen_df, FilterRenkonen_exp_df, step="FilerRenkonen")
  
  ### FilerIndel
  FilterIndel_df <- FilterIndel(input_df, write_csv=F, outdir=outdir, sep=sep)
  FilterIndel_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_indel_out.csv", sep=""), sep=sep)
  comp_FilterIndel <- compare_df(FilterIndel_df, FilterIndel_exp_df, step="FilterIndel")
  
  ### FilerCodonStop
  genetic_code = 5
  FilterCodonStop_df <- FilterCodonStop(input_df, write_csv=F, outdir=outdir, genetic_code=genetic_code, sep=sep)
  FilterCodonStop_exp_df = read_asv_table(filename=paste(test_dir, "test/test_file_stop_out.csv", sep=""), sep=sep)
  comp_FilterCodonStop <- compare_df(FilterCodonStop_df, FilterCodonStop_exp_df, step="FilterCodonStop")
  
  ### PoolReplicates
  digits = 0
  PoolReplicates_df <- PoolReplicates(input_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)
  PoolReplicates_exp_df <- read_asv_table_sample(filename=paste(test_dir, "test/test_file_pool_replicate_out.csv", sep=""), sep=sep)
  comp_PoolReplicates <- compare_df_sample(PoolReplicates_df, PoolReplicates_exp_df, step="PoolReplicates")
  
  
}


