#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first if_else
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym :=
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL

#' Test_MergeSortReads
#' 
#' Compare the Merge and SortReads output (using default values of vtam) of vtamR 
#' to pre-computed files obtained by vtam
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param vsearch_path path to vsearch executables
#' @param cutadapt_path path to cutadapt executables
#' @param delete_tmp Boolean.  Delete output folder
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @param sep Field separator character in input and output csv files.
#' @export
#'

Test_MergeSortReads <- function(test_dir="vtamR_test/", 
                                vsearch_path="vsearch", 
                                cutadapt_path="cutadapt", 
                                delete_tmp=T, 
                                sep=',', 
                                quiet=T
                                ){
  
  merge_pass <- F
  sortreads_pass <- F
  
  backup_wd <- getwd()
  setwd(test_dir)
  
  fastqinfo_df <- read.csv("data/fastqinfo_test.csv", header=T, sep=sep)
  fastq_dir <- "data/"
  outdir <-paste('out_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  check_dir(outdir)
  
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
  sep <- ","
  compress <- T
  merged_dir <- file.path(outdir, "merged")
  print("Runnig Merge")
  fastainfo_df <- Merge(fastqinfo=fastqinfo_df, 
                        fastq_dir=fastq_dir, 
                        vsearch_path=vsearch_path, 
                        outdir=merged_dir, 
                        fastq_ascii=fastq_ascii, 
                        fastq_maxdiffs=fastq_maxdiffs, 
                        fastq_maxee=fastq_maxee, 
                        fastq_minlen=fastq_minlen, 
                        fastq_maxlen=fastq_maxlen, 
                        fastq_minmergelen=fastq_minmergelen, 
                        fastq_maxmergelen=fastq_maxmergelen, 
                        fastq_maxns=fastq_maxns, 
                        fastq_truncqual=fastq_truncqual, 
                        fastq_minovlen=fastq_minovlen, 
                        fastq_allowmergestagger=fastq_allowmergestagger, 
                        sep=sep, compress=compress, 
                        quiet=quiet
                        )
  
  ### compare results to precomputed files by vtam
  vtam_out <- "vtam/merged/"
  vtamfiles <- list.files(path = vtam_out, pattern = "\\.fasta.gz$")
  vtamRfiles <- list.files(path = merged_dir, pattern = "\\.fasta.gz$")
  
  
  if(length(vtamfiles) == length(vtamRfiles)){# same number of files
    
    for(vtamf in vtamfiles){ # go through all files in vtamf
      vtamRf <- file.path(merged_dir, vtamf)
      vtamf <- file.path(vtam_out, vtamf)
      if(file.exists(vtamRf)){ # corresponding file exists for vtamR
        
        vtamRseq <- read_fasta_seq(filename=vtamRf, dereplicate = T)
        vtamRseq <- vtamRseq %>% arrange(asv)
        vtamseq <- read_fasta_seq(filename=vtamf, dereplicate = T)
        vtamseq <- vtamseq %>% arrange(asv)
        
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
  sorted_dir <- file.path(outdir, "sorted")
  check_reverse <- T
  tag_to_end <- F
  primer_to_end <-F
  cutadapt_error_rate <- 0.1 # -e in cutadapt
  cutadapt_minimum_length <- 50 # -m in cutadapt
  cutadapt_maximum_length <- 500 # -M in cutadapt
  compress <- T
  print("Runnig SortReads")
  sortedinfo_df <- SortReads(fastainfo=fastainfo_df, 
                             fasta_dir=merged_dir, 
                             outdir=sorted_dir, 
                             cutadapt_path=cutadapt_path, 
                             vsearch_path=vsearch_path, 
                             check_reverse=check_reverse, 
                             tag_to_end=tag_to_end, 
                             primer_to_end=primer_to_end, 
                             cutadapt_error_rate=cutadapt_error_rate, 
                             cutadapt_minimum_length=cutadapt_minimum_length, 
                             cutadapt_maximum_length=cutadapt_maximum_length,
                             sep=sep, 
                             compress=compress, 
                             quiet=quiet
                             )
  vtamR_csv <-  file.path(sorted_dir, "fastainfo.csv")
  ### compare output
  vtam_out <-  "vtam/sorted/"
  vtam_csv <-  "vtam/sorted/sortedinfo.tsv"
  fastainfo_vtam_df <- read.csv(vtam_csv, header=T, sep="\t")
  fastainfo_vtam_df <- fastainfo_vtam_df %>% select(-run, -marker)
  
  df <- full_join(fastainfo_vtam_df, sortedinfo_df, by=c("sample", "replicate"))
  
  # check if all output files are present
  if(any(is.na(df$sortedfasta)) | any(is.na(df$fasta))){
    setwd(backup_wd)
    vtam_missing <- df %>%
      filter(is.na(sortedfasta))
    print(vtam_missing)
    
    vtamR_missing <- df %>%
      filter(is.na(fasta))
    print(vtamR_missing)
    stop("Some output files are missing")
  }
  

  for(i in 1:nrow(df)){
    vtamRf <- file.path(sorted_dir, df$fasta[i])
    vtamf <- file.path(vtam_out, df$sortedfasta[i])
    
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
  
  if(delete_tmp){
    unlink(outdir, recursive=TRUE)
  }
  
  if(merge_pass){
    print("PASS: The results of merge are identical for vtam and vtamR")
  }
  if(sortreads_pass){
    print("PASS: The results of sortreads are identical for vtam and vtamR")
  }
  setwd(backup_wd)
}

#' Test_Filters
#' 
#' Run different filtering steps on a test input and compare results to expected 
#' (pre-computed) output
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param vsearch_path path to vsearch executables
#' @param swarm_path Character string: path to swarm executables. 
#' @param sep separator of the csv files
#' @param delete_tmp Boolean. Delete output folder
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @export
#'

Test_Filters <- function(test_dir="vtamR_test/", 
                         vsearch_path="vsearch", 
                         swarm_path="swarm", 
                         sep=",", 
                         delete_tmp=T, 
                         quiet=T
                         ){
  
  check_dir(test_dir)
  outdir <- paste('out_', 
  trunc(as.numeric(Sys.time())), 
  sample(1:100, 1), 
  sep=''
  )
  outdir <- file.path(tempdir(), outdir)
  check_dir(outdir)

  
#  outdir <- file.path(test_dir, "out")
  test_input_file <- file.path(test_dir, "test/test_file.csv")
  
  ### make input df
  input_df <- read_asv_table(filename=test_input_file, sep=sep) %>% 
    rename("asv_id"=seq_id)
  
  #Swarm
  test_input_file_swarm <- file.path(test_dir, "test/test_file_asv_id.csv")
  input_df_swarm <- read_asv_table(filename=test_input_file_swarm, sep=sep) %>% 
    rename("asv_id"=seq_id)
  swarm_out_df <- Swarm(input_df_swarm, swarm_path=swarm_path, by_sample=T, quiet=quiet)
  swarm_out_df$replicate <- as.integer(swarm_out_df$replicate)
  swarm_exp_df = read.csv(file=file.path(test_dir, "test/test_file_asv_id_swarm_out.csv"), sep=",")
  comp_swarm <- compare_df(swarm_out_df, swarm_exp_df, step="Swarm")
  
  
  #LFNglobalReadCount
  global_read_count_cutoff = 50
  global_read_count_cutoff_df <- LFNglobalReadCount(input_df, global_read_count_cutoff, sep=sep)
  global_read_count_cutoff_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_global_read_count50_out.csv"), 
    sep=sep
    )
  comp_LFN_global_read_count <- compare_df(global_read_count_cutoff_df, 
                                           global_read_count_cutoff_exp_df, 
                                           step="LFNglobalReadCount")
  
  
  input_df_tmp <- input_df %>%
    select(-asv_id)
  ### LFN_filters
  # LFNreadCount
  lfn_read_count_cutoff <- 10
  lfn_read_count_df <- LFNreadCount(input_df_tmp, cutoff=lfn_read_count_cutoff, sep=sep)
  lfn_read_count_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_min_read_count_out.csv"), 
    sep=sep
    )
  comp_LFN_read_count <- compare_df(lfn_read_count_df,
                                    lfn_read_count_exp_df, 
                                    step="LFNreadCount")
  
  # LFNsampleReplicate (by column)
  lfn_sample_replicate_cutoff <- 0.001
  lnf_sample_replicate_df <- LFNsampleReplicate(input_df_tmp, 
                                                cutoff=lfn_sample_replicate_cutoff, 
                                                sep=sep
                                                )
  lnf_sample_replicate_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_sample_replicate_out.csv"), 
    sep=sep
    )
  comp_LFN_sample_replicate <- compare_df(lnf_sample_replicate_df, 
                                          lnf_sample_replicate_exp_df, 
                                          step="LFNsampleReplicate"
                                          )
  
  # LFNvariant_replicate (by line)
  lnf_variant_cutoff = 0.002
  by_replicate = TRUE
  lnf_variant_replicate_df <- LFNvariant(input_df, 
                                         cutoff=lnf_variant_cutoff, 
                                         by_replicate=by_replicate, 
                                         sep=sep, 
                                         min_read_count_prop=0.7
                                         )
  lnf_variant_replicate_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_variant_replicate002_out.csv"),
    sep=sep
    )
  comp_LFN_variant_replicate <- compare_df(lnf_variant_replicate_df, 
                                           lnf_variant_replicate_exp_df,
                                           step="LFNvariant_replicate"
                                           )
  
  # LFNvariant (by line)
  lnf_variant_cutoff = 0.002
  by_replicate = FALSE
  lnf_variant_df <- LFNvariant(input_df, 
                               cutoff=lnf_variant_cutoff, 
                               by_replicate=by_replicate, 
                               sep=sep, 
                               min_read_count_prop=0.7
                               )
  lnf_variant_exp_df = read_asv_table(filename=file.path(
    test_dir, "test/test_file_variant002_out.csv"), 
    sep=sep
    )
  comp_LFN_variant <- compare_df(lnf_variant_df, 
                                 lnf_variant_exp_df, 
                                 step="LFNvariant"
                                 )
  
  
  # pool the results of the different filterLFN to one data frame; keep only occurrences that passed all filters
  lfn_pool_df <- PoolFilters(lfn_read_count_df, 
                             lnf_sample_replicate_df, 
                             lnf_variant_df, 
                             sep=sep
                             )
  lnf_pool_exp_df = read_asv_table(filename=file.path(
    test_dir, "test/test_file_pool_LFN_out.csv"), 
    sep=sep
    )
  comp_LFN_variant <- compare_df(lfn_pool_df, 
                                 lnf_pool_exp_df, 
                                 step="PoolFilters"
                                 )
  
  
  ### keep repeatable occurrences
  min_replicate_number <- 2
  FilterMinReplicateNumber_df <- FilterMinReplicate(input_df, 
                                                    min_replicate_number, 
                                                    sep=sep
                                                    )
  FilterMinReplicateNumber_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_repeat_out.csv"), 
    sep=sep
    )
  comp_FilterMinReplicateNumber <- compare_df(FilterMinReplicateNumber_df, 
                                              FilterMinReplicateNumber_exp_df, 
                                              step="FilterMinReplicate"
                                              )
  
  
  ### FilterPCRerror
  pcr_error_var_prop <- 0.1
  max_mismatch <- 1
  by_sample <- T
  sample_prop <- 0.3
  FilterPCRerror_df1 <- FilterPCRerror(input_df, 
                                       vsearch_path=vsearch_path, 
                                       pcr_error_var_prop=pcr_error_var_prop, 
                                       max_mismatch=max_mismatch, 
                                       by_sample=by_sample, 
                                       sample_prop=sample_prop, 
                                       sep=sep
                                       )
  FilterPCRerror_exp_df1 = read_asv_table(
    filename=file.path(test_dir, "test/test_file_pcr1_03_out.csv"), 
    sep=sep
    )
  comp_FilterPCRerror1 <- compare_df(FilterPCRerror_df1, 
                                     FilterPCRerror_exp_df1, 
                                     step="FilterPCRerror"
                                     )
  
  pcr_error_var_prop <- 0.1
  max_mismatch <- 1
  by_sample <- T
  sample_prop <- 0.6
  FilterPCRerror_df1_6 <- FilterPCRerror(input_df, 
                                         vsearch_path=vsearch_path, 
                                         pcr_error_var_prop=pcr_error_var_prop, 
                                         max_mismatch=max_mismatch, 
                                         by_sample=by_sample, 
                                         sample_prop=sample_prop, 
                                         sep=sep
                                         )
  FilterPCRerror_exp_df1_6 = read_asv_table(
    filename=file.path(test_dir, "test/test_file_pcr1_06_out.csv"), 
    sep=sep
    )
  comp_FilterPCRerror1_6 <- compare_df(FilterPCRerror_df1_6, 
                                       FilterPCRerror_exp_df1_6,
                                       step="FilterPCRerror"
                                       )
  
  
  pcr_error_var_prop <- 0.1
  max_mismatch <- 2
  by_sample <- T
  sample_prop <- 0.3
  FilterPCRerror_df2 <- FilterPCRerror(input_df, 
                                       vsearch_path=vsearch_path,
                                       pcr_error_var_prop=pcr_error_var_prop, 
                                       max_mismatch=max_mismatch, 
                                       by_sample=by_sample, 
                                       sample_prop=sample_prop,
                                       sep=sep
                                       )
  FilterPCRerror_exp_df2 = read_asv_table(
    filename=file.path(test_dir, "test/test_file_pcr2_03_out.csv"), 
    sep=sep
    )
  comp_FilterPCRerror_2 <- compare_df(FilterPCRerror_df2,
                                      FilterPCRerror_exp_df2, 
                                      step="FilterPCRerror"
                                      )
  
  ### FilterChimera
  test_input_file_cim <- file.path(test_dir, "test/test_file2.csv")
  input_df_chim <- read_asv_table(filename=test_input_file_cim, sep=sep)
  abskew=10
  by_sample = T
  sample_prop = 0.3
  FilterChimera_10_03_df <- FilterChimera(input_df_chim,
                                          vsearch_path=vsearch_path, 
                                          by_sample=by_sample, 
                                          sample_prop=sample_prop, 
                                          abskew=abskew, sep=sep)
  FilterChimera_10_03_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_chimera_out_03_10.csv"), 
    sep=sep
    )
  comp_FilterChimera_10_03 <- compare_df(FilterChimera_10_03_df, 
                                         FilterChimera_10_03_exp_df, 
                                         step="FilterChimera"
                                         )
  
  abskew=10
  by_sample = T
  sample_prop = 0.6
  FilterChimera_10_06_df <- FilterChimera(input_df_chim, 
                                          vsearch_path=vsearch_path, 
                                          by_sample=by_sample, 
                                          sample_prop=sample_prop, 
                                          abskew=abskew, 
                                          sep=sep
                                          )
  FilterChimera_10_06_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_chimera_out_06_10.csv"), 
    sep=sep
    )
  comp_FilterChimera_10_06 <- compare_df(FilterChimera_10_06_df, 
                                         FilterChimera_10_06_exp_df, 
                                         step="FilterChimera"
                                         )
  
  ### FilterRenkonen
  renkonen_distance_quantile = 0.9
  FilterRenkonen_df <- FilterRenkonen(input_df,
                                      renkonen_distance_quantile=renkonen_distance_quantile, 
                                      sep=sep
                                      )
  FilterRenkonen_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_FilterRenkonen_09_out.csv"), 
    sep=sep
    )
  comp_FilterRenkonen <- compare_df(FilterRenkonen_df, 
                                    FilterRenkonen_exp_df, 
                                    step="FilerRenkonen"
                                    )
  
  renkonen_distance_quantile = 0.8
  FilterRenkonen_df <- FilterRenkonen(input_df, 
                                      renkonen_distance_quantile=renkonen_distance_quantile, 
                                      sep=sep
                                      )
  FilterRenkonen_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_FilterRenkonen_08_out.csv"), 
    sep=sep
    )
  comp_FilterRenkonen <- compare_df(FilterRenkonen_df, 
                                    FilterRenkonen_exp_df, 
                                    step="FilerRenkonen"
                                    )
  
  ### FilerIndel
  FilterIndel_df <- FilterIndel(input_df, sep=sep)
  FilterIndel_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_indel_out.csv"), 
    sep=sep
    )
  comp_FilterIndel <- compare_df(FilterIndel_df, 
                                 FilterIndel_exp_df, 
                                 step="FilterIndel"
                                 )
  
  ### FilerCodonStop
  genetic_code = 5
  FilterCodonStop_df <- FilterCodonStop(input_df, 
                                        genetic_code=genetic_code, 
                                        sep=sep
                                        )
  FilterCodonStop_exp_df = read_asv_table(
    filename=file.path(test_dir, "test/test_file_stop_out.csv"),
    sep=sep
    )
  comp_FilterCodonStop <- compare_df(FilterCodonStop_df, 
                                     FilterCodonStop_exp_df, 
                                     step="FilterCodonStop"
                                     )
  
  ### PoolReplicates
  digits = 0
  PoolReplicates_df <- PoolReplicates(input_df, digits=digits, sep=sep)
  PoolReplicates_exp_df <- read_asv_table_sample(
    filename=file.path(test_dir, "test/test_file_pool_replicate_out.csv"), 
    sep=sep
    )
  PoolReplicates_exp_df <- PoolReplicates_exp_df %>% 
    rename("asv_id"=seq_id)
  comp_PoolReplicates <- compare_df_sample(PoolReplicates_df,
                                           PoolReplicates_exp_df, 
                                           step="PoolReplicates"
                                           )
  
  
  if(delete_tmp){
    unlink(outdir, recursive=TRUE)
  }
  
}


#' read_asv_table
#' 
#' Read asv table in wid format to a data frem in long format
#'  
#' @param filename name of the input file including full path; 
#' columns: asv, seq_id, plate.marker.sample.replicate columns containing read counts
#' @param sep separator; default ","
#' @export 
#' 
read_asv_table <- function(filename, sep=","){
  
  df <- read.csv(filename, sep=sep)
  long_df <- gather(df, key="sample.replicate", value ="read_count", -asv, -seq_id)  %>% 
    filter(read_count > 0)
  # separate column
  long_df <- separate(long_df, "sample.replicate", into=c("sample", "replicate"), sep="\\.")
  
  return(long_df)
}

#' compare_df
#' 
#' Compare two dataframes of and returns teh full join ob the two
#'  
#' @param df1 dataframe with columns: "asv", "sample","replicate","read_count"
#' @param df2 dataframe with columns: "asv", "sample","replicate","read_count"
#' @param step string to include in the FAIL or PASS message 
#' @export 
#' 
compare_df<- function(df1, df2, step=""){
  
  df1 <- df1 %>%
    select(asv, sample,replicate,"read_count_vtamR"="read_count")
  
  df1 <- full_join(df1, df2, by=c("sample", "replicate", "asv"))
  comp <- df1$read_count == df1$read_count_vtamR
  comp[(is.na(comp))] <- FALSE
  
  if(any(!comp)){
    result <- paste(step, ": FAIL", sep="")
  }else{
    result <- paste(step, ": PASS", sep="")
  }
  print(result)
  return(df1)
}

#' read_asv_table_sample
#' 
#' Read asv table in wid format to a data frem in long format
#'  
#' @param filename name of the input file including full path; 
#' columns: asv, seq_id, plate.marker.sample columns containing read counts
#' @param sep separator; default ","
#' @export 
#' 
read_asv_table_sample <- function(filename, sep=","){
  
  df <- read.csv(filename, sep=sep)
  long_df <- gather(df, key="sample", value ="read_count", -asv, -seq_id)  %>% 
    filter(read_count > 0)

  return(long_df)
}

#' compare_df_sample
#' 
#' Compare two dataframes of and returns teh full join ob the two
#'  
#' @param df1 dataframe with columns: "asv", "plate","marker", "sample","read_count"
#' @param df2 dataframe with columns: "asv", "plate","marker", "sample","read_count"
#' @param step string to include in the FAIL or PASS message 
#' @export 
#' 
compare_df_sample<- function(df1, df2, step=""){
  
  df1 <- df1 %>%
    select("asv", "sample","read_count_vtamR"="read_count", "asv_id")
  
  df1 <- full_join(df1, df2, by=c("sample", "asv", "asv_id"))
  comp <- df1$read_count == df1$read_count_vtamR
  comp[(is.na(comp))] <- FALSE
  
  if(any(!comp)){
    result <- paste(step, ": FAIL", sep="")
  }else{
    result <- paste(step, ": PASS", sep="")
  }
  print(result)
  return(df1)
}

#' Test_TaxAssign
#' 
#' Run TaxAssign on a test file and compare the output to mkLTG results
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param sep separator in csv files
#' @param taxonomy file containing the following 
#' columns: tax_id,parent_tax_id,rank,name_txt,old_tax_id (has been merged to another tax_id),
#' taxlevel 
#' (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 2: kingdom, 1: domain, 0: root)
#' @param blast_db BLAST database
#' @param blast_path path to BLAST executable
#' @param num_threads Number of CPUs
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @export
#'

Test_TaxAssign <- function(test_dir="vtamR_test/", 
                           sep=",", 
                           blast_path="blastn", 
                           blast_db="vtamR_test/test/db_test/COInr_reduced", 
                           taxonomy="vtamR_test/test/db_test/taxonomy_reduced.tsv", 
                           num_threads=1, 
                           quiet=T
                           ){
  
  check_dir(test_dir)
  input <- file.path(test_dir, "test/input_taxassign.csv")
  expeted_output <- file.path(test_dir, "test/test_taxassign_out.tsv")

  input_df <- read.csv(input) 
  input_df <- input_df %>%
    select(asv_id, asv)
  
  asv_tax <- TaxAssign(input_df, 
                       taxonomy=taxonomy,
                       blast_db=blast_db, 
                       blast_path=blast_path,
                       num_threads=num_threads, 
                       quiet=quiet
                       )
  
  expected_asv_tax <- read.table(expeted_output, sep="\t", header=T)
  expected_asv_tax <- expected_asv_tax %>%
    select(pid_exp=pid,
           ltg_taxid_exp=ltg_taxid, 
           ltg_name_exp=ltg_name, 
           ltg_rank_exp=ltg_rank,
           superkingdom_exp=superkingdom, 
           kingdom_exp=kingdom, phylum_exp=phylum,
           class_exp=class, order_exp=order, 
           family_exp=family, 
           genus_exp=genus,
           species_exp=species, 
           asv=sequence)
  
  asv_taxassign <- left_join(asv_tax, expected_asv_tax, by="asv")
  
  asv_taxassign$ltg_taxid[is.na(asv_taxassign$ltg_taxid)] <- 0
  asv_taxassign$ltg_taxid_exp[is.na(asv_taxassign$ltg_taxid_exp)] <- 0
  
  tmp_KO <- asv_taxassign %>%
    filter(ltg_taxid!=ltg_taxid_exp)
  
  if(nrow(tmp_KO)==0){
    cat("TaxAssign: PASS")
  }else{
    cat("TaxAssign: FAIL")
  }
}

#' Test_MakeKnownOccurrences
#' 
#' Run MakeKnownOccurrences on a test file and compare the output to expected results
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param sep separator in csv files
#' @param delete_tmp Boolean.  Delete output folder
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @export
#'

Test_MakeKnownOccurrences <- function(test_dir="vtamR_test/", 
                                      sep=",", 
                                      delete_tmp=T, 
                                      quiet=T
                                      ){
  # input dirs and files
  check_dir(test_dir)
  mock_composition <- file.path(test_dir, "test/mock_composition_test.csv")
  sortedinfo <- file.path(test_dir, "test/sortedinfo.csv")
  input_mock_composition <- file.path(test_dir, "test/input_test_known_occurrences.csv")
  # read data for data frame
  read_count_samples_df <- read.csv(input_mock_composition, sep=sep)
  ### Add asv_id to read_count_samples_df
  asv_unique <- read_count_samples_df %>%
    select(asv) %>%
    unique
  asv_unique$asv_id <- rownames(asv_unique)
  read_count_samples_df <- left_join(read_count_samples_df, asv_unique, by="asv")
  ###
  
  # output dirs and filenames
  outdir <-paste(test_dir, '/out_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  check_dir(outdir)

  known_occurrences <- file.path(outdir, "known_occurrences.csv")
  missing_occurrences <- file.path(outdir, "missing_occurrences.csv")
  # params
  # for each asv, if the proportion of reads in a habitat is below this cutoff, 
  # is is considered as an artifact in all samples of the habitat
  habitat_proportion= 0.5
  # run MakeKnownOccurrences
  TP_df <- MakeKnownOccurrences(read_count_samples_df, 
                                sortedinfo=sortedinfo, 
                                mock_composition=mock_composition, 
                                sep=sep, 
                                known_occurrences=known_occurrences, 
                                missing_occurrences=missing_occurrences, 
                                habitat_proportion=habitat_proportion
                                )
  
  # expected results
  expected_known_occurrences <- file.path(test_dir, "test/test_known_occurrences_out.csv")
  expected_known_occurrences_df = read.csv(expected_known_occurrences, sep=sep)
  expected_known_occurrences_df <- expected_known_occurrences_df %>%
    arrange(sample,action,asv)
  
  expected_missing_occurrences <- file.path(test_dir, "test/test_missing_occurrences_out.csv")
  expected_missing_occurrences_df = read.csv(expected_missing_occurrences, sep=sep)
  expected_missing_occurrences_df <- expected_missing_occurrences_df %>%
    arrange(sample,action,asv)
  
  # read and arrange output
  output_missing_occurrences_df = read.csv(missing_occurrences, sep=sep)
  output_missing_occurrences_df <- output_missing_occurrences_df %>%
    select(-asv_id) %>%
    arrange(sample,action,asv)
  
  output_known_occurrences_df = read.csv(known_occurrences, sep=sep)
  output_known_occurrences_df <- output_known_occurrences_df %>%
    select(-asv_id) %>%
    arrange(sample,action,asv)
  
  
  missing <- identical(output_missing_occurrences_df, expected_missing_occurrences_df)
  #   cat("\nOutput missing occurrences correspond to expected:")
  #  print(identical(output_missing_occurrences_df, expected_missing_occurrences_df))
  
  known <- identical(output_known_occurrences_df, expected_known_occurrences_df)
  #  cat("\nOutput known occurrences correspond to expected:")
  #  print(identical(output_known_occurrences_df, expected_known_occurrences_df))
  
  if(delete_tmp){
    unlink(outdir, recursive=TRUE)
  }
  
  if(missing & known){
    print("MakeKnownOccurrences: PASS")
  }else{
    print("MakeKnownOccurrences: FAIL")
  }
}

#' Test_Optimize
#' 
#' Run OptimizePCRerror, OptimizeLFNsampleReplicate and 
#' OptimizeLFNreadCountLFNvariant on a test file and compare the output to expected results
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param vsearch_path path to vsearch executable
#' @param delete_tmp Boolean.  Delete output folder
#' @param sep Field separator character in input and output csv files.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @export
#'
Test_Optimize <- function(test_dir="vtamR_test/", 
                          vsearch_path="vsearch", 
                          delete_tmp=T, 
                          sep=",", 
                          quiet=T
                          ){
  
  check_dir(test_dir)
  outdir <-paste('/out_', 
                 trunc(as.numeric(Sys.time())), 
                 sample(1:100, 1), sep='')
  outdir <- file.path(tempdir(), outdir)
  check_dir(outdir)
  # Attention, if optimize function is modified seriously, 
  # the expected output should be checked
  expected_OptimizePCRError <- file.path(test_dir, "test/OptimizePCRError.csv")
  expected_OptimizePCRError_df <- read.csv(
    expected_OptimizePCRError, 
    sep=",", 
    header=TRUE) %>%
    arrange(expected_asv, unexpected_asv)
  
  expected_OptimizeLFNsampleReplicate <- file.path(
    test_dir, 
    "test/OptimizeLFNsampleReplicate.csv")
  expected_OptimizeLFNsampleReplicate_df <- read.csv(
    expected_OptimizeLFNsampleReplicate, 
    sep=",", 
    header=TRUE
    ) %>%
    select(sample, replicate, action, read_count, 
           read_count_sample_replicate, lfn_sample_replicate_cutoff, asv, taxon) %>%
    arrange(asv, sample, replicate)
  
  expected_OptimizeLFNReadCountAndLFNvariant <- file.path(
    test_dir, "test/OptimizeLFNReadCountAndLFNvariant.csv"
    )
  expected_OptimizeLFNReadCountAndLFNvariant_df <- read.csv(
    expected_OptimizeLFNReadCountAndLFNvariant, 
    sep=",", 
    header=TRUE
    ) %>%
    select(lfn_sample_replicate_cutoff,
           pcr_error_var_prop,
           lfn_read_count_cutoff,
           lnf_variant_cutoff,
           FN,TP,FP
           ) %>%
    arrange(lfn_read_count_cutoff, lnf_variant_cutoff)
  
  # read input info
  input_test_optimize <- file.path(test_dir, "test/input_test_optimize.csv")
  sortedinfo <- file.path(test_dir, "test/sortedinfo.csv")
  sortedinfo_df <- read.csv(sortedinfo, sep=sep)
  mock_composition <- file.path(test_dir, "test/mock_composition_test.csv")
  known_occurrences <- file.path(test_dir, "test/known_occurrences.csv")
  read_count_df <- read.table(input_test_optimize, sep=sep, header=T)
  ### Add asv_id to read_count_df
  asv_unique <- read_count_df %>%
    select(asv) %>%
    unique
  asv_unique$asv_id <- rownames(asv_unique)
  read_count_df <- left_join(read_count_df, asv_unique, by="asv")
  ###
  
  ### PCRerror
  OptimizePCRError_df <- OptimizePCRerror(read_count_df, 
                                          mock_composition=mock_composition,
                                          max_mismatch=1,
                                          min_read_count=10, 
                                          vsearch_path=vsearch_path
                                          )
  OptimizePCRError_df <- OptimizePCRError_df %>%
    select(-expected_asv_id, -unexpected_asv_id) %>%
    arrange(expected_asv, unexpected_asv)
  # compare expected and observed results 
  diff_df <- full_join(OptimizePCRError_df, 
                       expected_OptimizePCRError_df,
                       by = names(OptimizePCRError_df),
                       suffix=c("_obs", "_exp")
                       ) %>%
    filter(rowSums(is.na(.))>0)
  
  if(nrow(diff_df)==0){
    print("OptimizePCRerror: PASS")
  }else{
    print("OptimizePCRerror: FAIL")
    print(diff_df)
  }
  
  ### OptimizeLFNsampleReplicate
  OptimizeLFNsampleReplicate_df <- OptimizeLFNsampleReplicate(
    read_count_df, 
    mock_composition=mock_composition
    )
  OptimizeLFNsampleReplicate_df <- OptimizeLFNsampleReplicate_df %>%
    select(sample, replicate, action, read_count, read_count_sample_replicate, 
           lfn_sample_replicate_cutoff, asv, taxon
           )  %>%
    arrange(asv, sample, replicate)
  # cannot use the full_join then filter for NA, since sme of the read_count can be NA in each of the input df => use identical
  if(identical(OptimizeLFNsampleReplicate_df, expected_OptimizeLFNsampleReplicate_df)){
    print("OptimizeLFNsampleReplicate: PASS")
  }else{
    print("OptimizeLFNsampleReplicate: FAIL")
  }
  
  ### OptimizeLFNreadCountLFNvariant
  OptimizeLFNReadCountAndLFNvariant_df <- OptimizeLFNreadCountLFNvariant(
    read_count_df, 
    known_occurrences=known_occurrences, 
    min_lfn_read_count_cutoff=10, 
    max_lfn_read_count_cutoff=50, 
    increment_lfn_read_count_cutoff=10,
    min_lnf_variant_cutoff=0.001,
    max_lnf_variant_cutoff=0.02,
    increment_lnf_variant_cutoff=0.005,
    by_replicate=TRUE,
    min_replicate_number=2,
    quiet=quiet
    )
  OptimizeLFNReadCountAndLFNvariant_df <- OptimizeLFNReadCountAndLFNvariant_df %>%
    select(lfn_read_count_cutoff,lnf_variant_cutoff,FN,TP,FP) %>%
    arrange(lfn_read_count_cutoff, lnf_variant_cutoff)
  
  # compare expected and observed results 
  diff_df <- full_join(OptimizeLFNReadCountAndLFNvariant_df, 
                       expected_OptimizeLFNReadCountAndLFNvariant_df, 
                       by = names(OptimizeLFNReadCountAndLFNvariant_df), 
                       suffix=c("_obs", "_exp")
                       ) %>%
    filter(rowSums(is.na(.))>0)
  
  if(delete_tmp){
    unlink(outdir, recursive=TRUE)
  }
  
  if(nrow(diff_df)==0){
    print("OptimizeLFNreadCountLFNvariant: PASS")
  }else{
    print("OptimizeLFNreadCountLFNvariant: FAIL")
    print(diff_df)
  }
}
