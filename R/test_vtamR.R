#' test_merge_and_sortreads
#' 
#' Compare the Merge and SortReads output (using default values of vtam) of vtamR to pre-computed files obtained by vtam
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param vsearch_path path to vsearch executables
#' @param cutadapt_path path to cutadapt executables
#' @export
#'

test_merge_and_sortreads <- function(test_dir="~/vtamR/vtamR_test/", vsearch_path="", cutadapt_path=""){
  
  merge_pass <- F
  sortreads_pass <- F
  
  backup_wd <- getwd()
  setwd(test_dir)

  fastqinfo_df <- read.csv("data/fastqinfo_mfzr_gz.csv", header=T, sep=sep)
  fastqdir <- "data/"
  outdir <- "out/"
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
  compress <- F
  merged_dir <- paste(outdir, "merged/", sep="")
  fastainfo_df <- Merge(fastqinfo_df=fastqinfo_df, fastqdir=fastqdir, vsearch_path=vsearch_path, outdir=merged_dir, fastq_ascii=fastq_ascii, fastq_maxdiffs=fastq_maxdiffs, fastq_maxee=fastq_maxee, fastq_minlen=fastq_minlen, fastq_maxlen=fastq_maxlen, fastq_minmergelen=fastq_minmergelen, fastq_maxmergelen=fastq_maxmergelen, fastq_maxns=fastq_maxns, fastq_truncqual=fastq_truncqual, fastq_minovlen=fastq_minovlen, fastq_allowmergestagger=fastq_allowmergestagger, sep=sep, compress=compress)
  
  ### compare results to precomputed files by vtam
  vtam_out <- "vtam/merged/"
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
  compress <- F
  
  fileinfo_df <- SortReads(fastainfo_df=fastainfo_df, fastadir=merged_dir, outdir=sorted_dir, cutadapt_path=cutadapt_path, vsearch_path=vsearch_path, check_reverse=check_reverse, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)
  vtamR_csv <-  paste(sorted_dir, "fastainfo.csv", sep="")
  ### compare output
  vtam_out <-  "vtam/sorted/"
  vtam_csv <-  "vtam/sorted/sortedinfo.tsv"
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
#' Run different filtering steps on a test input and compare results to expected (pre-computed) output
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param vsearch_path path to vsearch executables
#' @param sep separator of the csv files
#' @export
#'


test_filters <- function(test_dir="~/vtamR/vtamR_test/", vsearch_path="", sep=","){
  
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
  PoolReplicates_df <- PoolReplicates(input_df, digits=digits, write_csv=F, outdir=outdir, sep=sep)
  PoolReplicates_exp_df <- read_asv_table_sample(filename=paste(test_dir, "test/test_file_pool_replicate_out.csv", sep=""), sep=sep)
  comp_PoolReplicates <- compare_df_sample(PoolReplicates_df, PoolReplicates_exp_df, step="PoolReplicates")
  
  
}

#' test_make_known_occurrences
#' 
#' Run make_known_occurrences on a test file and compare the output to expected results
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param sep separator in csv files
#' @export
#'

test_make_known_occurrences <- function(test_dir="~/vtamR/vtamR_test/", sep=";"){
  # input dirs and files
  test_dir <- check_dir(test_dir)
  mock_composition <- paste(test_dir, "test/mock_composition_test.csv", sep="")
  fileinfo <- paste(test_dir, "test/fileinfo.csv", sep= "")
  input_mock_composition <- paste(test_dir, "test/input_test_known_occurrences.csv", sep="")
  # read data for data frame
  read_count_samples_df <- read.csv(input_mock_composition, sep=sep)
  
  # output dirs and filenames
  outdir <- paste(test_dir, "out/", sep="")
  outdir <- check_dir(outdir)
  known_occurrences <- paste(outdir, "known_occurrences.csv", sep= "")
  missing_occurrences <- paste(outdir, "missing_occurrences.csv", sep= "")
  # params
  habitat_proportion= 0.5 # for each asv, if the proportion of reads in a habitat is below this cutoff, is is considered as an artifact in all samples of the habitat
  # run make_known_occurrences
  make_known_occurrences(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=known_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=habitat_proportion)
  
  # expected results
  expected_known_occurrences <- paste(test_dir, "test/test_known_occurrences_out.csv", sep="")
  expected_known_occurrences_df = read.csv(expected_known_occurrences, sep=sep)
  expected_known_occurrences_df <- expected_known_occurrences_df %>%
    arrange(plate,marker,sample,action,asv)
  
  expected_missing_occurrences <- paste(test_dir, "test/test_missing_occurrences_out.csv", sep="")
  expected_missing_occurrences_df = read.csv(expected_missing_occurrences, sep=sep)
  expected_missing_occurrences_df <- expected_missing_occurrences_df %>%
    arrange(plate,marker,sample,action,asv)
  
  # read and arrange output
  output_missing_occurrences_df = read.csv(missing_occurrences, sep=sep)
  output_missing_occurrences_df <- output_missing_occurrences_df %>%
    arrange(plate,marker,sample,action,asv)
  
  output_known_occurrences_df = read.csv(known_occurrences, sep=sep)
  output_known_occurrences_df <- output_known_occurrences_df %>%
    arrange(plate,marker,sample,action,asv)
  

  missing <- identical(output_missing_occurrences_df, expected_missing_occurrences_df)
#   cat("\nOutput missing occurrences correspond to expected:")
#  print(identical(output_missing_occurrences_df, expected_missing_occurrences_df))
  
  knwon <- identical(output_known_occurrences_df, expected_known_occurrences_df)
#  cat("\nOutput known occurrences correspond to expected:")
#  print(identical(output_known_occurrences_df, expected_known_occurrences_df))
  
  if(missing & knwon){
    print("make_known_occurrences: PASS")
  }else{
    print("make_known_occurrences: FAIL")
  }
}

#' test_taxassign
#' 
#' Run TaxAssign on a test file and compare the output to mkLTG results
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param sep separator in csv files
#' @param taxonomy file containing the following columns: tax_id,parent_tax_id,rank,name_txt,old_tax_id(has been mered to tax_id),taxlevel (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 2: kingdom, 1: superkingdom, 0: root)
#' @param blast_db BLAST database
#' @param blast_path path to BLAST executable
#' @export
#'

test_taxassign <- function(test_dir="~/vtamR/vtamR_test/", sep=sep, blast_path=blast_path, blast_db=blast_db, taxonomy=taxonomy, num_threads=1){
  test_dir <- check_dir(test_dir)
  input <- paste(test_dir, "test/input_taxassign.csv", sep="")
  expeted_output <- paste(test_dir, "test/test_taxassign_out.tsv", sep="")
  outdir <- paste(test_dir, "out/taxassign/", sep="")
  outdir <- check_dir(outdir)
  outfile <- paste(outdir, "out_taxassign.csv", sep="")
  
  input_df <- read.csv(input) 
  input_df <- input_df %>%
    select(asv=sequence)
  
  asv_tax <- TaxAssign(df=input_df, ltg_params_df=ltg_params_df, taxonomy=taxonomy, blast_db=blast_db, blast_path=blast_path, outdir=outdir, num_threads=num_threads)
  
  expected_asv_tax <- read.table(expeted_output, sep="\t", header=T)
  expected_asv_tax <- expected_asv_tax %>%
    select(pid_exp=pid,ltg_taxid_exp=ltg_taxid, ltg_name_exp=ltg_name, ltg_rank_exp=ltg_rank, superkingdom_exp=superkingdom, kingdom_exp=kingdom, phylum_exp=phylum, class_exp=class, order_exp=order, family_exp=family, genus_exp=genus, species_exp=species, asv=sequence)
  
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
  return(asv_taxassign)
}

#' test_optimize
#' 
#' Run OptimizePCRError, OptimizeLFNsampleReplicate and OptimizeLFNReadCountAndLFNvariant on a test file and compare the output to expected results
#'  
#' @param test_dir directory of the test files (Default "~/vtamR/vtamR_test/")
#' @param sep separator in csv files
#' @param vsearch_path path to vsearch executable
#' @export
#'
test_optimize <- function(test_dir="vtamR_test/", vsearch_path=vsearch_path, sep=sep){
  
  test_dir <- check_dir(test_dir)
  outdir <- paste(test_dir, "out/optimize/", sep="")
  outdir <- check_dir(outdir)
  # Attention, if optimize function is modified seriously, the expected output should be checked
  expected_OptimizePCRError <- paste(test_dir, "test/OptimizePCRError.csv", sep="")
  expected_OptimizePCRError_df <- read.csv(expected_OptimizePCRError, sep=";", header=TRUE) %>%
    arrange(expected_asv, unexpected_asv)
  
  expected_OptimizeLFNsampleReplicate <- paste(test_dir, "test/OptimizeLFNsampleReplicate.csv", sep="")
  expected_OptimizeLFNsampleReplicate_df <- read.csv(expected_OptimizeLFNsampleReplicate, sep=";", header=TRUE) %>%
    arrange(asv, plate, marker, sample, replicate)
  
  expected_OptimizeLFNReadCountAndLFNvariant <- paste(test_dir, "test/OptimizeLFNReadCountAndLFNvariant.csv", sep="")
  expected_OptimizeLFNReadCountAndLFNvariant_df <- read.csv(expected_OptimizeLFNReadCountAndLFNvariant, sep=";", header=TRUE) %>%
    arrange(lfn_read_count_cutoff, lnf_variant_cutoff)
  
  # read input info
  input_test_optimize <- paste(test_dir, "test/input_test_optimize.csv", sep="")
  fileinfo <- paste(test_dir, "test/fileinfo.csv", sep="")
  fileinfo_df <- read.csv(fileinfo, sep=sep)
  mock_composition <- paste(test_dir, "test/mock_composition_test.csv", sep="")
  known_occurrences <- paste(test_dir, "test/known_occurrences.csv", sep="")
  read_count_df <- read.table(input_test_optimize, sep=sep, header=T)
  
  ### PCRerror
  OptimizePCRError(read_count_df, mock_composition=mock_composition, sep=sep, outdir=outdir, max_mismatch=1, min_read_count=10)
  output_OptimizePCRError <-  paste(outdir, "OptimizePCRError.csv", sep="")
  OptimizePCRError_df <- read.csv(output_OptimizePCRError, sep=";", header=TRUE) %>%
    arrange(expected_asv, unexpected_asv)
  if(identical(OptimizePCRError_df, expected_OptimizePCRError_df)){
    print("OptimizePCRError: PASS")
  }else{
    print("OptimizePCRError: FAIL")
  }
  ### OptimizeLFNsampleReplicate
  OptimizeLFNsampleReplicate(read_count_df, mock_composition=mock_composition, sep=sep, outdir=outdir)
  output_OptimizeLFNsampleReplicate <-  paste(outdir, "OptimizeLFNsampleReplicate.csv", sep="")
  OptimizeLFNsampleReplicate_df <- read.csv(output_OptimizeLFNsampleReplicate, sep=";", header=TRUE) %>%
    arrange(asv, plate, marker, sample, replicate)
  if(identical(OptimizeLFNsampleReplicate_df, expected_OptimizeLFNsampleReplicate_df)){
    print("OptimizeLFNsampleReplicate: PASS")
  }else{
    print("OptimizeLFNsampleReplicate: FAIL")
  }
  ### OptimizeLFNReadCountAndLFNvariant
  OptimizeLFNReadCountAndLFNvariant(read_count_df, known_occurrences=known_occurrences, sep=sep, outdir=outdir, min_lfn_read_count_cutoff=10, max_lfn_read_count_cutoff=50, increment_lfn_read_count_cutoff=10, min_lnf_variant_cutoff=0.001, max_lnf_variant_cutoff=0.02, increment_lnf_variant_cutoff=0.005, by_replicate=TRUE, lfn_sample_replicate_cutoff=0.001, pcr_error_var_prop=0.1, vsearch_path=vsearch_path, max_mismatch=1, by_sample=TRUE, sample_prop=0.8, min_replicate_number=2)
  output_OptimizeLFNReadCountAndLFNvariant <-  paste(outdir, "OptimizeLFNReadCountAndLFNvariant.csv", sep="")
  OptimizeLFNReadCountAndLFNvariant_df <- read.csv(output_OptimizeLFNReadCountAndLFNvariant, sep=";", header=TRUE)%>%
    arrange(lfn_read_count_cutoff, lnf_variant_cutoff)
  if(identical(OptimizeLFNReadCountAndLFNvariant_df, expected_OptimizeLFNReadCountAndLFNvariant_df)){
    print("OptimizeLFNReadCountAndLFNvariant: PASS")
  }else{
    print("OptimizeLFNReadCountAndLFNvariant: FAIL")
  }
}
