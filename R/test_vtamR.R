#' test_merge_and_sortreads
#' 
#' Compare the Merge and SortReads (using default values of vtam) of vtamR to precoputed file obtained by vtam
#'  
#' @param vsearch_path path to vsearch executables
#' @param cutadapt_path path to cutadapt executables
#' @export
#'

test_merge_and_sortreads <- function(vsearch_path="", cutadapt_path=""){
  
  merge_pass <- F
  sortreads_pass <- F
  
  backup_wd <- getwd()
  #  check_dir(test_dir)
  setwd("~/vtamR")
  ###
  # run Merge using the same parameters as vtam
  fastqinfo_df <- read.csv("vtamR_test/test/user_input/fastqinfo_mfzr_eu_unzipped.csv", header=T, sep=sep)
  fastqdir <- "local/fastq/"
  outdir <- "vtamR_test/test/merged/"
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
  fastainfo_df <- Merge(fastqinfo_df=fastqinfo_df, fastqdir=fastqdir, vsearch_path=vsearch_path, outdir=outdir, fastq_ascii=fastq_ascii, fastq_maxdiffs=fastq_maxdiffs, fastq_maxee=fastq_maxee, fastq_minlen=fastq_minlen, fastq_maxlen=fastq_maxlen, fastq_minmergelen=fastq_minmergelen, fastq_maxmergelen=fastq_maxmergelen, fastq_maxns=fastq_maxns, fastq_truncqual=fastq_truncqual, fastq_minovlen=fastq_minovlen, fastq_allowmergestagger=fastq_allowmergestagger, sep=sep, compress=compress)
  
  ### compare results to precomputed files by vtam
  vtam_out <- "vtamR_test/test/merged/"
  vtamfiles <- list.files(path = vtam_out, pattern = "\\.fasta$")
  vtamRfiles <- list.files(path = outdir, pattern = "\\.fasta$")
  
  
  if(length(vtamfiles) == length(vtamRfiles)){# same number of files
    
    for(vtamf in vtamfiles){ # go through all files in vtamf
      vtamRf <- paste(outdir, vtamf, sep="")
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
  fastainfo_df <- read.csv("vtamR_test/test/merged/fastainfo.csv", header=T, sep=sep)
  sorted_dir <- "vtamR_test/test/sorted/"
  check_reverse <- T
  tag_to_end <- F
  primer_to_end <-F
  cutadapt_error_rate <- 0.1 # -e in cutadapt
  cutadapt_minimum_length <- 50 # -m in cutadapt
  cutadapt_maximum_length <- 500 # -M in cutadapt
  compress <- "0"
  
  fileinfo_df <- SortReads(fastainfo_df=fastainfo_df, fastadir=outdir, outdir=sorted_dir, cutadapt_path=cutadapt_path, vsearch_path=vsearch_path, check_reverse=check_reverse, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)
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
    vtamRseq <- read_fasta_local(filename=vtamRf, use_id=F, case="uc", dereplicate=T)
    vtamRseq <- vtamRseq %>% arrange(sequence)
    vtamseq <- read_fasta_local(filename=vtamf, use_id=F, case="uc", dereplicate=T)
    vtamseq <- vtamseq %>% arrange(sequence)
    
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
}





