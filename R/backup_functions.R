#' FilerPCRerror_NW
#' 
#' Filter out all ASVs if they very similar (max_mismatch) to another more frequent ASV (pcr_error_var_prop)
#' The whole plate can be analysed et once (by_sampe=F)
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @param pcr_error_var_prop if the proportion of read_counts of two similar ASVs is bellow pcr_error_var_prop, the less abundant is flagged as a PCR error
#' @param max_mismatch maximum number of mismatches (gaps included) to consider two ASVs as similar
#' @param by_sample T/F, if T ASVs are flagged as an PCR error separately for each sample
#' @param sample_prop if by_sample=T, the ASV must be flagged as an PCRerror in sample_prop of the cases to be eliminated
#' @export
#' 
FilterPCRerror_NW <- function(read_count_df, by_sample=T, write_csv=F, outdir=outdir, pcr_error_var_prop=0.1, max_mismatch=1, sample_prop=0.8){
  
  # get unique list of ASVs with their total read_count in the run
  unique_asv_df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count = sum(read_count)) %>%
    arrange(desc(read_count))
  
  if(by_sample){ # sample by sample
    sample_list <- unique(read_count_df$sample)
    # loop over samples
    for(sample_loc in sample_list){
      # get unique list of ASVs with their total read_count in the sample
      unique_asv_df_sample <- read_count_df %>%
        filter(sample == sample_loc) %>%
        group_by(asv)%>%
        summarize(read_count = sum(read_count))%>%
        arrange(desc(read_count))
      # flag PCR errors
      unique_asv_df_sample <- flagPCRerrors(unique_asv_df_sample, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch)
      # remove read_count column
      unique_asv_df_sample <- unique_asv_df_sample[, -which(names(unique_asv_df_sample) == "read_count")]
      # add a column for for each sample to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df <- left_join(unique_asv_df, unique_asv_df_sample, by = "asv")
    }
  }
  else{ # whole dataset
    # add a column for for each sample to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
    unique_asv_df <- flagPCRerrors(unique_asv_df_sample, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch)
  }
  
  # count the number of times each ASV has been flagged and when it has not. Ignore NA, when the ASV is not present in the sample
  unique_asv_df$yes <- rowSums(unique_asv_df[3:ncol(unique_asv_df)] == 1, na.rm = TRUE)
  unique_asv_df$no <- rowSums(unique_asv_df[3:(ncol(unique_asv_df)-1)] == 0, na.rm = TRUE)
  # keep only ASVs, that are not flagged in sample_prop proportion of the samples where they are present  
  unique_asv_df <- unique_asv_df %>%
    filter(no/(yes+no) >= (1-sample_prop))
  
  # eliminate potential PCRerrors from read_count_df
  read_count_df <- read_count_df %>%
    filter(asv %in% unique_asv_df$asv)
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "Input.csv", sep=""))
  }
  
  return(read_count_df)
}

#' flagPCRerrors
#' 
#' Flag unique ASVs if they can be a PCR error
#' 
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param pcr_error_var_prop if the proportion of read_counts of two similar ASVs is bellow pcr_error_var_prop, the less abundant is flagged as a PCR error
#' @param max_mismatch maximum number of mismatches (gaps included) to consider two ASVs as similar
#' @export
#' 
flagPCRerrors <- function(unique_asv_df, pcr_error_var_prop=0.1, max_mismatch=1){
  
  # define DNAfull substitution matrix
  # https://rosalind.info/glossary/dnafull/
  scores <- matrix(c(5,-4,-4,-4,-4, 1, 1,-4,-4, 1,-4,-1,-1,-1,-2,
                     -4, 5,-4,-4,-4, 1,-4, 1, 1,-4,-1,-4,-1,-1,-2,
                     -4,-4, 5,-4, 1,-4, 1,-4, 1,-4,-1,-1,-4,-1,-2,
                     -4,-4,-4, 5, 1,-4,-4, 1,-4, 1,-1,-1,-1,-4,-2,
                     -4,-4, 1, 1,-1,-4,-2,-2,-2,-2,-1,-1,-3,-3,-1,
                     1, 1,-4,-4,-4,-1,-2,-2,-2,-2,-3,-3,-1,-1,-1,
                     1,-4, 1,-4,-2,-2,-1,-4,-2,-2,-3,-1,-3,-1,-1,
                     -4, 1,-4, 1,-2,-2,-4,-1,-2,-2,-1,-3,-1,-3,-1,
                     -4, 1, 1,-4,-2,-2,-2,-2,-1,-4,-1,-3,-3,-1,-1,
                     1,-4,-4, 1,-2,-2,-2,-2,-4,-1,-3,-1,-1,-3,-1,
                     -4,-1,-1,-1,-1,-3,-3,-1,-1,-3,-1,-2,-2,-2,-1,
                     -1,-4,-1,-1,-1,-3,-1,-3,-3,-1,-2,-1,-2,-2,-1,
                     -1,-1,-4,-1,-3,-1,-3,-1,-3,-1,-2,-2,-1,-2,-1 ,
                     -1,-1,-1,-4,-3,-1,-1,-3,-1,-3,-2,-2,-2,-1,-1,
                     -2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1), nrow=15, ncol=15)
  nt_list <- c("A" ,  "T" ,  "G" ,  "C"  , "S" ,  "W"  , "R"   ,"Y" , "K" ,  "M" ,  "B" ,  "V"  , "H" ,  "D"  , "N")
  DNAfull <-  matrix(scores, nrow = 15, ncol = 15, dimnames = list(nt_list, nt_list))
  
  # add an error column
  unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
  
  # go through all pairs of ASVs, starting with the 2nd most abundant
  for(fille in 2:length(unique_asv_df$asv)){ # potential daughter asvs
    for(parent in 1:(fille-1)){ # potential parents (higher of equal read_counts; then daughter)
      if((unique_asv_df$read_count[parent] * pcr_error_var_prop) >= unique_asv_df$read_count[fille] ){ # read_count of parent is high compared to daughter
        # align the two sequences
        #methods(class = class(alignment))
        alignment <- pairwiseAlignment(pattern = DNAString(unique_asv_df$asv[parent]), subject = DNAString(unique_asv_df$asv[fille]), substitutionMatrix = DNAfull, type = "global")
        # get he number of mismatches (gaps included)
        nmismatch <- (width(alignedPattern(alignment)) - nmatch(alignment))
        if(nmismatch <= max_mismatch){ # mark the daughter sequence as a probable error
          unique_asv_df$PCRerror[fille] <- 1
          break # we found a parent, stop searching
        }
      }
    }
  }
  return(unique_asv_df)
}

#' Read all fasta files in fileinfo file to a data frame
#' 
#' Read all fasta file in fileinfo_df data frame (plate	marker	sample	replicate	filename).
#' Dereplicate reads to ASVs.
#' Count the number of reads of each ASVs in each plate-marker-sample-replicate.
#' Returns a df ("asv", "plate", "marker", "sample", "replicate", "read_count").
#' 
#' @param fileinfo_df data frame with columns:  plate, marker, sample,  replicate, filename, (optional: sample_type(mock/negative/real), habitat)
#' @param dir name of the directory with fasta file 
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @export
read_fastas_from_fileinfo_backup <- function (fileinfo_df, dir="", write_csv=F, outdir=NA, sep=",") {
  # read all fasta files in fileinfo to a read_count_df
  if(nchar(dir)>0){
    dir <- check_dir(dir)
  }
  # define empty read_count_df to pool the results of variables
  read_count_df <- data.frame(asv=character(),
                              read_count=integer(),
                              plate=character(),
                              marker=character(),
                              sample=character(),
                              replicate=character())
  # read all fasta files in fileinfo and count the reads
  for(i in 1:length(fileinfo_df$filename)){
    fas <- paste(dir, fileinfo_df$filename[i], sep = "")
    print(fas)
    if(file.size(fas) < 50){ #' an empty file gzipped has 42 size => skip these files, An unzipped fasta with 80 nt is bigger than 50
      next
    }
    
    read_count_df_tmp <- read_fasta_count_reads(fas)
    read_count_df_tmp$plate <- fileinfo_df[i,"plate"]
    read_count_df_tmp$marker <- fileinfo_df[i,"marker"]
    read_count_df_tmp$sample <- fileinfo_df[i,"sample"]
    read_count_df_tmp$replicate <- fileinfo_df[i,"replicate"]
    read_count_df <- rbind(read_count_df, read_count_df_tmp)
  }
  rm(read_count_df_tmp)
  # reorder columns
  read_count_df <- read_count_df[, c("asv", "plate", "marker", "sample", "replicate", "read_count")]
  
  if(write_csv){
    outdir <- check_dir(outdir)
    write.table(read_count_df, file = paste(outdir, "Input.csv", sep=""),  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' Read a fasta file, demultiplex reads and count them
#' 
#' Read only sequences from input fasta file.
#' Transform all nt to upper case.
#' Demultiplex reads to ASVs.
#' Count the number of reads of each ASVs.
#' Return a tibble with asv and nb_reads.
#' 
#' @param file fasta file
#' @export
read_fasta_count_reads_backup <- function (file) {
  fas <- read.fasta(file, seqonly = T)
  # transform list to matrix (1 column)
  fas <-as.matrix(fas)
  # all nt to upper case letters
  fas <-toupper(fas)
  # transform matrix to df. If I want to do it directly from list, it will put each sequence to one separate column
  fas <- as.data.frame(fas)
  colnames(fas) <- c("asv")
  
  # make unique list of sequences and count them
  read_count_df <- fas %>%
    group_by(asv) %>%
    summarize(read_count=length(asv))
  return(read_count_df)
}

#' select_sequences_backup
#' 
#' Select sequences from the input file that correspond to the vector of integers
#'  
#' @param file input fasta file; can be  uncompressed o compressed in gz format only 
#' @param fasta_dir path to the input fasta files
#' @param outdir directory for the output files; same compression as the input file
#' @param random_integers vector of random integers between 1 and the number of sequences in the input file
#' @export
#'
select_sequences_backup <- function(fasta_dir="", file, outdir="", random_integers){
  
  fasta_dir<- check_dir(fasta_dir)
  outdir<- check_dir(outdir)
  # can deal with uncompressed files and gz compressed files. Zip files should be decompressed previously
  outfile <- paste(outdir, file, sep="")
  #  outfile <- sub(".gz", "", outfile)
  
  
  filename <- paste(fasta_dir, file, sep="")
  if(endsWith(filename, ".gz")){
    file_connection <- gzfile(filename, "rb") 
    outfile_connection <- gzfile(outfile, "wb")
  }else{
    file_connection <- file(filename, "r")
    outfile_connection <- file(outfile, "w")
  }
  
  random_integers <- sort(random_integers)
  count = 0
  bool = F
  while (length(line <- readLines(file_connection, n = 1)) > 0) {
    if(startsWith(line, ">")){
      bool = F
      if(length(random_integers)==0){
        break
      }
      count <- count + 1
      if(count == random_integers[1]){
        writeLines(line, con=outfile_connection)
        bool = T
        random_integers <- random_integers[-1]
      }
    }else if(bool){
      writeLines(line, con=outfile_connection)
    }
  }
  close(outfile_connection)
  close(file_connection)
  
  #  if(endsWith(file, ".gz")){
  #    outfile_gz <- paste(outfile_gz, ".gz", sep="")
  
  #  }
}

#' count_seq_backup
#' 
#' Count the number of sequences in the input fasta file. Input can be uncompressed, gz, or zip file
#'  
#' @param file fasta file, without it's path
#' @param dir directory that contains the input file
#' @export
#' 
count_seq_backup <- function(dir="", file){
  
  dir <- check_dir(dir)
  file <- paste(dir, file, sep="")
  
  count = 0
  # open the file for reading in function of its type
  if(endsWith(file, ".gz")){
    file_connection <- gzfile(file, "rb") 
  }else{
    file_connection <- file(file, "r")
  }
  
  # count sequences
  while (length(line <- readLines(file_connection, n = 1)) > 0) {
    if(startsWith(line, ">")){
      count <- count + 1
    }
  }
  # close file
  close(file_connection)
  return(count)
}


#' RandomSeq_backup
#' 
#' Random select n sequences from each input fasta file. The output is the same compression type (if any) as the input
#'  
#' @param fastainfo_df data frame with a 'fasta' column containing input file names; files can be compressed in gz and zip format
#' @param n integer; the number of randomly selected sequences 
#' @param fasta_dir directory that contains the input fasta files
#' @param outdir directory for the output files
#' @export
#' 
RandomSeq_backup <- function(fastainfo_df, fasta_dir="", outdir="", n){
  
  fasta_dir<- check_dir(fasta_dir)
  outdir<- check_dir(outdir)
  
  unique_fasta <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(unique_fasta)){ # go throgh all fasta files
    input_fasta <- unique_fasta[i]
    # Zip files should be unzipped once for sequence count and select_seq, than results re-zipped, wd is changed several times
    if(endsWith(input_fasta, ".zip")){
      # change to fasta_dir
      backup_wd <- getwd()
      setwd(fasta_dir)
      unzip(input_fasta)
      unzipped_file <- sub(".zip", "", input_fasta)
      # count the number of sequences in the input file
      seq_n <- count_seq(dir="", file=unzipped_file)
      # Generate k random integers between 1 and n
      set.seed(Sys.time())
      random_integers <- sample(1:seq_n, size = n, replace = FALSE)
      # write sequences corresponding to the random numbers to an output file
      setwd(backup_wd) # get back to the original wd, so input and output filepath can be habled more easily
      select_sequences(fasta_dir=fasta_dir, file=unzipped_file, outdir=outdir, random_integers)
      # delete unzipped input file
      setwd(fasta_dir)
      file.remove(unzipped_file)  # remove unzipped input file
      # zip output file
      setwd(backup_wd)
      setwd(outdir)
      zip(input_fasta, unzipped_file) # zip outfile
      file.remove(unzipped_file) # rm unzipped outfile 
      setwd(backup_wd)
    }else{
      seq_n <- count_seq(dir=fasta_dir, file=input_fasta)
      # Generate k random integers between 1 and n
      set.seed(Sys.time())
      random_integers <- sample(1:seq_n, size = n, replace = FALSE)
      select_sequences(fasta_dir=fasta_dir, file=input_fasta, outdir=outdir, random_integers)
    }
    
  }
  
}


RandomSeq_with_zip <- function(fastainfo_df, fasta_dir="", outdir="", n, randseed=0){
  # zipping and unzipping is quite long. It is probably better to work on uncompressed files on windows
  # quite fast for uncompressed ad gz files
  fasta_dir<- check_dir(fasta_dir)
  outdir<- check_dir(outdir)
  
  unique_fasta <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(unique_fasta)){ # go through all fasta files
    input_fasta <- unique_fasta[i]
    # Zip files should be unzipped once for sequence count and select_seq, than results re-zipped, wd is changed several times
    if(F && endsWith(input_fasta, ".zip")){
      # change to fasta_dir
      backup_wd <- getwd()
      setwd(fasta_dir)
      unzip(input_fasta)
      setwd(backup_wd) # get back to the original wd, so input and output filepath can be handled more easily
      unzipped_file <- sub(".zip", ".fasta", input_fasta)
      input_fasta_p <- paste(fasta_dir, unzipped_file, sep="")
      output_fasta_p <- paste(outdir, unzipped_file, sep="")
      # count the number of sequences in the input file
      seq_n <- count_seq(filename=input_fasta_p)
      #random sample with vsearch
      options(scipen=100)
      vsearch_cmd <- paste("vsearch --fastx_subsample ", input_fasta_p, " --fastaout ", output_fasta_p, " --sample_size ", n, " --randseed ", randseed, sep="")
      print(vsearch_cmd)
      system(vsearch_cmd)
      options(scipen=0)
      
      # delete unzipped input file
      setwd(fasta_dir)
      file.remove(unzipped_file)  # remove unzipped input file
      # zip output file
      setwd(backup_wd)
      setwd(outdir)
      zip(input_fasta, unzipped_file) # zip outfile
      file.remove(unzipped_file) # rm unzipped outfile 
      setwd(backup_wd)
    }else{ # uncompressed or gz
      input_fasta_p <- paste(fasta_dir, input_fasta, sep="")
      output_fasta_p <- paste(outdir, input_fasta, sep="")
      seq_n <- count_seq(filename=input_fasta_p)
      if(n > seq_n ){ # not enough seq
        file.copy(input_fasta_p, output_fasta_p)
        msg <- paste("WARNING:", input_fasta_p,"has",seq_n,"sequences, which is lower than", n,". The file is copied to the",outdir,"directory without subsampling", sep=" ")
        print(msg)
      }else{ # enough seq => resample
        options(scipen=100)
        output_fasta_p <- gsub(".gz", "", output_fasta_p)
        vsearch_cmd <- paste("vsearch --fastx_subsample ", input_fasta_p, " --fastaout ", output_fasta_p, " --sample_size ", n, " --randseed ", randseed, sep="")
        print(vsearch_cmd)
        system(vsearch_cmd)
        options(scipen=0)
        # compress the output file if the input was gz
        if(endsWith(input_fasta_p, ".gz")){
          # Specify the path for the gzipped output file
          outfile_gz <- paste(output_fasta_p, ".gz", sep="")
          # Open the existing uncompressed file for reading
          file_content <- readBin(output_fasta_p, "raw", file.info(output_fasta_p)$size)
          # Create a gzipped copy of the file
          gz <- gzfile(outfile_gz, "wb")
          writeBin(file_content, gz)
          close(gz)
          file.remove(output_fasta_p)
        }
      }
    }  # uncompressed or gz
  }# all files
}


