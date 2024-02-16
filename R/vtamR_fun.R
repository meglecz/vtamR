#' Check directory
#' 
#' Create dir if does not exists.
#' Add slash to the end of the directory name.
#' 
#' @param dir directory name
#' @export
check_dir <- function(dir){
  
  if(dir == ""){# present dir => do not add /
    return(dir)
  }else{
    if(!(endsWith(dir, "/"))){
      dir <- paste(dir, '/', sep="")
    }
    if(!dir.exists(dir)){
      dir.create(dir, recursive =TRUE)
    }
    return(dir)
  }
}


#' Get read, variant, sample and replicate counts.
#' Complete the stat_df with the above statistics.
#' 
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param stat_df data frame with the following variables: parameters, asv_count, read_count, read_count, sample_count, sample_replicate_count
#' @param stage name of the filtering step; it is used for labeling the line in the stat_df
#' @param params parameter value (or their concatenation if more than one) used for the filtering step
#' @export
#' 
get_stat <- function(read_count_df, stat_df, stage, params=NA){
  #define a temporary data frame
  df <- data.frame(parameters=character(),
                   asv_count=integer(),
                   read_count=integer(),
                   sample_count=integer(),
                   sample_replicate_count=integer())
  # get 4 different counts and place then into a data fram
  sample_repl <- paste(read_count_df$sample, read_count_df$replicate, sep="-")
  df[1,"parameters"] <-params
  df[1,"asv_count"] <-length(unique(read_count_df$asv))
  df[1,"read_count"] <-  sum(read_count_df$read_count)
  df[1,"sample_count"] <-length(unique(read_count_df$sample))
  df[1,"sample_replicate_count"] <-length(unique(sample_repl))
  # add rowname
  rownames(df) <- c(stage)
  # add new data to stat_df
  rbind(stat_df, df)
}

#' Merge
#' 
#' Merge forward are reverse fastq read pairs to fasta.
#' Output fasta files can be compressed if compress option is used.
#' The output fastainfo.csv file is similar to the fastqinfo file, but the fastq columns are replaced by a fasta column with the name of the output files
#' Returns data frame corresponding to the output csv file
#'  
#' @param fastqinfo_df data frame with column: tag_fw,primer_fw,tag_rv,primer_rv,sample,sample_type(mock/negative/real),habitat(optional),replicate,fastq_fw,fastq_rv  
#' @param fastqdir directory with input fastq files (listed in fastqinfo_df$fastq_fw and fastqinfo_df$fastq_rv)
#' @param vsearch_path path to vsearch executables
#' @param outdir output directory
#' @param fastq_ascii [33/64] ASCII character number used as the basis for the FASTQ quality score; default: 33 
#' @param fastq_maxdiffs maximum number of non-matching nucleotides allowed in the overlap region; default:10
#' @param fastq_maxee discard sequences with more than the specified number of expected errors; default:
#' @param fastq_minlen discard sequences with less than fastq_minlen bases; default: 50
#' @param fastq_maxlen discard sequences with more than fastq_maxlen bases; default: 500
#' @param fastq_minmergelen minimum length of the merged sequence; default: 50
#' @param fastq_maxmergelen maximum length of the merged sequence; default: 1000
#' @param fastq_maxns discard sequences with more than fastq_maxns of N’s; default: 0
#' @param fastq_truncqual truncate sequences starting from the first base with fastq_truncqual base quality score value or lower; default: 10
#' @param fastq_minovlen minimum overlap between the merged reads; default: 50
#' @param fastq_allowmergestagger [T/F] allow to merge staggered read pairs. Staggered pairs are pairs where the 3’ end of the reverse read has an overhang to the left of the 5’ end of the forward read; default: F
#' @param sep separator in input and output csv files; default: ","
#' @param compress [T/F]; Compress output file to gzip format; Default=F
#' @export
#'

Merge <- function(fastqinfo_df, fastqdir, vsearch_path="", outdir="", fastq_ascii=33, fastq_maxdiffs=10, fastq_maxee=1, fastq_minlen=50, fastq_maxlen=500, fastq_minmergelen=50, fastq_maxmergelen=1000, fastq_maxns=0, fastq_truncqual=10, fastq_minovlen=50, fastq_allowmergestagger=F, sep=",", compress=F){
  
  
  vsearch_path<- check_dir(vsearch_path)
  outdir<- check_dir(outdir)
  # get unique list of fastq file pairs
  tmp <- fastqinfo_df %>%
    select(fastq_fw, fastq_rv)
  tmp <- unique(tmp)
  tmp$fasta <- NA
  
  for(i in 1:nrow(tmp)){# for each file pairs
    
    # use the name of the fw fastq file and replace extension by fasta (uncompressed)
    outfile <- sub("\\..*", ".fasta", tmp[i,1])
    tmp$fasta[i] <- outfile
    outfile <- paste(outdir, outfile, sep="")
    # add path to input filenames
    fw_fastq <- paste(fastqdir, tmp[i,1], sep="")
    rv_fastq <- paste(fastqdir, tmp[i,2], sep="")
    
    if(!is_linux() && endsWith(fw_fastq, ".gz")){ #Decompress input files, since they are cannot be treated directly by vsearch on the OS
      fw_fastq <- decompress_file(fw_fastq, remove_input=F)
      rv_fastq <- decompress_file(rv_fastq, remove_input=F)
    }
    
    if(fw_fastq == outfile){ # stop the run if input and output files have the same name
      stop("ERROR: Input and output directories for fastq and fasta files are indentical. Please, give a different output directory")
    }
    
    # vsearch can accept gz files in linux, but not on windows, the outfile is always uncompressed
    vsearch <- paste(vsearch_path, "vsearch --fastq_mergepairs ", fw_fastq, " --reverse ", rv_fastq ," --fastaout ",outfile," --quiet --fastq_ascii ",fastq_ascii," --fastq_maxdiffs ", fastq_maxdiffs, " --fastq_maxee ", fastq_maxee, " --fastq_minlen ", fastq_minlen, " --fastq_maxlen ",fastq_maxlen, " --fastq_minmergelen ",fastq_minmergelen," --fastq_maxmergelen ",fastq_maxmergelen," --fastq_maxns ", fastq_maxns, " --fastq_truncqual ", fastq_truncqual, " --fastq_minovlen ", fastq_minovlen, sep="")
    if(fastq_allowmergestagger){ # if reads are longer than the amplicon
      paste(vsearch, " --fastq_allowmergestagger", sep="")
    }
    print(vsearch)
    system(vsearch)
    
    # vsearch produces uncompressed files even if input is compressed => compress output file
    if(compress){ 
      out <- compress_file(filename=outfile, remove_input=T)
      if(!endsWith(tmp$fasta[i], ".gz")){ # correct output filename in fastainfo if necessary
        tmp$fasta[i] <- paste(tmp$fasta[i], ".gz", sep="")
      }
    }
    
    original_fw_fastq <- paste(fastqdir, tmp[i,1], sep="")
    if( original_fw_fastq != fw_fastq){# the input fastq has been unzipped for vsearch => rm unzipped file to free space
      file.remove(fw_fastq)
      file.remove(rv_fastq)
    }
    
    if(F){ # compress == "zip" Keep this part as a backup, but finally it should be deleted
      # when zipping a file using a full path, the zip file will contain the embedded directory structure
      # To avoid this, change the wd to the output dir, zip the file and change back to the orignal wd
      backup_dir <- getwd()
      setwd(outdir)
      fasta <- tmp$fasta[i]
      fasta_zip <- paste(fasta, ".zip", sep="")
      zip(fasta_zip, fasta)
      file.remove(fasta)
      setwd(backup_dir)
      tmp$fasta[i] <- fasta_zip
    }
  }
  # make fastainfo file
  fastainfo_df <- left_join(fastqinfo_df, tmp, by=c("fastq_fw", "fastq_rv")) %>%
    select(-fastq_fw, -fastq_rv)
  write.table(fastainfo_df, file = paste(outdir, "fastainfo.csv", sep=""),  row.names = F, sep=sep)
  
  return(fastainfo_df)
  
}

#' is_linux
#' 
#' Returns TRUE if operating system is linux or macOS, FALSE otherwise
#'  
#' 
#' @export 
#
is_linux <- function(){
  
  system_info <- Sys.info()
  
  # Check the operating system
  if (startsWith(system_info["sysname"], "Windows")) {
    return(FALSE)
  } else if (startsWith(system_info["sysname"], "Darwin")) {
    return(TRUE)
  } else if (startsWith(system_info["sysname"], "Linux")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' decompress_file
#' 
#' Decompress the input gzipped file. Returns the output filename
#' 
#' @param filename name of the gzipped file, including path
#' @param remove_input [T/F] If TRUE the input file is deleted
#' @export 
#
decompress_file <- function(filename="", remove_input=F){
  
  # make output filename
  outfile <- gsub(".gz", "", filename)
  
  if(outfile == filename){
    stop("The input file must have .gz extention")
  }
  # read compressed file
  compressed_con <- gzfile(filename, "rb")
  text_content <- readLines(compressed_con)
  close(compressed_con)
  # writ uncompressed file
  writeLines(text_content, outfile)
  if(remove_input){
    file.remove(filename)
  }
  return(outfile)
}
#' compress_file
#' 
#' Compress input file to gzip format. Returns the output (compressed) filename
#' 
#' @param filename name of the gzipped file, including path
#' @param remove_input [T/F] If TRUE the input file is deleted
#' @export 
#
compress_file <- function(filename="", remove_input=F){
  
  # Specify the path for the gzipped output file
  outfile_gz <- paste(filename, ".gz", sep="")
  # Open the existing uncompressed file for reading
  file_content <- readBin(filename, "raw", file.info(filename)$size)
  # Create a gzipped copy of the file
  gz <- gzfile(outfile_gz, "wb")
  writeBin(file_content, gz)
  close(gz)
  
  if(remove_input){
    file.remove(filename)
  }
  return(outfile_gz)
}

#' RandomSeq
#' 
#' Random select n sequences from each input fasta file. The output is the same compression type (if any) as the input
#'  
#' @param fastainfo_df data frame with a 'fasta' column containing input file names; files can be compressed in gz and zip format
#' @param n integer; the number of randomly selected sequences 
#' @param fasta_dir directory that contains the input fasta files
#' @param outdir directory for the output files
#' @param vsearch_path path to vsearch
#' @param randseed positive integer; seed for random sampling; 0 by default means to use a pseudo-random seed, a given non-zero seed produce always the same result
#' @param compress [T/F]; If TRUE, output files are compressed in the same format as the input. Otherwise the output is uncompressed;
#' @export
#' 
RandomSeq <- function(fastainfo_df, fasta_dir="", outdir="", vsearch_path="", n, randseed=0, compress=F){
  # quite fast for uncompressed and gz files
  fasta_dir<- check_dir(fasta_dir)
  vsearch_path<- check_dir(vsearch_path)
  outdir<- check_dir(outdir)
  
  unique_fasta <- unique(fastainfo_df$fasta)
  
  for(i in 1:length(unique_fasta)){ # go through all fasta files
    input_fasta <- unique_fasta[i]
    
    if(endsWith(input_fasta, ".zip")){
      stop("Zip files are not supported")
    }
    
    input_fasta_p <- paste(fasta_dir, input_fasta, sep="")
    output_fasta_p <- paste(outdir, input_fasta, sep="")
    
    if(!is_linux() && endsWith(input_fasta_p, ".gz")){ #Decompress input files, since they are cannot be treated directly by vsearch on the OS
      input_fasta_p <- decompress_file(input_fasta_p, remove_input=F)
    }
    
    seq_n <- count_seq(filename=input_fasta_p)
    if(n > seq_n ){ # not enough seq
      file.copy(input_fasta_p, output_fasta_p)
      msg <- paste("WARNING:", input_fasta_p,"has",seq_n,"sequences, which is lower than", n,". The file is copied to the",outdir,"directory without subsampling", sep=" ")
      print(msg)
    }else{ # enough seq => resample
      options(scipen=100) # do not transform large numbers to scentific forms, since it would lead to an error in vsearch
      output_fasta_p <- gsub(".gz", "", output_fasta_p)
      vsearch_cmd <- paste(vsearch_path, "vsearch --fastx_subsample ", input_fasta_p, " --fastaout ", output_fasta_p, " --sample_size ", n, " --randseed ", randseed, sep="")
      print(vsearch_cmd)
      system(vsearch_cmd)
      options(scipen=0)
    }
    
    if(compress){ # compress the output file 
      outfile_gz <- compress_file(filename=output_fasta_p, remove_input=T)
    }
    
    if(!is_linux() && endsWith(input_fasta, ".gz")){ # delete decompressed input files
      file.remove(input_fasta_p)
    }
    
  }# all files
}

#' count_seq
#' 
#' Count the number of sequences in the input fasta file. Input can be uncompressed of gz file, but zip file are not handled
#'  
#' @param file fasta file including path
#' @export
count_seq <- function(filename=filename){
  # work on uncompressed and gz files
  if(endsWith(filename, "gz")){
    file_connection <- gzfile(filename, "rb")
  }else{
    file_connection <- file(filename, "r")
  }
  data <- as.data.frame(readLines(file_connection, n = -1))
  close(file_connection)
  
  colnames(data) <- c("read")
  data <- data %>%
    filter(startsWith(read, ">"))
  
  count = nrow(data)
  rm(data)
  
  return(count)
}


#' SortReads
#' 
#' Demultiplex each input fasta file using the tag combinations at the extremities of the merged reads.
#' Trim primers from demultiplexed reads.
#' Output fasta files can be compressed if compress option is used.
#' The output fileinfo.csv file is similar to the fastainfo file, but the but do not have tag and primer columns.
#' Returns data frame of corresponding to the output csv file.
#'  
#' @param fastainfo_df data frame with column: tag_fw,primer_fw,tag_rv,primer_rv,sample,sample_type(mock/negative/real),habitat(optional),replicate,fasta  
#' @param fastadir directory with input fasta files (listed in fastainfo_df$fasta)
#' @param vsearch_path path to vsearch executables
#' @param cutadapt_path path to cutadapt executables
#' @param outdir output directory
#' @param check_reverse [T/F] if TRUE, ckeck the reverse comlementary sequences of the input fasta as well; default: F
#' @param tag_to_end tags are at the extremity of the reads (starting at the first base); default: T
#' @param primer_to_end primers follow directly the tags (no heterogeneity spacer); default: T
#' @param cutadapt_error_rate maximum proportion of errors between primers and reads (for tags, exact match is required); default: 0.1
#' @param cutadapt_minimum_length minimum length of the trimmed sequence; default: 50
#' @param cutadapt_maximum_length maximum length of the merged sequence; default: 500
#' @param sep separator in input and output csv files; default: ","
#' @param compress [T/F]; compress output to gzip format; Deault=F
#' @export

SortReads <- function(fastainfo_df, fastadir, outdir="", cutadapt_path="" ,vsearch_path="", check_reverse=F, tag_to_end=T, primer_to_end=T, cutadapt_error_rate=0.1,cutadapt_minimum_length=50,cutadapt_maximum_length=500, sep=",",  compress=F){
  #########
  # SortReads_no_reverse does the whole demultilexing, trimming and compress on the + strand
  # If sequences are not oriented, the -strand should be checked => 
  # run SortReads_no_reverse of plus strand and on - strand after swapping fw and rev tags and primers,
  # take the reverse complement of the -strand results (vsearch)
  # pool the results of the 2 strands
  # compress if necessary
  
  # run on strand +
  if(check_reverse){
    #### use +strand, output to sorted_dir, uncompressed
    sortedinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df, fastadir=fastadir, outdir=outdir, cutadapt_path=cutadapt_path, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=F)
    
    #### use - strand
    # swap fw and rv tags and primers
    fastainfo_df_tmp <- fastainfo_df %>%
      select(tag_fw_tmp = tag_rv, tag_rv_tmp = tag_fw, primer_fw_tmp = primer_rv, primer_rv_tmp = primer_fw, sample, sample_type,habitat, replicate, fasta) %>%
      select(tag_fw = tag_fw_tmp, tag_rv = tag_rv_tmp, primer_fw = primer_fw_tmp, primer_rv = primer_rv_tmp, sample, sample_type,habitat, replicate, fasta)
    # make temp dir 
    outdir <- check_dir(outdir)
    rc_dir <-paste(outdir, 'rc_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    rc_dir <- check_dir(rc_dir)
    # run sortreads on for reverse strand
    sortedinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df_tmp, fastadir=fastadir, outdir=rc_dir, cutadapt_path=cutadapt_path, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=F)
    
    ### reverse complement and pool
    # get list of files demultiplexed on - strand
    files <- list.files(path = rc_dir, pattern=".fasta")
    # Filter the files based on the motif using regular expressions
    # reverse complement sequences on the minus stand, and append info to the plus strand output
    files <- grep(pattern = "\\.fasta", x = files, value = TRUE)
    for(i in 1:length(files)){
      plus <- paste(outdir, files[i], sep="")
      minus <- paste(rc_dir, files[i], sep="")
      minus_rc <- paste(rc_dir, "rc_", files[i], sep="")
      # reverse complement sequences in minus file
      rev_comp_cmd <- paste(vsearch_path, "vsearch --fastx_revcomp ", minus, " --fastaout ", minus_rc, " --quiet", sep="")
      print(rev_comp_cmd)
      system(rev_comp_cmd)
      # append content of minus_rc to plus file
      file.append(plus, minus_rc)
    }
    
    # delete temporary reverse_comp dir
    unlink(rc_dir, recursive = TRUE)
    
    ### compress
    if(compress){
      
      for(i in 1:nrow(sortedinfo_df)){
        
        file <- sortedinfo_df$filename[i]
        sortedinfo_df$filename[i] <- paste(file, ".gz", sep="") # correct output filename
        file <- paste(outdir, file, sep="") # add path
        file_gz <- compress_file(file, remove_input=T) # compress file
      }
      write.table(sortedinfo_df, file = paste(outdir, "sortedinfo.csv", sep=""),  row.names = F, sep=sep) 
    }
  }
  else{
    # check only + strand
    sortedinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df, fastadir=fastadir, outdir=outdir, cutadapt_path=cutadapt_path, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)
  }
  return(sortedinfo_df)
  
}

#' SortReads_no_reverse
#' 
#' Same as SortReads, but do not check the reverse complement of the sequences.
#' Demultiplex each input fasta file using the tag combinations at the extremities of the merged reads.
#' Trim primers from demultiplexed reads.
#' Output fasta files can be compressed if compress option is used.
#' The output fileinfo.csv file is similar to the fastainfo file, but the but do not have tag and primer columns.
#' Returns data frame of corresponding to the output csv file.
#'  
#' @param fastainfo_df data frame with column: tag_fw,primer_fw,tag_rv,primer_rv,sample,sample_type(mock/negative/real),habitat(optional),replicate,fasta  
#' @param fastadir directory with input fasta files (listed in fastainfo_df$fasta)
#' @param cutadapt_path path to cutadapt executables
#' @param outdir output directory
#' @param tag_to_end tags are at the extremity of the reads (starting at the first base); default: T
#' @param primer_to_end primers follow directly the tags (no heterogeneity spacer); default: T
#' @param cutadapt_error_rate maximum proportion of errors between primers and reads (for tags, exact match is required); default: 0.1
#' @param cutadapt_minimum_length minimum length of the trimmed sequence; default: 50
#' @param cutadapt_maximum_length maximum length of the merged sequence; default: 500
#' @param sep separator in input and output csv files; default: ","
#' @param compress [T/F]; compress output to gzip files.
#' @export
SortReads_no_reverse <- function(fastainfo_df, fastadir, outdir="", cutadapt_path="", tag_to_end=T, primer_to_end=T, cutadapt_error_rate=0.1,cutadapt_minimum_length=50,cutadapt_maximum_length=500, sep=",",  compress=0){
  # do the complete job of demultiplexing and trimming of input file without checking the reverse sequences
  
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
  
  for(i in 1:length(fastas)){ # for each input fasta
    # make temp dir for tag file and tag-trimmed files
    tmp_dir <-paste(outdir, 'tmp_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
    tmp_dir <- check_dir(tmp_dir)
    
    # select lines in fastainfo_df that corresponds to a given input fasta file
    fasta_file <- fastas[i]
    df <- fastainfo_df %>%
      filter(fasta==fasta_file)
    
    # make a tags.fasta file with all tag combinations of the fasta to be demultiplexed
    tag_file <- make_adapter_fasta(fastainfo_df, fasta_file=fasta_file, tag_to_end=tag_to_end, outdir=tmp_dir)
    # add path
    fasta_file <- paste(fastadir, fasta_file, sep="")
    # demultiplex fasta, write output to tmp file
    demultiplex_cmd = paste(cutadapt_path, "cutadapt --cores=0 -e 0 --no-indels --trimmed-only -g file:", tag_file," -o ", tmp_dir, "tagtrimmed-{name}.fasta ", fasta_file, sep="")
    print(demultiplex_cmd)
    system(demultiplex_cmd)
    
      # for a given marker, there is only one primer combination
      primer_fwl <- df[1,"primer_fw"]
      primer_rvl <- df[1,"primer_rv"]
      primer_rvl_rc <- reverse_complement(primer_rvl)
      
      for(f in 1:nrow(df)){# go through each de-multiplexed, tag-trimmed file and trim primers
        outfilename <- paste(df[f,"sample"], df[f,"replicate"], sep="-")
        outfilename <- paste(outfilename, ".fasta", sep="")
        if(compress){
          outfilename <- paste(outfilename, ".gz", sep="")
        }
        # complete fastainfo_df with output fasta name
        fastainfo_df$filename[which(fastainfo_df$sample==df[f,"sample"] & fastainfo_df$replicate==df[f,"replicate"])] <- outfilename
        # add path to output file
        primer_trimmed_file <- paste(outdir, outfilename, sep="")
        tag_trimmed_file <- paste(tmp_dir, "tagtrimmed-", df[f,"tag_fw"], "-", df[f,"tag_rv"], ".fasta", sep="")
        if(primer_to_end){
          primer_trim_cmd <- paste(cutadapt_path, "cutadapt --cores=0 -e ",cutadapt_error_rate ," --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length ," --maximum-length ", cutadapt_maximum_length, " -g ^", primer_fwl, "...", primer_rvl_rc, "$ --output ", primer_trimmed_file, " ", tag_trimmed_file, sep="")
        }
        else{
          primer_trim_cmd <- paste(cutadapt_path, "cutadapt --cores=0 -e ",cutadapt_error_rate ," --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length ," --maximum-length ", cutadapt_maximum_length, ' -g "', primer_fwl, ';min_overlap=',nchar(primer_fwl),'...', primer_rvl_rc,  ';min_overlap=',nchar(primer_rvl_rc),'" --output ', primer_trimmed_file, " ", tag_trimmed_file, sep="")
        }
        print(primer_trim_cmd)
        system(primer_trim_cmd)
      } # end tag-trimmed 
    # delete the tmp dir wit the tag-trimmed files
    unlink(tmp_dir, recursive = TRUE)
  }# end fasta
  
  # make sortedinfo file
  fastainfo_df <- fastainfo_df %>%
    select(-tag_fw, -primer_fw, -tag_rv, -primer_rv, -fasta)
  write.table(fastainfo_df, file = paste(outdir, "sortedinfo.csv", sep=""),  row.names = F, sep=sep) 
  return(fastainfo_df)
}

#' make_adapter_fasta
#' 
#' Make a fasta file with tag combinations
#' 
#' @param fastainfo_df data frame with column: tag_fw,tag_rv,fasta  
#' @param fasta_file name of the fasta file for wich tag-combination should be written to a tags.fasta file
#' @param tag_to_end [T/F] if T use anchored search (tags are at the extremity of the sequence)
#' @param outdir output directory
#' @export
#'
make_adapter_fasta <- function(fastainfo_df, fasta_file, tag_to_end, outdir){
  # select tag combinations pour the fasta file
  tags <- fastainfo_df %>%
    filter(fasta==fasta_file) %>%
    select(tag_fw, tag_rv)
  # make unique tag combinations and add necessary columns
  tags <- unique(tags)
  tags$tag_fw <- toupper(tags$tag_fw)
  tags$tag_rv <- toupper(tags$tag_rv)
  tags$tag_rv_rc <- lapply(tags$tag_rv, reverse_complement)
  tags$tag_fwl <- lapply(tags$tag_fw, nchar)
  tags$tag_rvl <- lapply(tags$tag_rv, nchar)
  
  # Specify the file path
  tag_file <- paste(outdir, "tags.fasta", sep="")
  # initialise te content of the tag_file
  text <- c()
  if(tag_to_end){
    #>tag_fw-tag_rv
    #^tcgatcacgatgt...gctgtagatcgaca$
    for(j in 1:nrow(tags)){
      title <- paste(">", tags[j,"tag_fw"], "-",   tags[j,"tag_rv"], sep="")
      seq <- paste("^", tags[j,"tag_fw"], "...",   tags[j,"tag_rv_rc"], "$",sep="")
      text <- c(text, c(title,seq))
    }
  }else{
    #>tag_fw-tag_rv
    #tcgatcacgatgt;min_overlap=13...gctgtagatcgaca;min_overlap=14
    for(j in 1:nrow(tags)){
      title <- paste(">", tags[j,"tag_fw"], "-",   tags[j,"tag_rv"], sep="")
      seq <- paste(tags[j,"tag_fw"], ";min_overlap=", tags[j,"tag_fwl"], "...",  tags[j,"tag_rv_rc"], ";min_overlap=", tags[j,"tag_rvl"], sep="")
      text <- c(text, c(title,seq))
    }
  }
  # write file
  writeLines(text, tag_file)
  
  return(tag_file)
}


#' reverse_complement
#' 
#' Reverse complement a DNA sequence. IUPAC codes are accepted
#' 
#' @param sequence DNA sequence
#' @export
#'

reverse_complement <- function(sequence){
  # define complementary bases
  comp <- data.frame(orig=c("A","T","C","G","R","Y","W","S","M","K","B","H","D","V","N","a","t","c","g","r","y","w","s","m","k","b","h","d","v","n"),
                     complement=c("T","A","G","C","Y","R","W","S","K","M","V","D","H","B","N","t","a","g","c","y","r","w","s","k","m","v","d","h","b","n")
  )
  # revers, split sequneces and make dataframe
  sequence_df <- data.frame(reverse=rev(strsplit(sequence, NULL)[[1]]))
  # add complimetary nt
  sequence_df <- left_join(sequence_df, comp, by=c("reverse"="orig"))
  if(any(is.na(sequence_df$complement))){
    print(sequence)
    stop("ERROR: Sequence contains non-IUPAC character")
  }
  # collaps vector to string
  reverse_comp <- paste(sequence_df$complement, collapse = "")
  
  return(reverse_comp)
}

#' Read all fasta files in fileinfo file to a data frame
#' 
#' Read all fasta file in sortedinfo_df data frame (sample	replicate	file name).
#' Dereplicate reads to ASVs.
#' Count the number of reads of each ASVs in each sample-replicate.
#' Returns a df ("asv", "sample", "replicate", "read_count")
#' 
#' @param sortedinfo_df data frame with columns:  sample,  replicate, file name, (optional: sample_type(mock/negative/real), habitat)
#' @param dir name of the directory with fasta files 
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param sep separator in csv files; default: ","
#' @param asv_list name of the file, containing asvs and asv_ids from earlier analyses. Optional. It is used to homogenize asv_ids between different data sets
#' @param updated_asv_list name of the output file, containing the updated ASVs. Optional.
#' @export
#' 
read_fastas_from_fileinfo <- function (sortedinfo_df, dir="", outfile="", sep=",", asv_list="", updated_asv_list="") {
  # read all fasta files in fileinfo to a read_count_df
  if(nchar(dir)>0){
    dir <- check_dir(dir)
  }
  # define empty read_count_df to pool the results of variables
  read_count_df <- data.frame(asv=character(),
                              read_count=integer(),
                              sample=character(),
                              replicate=character())
  # read all fasta files in fileinfo and count the reads
  for(i in 1:length(sortedinfo_df$filename)){
    fas <- paste(dir, sortedinfo_df$filename[i], sep = "")
    print(fas)
    if(file.size(fas) < 50){ #' an empty file gzipped has 42 size => skip these files, An unzipped fasta with 80 nt is bigger than 50
      next
    }
    
    read_count_df_tmp <- read_fasta_seq(filename=fas, dereplicate=T) # returns data frame with asv and read_cunt columns
    read_count_df_tmp$sample <- sortedinfo_df[i,"sample"]
    read_count_df_tmp$replicate <- sortedinfo_df[i,"replicate"]
    read_count_df <- rbind(read_count_df, read_count_df_tmp)
  }
  rm(read_count_df_tmp)
  # reorder columns
  read_count_df <- read_count_df[, c("asv", "sample", "replicate", "read_count")]
  
  # add asv_id, but do not add new asv s to asv_list; it is probably too early at this stage, better to do it after swarm
  read_count_df <- add_ids(read_count_df, asv_list=asv_list, sep=sep,  updated_asv_list=updated_asv_list)
  
  # write read_count table
  if(outfile != ""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' add_ids
#' 
#' Add asv_ids to a data frame that have an asv column. Can take into account already existing asv and asv_id (from earlier dataset)
#' If update_asv_list, the input asv_list id completed by new asvs and asv_ids
#' Returns a data frame with as asv_id column on top of the input.
#' 
#' @param read_count_df data frame with columns:  asv, sample,  replicate, read_count
#' @param asv_list name of the file, containing asvs and asv_ids from earlier analyses. Optional. It is used to homogenize asv_ids between different data sets
#' @param updated_asv_list name of the output file, containing the updated ASVs. Optional.
#' @param sep separator in csv files; default: ","
#' @export
#' 
add_ids <- function(read_count_df, asv_list=asv_list, updated_asv_list="", sep=","){
  
  if(asv_list != ""){  # read already existing asvs, if the file is given
    asv_df <- read.csv(asv_list, sep=sep, header=TRUE)
    check_one_to_one_relationship(asv_df)
  }else{
    asv_df <- data.frame("asv_id"=integer(),
                         "asv"=as.character()
                         )
  }
  # list of unique asvs
  asv_uniq <- unique(read_count_df$asv)
  # list of unique asvs, not in the asv_df
  new_asvs <- asv_uniq[!asv_uniq %in% asv_df$asv]
  max_id <- max(asv_df$asv_id)
  new_ids <- seq(from =max_id+1, to = (max_id + length(new_asvs)), by=1)
  new_asvs_df <- data.frame(
    "asv_id" = new_ids, "asv"=new_asvs)
  # add new asvs to asv_df
  asv_df <- rbind(asv_df, new_asvs_df)
  # add asv_id to read_count_df
  read_count_df <- left_join(read_count_df, asv_df, by="asv") %>%
    select(asv_id, sample, replicate, read_count, asv)
  
  # if asv_list should be updated, write it to a new file
  if(updated_asv_list != ""){
#    str <- paste("_", trunc(as.numeric(trunc(Sys.time()))), ".", sep="")
#    new_file <- sub("\\.", str, asv_list)
    write.table(asv_df, file=updated_asv_list, row.names = FALSE, sep=sep)
  }
  
  return(read_count_df)
}

#' read_fasta_seq
#' 
#' Read sequences from a fasta file. Fasta can be gz compressed or uncompressed
#' Returns a either a data frame with read read in a line, or a data frame with asv and read counts (if dereplicate==T) 
#' 
#' @param filename name of the input file including full path
#' @param dereplicate [T/F] if T, return asvs with read counts
#' @export 
#' 

read_fasta_seq <- function(filename=filename, dereplicate=F){
  # can deal with sequences in multiple lines
  # only R-base
  # quicker than read.fasta from seqinr
  if(endsWith(filename, ".gz")){
    file_connection <- gzfile(filename, "rb")
  }else{
    file_connection <- file(filename, "r")
  }
  data <- readLines(file_connection, n = -1)
  close(file_connection)
  
  data <- gsub(" ", "_", data)
  data <- gsub(">[^ ]+", ">", data, fixed=F, perl=T)
  data <- do.call(paste, c(as.list(data), sep = ""))
  data <- as.data.frame(strsplit(data, ">"))
  colnames(data) <- c("read")
  data <- data %>%
    filter(!(read==""))
  
  data$read <-toupper(data$read)
  
  if(dereplicate){
    data <- data %>%
      group_by(read) %>%
      summarize(read_count = length(read)) %>%
      select(asv=read, read_count)
  }
  return(data)
}


#' swarm
#' 
#' Run swarm (https://github.com/torognes/swarm) on input read_count_df data frame, pool variants of the same cluster sum reads of the underlying ASVs
#' Return a data frame with the same structure as the input
#' Swarm can be run sample by sample (by_sample=T) or for the whole data set in ine go
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param outfile name of the output file with asv_id, sample, replicate, read_count and asv columns; Optional; If empty the results are not written to a file
#' @param swarm_path Path to th swarm executable (Default: TRUE)
#' @param by_sample [T/F], if TRUE, swarm run separately for each sample
#' @param num_threads Number of CPUs
#' @param swarm_d positive integer, d parameter for swarm (1 by default); maximum number of differences allowed between two amplicons, meaning that two amplicons will be grouped if they have d (or less) differences.
#' @param fastidious [T/F] when working with d = 1, perform a second clustering pass to reduce the number of small clusters (Default: TRUE)
#' @param write_csv [T/F]; write read_counts to csv file; default=FALSE
#' @param sep separator for the output file
#' @export
#' 

swarm <- function(read_count_df, outfile="", swarm_path="", num_threads=1, swarm_d=1, fastidious=T, write_csv=F, sep=",", by_sample=T){
  
  if(by_sample){
    # make an empty output df, with the same columns and variable types as read_count_df
    out_df <- read_count_df %>%
      filter(asv=="")
    
    # get list of samples 
    sample_list <- unique(read_count_df$sample)
    
    # run swarm for each sample
    for(s in sample_list){
      print(s)
      # select occurrences for sample
      df_sample <- read_count_df %>%
        filter(sample==s)
      # run swarm
      df_sample <- run_swarm(df_sample, swarm_path=swarm_path, num_threads=num_threads, swarm_d=swarm_d, fastidious=fastidious)
      # add output of the sample to the total data frame
      out_df <- rbind(out_df, df_sample)
    }
  }else{ # run swarm for all samples together
    out_df <- run_swarm(read_count_df, swarm_path=swarm_path, num_threads=num_threads, swarm_d=swarm_d, fastidious=fastidious)
  }
  
  if(outfile != ""){
    write.table(out_df, file = outfile,  row.names = F, sep=sep)
  }
  return(out_df)
  
}

#' run_swarm
#' 
#' Run swarm (https://github.com/torognes/swarm) on input read_count_df data frame, pool variants of the same cluster sum reads of the underlying ASVs
#' Return a data frame with the same structure as the input
#' 
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param swarm_path Path to th swarm executable (Default: TRUE)
#' @param num_threads Number of CPUs
#' @param swarm_d positive integer, d parameter for swarm (1 by default); maximum number of differences allowed between two amplicons, meaning that two amplicons will be grouped if they have d (or less) differences.
#' @param fastidious [T/F] when working with d = 1, perform a second clustering pass to reduce the number of small clusters (Default: TRUE)
#' @export
#' 

run_swarm <- function(read_count_df, swarm_path="", num_threads=1, swarm_d=1, fastidious=T){
  
  tmp_dir <-paste('tmp_swarm_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  tmp_dir <- check_dir(tmp_dir)
  swarm_path <- check_dir(swarm_path)
  
  ### make df with unique asv and read_count
  df_unique <- read_count_df %>%
    group_by(asv, asv_id) %>%
    summarize(sum_read_count = sum(read_count), .groups="drop_last") %>%
    ungroup() 

  ### make a fasta with dereplicated sequences  
  input_swarm <- paste(tmp_dir, "swarm_input.fasta", sep="")
  writeLines(paste(">", df_unique$asv_id, "_", df_unique$sum_read_count, "\n", df_unique$asv, sep="" ), input_swarm)
  
  df_unique <- df_unique %>%
    select(-sum_read_count)
  ### run swarm
  #  representatives <- paste(tmp_dir, "representatives.fasta", sep="")
  clusters <- paste(tmp_dir, "clusters.txt", sep="")
  swarm <- paste(swarm_path, "swarm -d ",swarm_d," -t ", num_threads, " -o ", clusters, sep="")
  if(fastidious){
    swarm <- paste(swarm, "-f", sep=" ")
  }
  swarm <- paste(swarm, input_swarm, sep=" ")
  print(swarm)
  system(swarm)
  
  ###
  # pool clusters in read_count_df
  ###
  # make a data frame with representative and clustered columns, where clustered has all swarm input sequences id, and  representative is the name of the cluster they belong to
  cluster_df <- read.table(clusters, fill =TRUE, strip.white=TRUE, header = FALSE)
  cluster_df <- data.frame(representative = rep(cluster_df$V1, each = ncol(cluster_df)),
                           clustered = as.vector(t(cluster_df[,])))
  # delete line with no values in clustered
  cluster_df <- cluster_df %>%
    filter(clustered != "")
  # delete read counts from id
  cluster_df$representative <- sub("_[0-9]+", "", cluster_df$representative )
  cluster_df$representative <- as.numeric(cluster_df$representative)
  cluster_df$clustered <- sub("_[0-9]+", "", cluster_df$clustered )
  cluster_df$clustered <- as.numeric(cluster_df$clustered)
  # add representative asv
  cluster_df <- left_join(cluster_df, df_unique, by= c("representative" = "asv_id")) %>%
    select("representative_id"=representative, representative_asv=asv, "clustered_id"=clustered)

  # free space
  remove(df_unique)
  unlink(input_swarm)
  unlink(clusters)
  unlink(tmp_dir, recursive=TRUE)
  
  # replace asv by representative sequences in read_count_df
  read_count_df <- left_join(read_count_df, cluster_df,  by= c("asv_id" = "clustered_id"))
  read_count_df <- read_count_df %>%
    select(-asv, -asv_id) %>%
    group_by(representative_asv, representative_id, sample, replicate) %>%
    summarize(read_count_cluster=sum(read_count), .groups="drop_last") %>%
    rename("asv" = representative_asv, "read_count"=read_count_cluster, "asv_id"=representative_id) %>%
    ungroup()
  
  
  return(read_count_df)
}

#' check_one_to_one_relationship
#' 
#' Check if in th input data frame there is a one to one relationship between unique asvs and unisq asv_ids. If yes, returns TRUE, otherwise quit the run
#' The same asv - asv_id combination can appear more than once in the data frame
#' 
#' @param df data frame with the following variables: asv_id, asv (can have other columns as well)
#' @export
#' 
check_one_to_one_relationship <- function(df){
  
  # make unique asv-asv_id combinations
  df <- df %>%
    select(asv_id, asv) %>%
    distinct()
  
  # check if more then one asv per asv_id
  unique_asv_id <- df %>%
    group_by(asv_id) %>%
    summarize(count= length(asv)) %>%
    filter(count>1)
  if(nrow(unique_asv_id) > 0 ){
    print(unique_asv_id)
    stop("Some of the the asv_ids belong to multile asv")
  }
  
  # check if more then one asv_id per asv
  unique_asv <- df %>%
    group_by(asv) %>%
    summarize(count= length(asv_id)) %>%
    filter(count>1)
  if(nrow(unique_asv) > 0 ){
    print(unique_asv)
    stop("Some of the asv has multiple asv_ids")
  }
  
  return(TRUE)
}

#' update_asv_list
#' 
#' Pools unique asv - asv_id combinations in the input data frame and asv - asv_id combinations in the input file
#' The input file is typically a csv file containing asv seen in earlier runs with their asv_id.
#' If there is a conflict within or between the input data quits with a error message. Otherwise write the updated asv list with their asv_id to the outfile
#' 
#' @param read_count_df data frame with columns:  asv_id, sample,  replicate, read_count, asv
#' @param asv_list name of the file, containing asvs and asv_ids from earlier analyses. Optional. It is used to homogenize asv_ids between different data sets
#' @param outfile Name of the output file; if empty, write a new file using the name asv_list, completed by the number of seconds from 01/01/1970
#' @param sep separator in csv files; default: ","
#' @export
#' 
update_asv_list <- function(read_count_df, asv_list=asv_list, outfile="", sep=","){
  
  # read earlier ASV list
  if(asv_list != ""){  # read already existing asvs, if the file is given
    asv_df <- read.csv(asv_list, sep=sep, header=TRUE)
  }else{
    asv_df <- data.frame("asv_id"=integer(),
                         "asv"=as.character()
    )
  }
  
  #  asv_df[2,1] <- 1
  if(check_one_to_one_relationship(asv_df)){
    print("One to one relationship between asv_id and asv in input asv_list")
  }
  # make a dataframe with the unique combinations of asv_id-asv
  new_df <- read_count_df %>%
    select(asv_id, asv) %>%
    distinct()
  if(check_one_to_one_relationship(new_df)){
    print("One to one relationship between asv_id and asv in input read_count_df")
  }
  # pool earlier asvs and new ones and avoid redundancy
  asv_df <- rbind(asv_df, new_df) %>%
    distinct() %>%
    arrange(asv_id)
  if(check_one_to_one_relationship(asv_df)){
    print("One to one relationship between asv_id and asv in output asv_list")
  }
  
  if(outfile == ""){
    str <- paste("_", trunc(as.numeric(trunc(Sys.time()))), ".", sep="")
    outfile <- sub("\\.", str, asv_list)
  }
  write.table(asv_df, file=outfile, row.names = FALSE, sep=sep)
  
  
}

#' LFN_global_read_count
#' 
#' Eliminate ASVs with less than cutoff reads in the dataset.
#' Returns the filtered read_count_df data frame.
#' 
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param cutoff minimum number of reads for an ASV in the whole dataset; default=10
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param sep separator in csv files; default: ","
#' @export
LFN_global_read_count <- function (read_count_df, cutoff=10, outfile="", sep=",") {
  df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count=sum(read_count)) %>%
    filter(read_count > cutoff)
  read_count_df <- filter(read_count_df, (asv %in% df$asv))
  
  if(outfile != ""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' LFN_read_count
#' 
#' Eliminate occurrences with less than cutoff reads.
#' Returns the filtered read_count_df data frame.
#' 
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param cutoff minimum number of reads for an occurrence; default=10
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param sep separator in csv files; default: ","
#' @export
LFN_read_count <- function (read_count_df, cutoff=10, outfile="", sep=",") {
  read_count_df <- filter(read_count_df,  (read_count >= cutoff))
  if(outfile != ""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' LFN_sample_replicate
#' 
#' Eliminate occurrences where the read_count/sum(read_count of the sample-replicate) is less than cutoff.
#' Returns the filtered read_count_df data frame.
#' 
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param cutoff minimum proportion for an occurrence within a sample-replicate; default=0.001
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param sep separator in csv files; default: ","
#' @export
LFN_sample_replicate <- function (read_count_df, cutoff=0.001, outfile="", sep=",") {
  
  sum_by_column_df <- read_count_df %>%
    group_by(sample,replicate) %>%
    summarize(sr_sum = sum(read_count), .groups="drop_last")
  
  read_count_df <- left_join(read_count_df, sum_by_column_df, by=c("sample","replicate")) %>%
    filter(read_count/sr_sum >= cutoff) %>%
    select(-sr_sum) %>%
    ungroup()
  
  if(outfile != ""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' LFN_variant
#' 
#' If by_replicate=F: Eliminate occurrences where the read_count/sum(read_count of the asv) is less than cutoff.
#' If by_replicate=T: Eliminate occurrences where the read_count/sum(read_count of the asv in its replicate) is less than cutoff.
#' Returns the filtered read_count_df data frame.
#' 
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param cutoff minimum proportion for an occurrence within an asv or an asv-replicate; default=0.001
#' @param by_replicate T/F; default=FALSE
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param sep separator in csv files; default: ","
#' @export
LFN_variant <- function (read_count_df, cutoff=0.001, by_replicate=FALSE, outfile="", sep=",") {
  
  if(by_replicate){
    sum_by_asv <- read_count_df %>%
      group_by(asv,replicate) %>%
      summarize(asv_sum = sum(read_count), .groups="drop_last")
    read_count_df <- left_join(read_count_df, sum_by_asv, by=c("asv", "replicate")) %>%
      filter(read_count/asv_sum >= cutoff) %>%
      select(-asv_sum)
    
  }else{
    sum_by_asv <- read_count_df %>%
      group_by(asv) %>%
      summarize(asv_sum = sum(read_count))
    read_count_df <- left_join(read_count_df, sum_by_asv, by="asv") %>%
      filter(read_count/asv_sum >= cutoff) %>%
      select(-asv_sum)
  }
  
  if(outfile != ""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' pool_LFN
#' 
#' pool all count_read_df data frames, and keep only occurrences present in all filters
#' 
#' @param ... a list of data frames
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param sep separator in csv files; default: ","
#' @export
#' 
pool_LFN <- function (... , outfile="", sep=",") {
  df_list <- list(...)
  merged <-  df_list[[1]]
  for(i in 2:length(df_list)){
    merged <- inner_join(merged, df_list[[i]])
  }
  
  if(outfile != ""){
    write.table(merged, file = outfile,  row.names = F, sep=sep)
  }
  return(merged)
}

#' FilterMinReplicateNumber
#' 
#' Filter out all occurrences where the asv in not present in at least cutoff number of replicates.
#' Returns the filtered read_count_df data frame.
#'  
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param cutoff Minimum number of replicates; default=2
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @export
#'
FilterMinReplicateNumber <- function(read_count_df, cutoff=2, outfile="", sep=","){
  # read_count_df <- df
  # add a temporary column with asv and sample concatenated
  read_count_df$tmp <- paste(read_count_df$asv,  read_count_df$sample, sep="-")
  # make a df_tmp containing the number of replicates for each asv-sample combination
  df_tmp <- read_count_df  %>%
    group_by(tmp) %>%
    summarize(repl_number=length(tmp))  %>%
    filter(repl_number >= cutoff)
  # keep only asv-sample if present at least in min_replicate_number replicates
  read_count_df <- filter(read_count_df, (read_count_df$tmp %in% df_tmp$tmp))
  read_count_df$tmp <- NULL
  
  if(outfile !=""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' FilterIndel
#' 
#' Filter out all ASVs, if the modulo 3 of their length is not the same as that of the majority of the ASVs
#'  
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @export
#' 
FilterIndel <- function(read_count_df, outfile="", sep=","){
  # add a column with the modulo 3 of the length of the asvs
  read_count_df$mod3 <- nchar(read_count_df$asv) %% 3
  # make a tibble with modulo3 of the length of the ASVs and their count 
  # ASVs are counted as many times as they occur, so the most frequent ASV have higher weight
  # read_counts are not taken into account
  tmp <- read_count_df %>%
    group_by(mod3) %>%
    summarize(length_modulo=length(mod3)) %>%
    arrange(desc(length_modulo))
  # get the modulo 3 the most frequent
  my_modulo3 <- as.integer(tmp[1,"mod3"])
  # select only the lines with asv length compatible with the most frequent modulo3
  read_count_df <- read_count_df %>%
    filter(mod3 == my_modulo3)
  
  # delete the temporary column
  read_count_df$mod3 <- NULL
  
  if(outfile !=""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' codon_stops_from_genetic_code
#' 
#' Returns a vector of codon stops corresponding the the genetic code number in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
#'  
#' @param genetic_code genetic code number from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
#' @export
#' 
codon_stops_from_genetic_code <- function(genetic_code=5){
  if(genetic_code == 1){
    return(c("TAA","TAG","TGA"))
  }
  else if(genetic_code == 2){
    return(c("TAA","TAG","AGA", "AGG"))
  }
  else if(genetic_code == 3 || genetic_code == 4 || genetic_code == 5 || genetic_code == 9 || genetic_code == 10 || genetic_code == 13 || genetic_code == 21 || genetic_code == 24 || genetic_code == 25 || genetic_code == 31){
    return(c("TAA","TAG"))
  }
  else if(genetic_code == 6 || genetic_code == 14 || genetic_code == 33){
    return(c("TAG"))
  }
  else if(genetic_code == 16){
    return(c("TAA","TGA"))
  }
  else if(genetic_code == 11 || genetic_code == 12 || genetic_code == 26 || genetic_code == 28 ){
    return(c("TAA","TAG","TGA"))
  }
  else if(genetic_code == 22){
    return(c("TCA","TAA", "TGA"))
  }
  else if(genetic_code == 23){
    return(c("TTA","TAA", "TAG", "TGA"))
  }
  else if(genetic_code == 27 || genetic_code == 29 || genetic_code == 30){
    return( c("TGA"))
  }
  else{
    return(c())
  }
}


#' FilterCodonStop
#' 
#' Filter out all ASVs, if there is a codon stop in all three frames of the direct strand.
#' Returns the filtered read_count_df data frame.
#'  
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param genetic_code genetic code number from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
#' @export
#' 
FilterCodonStop <- function (read_count_df, outfile="", genetic_code=5, sep=","){
  
  codon_stops <- codon_stops_from_genetic_code(genetic_code=genetic_code)
  if(length(codon_stops) == 0){
    print("WARNING: The Genetic Code Number provided does not correspond to any code in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes\nThe FilterCodonStop step is skipped")
    return(read_count_df)
  }
  
  # get unique list of variants to a df
  unique_asv_df <- data.frame(asv = unique(read_count_df$asv))
  
  # define a function to transform each sequence to a list of nucleotides
  seq_to_list_of_nt <- function(seq){
    return(strsplit(seq, "")[[1]])
  }
  # define a function to group a list of nucleotides to pieces of three (codons) starting at start_pos
  pool_nt_to_codons <- function(seq, start_pos){
    return(splitseq(seq, frame = start_pos, word = 3))
  }
  # check if a codon is present among a list of codons
  check_codon_stops <- function(codon_list, codon_stops=codon_stops){
    if(any(codon_stops %in% codon_list)){
      return(as.numeric(1))
    }
    else{
      return(as.numeric(0))
    }
  }
  # apply the function seq_to_list_of_nt to each unique sequence
  seqs1 <- lapply(unique_asv_df$asv, seq_to_list_of_nt)
  
  # Go through the 3 reading frames and check if there is a codon sotp in each of them
  for(j in 0:2){
    # get codons in reading frame j
    seqs2 <- lapply(seqs1, pool_nt_to_codons, j)
    # check if there is at least one codon stop among the codons; retiurn 1 if CodonStop, 0 if not
    unique_asv_df[[j+2]] <- as.vector(lapply(seqs2, check_codon_stops, codon_stops))
  }
  
  # Convert columns 2 to 4 to numeric and make a column with the sum of the 3 frames
  unique_asv_df[, 2:4] <- sapply(unique_asv_df[, 2:4], as.numeric)
  unique_asv_df$CodonStop <- rowSums(unique_asv_df[,2:4])
  
  # Keep only sequences where there is at least one frame without codon stop
  unique_asv_df <-unique_asv_df %>%
    filter(CodonStop < 3)
  # filter out ASV from read_count_df, where here is a codon stop in all reading frame
  read_count_df <- filter(read_count_df, (asv %in% unique_asv_df$asv))
  
  if(outfile !=""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' write_fasta
#' 
#' Write a fasta file using the sequences as ID
#'  
#' @param sequences list of sequences
#' @param filename output file name, including path
#' @param seq_as_id T/F; It TRUE, use sequences as seqID (FALSE by default)
#' @export
#' 
write_fasta <- function(sequences, filename, seq_as_id=F) {
  # Open the file for writing
  file <- file(filename, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(sequences)) {
    if(seq_as_id){
      header <- paste0(">", sequences[[i]])
    }else{
      header <- paste0(">", i)
    }
    writeLines(c(header, sequences[[i]], ""), file)
  }
  # Close the file
  close(file)
}

#' flagPCRerror_vsearch
#' 
#' Identify potential PCRerrors: ASVs very similar (max_mismatch) to another more frequent ASV (pcr_error_var_prop) in a data frame.
#' Adds a column to the input data frame, with 1 if sequence is a probable PCR error, and 0 otherwise
#'  
#' @param unique_asv_df data frame with the following variables: asv, read_count; ASVs should be unique
#' @param pcr_error_var_prop if the proportion of read_counts of two similar ASVs is bellow pcr_error_var_prop, the less abundant is flagged as a PCR error
#' @param max_mismatch maximum number of mismatches (gaps included) to consider two ASVs as similar
#' @param vsearch_path path to vsearch executable; can be empty if vsearch in the the PATH
#' @export
#' 
flagPCRerror_vsearch <- function(unique_asv_df, vsearch_path="", pcr_error_var_prop=0.1, max_mismatch=1){
  
  vsearch_path <- check_dir(vsearch_path)
  
  # no ASV in the unique_asv_df => return a dataframe with 0 for all ASVs in PCRerror column
  if(length(unique_asv_df$asv) == 0){ 
    unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
    return(unique_asv_df)
  }
  
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste('tmp_PCRerror_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  outdir_tmp <- check_dir(outdir_tmp)
  
  # make fasta file with unique reads; use sequences as ids
  fas <- paste(outdir_tmp, 'unique.fas', sep="")
  write_fasta(unique_asv_df$asv, fas, seq_as_id=T)
  # vsearch --usearch_global to find highly similar sequence pairs
  vsearch_out <- paste(outdir_tmp, 'unique_vsearch_out.out', sep="")
  #  vsearch <- paste(vsearch_path, "vsearch --usearch_global ", fas, " --db ", fas, " --quiet --iddef 1 --self --id 0.90 --maxaccepts 0 --maxrejects 0 --userfields 'query+target+ids+aln' --userout ", vsearch_out, sep="")
  vsearch <- paste(vsearch_path, "vsearch --usearch_global ", fas, " --db ", fas, ' --quiet --iddef 1 --self --id 0.90 --maxaccepts 0 --maxrejects 0 --userfields "query+target+ids+aln" --userout ', vsearch_out, sep="")
  #https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/system
  system(vsearch)
  
  # no vsearch hit => return unique_asv_df completed with a PCRerror, with 0 for all ASVs
  if(file.size(vsearch_out) == 0){
    unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
    # Delete the temp directory
    unlink(outdir_tmp, recursive = TRUE)
    return(unique_asv_df)
  }
  
  # read vsearch results
  results_vsearch<- read.csv(vsearch_out, header = FALSE, sep="\t")
  colnames(results_vsearch) <- c("query","target","nb_ids","aln")
  # none of the values easily outputted by vsearch take into the extenal gaps as a diff => correct this, based on the alnlen and the number of identities
  results_vsearch$nb_diff <- nchar(results_vsearch$aln) - results_vsearch$nb_ids
  # delete unnecessary columns
  results_vsearch <- select(results_vsearch, -c(nb_ids, aln))
  # keep only pairs with max_mismatch differences 
  results_vsearch <- results_vsearch %>%
    filter(nb_diff <= max_mismatch)
  
  # rename columns and add read counts to query and target ASVs in results_vsearch from unique_asv_df
  results_vsearch <- rename(results_vsearch, asv = query)
  results_vsearch <- left_join(results_vsearch, unique_asv_df, by="asv")
  results_vsearch <- rename(results_vsearch, qasv = asv)
  results_vsearch <- rename(results_vsearch, asv = target)
  results_vsearch <- rename(results_vsearch, qread_count = read_count)
  results_vsearch <- left_join(results_vsearch, unique_asv_df, by="asv")  
  results_vsearch <- rename(results_vsearch, tread_count = read_count)
  results_vsearch <- rename(results_vsearch, tasv = asv)
  
  # flag target ASV as a PCR error, if low read_count compared to query ASV
  results_vsearch$PCRerror_target <- ((results_vsearch$qread_count * pcr_error_var_prop) >= results_vsearch$tread_count)
  # keep only one column (tasv) with unique ASVs, that were flagged as PCRerror
  results_vsearch <- results_vsearch %>%
    filter(PCRerror_target==TRUE) %>%
    group_by(tasv) %>%
    select(tasv)
  # complete unique_asv_df with a PCRerror column
  unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
  unique_asv_df$PCRerror[unique_asv_df$asv %in% results_vsearch$tasv] <- 1
  
  # Delete the temp directory
  unlink(outdir_tmp, recursive = TRUE)
  
  return(unique_asv_df)
}

#' FilerPCRerror
#' 
#' Filter out an ASVs if it is very similar (max_mismatch) to another more frequent ASV (pcr_error_var_prop).
#' The whole data set can be analyzed at once (by_sample=F) or sample by sample.
#' Returns the filtered read_count_df data frame.
#'  
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param pcr_error_var_prop if the proportion of read_counts of two similar ASVs is less or equal to pcr_error_var_prop, the less abundant is flagged as a PCR error
#' @param max_mismatch maximum number of mismatches (gaps included) to consider two ASVs as similar
#' @param by_sample T/F, if T ASVs are flagged as an PCR error separately for each sample
#' @param sample_prop if by_sample=T, the ASV must be flagged as an PCRerror in sample_prop of the cases to be eliminated
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param vsearch_path path to vsearch executable; can be empty if vsearch in the the PATH
#' @export
#' 
FilterPCRerror <- function(read_count_df, outfile="", vsearch_path="", pcr_error_var_prop=0.1, max_mismatch=1, by_sample=T, sample_prop=0.8, sep=","){
  
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
      
      # flag PCR errors; add one column to unique_asv_df for each sample with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df_sample <- flagPCRerror_vsearch(unique_asv_df_sample, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch)
      
      # remove read_count column
      unique_asv_df_sample$read_count <- NULL
      #unique_asv_df_sample <- unique_asv_df_sample[, -which(names(unique_asv_df_sample) == "read_count")]
      # add a column for for each sample to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df <- left_join(unique_asv_df, unique_asv_df_sample, by = "asv")
    }
  }
  else{ # whole dataset
    # add a column to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
    unique_asv_df <- flagPCRerror_vsearch(unique_asv_df, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch)
  }
  
  # count the number of times each ASV has been flagged and when it has not. Ignore NA, when the ASV is not present in the sample
  unique_asv_df$yes <- rowSums(unique_asv_df[3:ncol(unique_asv_df)] == 1, na.rm = TRUE)
  unique_asv_df$no <- rowSums(unique_asv_df[3:(ncol(unique_asv_df)-1)] == 0, na.rm = TRUE)
  # keep only ASVs, that are not flagged in sample_prop proportion of the samples where they are present  
  unique_asv_df <- unique_asv_df %>%
    filter(yes/(yes+no) >= sample_prop)
  
  # eliminate potential PCRerrors from read_count_df
  read_count_df <- read_count_df %>%
    filter(!asv %in% unique_asv_df$asv)
  
  if(outfile !=""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}


#' flagChimera
#' 
#' Select chimeras in a dataframe of unique ASVs.
#' Add a column to the input data frame if the ASV is a probable chimera and 0 otherwise.
#'  
#' @param unique_asv_df data frame with the following variables: asv, read_count
#' @param abskew A chimera must be at least abskew times less frequent that the parental ASVs 
#' @param vsearch_path path to vsearch executable; can be empty if vsearch in the the PATH
#' @export
#'

flagChimera <- function(unique_asv_df, vsearch_path="", abskew=2){
  
  vsearch_path <- check_dir(vsearch_path)
  
  # no ASV in the unique_asv_df => return a data frame with 0 for all ASVs in Chimera column
  if(length(unique_asv_df$asv) == 0){ 
    unique_asv_df$chimera <- rep(0, length(unique_asv_df$asv))
    return(unique_asv_df)
  }
  
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste('tmp_FilterChimeara_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  outdir_tmp <- check_dir(outdir_tmp)
  
  # make fasta file with unique reads; use sequences as ids
  fas <- paste(outdir_tmp, 'unique.fas', sep="")
  # Open the file for writing
  file <- file(fas, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(unique_asv_df$asv)) {
    header <- paste0(">",unique_asv_df$asv[i], ";size=", unique_asv_df$read_count[i], sep="")
    writeLines(c(header, unique_asv_df$asv[i], ""), file)
  }
  close(file)
  
  # vsearch --usearch_global to find highly similar sequence pairs
  vsearch_out <- paste(outdir_tmp, 'unique_vsearch_out.out', sep="")
  vsearch <- paste(vsearch_path, "vsearch --uchime3_denovo ", fas, " --quiet --abskew ", abskew ," --uchimeout  ", vsearch_out, sep="")
  system(vsearch)
  
  # no vsearch hit => return unique_asv_df completed with a PCRerror, with 0 for all ASVs
  if(file.size(vsearch_out) == 0){
    unique_asv_df$chimera <- rep(0, length(unique_asv_df$asv))
    # Delete the temp directory
    unlink(outdir_tmp, recursive = TRUE)
    return(unique_asv_df)
  }
  
  # read vsearch results
  results_vsearch<- read.csv(vsearch_out, header = FALSE, sep="\t")
  # keep only pertinent columns
  results_vsearch <- select(results_vsearch, c(2, ncol(results_vsearch)))
  colnames(results_vsearch) <- c("asv", "chimera")
  results_vsearch$asv <- gsub(";size=[0-9]+", "", results_vsearch$asv)
  # keep only chimeras
  results_vsearch <- results_vsearch %>%
    filter(chimera == "Y")
  
  # complete unique_asv_df with chimara info
  unique_asv_df$chimera <- rep(0, length(unique_asv_df$asv))
  unique_asv_df$chimera[unique_asv_df$asv %in% results_vsearch$asv] <- 1
  
  # Delete the temp directory
  unlink(outdir_tmp, recursive = TRUE)
  return(unique_asv_df)
}

#' FilterChimera
#' 
#' Filter out Chimeras.
#' Returns the filtered read_count_df data frame
#'  
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param abskew A chimera must be at least abskew times less frequent that the parental ASVs 
#' @param by_sample T/F, if T ASVs are flagged as chimera separately for each sample
#' @param sample_prop if by_sample=T, the ASV deleted if they are flagged as chimera in at least sample_prop of the samples among the sample they are present
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @param vsearch_path path to vsearch executable; can be empty if vsearch in the the PATH
#' @export
#' 
FilterChimera <- function(read_count_df, outfile="", vsearch_path="", by_sample=T, sample_prop=0.8, abskew=abskew, sep=","){
  
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
      
      # flag chimeras; add one column to unique_asv_df for each sample with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df_sample <- flagChimera(unique_asv_df_sample, vsearch_path=vsearch_path, abskew=abskew)
      
      # remove read_count column
      unique_asv_df_sample <- select(unique_asv_df_sample, -c("read_count"))
      # add a column for for each sample to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df <- left_join(unique_asv_df, unique_asv_df_sample, by = "asv")
    }
  }else{ # whole dataset
    # add a column to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
    unique_asv_df <- flagChimera(unique_asv_df, vsearch_path=vsearch_path, abskew=abskew)
  }
  
  # count the number of times each ASV has been flagged and when it has not. Ignore NA, when the ASV is not present in the sample
  unique_asv_df$yes <- rowSums(unique_asv_df[3:ncol(unique_asv_df)] == 1, na.rm = TRUE)
  unique_asv_df$no <- rowSums(unique_asv_df[3:(ncol(unique_asv_df)-1)] == 0, na.rm = TRUE)
  # keep only ASVs, that are flagged in sample_prop proportion of the samples where they are present  
  unique_asv_df <- unique_asv_df %>%
    filter(yes/(yes+no) >= sample_prop)
  
  # eliminate potential Chimeras from read_count_df
  read_count_df <- read_count_df %>%
    filter(!asv %in% unique_asv_df$asv)
  
  if(outfile !=""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)
}

#' make_renkonen_df
#' 
#' Calculate renkonen distances between all pairs of replicates within samples.
#' Return a data frame with the following columns: sample, replicate1, replicate2, renkonen_d
#'  
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @export
#'
make_renkonen_df <- function(read_count_df){
  
  # list of samples
  sample_list <- unique(read_count_df$sample)
  # fine empty dataframe
  renkonen_df <- data.frame("sample" = character(),
                            "replicate1"  = character(),
                            "replicate2"  = character(),
                            "renkonen_d" = numeric())
  # loop over samples
  for(samp in sample_list){
    # get data for a given samples
    sample_df <- read_count_df %>%
      filter(sample == samp)
    # list of replicates for the sample
    replicate_list <- unique(sample_df$replicate)
    
    # loop over all pairs of replicates within sample
    for(i in 1:(length(replicate_list)-1)){
      # data frame for samp replicate i
      dfi <- filter(sample_df, replicate == replicate_list[[i]])
      for(j in ((i+1):length(replicate_list))){
        # data frame for samp replicate j
        dfj <- filter(sample_df, replicate == replicate_list[[j]])
        # make a data frame with all ASVs in at least one of the 2 replicates
        df <- full_join(dfi, dfj, by="asv")
        # replace NA by 0
        df <- df %>%
          mutate(read_count.x = ifelse(is.na(read_count.x), 0, read_count.x)) %>%
          mutate(read_count.y = ifelse(is.na(read_count.y), 0, read_count.y))
        # calculate  number of reads for variant x in replicate i / number of reads in replicate i
        df$read_count.x <- df$read_count.x/ sum(df$read_count.x)
        df$read_count.y <- df$read_count.y/ sum(df$read_count.y)
        # minimum of the above proportion between the 2 replicates
        df$min <- pmin(df$read_count.x, df$read_count.y)
        rdist <- 1- sum(df$min)
        # add line to renkonen_df
        new_line <- data.frame(sample = samp, replicate1 = as.character(replicate_list[[i]]), replicate2 = as.character(replicate_list[[j]]), renkonen_d = rdist)
        renkonen_df <- rbind(renkonen_df, new_line)
      }
    }
  }
  return(renkonen_df)
}

#' FilterRenkonen
#' 
#' Filter out all replicates that have renkonen distances above cutoff to most other replicates of the sample.
#' Returns the filtered read_count_df data frame
#'  
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param renkonen_distance_quantile quantile renkonen distance to determine cutoff value (among all distances, the lowest renkonen_distance_quantile proportion of values are considered as bellow cutoff)
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @export
#' 
FilterRenkonen <- function(read_count_df, outfile="", renkonen_distance_quantile=0.9, sep=","){
  # calculate Renkonen distances between all pairs of replicates of within sample
  renkonen_df <- make_renkonen_df(read_count_df)
  # determine the cut off renkonen distance; values > cutoff are considered as high
  renkonen_df <- renkonen_df %>%
    arrange(renkonen_d)
  last_row <- floor(length(renkonen_df$renkonen_d) * renkonen_distance_quantile)
  cutoff <- renkonen_df$renkonen_d[last_row]
  
  # get list of samples
  sample_list <- unique(renkonen_df$sample)
  # filter out replicates sample by smaple
  for(samp in sample_list){
    sample_df <- renkonen_df %>%
      filter(sample == samp)
    # make a df with the replicate columns exchanged
    sample_tmp <- data.frame("sample" = sample_df$sample,
                             "replicate1"  = sample_df$replicate2,
                             "replicate2"  = sample_df$replicate1,
                             "renkonen_d" = sample_df$renkonen_d)
    # complete sample_df to include distances between repl X and Y and also Y and X
    sample_df <- rbind(sample_df, sample_tmp)
    # get unique list of replicates
    replicate_list <-  unique(sample_df$replicate1)
    # the minimum number of distances to be bellow cutoff, to keep the replicate
    min_number_of_distances_bellow_cutoff <- (length(replicate_list) -1) / 2
    
    # keep only distances above cutoff in sample_df
    sample_df <- sample_df %>%
      filter(renkonen_d > cutoff)
    
    # count for each replicate the number of distances above ctuoff 
    sample_df <- sample_df %>%
      group_by(replicate1) %>%
      summarize(n_dist=length(renkonen_d)) %>%
      filter(n_dist > min_number_of_distances_bellow_cutoff)
    
    # eliminate replicates with too many distances above cutoff
    read_count_df <- read_count_df %>%
      filter(!(sample == samp & replicate %in% sample_df$replicate1))
  }
  
  if(outfile !=""){
    write.table(read_count_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_df)  
}

#' Pool replicates sample by sample
#' 
#' Take the mean non-zero read counts over replicates for each sample and asv.
#' Returns a data frame with the following columns: asv, plate, marker, sample, mean_read_count (over replicates)
#'  
#' @param read_count_df data frame with the following variables: asv_id, sample, replicate, read_count, asv
#' @param digits round the mean read counts to digits
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if no file name provided, only a data frame is returned
#' @export
#'
PoolReplicates <- function(read_count_df, digits=0, outfile="", sep=","){
  
  read_count_samples_df <- read_count_df %>%
    group_by(asv_id,sample,asv) %>%
    summarize(mean_read_count = mean(read_count), .groups="drop_last") %>%
    select(asv_id, sample, mean_read_count, asv)
  
  read_count_samples_df$mean_read_count <- round(read_count_samples_df$mean_read_count, digits =digits)
  
  if(outfile !=""){
    write.table(read_count_samples_df, file = outfile,  row.names = F, sep=sep)
  }
  return(read_count_samples_df)
}
