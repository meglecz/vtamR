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
#' @param outfile Name of the output csv file with asv_id, sample, replicate, read_count and asv as columns; if NULL, only a data frame is returned
#' @param sep separator in csv files; default: ","
#' @param asv_list name of the file, containing asvs and asv_ids from earlier analyses. Optional. It is used to homogenize asv_ids between different data sets
#' @param update_asv_list [T/F]; If TRUE, write a new file (columns:asv_id, asv) completed with new ASVs
#' @export
#' 
read_fastas_from_fileinfo <- function (sortedinfo_df, dir="", outfile="", sep=",", asv_list="", update_asv_list=F) {
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
    
    read_count_df_tmp <- read_fasta_seq(filename=fas, dereplicate=T)
    read_count_df_tmp$sample <- sortedinfo_df[i,"sample"]
    read_count_df_tmp$replicate <- sortedinfo_df[i,"replicate"]
    read_count_df <- rbind(read_count_df, read_count_df_tmp)
  }
  rm(read_count_df_tmp)
  # reorder columns
  read_count_df <- read_count_df[, c("asv", "sample", "replicate", "read_count")]
  
  # add asv_id, but do not add new asv s to asv_list; it is probably too early at this stage, better to do it after swarm
  read_count_df <- add_ids(read_count_df, asv_list=asv_list, sep=sep,  update_asv_list=update_asv_list)
  
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
#' @param update_asv_list [T/F]; If TRUE, write a new file (columns:asv_id, asv) completed with new ASVs
#' @param sep separator in csv files; default: ","
#' @export
#' 
add_ids <- function(read_count_df, asv_list=asv_list, sep=sep, update_asv_list=T){
  
 
  if(asv_list != ""){  # read already existing asvs, if the file is given
    asv_df <- read.csv(asv_list, sep=sep, header=TRUE)
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
  if(update_asv_list){
    str <- paste("_", trunc(as.numeric(trunc(Sys.time()))), ".", sep="")
    new_file <- sub("\\.", str, asv_list)
    write.table(asv_df, file=new_file, row.names = FALSE, sep=sep)
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

