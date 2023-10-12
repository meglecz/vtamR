#' Merge
#' 
#' Merge forward are reverse fastq read pairs to fasta.
#' Output fasta files can be compressed if compress option is used.
#' The output fastainfo.csv file is similar to the fastqinfo file, but the fastq columns are replaced by a fasta column with the name of the output files
#' Returns data frame corresponding to the output csv file
#'  
#' @param fastqinfo_df data frame with column: tag_fw,primer_fw,tag_rv,primer_rv,plate,marker,sample,sample_type(mock/negative/real),habitat(optional),replicate,fastq_fw,fastq_rv  
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
#' @param sep separator in input and output csv files ;default: ","
#' @param compress [gz/zip]; compress output to gz or zip files. Do not compress by default;
#' @export
#'

Merge <- function(fastqinfo_df, fastqdir, vsearch_path="", outdir="", fastq_ascii=33, fastq_maxdiffs=10, fastq_maxee=1, fastq_minlen=50, fastq_maxlen=500, fastq_minmergelen=50, fastq_maxmergelen=1000, fastq_maxns=0, fastq_truncqual=10, fastq_minovlen=50, fastq_allowmergestagger=F, sep=",", compress=0){

  outdir<- check_dir(outdir)
  # get unique list of fastq file pairs
  tmp <- fastqinfo_df %>%
    select(fastq_fw, fastq_rv)
  tmp <- unique(tmp)
  tmp$fasta <- NA
  
  for(i in 1:nrow(tmp)){# for each file pairs
    # use the name of the fw file and replace extension by fasta
    outfile <- sub("\\..*", ".fasta", tmp[i,1])
    tmp$fasta[i] <- outfile
    outfile <- paste(outdir, outfile, sep="")
    # add path to input filenames
    fw_fastq <- paste(fastqdir, tmp[i,1], sep="")
    rv_fastq <- paste(fastqdir, tmp[i,2], sep="")
    
    if(fw_fastq == outfile){ # stop the run if input and output files have the same name
      stop("ERROR: Input and output directories for fastq and fasta files are indentical. Please, give a different output directory")
    }
    
    vsearch <- paste(vsearch_path, "vsearch --fastq_mergepairs ", fw_fastq, " --reverse ", rv_fastq ," --fastaout ",outfile," --quiet --fastq_ascii ",fastq_ascii," --fastq_maxdiffs ", fastq_maxdiffs, " --fastq_maxee ", fastq_maxee, " --fastq_minlen ", fastq_minlen, " --fastq_maxlen ",fastq_maxlen, " --fastq_minmergelen ",fastq_minmergelen," --fastq_maxmergelen ",fastq_maxmergelen," --fastq_maxns ", fastq_maxns, " --fastq_truncqual ", fastq_truncqual, " --fastq_minovlen ", fastq_minovlen, sep="")
    if(fastq_allowmergestagger){ # if reads are longer than the amplicon
      paste(vsearch, " --fastq_allowmergestagger", sep="")
    }
        print(vsearch)
    system(vsearch)
    
    if(compress == "gz"){
      # Specify the path for the gzipped output file
      outfile_gz <- paste(outfile, ".gz", sep="")
      # Open the existing uncompressed file for reading
      file_content <- readBin(outfile, "raw", file.info(outfile)$size)
      # Create a gzipped copy of the file
      gz <- gzfile(outfile_gz, "wb")
      writeBin(file_content, gz)
      close(gz)
      file.remove(outfile)
      tmp$fasta[i] <- paste(tmp$fasta[i], ".gz", sep="")
    }
    if(compress == "zip"){
      # when zipping a file using a full path, the zip file will contain the embedded directory structure
      # To avoid this, change the wd to the output dir zip, the file and change back to the orignal wd
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


#' make_adapater_fasta
#' 
#' Make a fasta file with tag combinations
#' 
#' @param fastainfo_df data frame with column: tag_fw,tag_rv,fasta  
#' @param fasta_file name of the fasta file for wich tag-combination should be written to a tags.fasta file
#' @param tag_to_end [T/F] if T use anchored search (tags are at the extremity of the sequence)
#' @param outdir output directory
#' @export
#'
make_adapater_fasta <- function(fastainfo_df, fasta_file, tag_to_end, outdir){
  # select tag combinations pour le fasta file
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

#' SortReads
#' 
#' Demultiplex each input fasta file using the tag combinations at the extremities of the merged reads.
#' Trim primers from demultiplexed reads.
#' Output fasta files can be compressed if compress option is used.
#' The output fileinfo.csv file is similar to the fastainfo file, but the but do not have tag and primer columns.
#' Returns data frame of corresponding to the output csv file.
#'  
#' @param fastainfo_df data frame with column: tag_fw,primer_fw,tag_rv,primer_rv,plate,marker,sample,sample_type(mock/negative/real),habitat(optional),replicate,fasta  
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
#' @param compress [gz/zip]; compress output to gz or zip files. Do not compress by default;
#' @export

SortReads <- function(fastainfo_df, fastadir, outdir="", cutadapt_path="" ,vsearch_path="", check_reverse=F, tag_to_end=T, primer_to_end=T, cutadapt_error_rate=0.1,cutadapt_minimum_length=50,cutadapt_maximum_length=500, sep=",",  compress=0){
  #########
  # SortReads_no_reverse dos the wole demultilexin, trimming and compress on the + strane
  # If seqences are notorienetd, the -strand sould be checked => 
  # run SortReads_no_reverse of plus starnd and on - starnd after swapping fw and rev tags and primers,
  # take the reverse complement of the -stard results
  # pool the results of the 2 strands
  # compress if necessary
  #########
  #!!!!!!!!!!!!!!! check latter the zip option for compress, especially for windows
  
  # run on strand +
  if(check_reverse){
    #### use +strand, output to sorted_dir, uncompressed
    fileinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df, fastadir=fastadir, outdir=outdir, cutadapt_path=cutadapt_path, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=0)
    
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
    fileinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df_tmp, fastadir=fastadir, outdir=rc_dir, cutadapt_path=cutadapt_path, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=0)
    
    ### reverse complment and pool
    # get list of files demultiplexed on - strand
    files <- list.files(path = rc_dir)
    # Filter the files based on the motif using regular expressions
    # reverse complement sequences on the minus stand, and append info to the plus strand output
    files <- grep(pattern = "\\.fasta", x = files, value = TRUE)
    for(i in 1:length(files)){
      plus <- paste(outdir, files[i], sep="")
      minus <- paste(rc_dir, files[i], sep="")
      minus_rc <- paste(rc_dir, "rc_", files[i], sep="")
      # reverse complement sequences in minus file
      rev_comp_cmd <- paste("vsearch --fastx_revcomp ", minus, " --fastaout ", minus_rc, " --quiet", sep="")
      print(rev_comp_cmd)
      system(rev_comp_cmd)
      # append content of minus_rc to plus file
      file.append(plus, minus_rc)
    }
    
    # delete temporary reverse_comp dir
    unlink(rc_dir, recursive = TRUE)
    
    ### compress
    if(compress == "gz"){
      
      for(i in 1:nrow(fileinfo_df)){
        
        file <- fileinfo_df$filename[i]
        
        # Specify the path for the gzipped output file
        file <- paste(outdir, file, sep="")
        file_gz <- paste(file, ".gz", sep="")
        
        # Open the existing uncompressed file for reading
        file_content <- readBin(file, "raw", file.info(file)$size)
        # Create a gzipped copy of the file
        gz <- gzfile(file_gz, "wb")
        writeBin(file_content, gz)
        close(gz)
        file.remove(file)
        fileinfo_df$filename[i] <- paste(fileinfo_df$filename[i], ".gz", sep="")
      }
      write.table(fileinfo_df, file = paste(outdir, "fileinfo.csv", sep=""),  row.names = F, sep=sep) 
    }
    if(compress == "zip"){
      
      for(i in 1:nrow(fileinfo_df)){
        
        file <- fileinfo_df$filename[i]
        
        file_zip <- paste(file, ".zip", sep="")
        zip(file_zip, file)
        file.remove(outfile)
        fileinfo_df$filename[i] <- paste(fileinfo_df$filename[i], ".gz", sep="")
      }
      # overwrite fileinfo file using zipped filenames
      write.table(fileinfo_df, file = paste(outdir, "fileinfo.csv", sep=""),  row.names = F, sep=sep) 
    }
  }
  else{
    # check only + strand
    fileinfo_df <- SortReads_no_reverse(fastainfo_df=fastainfo_df, fastadir=fastadir, outdir=outdir, cutadapt_path=cutadapt_path, tag_to_end=tag_to_end, primer_to_end=primer_to_end, cutadapt_error_rate=cutadapt_error_rate, cutadapt_minimum_length=cutadapt_minimum_length, cutadapt_maximum_length=cutadapt_maximum_length, sep=sep, compress=compress)
  }
  return(fileinfo_df)
  
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
#' @param fastainfo_df data frame with column: tag_fw,primer_fw,tag_rv,primer_rv,plate,marker,sample,sample_type(mock/negative/real),habitat(optional),replicate,fasta  
#' @param fastadir directory with input fasta files (listed in fastainfo_df$fasta)
#' @param cutadapt_path path to cutadapt executables
#' @param outdir output directory
#' @param tag_to_end tags are at the extremity of the reads (starting at the first base); default: T
#' @param primer_to_end primers follow directly the tags (no heterogeneity spacer); default: T
#' @param cutadapt_error_rate maximum proportion of errors between primers and reads (for tags, exact match is required); default: 0.1
#' @param cutadapt_minimum_length minimum length of the trimmed sequence; default: 50
#' @param cutadapt_maximum_length maximum length of the merged sequence; default: 500
#' @param sep separator in input and output csv files; default: ","
#' @param compress [gz/zip]; compress output to gz or zip files. Do not compress by default;
#' @export
SortReads_no_reverse <- function(fastainfo_df, fastadir, outdir="", cutadapt_path="", tag_to_end=T, primer_to_end=T, cutadapt_error_rate=0.1,cutadapt_minimum_length=50,cutadapt_maximum_length=500, sep=",",  compress=0){
  # do the complete job of demultiplexing and trimming of input file without checking the reverse sequences
  # !!!!! Makes gz files using cutadapt, but adapt this to windows, later
  
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
#          primer_trim_cmd <- paste(cutadapt_path, "cutadapt --cores=0 -e ",cutadapt_error_rate ," --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length ," --maximum-length ", cutadapt_maximum_length, " -g '", primer_fwl, ";min_overlap=",nchar(primer_fwl),"...", primer_rvl_rc,  ";min_overlap=",nchar(primer_rvl_rc),"' --output ", primer_trimmed_file, " ", tag_trimmed_file, sep="")
          primer_trim_cmd <- paste(cutadapt_path, "cutadapt --cores=0 -e ",cutadapt_error_rate ," --no-indels --trimmed-only --minimum-length ", cutadapt_minimum_length ," --maximum-length ", cutadapt_maximum_length, ' -g "', primer_fwl, ';min_overlap=',nchar(primer_fwl),'...', primer_rvl_rc,  ';min_overlap=',nchar(primer_rvl_rc),'" --output ', primer_trimmed_file, " ", tag_trimmed_file, sep="")
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

#' count_seq
#' 
#' Count the number of sequences in the input fasta file. Input can be uncompressed, gz, or zip file
#'  
#' @param file fasta file, without it's path
#' @param dir directory that contains the input file
#' @export

count_seq <- function(dir="", file){
  
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



#' RandomSeq
#' 
#' Random select n sequences from each input fasta file. The output is the same compression type (if any) as the input
#'  
#' @param fastainfo_df data frame with a 'fasta' column containing input file names; files can be compressed in gz and zip format
#' @param n integer; the number of randomly selected sequences 
#' @param fasta_dir directory that contains the input fasta files
#' @param outdir directory for the output files
#' @export
#' 
RandomSeq <- function(fastainfo_df, fasta_dir="", outdir="", n){
  
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

#' select_sequences
#' 
#' Select sequences from the input file that correspond to the vector of integers
#'  
#' @param file input fasta file; can be  uncompressed o compressed in gz format only 
#' @param fasta_dir path to the input fasta files
#' @param outdir directory for the output files; same compression as the input file
#' @param random_integers vector of random integers between 1 and the number of sequences in the input file
#' @export
#'
select_sequences <- function(fasta_dir="", file, outdir="", random_integers){
  
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