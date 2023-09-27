#' Merge
#' 
#' Merge forward are reverse fastq read pairs to fasta.
#' Output fasta files can be compressed if compress option is used.
#' The output fastainfo file is similar to the fastqinfo file, but the fastq columns are replaced by a fasta column with the name of the output files
#' Returns data frame of from the outout csv file
#'  
#' @param fastqinfo_df data frame with column: tag_fw,primer_fw,tag_rv,primer_rv,plate,marker,sample,sample_type(mock/negative/real),habitat(optional),replicate,fastq_fw,fastq_rv  
#' @param fastqdir directory with input fastq files
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
      outfile_zip <- paste(outfile, ".zip", sep="")
      zip(outfile_zip, outfile)
      tmp$fasta[i] <- paste(tmp$fasta[i], ".zip", sep="")
      file.remove(outfile)
    }
    
  }
  # make fastinfo file
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
  # define complemetary bases
  comp <- data.frame(orig=c("A","T","C","G","R","Y","M","K","B","H","D","V","a","t","c","g","r","y","m","k","b","h","d","v"),
                     complement=c("T","A","G","C","Y","R","K","M","V","D","H","B","t","a","g","c","y","r","k","m","v","d","h","b")
  )
  # revers, split sequneces and make dataframe
  sequence_df <- data.frame(reverse=rev(strsplit(sequence, NULL)[[1]]))
  # add complimetary nt
  sequence_df <- left_join(sequence_df, comp, by=c("reverse"="orig"))
  # collaps vector to string
  reverse_comp <- paste(sequence_df$complement, collapse = "")
  
  return(reverse_comp)
}
