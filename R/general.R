#' Check directory
#' 
#' Create dir if does not exists.
#' Add slash to the end of the directory name.
#' 
#' @param dir directory name
#' @export
check_dir <- function(dir){
  if(!(endsWith(dir, "/"))){
    dir <- paste(dir, '/', sep="")
  }
  if(!dir.exists(dir)){
    dir.create(dir, recursive =TRUE)
  }
  return(dir)
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
