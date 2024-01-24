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


#' read_asv_table
#' 
#' Read asv table in wid format to a data frem in long format
#'  
#' @param filename name of the input file including full path; columns: asv, seq_id, plate.marker.sample.replicate columns containing read counts
#' @param sep separator; default ","
#' @export 
#' 
read_asv_table <- function(filename, sep=","){
  
  df <- read.csv(filename, sep=sep)
  long_df <- gather(df, key="plate.marker.sample.replicate", value ="read_count", -asv, -seq_id)  %>% 
    filter(read_count > 0)
  # separate column
  long_df <- separate(long_df, "plate.marker.sample.replicate", into=c("plate", "marker", "sample", "replicate"), sep="\\.")
  
  return(long_df)
}

#' read_asv_table_sample
#' 
#' Read asv table in wid format to a data frem in long format
#'  
#' @param filename name of the input file including full path; columns: asv, seq_id, plate.marker.sample columns containing read counts
#' @param sep separator; default ","
#' @export 
#' 
read_asv_table_sample <- function(filename, sep=","){
  
  df <- read.csv(filename, sep=sep)
  long_df <- gather(df, key="plate.marker.sample", value ="mean_read_count", -asv, -seq_id)  %>% 
    filter(mean_read_count > 0)
  # separate column
  long_df <- separate(long_df, "plate.marker.sample", into=c("plate", "marker", "sample"), sep="\\.")
  
  return(long_df)
}

#' compare_df
#' 
#' Compare two dataframes of and returns teh full join ob the two
#'  
#' @param df1 dataframe with columns: "asv", "plate","marker", "sample","replicate","read_count"
#' @param df2 dataframe with columns: "asv", "plate","marker", "sample","replicate","read_count"
#' @param step string to include in the FAIL or PASS message 
#' @export 
#' 
compare_df<- function(df1, df2, step=""){
  
  df1 <- df1 %>%
    select("asv", "plate","marker", "sample","replicate","read_count_vtamR"="read_count")
  
  df1 <- full_join(df1, df2, by=c("plate", "marker", "sample", "replicate", "asv"))
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
    select("asv", "plate","marker", "sample","mean_read_count_vtamR"="mean_read_count")
  
  df1 <- full_join(df1, df2, by=c("plate", "marker", "sample", "asv"))
  comp <- df1$mean_read_count == df1$mean_read_count_vtamR
  comp[(is.na(comp))] <- FALSE
  
  if(any(!comp)){
    result <- paste(step, ": FAIL", sep="")
  }else{
    result <- paste(step, ": PASS", sep="")
  }
  print(result)
  return(df1)
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

