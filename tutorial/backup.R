make_known_occurrences <- function(read_count_samples_df, fileinfo="", mock_composition="", sep=",", out=""){
  
  # read info on samples types and keep only relevant info
  fileinfo_df <- read.csv(fileinfo, header=T, sep=sep) %>%
    select(plate, marker, sample, sample_type, habitat)
  # make a pms (plate.marker.sample) column
  fileinfo_df$pms <- paste(fileinfo_df$plate, fileinfo_df$marker, fileinfo_df$sample, sep=".")
  # get unique lines to avoid replicates
  fileinfo_df <- unique(fileinfo_df)
  
  # define data frame for known occurrences
  occurrence_df <- read_count_samples_df
  occurrence_df$pms <- paste(occurrence_df$plate, read_count_samples_df$marker, read_count_samples_df$sample, sep=".")
  occurrence_df$pmsasv <- paste(occurrence_df$pms, read_count_samples_df$asv, sep=".")
  occurrence_df$action <- rep(NA, nrow(occurrence_df))
  
  # flag occurrences in negative control samples as delete
  occurrence_df <- flag_from_negative_controls(occurrence_df, fileinfo_df)
  
  # flag all expected occurrences in mock samples as "keep", NA for tolerate, and delete for all others
  occurrence_df <- flag_from_mock(occurrence_df, mock_composition, fileinfo_df, sep=sep)
  
  # flag occurrences as keep with low read count in habitat, compared to the others habitats
  occurrence_df <- flag_from_habitat(occurrence_df, fileinfo_df) 
  
  # keep nly relevant colums and line, sort data
  occurrence_df <- occurrence_df %>%
    select(plate,marker,sample,action,asv) %>%
    filter(!is.na(action)) %>%
    arrange(plate, marker, sample, action)
  # write to outfile
  write.table(occurrence_df, file=out, row.names = F, sep=sep)
}

flag_from_negative_controls <- function(occurrence_df, fileinfo_df){
  # keep only negative controls in fileinfo_df
  fileinfo_df <- fileinfo_df %>%
    filter(sample_type=="negative")
  # get unique list of pms with negative controls
  neg <- unique(fileinfo_df$pms)
  # set all occurrences to keep if in negative control
  occurrence_df$action <- ifelse(occurrence_df$pms %in% neg, "delete", occurrence_df$action)
  return(occurrence_df)
}

flag_from_mock <- function(occurrence_df, mock_composition, fileinfo_df, sep=","){
  
  mock_composition_df <- read.csv(mock_composition, header=T, sep=sep)
  mock_composition_df$pms <- paste(mock_composition_df$plate, mock_composition_df$marker, mock_composition_df$sample, sep=".")
  mock_composition_df$pmsasv <- paste(mock_composition_df$pms, mock_composition_df$asv, sep=".")
  
  # set temporarily all occurrences in mock as delete
  mock_samples <- unique(mock_composition_df$pms)
  occurrence_df$action <- ifelse(occurrence_df$pms %in% mock_samples, "delete", occurrence_df$action)
  
  # reset expected occurrences in mocks as keep
  keep_occurrences <- mock_composition_df %>%
    filter(action=="keep")
  occurrence_df$action <- ifelse(occurrence_df$pmsasv %in% keep_occurrences$pmsasv, "keep", occurrence_df$action)
  
  # reset tolerate occurrences in mocks as NA
  tolerate_occurrences <- mock_composition_df %>%
    filter(action=="tolerate")
  occurrence_df$action <- ifelse(occurrence_df$pmsasv %in% tolerate_occurrences$pmsasv, NA, occurrence_df$action)
  
  return(occurrence_df)
}


flag_from_habitat <- function(occurrence_df, fileinfo_df){
  # add habitat and sample_type to occurrence_df
  occurrence_df <- left_join(occurrence_df, fileinfo_df)
  # group by asv and habitat and count the total number of reads for each habitat-asv combination
  tmp <- occurrence_df %>%
    group_by(habitat, asv) %>%
    summarize(habitat_read_count=sum(mean_read_count)) %>%
    filter(!is.na(habitat))
  # count the number of habitats for each asv and keep only the ones present in at least two different habitats
  tmp2 <- tmp %>%
    group_by(asv) %>%
    summarize(nb_habitat=length(asv)) %>%
    filter(nb_habitat>1)
  # keep only selected asvs in tmp
  tmp <- tmp[tmp$asv %in% tmp2$asv, ]
  # get the total readcount for each asv in tmp
  tmp3 <- tmp %>%
    group_by(asv) %>%
    summarize(sum_read_count = sum(habitat_read_count))
  # add total readcount of asv to tmp and keep only lines where asv ih babitats where it is less frequent than in the others
  tmp <- left_join(tmp, tmp3, by="asv")
  tmp <- tmp[tmp$habitat_read_count/tmp$sum_read_count < 0.5, ]
  # keep only pertinent columns in tmp and add hab_action column with "delete"
  tmp <- tmp %>%
    select(habitat, asv)
  tmp$hab_action <- rep("delete", nrow(tmp))
  
  occurrence_df <- left_join(occurrence_df, tmp, by=c("habitat", "asv"))
  occurrence_df$action[which(occurrence_df$hab_action=="delete")] <- "delete"
  
  occurrence_df <- occurrence_df %>%
    select(-hab_action) %>%
    select(-sample_type) %>%
    select(-habitat)
  
  return(occurrence_df)
}





####################################################
####################################################
####################################################
####################################################

# info divers


decompress<- function(file=input_fasta, dir=fasta_dir){
  backup_wd <- getwd()
  setwd(dir)
  
  decompressed_file <- ""
  
  if(endsWith(file, ".zip")){
    unzip(file)
    decompressed_file <- sub(".zip", "", file)
  }else if(endsWith(file, ".gz")){
    # Open a connection to the gzipped file
    file_connection <- gzfile(file, "rb")  # "rb" stands for read binary mode
    
    # Read the contents of the gzipped file
    file_contents <- readLines(file_connection)
    
    # Close the file connection
    close(file_connection)
    
  }else{
    setwd(backup_wd)
    stop("Only gz and zip files are supported")
  }
  
  
  # reset wd
  setwd(backup_wd)
}


compressed_file <- compress_file("/home/meglecz/vtamR/local/out/sorted/plate1-MFZR-14ben01-1.fasta", compress="gz", remove_uncompressed=T)
compress_file <- function(file, compress="gz", remove_uncompressed=T){
  #### !!!!! zip is not working correctly in minux, since it produces and archive with the embedded directory structure in the path => correct it when using windows
  outfile <- file
  if(compress == "gz"){
    # Specify the path for the gzipped output file
    file_gz <- paste(file, ".gz", sep="")
    # Open the existing uncompressed file for reading
    file_content <- readBin(file, "raw", file.info(file)$size)
    # Create a gzipped copy of the file
    gz <- gzfile(file_gz, "wb")
    writeBin(file_content, gz)
    close(gz)
    if(remove_uncompressed){
      file.remove(file)
    }
    outfile <- file_gz
  }
  if(compress == "zip"){
    file_zip <- paste(file, ".zip", sep="")
    zip(file_zip, file, recursive = FALSE)
    if(remove_uncompressed){
      file.remove(file)
    }
    outfile <- file_zip
  }
  
  return(outfile)
}


#' write_fasta_seq_as_id
#' 
#' Write a fasta file using the sequences as ID
#'  
#' @param sequences list of sequences
#' @param filename output file name, including path
#' @export
#' 
write_fasta_seq_as_id <- function(sequences, filename) {
  # Open the file for writing
  file <- file(filename, "w")
  # Iterate over the sequences and write them to the file
  for (i in seq_along(sequences)) {
    header <- paste0(">", sequences[[i]])
    writeLines(c(header, sequences[[i]], ""), file)
  }
  # Close the file
  close(file)
}



### run blast with seqinr
library("seqinr")
# Load the query sequence from a file or define it directly
query_sequence <- read.fasta(file = "/home/meglecz/vtamR/local/small_test/small_test1_1.fas")
# Perform local sequence search using blast()
results <- blast(query = query_sequence, database = "/home/meglecz/vtamR/local/small_test/small_test1_1.fas", type = "DNA")
# Print the results
print(results)


if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")
library("DECIPHER")

# Load the query sequence from a file or define it directly
query_sequence <- readDNAStringSet("/home/meglecz/vtamR/local/small_test/small_test1_1.fas")

# Load the reference sequences from a file or define them directly
reference_sequences <- readDNAStringSet("database.fasta")

# Perform local sequence search using blastn()
results <- blastn(query_sequence, reference_sequences)

# Print the results
print(results)

### get the list of function in a package
functions_seqinr <- ls("package:seqinr", all.names = TRUE)
print(functions_seqinr)


### Define the command to run the third-party program
command <- "path/to/program"
args <- c("arg1", "arg2")
# Run the command using system2()
output <- system2(command, args, stdout = TRUE, stderr = TRUE)
# Print the output
print(output)


command <- "path/to/program arg1 arg2"
# Run the command using system()
system(command)



# start blast from R
myPipe <- pipe( "blastall -p blastp -i text.fasta -d data.fasta" )
results <- read.table( myPipe )
colnames( blastResults ) <- c( "QueryID",  "SubjectID", "Perc.Ident",
                               "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                               "S.start", "S.end", "E", "Bits" )


# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html



# Create a connection to an external command (e.g., "ls" command on Unix/Linux)
cmd <- pipe("ls", "r")  # Open the "ls" command for reading

# Read the output from the command
output <- readLines(cmd)

# Close the connection
close(cmd)

# Print the output
print(output)



makeblastdb -in local/out/small/tmp_1687945245/test.fas -dbtype nucl
blastn –db local/out/small/tmp_1687945245/test.fas –query local/out/small/tmp_1687945245/test.fas –outfmt 6 –out local/out/small/tmp_1687945245/test_blastout.out

pipe_vsearch <- pipe(vsearch)
results <- read.table( pipe_vsearch )

blast <- "blastn -db local/out/small/tmp_1687945245/test.fas -query local/out/small/tmp_1687945245/test.fas -outfmt 6"
myPipe <- pipe(blast)
results <- read.table( myPipe )
colnames( results ) <- c( "QueryID",  "SubjectID", "Perc.Ident",
                          "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                          "S.start", "S.end", "E", "Bits" )

