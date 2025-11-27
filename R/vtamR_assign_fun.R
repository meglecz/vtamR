#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first 
#' @importFrom dplyr relocate if_else across slice_head
#' @importFrom utils read.csv write.table read.table read.delim count.fields
NULL



#' TaxAssignRDP
#' 
#' Assign 16S bacteria/Archaea sequences to taxa using the rRDP package
#' 
#' This function needs the rRDP, rRDPData and Biostrings packages.
#' 
#' By default, the 16S Bacterial/Archaea RDP database is used from rRDPData. 
#' (http://sourceforge.net/projects/rdp-classifier/). 
#' Otherwise the directory name containing the a trained classifier object created 
#' with trainRDP(), can be provided.
#'  
#' @param asv Data frame or csv file containing an asv and asv_id columns.
#' @param dir Directory with the data for an existing classifier 
#' (created with trainRDP()). If NULL, default classifier trained with 
#' the data shipped with RDP is used. See RDP documentation (?rdp)
#' @param max_memory Positive integer: RAM available for rRDP in Gb.
#' @param confidence Numeric; Value between 0-1. Minimum confidence value, to accept
#' the taxonomic assignment
#' @param rm_chloroplast Logical; If true, reset the lineage to NA if class is
#' Chloroplast
#' @param sep Field separator character in input and output csv files.
#' @param outfile Character string: csv file name to print the output data 
#' frame if necessary. If empty, no file is written.
#' @param quiet logical: If TRUE, suppress informational messages and only 
#' show warnings or errors.
#' @returns Data frame with the following columns:
#' asv_id,domain,kingdom,phylum,class,order,family,genus
#' @examples
#' \dontrun{
#' taxa <- TaxAssignRDP(asv=read_count_df, confidence=0.7, max_memory=8, rm_chloroplast=FALSE)
#' }
#' @export
#'
#'
TaxAssignRDP <- function(
  asv,
  dir=NULL,
  max_memory=1,
  confidence=0.8,
  rm_chloroplast=TRUE,
  outfile="",
  quiet=TRUE,
  sep=","){
  
  ###### test if rRDP is installed
  if (!requireNamespace("rRDP", quietly = TRUE) || !requireNamespace("rRDPData", quietly = TRUE)) {
    stop(
      "Package 'rRDP' and 'rRDPData' are required for this function.\n",
      "Please install them with:\n",
      "  if (!requireNamespace('BiocManager', quietly = TRUE))\n",
      "    install.packages('BiocManager')\n",
      "  BiocManager::install('rRDP')\n",
      "  BiocManager::install('rRDPData')",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop("The Biostrings package is required for this function. Please install it via BiocManager::install('Biostrings').")
      stop(
        "The Biostrings package is required for this function.\n",
        "Please install it with:\n",
        "  if (!requireNamespace('BiocManager', quietly = TRUE))\n",
        "    install.packages('BiocManager')\n",
        "  BiocManager::install('Biostrings')",
        call. = FALSE
      )
    }

  
  
  # can accept df or file as an input
  if(is.character(asv)){
    # read known occurrences
    asv_df <- read.csv(asv, header=T, sep=sep)
  }else{
    asv_df <- asv
  }
  
  ### get unique ASVs
  asv_df <- asv_df %>%
    group_by(asv_id) %>%
    summarize(asv = head(asv, n=1))
  #### make a DNAStringSet
  dna_set <- Biostrings::DNAStringSet(asv_df$asv)
  names(dna_set) <- asv_df$asv_id
  
  #### assign 
  java_arg = paste("-Xmx", max_memory, "g", sep="")
  taxa <- predict(rRDP::rdp(dir=dir), 
                  confidence=confidence, 
                  java_args=java_arg,
                  dna_set)
  
  
  #### add asv 
  taxa$asv_id <- as.numeric(rownames(taxa))
  taxa <- left_join(taxa, asv_df, by="asv_id")
  
  #### Complete df with missing taxlevels
  # Define taxonomy columns
  tax_cols <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  # Add missing columns as NA
  missing_cols <- setdiff(tax_cols, names(taxa))
  if (length(missing_cols) > 0) {
    taxa[missing_cols] <- NA_character_
  }
  # re-order columns
  taxa <- taxa %>%
    relocate(any_of(tax_cols), .after = asv_id)
  
  #### remove chloroplast
  if (rm_chloroplast) {
    taxa <- taxa %>%
      mutate(
        across(c(domain, phylum, class, order, family, genus),
               ~ if_else(class == "Chloroplast", NA_character_, .x))
      )
  }
  #### print outfile
  if(outfile != ""){
    check_dir(outfile, is_file=TRUE)
    write.table(taxa, file = outfile,  row.names = F, sep=sep)
  }
  
  return(taxa)
}
