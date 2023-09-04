#' run_blast
#' 
#' Run BLAST
#'  
#' @param seqs vector containing the query sequences
#' @param blast_db BLAST DB incuding path
#' @param blast_path path to BLAST executables
#' @param outdir output directory
#' @param qcov_hsp_perc minimum query coverage
#' @param perc_identity minimum percentage of identity
#' @param num_threads number of threads
#' @export
#'
run_blast <- function(df, blast_db="", blast_path="", outdir="", qcov_hsp_perc=70, perc_identity=70, num_threads=1){
 
   # make fasta file with unique reads; use numbers as ids  
  seqs <- unique(df$asv)
  fas <- paste(outdir, 'unique.fas', sep="")
  write_fasta(seqs, fas, seq_as_id=F)
  # define the name of the output file
  blast_out <- paste(outdir, 'blast.out', sep="")
  
  task = "megablast"
  e = 1e-20
  dust = "yes"
  max_target_seqs=500
  
  
  blast <- paste(blast_path, "blastn -task ", task, " -db ",blast_db ," -query ",fas," -evalue ",e," -out ",blast_out," -outfmt '6 qseqid pident qcovhsp staxids' -dust ",dust," -qcov_hsp_perc ",qcov_hsp_perc," -perc_identity ",perc_identity," -num_threads ",num_threads," -max_target_seqs ",max_target_seqs, sep="")
  system(blast)
  
  # read BLAST results 
  blast_res <- read.delim(blast_out, header=F, sep="\t", fill=T, quote="")
  colnames(blast_res) <- c("qseqid","pident","qcovhsp","staxids") 
  return(blast_res)
}

#' delete_1_by_row
#' 
#' Delete all 1 from the beginning of a raw , shift the other values and replace de missing ones at the end by NA
#'  
#' @param row vector
#'
delete_1_by_row <- function(row) {
  
  n <- length(row) 
  # Remove all occurrences of 1
  row <- row[row != 1]
  
  # Create a new row with NA at the end
  new_row <- c(row, rep(NA, n - length(row)))
  
  return(new_row)
}
#' get_lineage_ids
#' 
#' Get the complete taxonomic lineage of each taxID in the input vector; return taxIDs of each taxa in the lineage
#'  
#' @param taxids vector taxIDs (taxonomic IDs)
#' @param tax_df data frame with the following columns: tax_id, parent_tax_id, rank, name_txt, taxlevel (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 2: kingdom, 1: superkingdom, 0: root)
#' @export
#'
get_lineage_ids <- function(taxids, tax_df){
  
  # taxids is a vector of taxids; there can be duplicated values
  lineages <- as.data.frame(taxids)
  colnames(lineages) <- c("tax_id")
  
  # if input is dataframe with blast output
  #lineages <- df_intern %>%
  #  select(staxids) %>%
  #  rename(tax_id=staxids)
  
  i <- 1
  while(i < 100){
    # use i as name instead of tax_id
    new_colname <- as.character(i)
    # add parent_tax_id and rename columns
    lineages <- left_join(lineages, tax_df, by="tax_id")%>%
      select(-rank, -name_txt, -taxlevel) %>%
      # !! = interprent the variable
      rename(!!new_colname :=tax_id, "tax_id"=parent_tax_id)
    
    i <- i+1
    # stop if all lines has the same value (usually 1)
    tid_list <- unique(lineages$tax_id)
    if(length(tid_list) == 1 && tid_list[1] ==1){
      break
    }
  }
  # delete the last column, where all values are 1
  lineages <- lineages %>%
    select(-tax_id)
  # reverse order of columns
  lineages <- lineages[, ncol(lineages):1]
  # Apply the function to each row of the lineages data frame: 
  # delete all 1, shift the remaining elements of each row to the beginning, 
  # and replace missing values at the end of the row bu NA
  lineages <- as.data.frame(t(apply(lineages, 1, delete_1_by_row)))
  # add as a first column the taxid, so thay can easily me accessed
  lineages <- cbind(taxids, lineages)
  
  return(lineages)
}
#' make_ltg
#' 
#' Determine the Lowest Taxonomic Group (LTG) that contains phit percentage of the input taxids
#'  
#' @param taxids vector taxIDs (taxonomic IDs); there can be duplicated values
#' @param lineages data frame: taxids in the first column followed by their lineages represented taxids (starting at the lowest resolution)
#' @export
#'
make_ltg <- function(taxids, lineages, phit=70){
  # taxids is a vector of taxids; there can be duplicated values
  # mak a data frame from the vector
  lin <- as.data.frame(taxids)
  colnames(lin) <- c("staxid")
  
  # add lineage to each taxid
  lin <- left_join(lin, lineages, by=c("staxid" = "taxids")) %>%
    select(-where(~all(is.na(.)))) # delete columns if all elements are NA
  
  ltg <- NA
  if(length(unique(lin$staxid)) == 1){ # only one taxid among hits; avoid loops
    ltg <-lin$staxid[1]
  }else{
    for(i in 2:ncol(lineages)){ # start from low resolution
      tmp <- as.data.frame(lineages[,i])
      colnames(tmp) <- c("tax_id")
      # get unique taxids, and their numbers in the i-th column
      tmp <- tmp %>%
        group_by(tax_id) %>%
        summarize(nhit=length(tax_id)) %>%
        arrange(desc(nhit))
      
      # if the taxid with the highest number of sequences does not contain at least phit percent of the hits, stop
      max_hitn <- as.numeric(tmp[1,"nhit"])
      if(max_hitn/sum(tmp[,"nhit"]) < phit/100){
        break
      }
      ltg <- as.numeric(tmp[1,"tax_id"])
    }
  }
  return(ltg)
}

#' update_taxids
#' 
#' Replaces old taxids by valid ones
#'  
#' @param df data frame with the following columns: qseqid,pident,qcovhsp,staxids
#' @param old_taxid data frame with the following columns:  tax_id,old_tax_id
#' @export
#'

update_taxids <- function(df, old_taxid){
  # df is a data frame with the following columns: qseqid,pident,qcovhsp,staxids
  # old_taxid is a data frame with the following columns:  tax_id,old_tax_id
  
  # replace old taxids (if any) in df by up to date ones 
  df <- left_join(df, old_taxid, by=c("staxids" = "old_tax_id"))
  df$staxids[which(!is.na(df$tax_id))] <- df$tax_id[which(!is.na(df$tax_id))]
  # delete tax_id column since the values (if non NA were used to replace staxids)
  df <- df %>%
    select(-tax_id)
  return(df)
}
