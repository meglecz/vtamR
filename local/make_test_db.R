library("devtools")
library("roxygen2")
library("seqinr") # splitseq for FilterCodonStop
library("dplyr")
library("tidyr") # gather for read_asv_table; pivot_wider in WriteAsVtable and stat_sample !!sym
#library("utils") # to handle zipped files
library("ggplot2") 

# input
blast_out <- "/home/meglecz/vtamR/local/blast.out"
taxonomy <- "~/mkLTG/COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam_taxonomy.tsv"
new_taxonomy_file <- "/home/meglecz/vtamR/vtamR_test/data/db_test/taxonomy_reduced.tsv"
COInr <- "/home/meglecz/mkCOInr/COInr/COInr_2023_05_03/COInr.tsv"
blast_path <- "~/ncbi-blast-2.11.0+/bin/"

# output
out_db_name <- "/home/meglecz/vtamR/vtamR_test/data/db_test/COInr_reduced"
taxids <-"/home/meglecz/vtamR/vtamR_test/data/db_test/COInr_reduced_taxids.tsv"
new_db_fas <-"/home/meglecz/vtamR/vtamR_test/data/db_test/COInr_reduced.fas"


#### Read taxonomy info 
# read taxonomy file; quote="" is important, since some of the taxon names have quotes and this should be ignored
tax_df <- read.delim(taxonomy, header=T, sep="\t", fill=T, quote="")
# make data frame with old taxids as line numbers and taxids in a columns
old_taxid <- tax_df %>%
  filter(!is.na(old_tax_id)) %>%
  select(tax_id, old_tax_id)
# delete old_tax_ids from tax_df and make taxids unique
tax_df <- tax_df %>%
  select(-old_tax_id)
tax_df <- unique(tax_df)
####

### read read results to data frame (blast_res columns: "qseqid","pident","qcovhsp","staxids"); takes into account multiple taxids in one cell
blast_res <- read_blast_res(file=blast_out)
# add update old taxids to valid ones
blast_res <- update_taxids(blast_res, old_taxid)
# get the list of the taxids of the sequences (WO the taxids of the lineages)
unique_lowest_resolution_taxids <- unique(blast_res$staxids)

### make a lineage for each taxid in blastres
lineages <- get_lineage_ids(unique(blast_res$staxids), tax_df)
# get all taxid in the blast res and all taxids in their lineages
unique_taxids <- as.data.frame(unique(na.omit(unlist(lineages))))
colnames(unique_taxids) <- c("tax_id")

new_tax_df <- left_join(unique_taxids, tax_df, by="tax_id")
new_tax_df$old_tax_id <- ""
write.table(new_tax_df, new_taxonomy_file, sep="\t", row.names = F)
# save memory
rm(df)
rm(tax_df)
rm(lineages)
rm(blast_res)
rm(old_taxid)
rm(new_tax_df)

### select sequneces of the taxa in unique_lowest_resolution_taxids
seq_ref <- read.table(COInr, sep="\t", header=T)
seq_ref <- seq_ref %>%
  filter(taxID %in% unique_lowest_resolution_taxids)
seq_ref$taxID <- as.numeric(seq_ref$taxID)

# tsv with seq id taxid
taxids_df <- seq_ref %>%
  select(seqID, taxID)
write.table(taxids_df, file=taxids, sep="\t", row.names = F, col.names=FALSE, quote=FALSE)

# fasta with seq id and seq
file_connection <- file(new_db_fas, "w")
seq_ref$seqID <- paste('>', seq_ref$seqID, sep="")
writeLines(paste(seq_ref$seqID, seq_ref$sequence, sep="\n"), con=file_connection, sep="\n")
close(file_connection)

# format blast DB
makeblastdb <- paste(blast_path, 'makeblastdb -dbtype nucl -in ', new_db_fas, ' -parse_seqids -taxid_map ', taxids,' -out ', out_db_name, sep="")
system(makeblastdb)


