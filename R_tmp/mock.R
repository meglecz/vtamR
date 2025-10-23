
################################################"
setwd("C:/Users/emese/vtamR")
blast_path <- "C:/Users/Public/blast-2.16.0+/bin/blastn"
num_threads <- 4
outdir <- "C:/Users/emese/mock"


library("dplyr")
library("devtools")
library("roxygen2")
load_all(".")
roxygenise()
usethis::use_roxygen_md()


 

####  make Taxid File by hand with columns: seqid  taxid
# j'attache le fichier taxID.

### make BLAST DB
# Run this just once, since it is the same DB for the 2 markers
############################
run_system2 <- function(path, args, quiet = FALSE) {
  
  # system2 cannot use ~ as a home
  path <- path.expand(path)
  
  if (!quiet) {
    # Show the full command that will be run
    cat("Running command:\n")
    cat(path, paste(shQuote(args), collapse = " "), "\n")
    
    system2(
      command = path,
      args = args,
      stdout = "",
      stderr = ""
    )
    
  } else {
    output <- suppressWarnings(system2(
      command = path,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    ))
    
    # Extract only error/warning/fail lines
    errors_only <- grep("error|warning|fail", output, ignore.case = TRUE, value = TRUE)
    
    # Remove lines containing "Mean expected error"
    errors_only <- grep("Mean expected error|Mean observed errors", errors_only, ignore.case = TRUE, value = TRUE, invert = TRUE)
    
    if (length(errors_only) > 0) {
      cat(errors_only, sep = "\n")
    }
  }
}
############################

taxids <- "C:/data/Diane/data/TaxID_blast.tsv"
fas <- "C:/data/Diane/data/Fasta_seq_Tpos_19sp_AgroBat1.txt"
blastdb_path = blast_path
blastdb_path <- sub("blastn$", "makeblastdb", blastdb_path)
bdmock <- "C:/data/Diane/data/db_mock/db"

args = c(
  "-dbtype", "nucl",
  "-in", fas,
  "-taxid_map", taxids,
  "-out", bdmock,
  "-parse_seqids"
)
run_system2(blastdb_path, args, quiet=FALSE)
####

### Get list of mocks
# Run this just once, since it is the same list for the 2 markers
mocks <- read.csv("C:/data/Diane/data/mock_csv_DZL_keep_FwR5_ab01.csv")
mocks <- mocks %>%
  select(sample) %>%
  distinct()

### Select sequences in mock samples
# From here, run everything once for each marker
#read_count_df <- read.csv("C:/data/Diane/data/7_FilterRenkonen_FwR5_ab01.csv")
read_count_df <- read.csv("C:/data/Diane/data/7_FilterRenkonen_zfra_ab01.csv")
# select mock samples
read_count_df <- read_count_df %>%
  filter(sample %in% mocks$sample)

#### taxassign
ltg_params_df = data.frame( pid=c(100,97,95),
                            pcov=c(70,70,70),
                            phit=c(70,70,70),
                            taxn=c(1,1,2),
                            seqn=c(1,1,2),
                            refres=c(8,8,8),
                            ltgres=c(8,8,8))

#outfile = "C:/data/Diane/FwR5.tsv"
outfile = "C:/data/Diane/zfra.tsv"
taxa <- TaxAssign(asv= read_count_df, 
          taxonomy = "C:/Users/Public/COInr_for_vtam_2025_05_23_dbV5/COInr_for_vtam_taxonomy.tsv", 
          blast_db = bdmock, 
          blast_path=blast_path, 
          ltg_params=ltg_params_df, 
          outfile=outfile, 
          quiet=F, 
          fill_lineage=TRUE)

### filter and organize taxa
# For each ASV, get the sum of the read counts in the mocks, and the number of mocks, where the ASV is present.
df <- read_count_df %>%
  group_by(asv_id) %>%
  summarize(total_read_count_mock = sum(read_count), number_mock=n_distinct(sample))
# complete taxa with read and sample count info
taxa <- left_join(taxa, df, by="asv_id") 
# select columns and sort the lines by species and decreasing read counts
taxa <- taxa %>%
  filter(!is.na(ltg_name)) %>%
  select(asv_id, total_read_count_mock, number_mock, species, genus, family, order, class, phylum, pid, asv) %>%
  arrange(species, desc(total_read_count_mock))
# get the first line for each species (the one with the highest number of reads.)
taxa_select <- taxa %>%
  group_by(species) %>%
  slice_head(n = 1) %>%
  ungroup()

#write.csv(taxa_select, file="C:/data/Diane/FwR5_mock.csv")
write.csv(taxa_select, file="C:/data/Diane/zfra_mock.csv")


