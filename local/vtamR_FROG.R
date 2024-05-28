install.packages("devtools")
install.packages("roxygen2")
install.packages("seqinr")
install.packages("tidyr")

library("devtools")
library("roxygen2")
library("seqinr")
library("dplyr")
library("tidyr")
library("utils") # to handle zipped files
library("stringr")
#library("Biostrings")


computer <- "Bombyx" # Bombyx/Endoume/Windows
if(computer == "Bombyx"){
  vtam_dir <- "~/vtamR"
  cutadapt_path="/home/meglecz/miniconda3/envs/vtam_2/bin/"
  vsearch_path = ""
  blast_path="~/ncbi-blast-2.11.0+/bin/" # bombyx
  db_path="~/mkLTG/COInr_for_vtam_2022_05_06_dbV5/"
  fastqdir <- "vtamR_test/data/"
  fastqinfo <- "vtamR_test/data/fastqinfo_mfzr_gz.csv"
  outdir <- "vtamR_test/out/"
#  fastqdir <- "/home/meglecz/vtamR_large_files/fastq/"
#  fastqinfo <- "/home/meglecz/vtamR_large_files/user_input/fastqinfo_mfzr.csv"
#  outdir <- "/home/meglecz/vtamR_large_files/out/"
  mock_composition <- "local/user_input/mock_composition_mfzr_eu.csv"
  num_threads=8
  compress = T
} else if (computer == "Endoume"){
  vtam_dir <- "~/vtamR"
  cutadapt_path="/home/emese/miniconda3/bin/"
  vsearch_path = "/home/emese/miniconda3/bin/"
  blast_path= "" # deactivate conda
  db_path= "/home/emese/mkCOInr/COInr/COInr_for_vtam_2023_05_03_dbV5/"
#  fastqdir <- "local/fastq/"
  fastqdir <- "vtamR_test/data/"
  fastqinfo <- "vtamR_test/data/fastqinfo_mfzr_gz.csv"
  outdir <- "vtamR_test/out/"
  num_threads=8
  compress = T
}else if (computer == "Windows"){
  vtam_dir <- "C:/Users/emese/vtamR/"
  cutadapt_path="C:/Users/Public/"
  vsearch_path = "C:/Users/Public/vsearch-2.23.0-win-x86_64/bin/"
  blast_path="C:/Users/Public/blast-2.14.1+/bin/"
  db_path="C:/Users/Public/COInr_for_vtam_2023_05_03_dbV5/"
#  fastqdir <- "C:/Users/emese/vtamR_private/fastq/"
  fastqdir <- "vtamR_test/data/"
  fastqinfo <- "vtamR_test/data/fastqinfo_mfzr_gz.csv"
  outdir <- "vtamR_test/out/"
  num_threads=4
  compress = F
}
sep=";"
setwd(vtam_dir)

taxonomy=paste(db_path, "COInr_for_vtam_taxonomy.tsv", sep="")
blast_db=paste(db_path, "COInr_for_vtam", sep="")


ltg_params_df = data.frame( pid=c(100,97,95,90,85,80),
                            pcov=c(70,70,70,70,70,70),
                            phit=c(70,70,70,70,70,70),
                            taxn=c(1,1,2,3,4,4),
                            seqn=c(1,1,2,3,4,4),
                            refres=c("species","species","species","genus","family","family"),
                            ltgres=c("species","species","species","species", "genus","genus")
)

ltg_params_df = data.frame( pid=c(100,97,95,90,85,80),
                            pcov=c(70,70,70,70,70,70),
                            phit=c(70,70,70,70,70,70),
                            taxn=c(1,1,2,3,4,4),
                            seqn=c(1,1,2,3,4,4),
                            refres=c(8,8,8,7,6,6),
                            ltgres=c(8,8,8,8,7,7)
)



#setwd("D:/vtamR")
# load local packages
load_all(".")
roxygenise() # Builds the help files
usethis::use_roxygen_md() # rebuild the help files ?

setwd("/home/meglecz/ITS_pipeline/FROGS/filter_vtamR/")
fileinfo<- "user_input/fileinfo_its21_bqt1.csv"
############################################################
filter_contaminants <- function(read_count_df, fileinfo="", sep=","){
  
  # select negative samples (fileinfo_df$sample)
  fileinfo_df <- read.csv(fileinfo, header=T, sep=sep) %>%
    select(plate, marker, sample, sample_type) %>%
    filter(sample_type == "negative")
  
  # get the list of asv present in at east one negative control, and their maximum read count among all neg control samples
  max_df <- read_count_df %>%
    group_by(asv) %>%
    summarize("max_read_count" = max(read_count))
  
  df_tmp <- left_join(read_count_df, max_df, by="asv") 
  df_tmp <- df_tmp %>%
    filter(max_read_count == read_count) %>%
    filter(sample %in% fileinfo_df$sample)
  
  
  
}

############################################################
# differs from make_known_occurrences: return a DF with FP, FN, TP
make_known_occurrences2 <- function(read_count_samples_df, fileinfo="", mock_composition="", sep=",", out="", missing_occurrences="", habitat_proportion=0.5){
  
  # read info on samples types and keep only relevant info
  fileinfo_df <- read.csv(fileinfo, header=T, sep=sep) %>%
    select(plate, marker, sample, sample_type, habitat)
  # get unique lines to avoid replicates
  fileinfo_df <- unique(fileinfo_df)
  
  # define data frame for known occurrences
  occurrence_df <- read_count_samples_df
  # add habitat and sample_type to occurrence_df
  occurrence_df <- left_join(occurrence_df, fileinfo_df, by=c("plate", "marker", "sample"))
  # add action column
  occurrence_df$action <- rep(NA, nrow(occurrence_df))
  
  # flag occurrences in negative control samples as delete
  occurrence_df <- flag_from_negative_controls(occurrence_df, fileinfo_df)
  # flag all expected occurrences in mock samples as "keep", NA for tolerate, and delete for all others
  occurrence_df <- flag_from_mock(occurrence_df, mock_composition, fileinfo_df, sep=sep)
  # flag occurrences as delete with low read count in habitat, compared to the others habitats
  occurrence_df <- flag_from_habitat(occurrence_df, fileinfo_df, habitat_proportion=habitat_proportion) 
  
  # keep only relevant columns and lines, sort data
  occurrence_df <- occurrence_df %>%
    select(plate,marker,sample,action,asv) %>%
    filter(!is.na(action)) %>%
    arrange(plate, marker, sample, action)
  # write to outfile
  write.table(occurrence_df, file=out, row.names = F, sep=sep)
  
  # count the number of FP and expected TP
  FP <- nrow(occurrence_df %>%
               filter(action=="delete"))
  TP <- nrow(occurrence_df %>%
               filter(action=="keep"))
  
  # count the number of FN and write missing_occurrences, if filename is defined
  FN <-make_missing_occurrences2(read_count_samples_df, mock_composition=mock_composition, sep=sep, out=missing_occurrences)
  # real TP is the expected occurrences - FN
  TP <- TP - FN
  count_df <- data.frame("TP" = c(TP),
                         "FP" = c(FP),
                         "FN" = c(FN))
  return(count_df)
}


# differs from make_missing_occurrences: return a DF with FP, write missing_occurrences onmy if out !=""
make_missing_occurrences2 <- function(read_count_samples_df, mock_composition="", sep=",", out=""){
  
  # read mock composition to a df
  mock_comp <- read.csv(mock_composition, header=T, sep=sep) %>%
    filter(action=="keep")
  # add mean_read_count to df from read_count_samples_df, and keep only if value is NA
  df <- left_join(mock_comp, read_count_samples_df,  by=c("plate", "marker", "sample", "asv")) %>%
    filter(is.na(mean_read_count)) %>%
    select(-"mean_read_count")
  
  # write to outfile
  if(out != ""){
    write.table(df, file=out, row.names = F, sep=sep)
  }
  
  FN <- nrow(df %>%
               filter(action=="keep"))
  return(FN)
}
############################################################



##############
# Read frogs abundance into and transform it to a long read_count_df format
#############
setwd("/home/meglecz/ITS_pipeline/FROGS/filter_vtamR/")
fileinfo<- "user_input/fileinfo_its21_bqt1.csv"

frogs_abundance_file <- "../Galaxy23-[FROGS_BIOM_to_TSV__abundance_WO_ITSx].tsv"
df <- read.csv(frogs_abundance_file, sep="\t")

read_count_df <- df %>%
  select(-"X.comment",-"rdp_tax_and_bootstrap",-"blast_taxonomy",-"blast_subject",-"blast_perc_identity",-"blast_perc_query_coverage",
         -"blast_evalue",-"blast_aln_length",-"observation_name",-"observation_sum", -seed_id)

read_count_df <- pivot_longer(read_count_df, cols=-"seed_sequence", names_to="samples", values_to="read_count")
read_count_df <- read_count_df %>%
  filter(read_count > 0)
read_count_df$plate <- "Bqt1"
read_count_df$marker <- "ITS21"
read_count_df$samples <- gsub("Bqt1_ITS21_","", read_count_df$samples)
read_count_df$sample <- NA
read_count_df$replicate <- NA

for(i in 1:nrow(read_count_df)){
#  print(read_count_df$samples[i])
  replicat <- str_extract(read_count_df$samples[i], "_[123]$")
  read_count_df$sample[i] <- sub(replicat, "", read_count_df$samples[i])
  read_count_df$replicate[i] <- sub('_', "", replicat)
}
read_count_df <- read_count_df %>%
  select("asv" = seed_sequence,"plate","marker","sample","replicate","read_count")

read_count_df_backup <- read_count_df

############
#Stat frogs (no ITSx) without further filtering
############
read_count_df <- read_count_df_backup
outdir <- "1_out_frogs_raw/"
outdir <- check_dir(dir=outdir)
mock_composition <- "user_input/mock_composition_its21.csv"
# define stat data frame that will be completed with counts after each step
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())

stat_df <- get_stat(read_count_df, stat_df, stage="Frog_no_filter", params="")

###
### PoolReplicates
###
digits = 0
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)

###
### Count FP, FN, TP
###
know_occurrences = paste(outdir, "known_occurrences.csv", sep="")
missing_occurrences = paste(outdir, "missing_occurrences.csv", sep="")
TP_count <- make_known_occurrences2(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=know_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=0.5)

###
### write ASV table and stat file
###
outfile=paste(outdir, "FROGS_asvtable.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)
# write sequence and variant counts after each step
write.csv(stat_df, file = paste(outdir, "count_stat.csv", sep=""))


############
#Stat frogs (no ITSx) with further filtering
############
read_count_df <- read_count_df_backup
outdir <- "2_out_frogs_filter_light/"
outdir <- check_dir(dir=outdir)
# define stat data frame that will be completed with counts after each step
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())
###
### LFN_global_read_count
###
# Eliminate variants with less than global_read_count_cutoff reads in the dataset
global_read_count_cutoff = 10
read_count_df <- LFN_global_read_count(read_count_df, global_read_count_cutoff, write_csv=T, outdir=outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="LFN_global_read_count", params=global_read_count_cutoff)

###
### LFN_filters
###
# LFN_read_count
lfn_read_count_cutoff <- 10
read_count_df <- LFN_read_count(read_count_df, cutoff=lfn_read_count_cutoff, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="LFN_read_count", params=lfn_read_count_cutoff)

###
### keep repeatable occurrences
###
min_replicate_number <- 2
read_count_df <- FilterMinReplicateNumber(read_count_df, min_replicate_number, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterMinReplicateNumber", params=min_replicate_number)

###
### PoolReplicates
###
digits = 0
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)

###
### write ASV table and stat file
###
outfile=paste(outdir, "FROGS_asvtable_min_read_count_FilterMinReplicateNumber.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)
# write sequence and variant counts after each step
write.csv(stat_df, file = paste(outdir, "count_stat.csv", sep=""))

###
### Count FP, FN, TP
###
know_occurrences = paste(outdir, "known_occurrences.csv")
missing_occurrences = paste(outdir, "missing_occurrences.csv")
TP_count <- make_known_occurrences2(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=know_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=0.5)



############
#frogs (no ITSx) find optimal params and stat
############

read_count_df <- read_count_df_backup
outdir <- "3_out_frogs_optimized/"
outdir <- check_dir(dir=outdir)
# define stat data frame that will be completed with counts after each step
stat_df <- data.frame(parameters=character(),
                      asv_count=integer(),
                      read_count=integer(),
                      sample_count=integer(),
                      sample_replicate_count=integer())


###
### PoolReplicates
###
digits = 0
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)

###
### make_known_occurrences2
###
known_occurrences = paste(outdir, "known_occurrences.csv", sep="")
missing_occurrences = paste(outdir, "missing_occurrences.csv", sep="")
TP_count <- make_known_occurrences2(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=known_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=0.5)


# start optimize from almost unfiltered data (filtered only by frogs)
optimize_read_count_df <- read_count_df
###
### OptimizePCRError
###
optimize_dir = paste(outdir, "optimize", sep="")
OptimizePCRError(optimize_read_count_df, mock_composition=mock_composition, sep=sep, outdir=optimize_dir, max_mismatch=2, min_read_count=10)

###
### OptimizeLFNsampleReplicate
###
OptimizeLFNsampleReplicate(optimize_read_count_df, mock_composition=mock_composition, sep=sep, outdir=optimize_dir)


###
### OptimizeLFNReadCountAndLFNvariant
###
lfn_read_count_cutoff=10
lnf_variant_cutoff=0.0001
by_replicate=T

lfn_sample_replicate_cutoff <- 0.0005
pcr_error_var_prop=0.1
max_mismatch=2
sample_prop=0.8
by_sample=T
min_replicate_number=2
OptimizeLFNReadCountAndLFNvariant(optimize_read_count_df, known_occurrences=known_occurrences, sep=sep, outdir=optimize_dir, min_lfn_read_count_cutoff=lfn_read_count_cutoff, min_lnf_variant_cutoff=lnf_variant_cutoff, by_replicate=by_replicate, lfn_sample_replicate_cutoff=lfn_sample_replicate_cutoff, pcr_error_var_prop=pcr_error_var_prop, vsearch_path=vsearch_path, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, min_replicate_number=min_replicate_number)


###
### LFN_filters
###
# LFN_read_count
read_count_df_LFN_read_count <- LFN_read_count(read_count_df, cutoff=lfn_read_count_cutoff, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df_LFN_read_count, stat_df, stage="LFN_read_count", params=lfn_read_count_cutoff)

# LFN_sample_replicate (by column)
read_count_df_lnf_sample_replicate <- LFN_sample_replicate(read_count_df, cutoff=lfn_sample_replicate_cutoff, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df_lnf_sample_replicate, stat_df, stage="LFN_sample_replicate", params=lfn_sample_replicate_cutoff)


# LFN_sample_variant (by line)
read_count_df_lnf_variant <- LFN_variant(read_count_df, cutoff=lnf_variant_cutoff, by_replicate, write_csv=T, outdir = outdir, sep=sep)
param_values <- paste(lnf_variant_cutoff, by_replicate, sep=";")
stat_df <- get_stat(read_count_df_lnf_variant, stat_df, stage="LFN_variant", params=param_values)


# pool the results of the different filterLFN to one data frame; keep only occurrences that passed all filters
read_count_df <- pool_LFN(read_count_df_LFN_read_count, read_count_df_lnf_variant, read_count_df_lnf_sample_replicate, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterLFN")
# delete temporary data frames
read_count_df_lfn_read_count <- NULL
read_count_df_lnf_variant <- NULL
read_count_df_lnf_sample_replicate <- NULL

###
### keep repeatable occurrences
###
min_replicate_number <- 2
read_count_df <- FilterMinReplicateNumber(read_count_df, min_replicate_number, write_csv=T, outdir = outdir, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilterMinReplicateNumber", params=min_replicate_number)

###
### FilerPCRerror
###
pcr_error_var_prop <- 0.1
max_mismatch <- 2
by_sample <- T
sample_prop <- 0.8
read_count_df <- FilterPCRerror(read_count_df, write_csv=T, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, sep=sep)
params <- paste(pcr_error_var_prop, max_mismatch, by_sample, sample_prop, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="FilerPCRerror", params=params)

###
### FilterChimera
###
abskew=2
by_sample = T
sample_prop = 0.8
read_count_df <- FilterChimera(read_count_df, write_csv=T, outdir=outdir, vsearch_path=vsearch_path, by_sample=by_sample, sample_prop=sample_prop, abskew=abskew, sep=sep)
params <- paste(abskew, by_sample, sample_prop, sep=";")
stat_df <- get_stat(read_count_df, stat_df, stage="FilterChimera", params=params)

###
### FilerRenkonen
###
# Renkonen index:
# PS = summ(min(p1i, p2i))
# p1i = number of reads for variant i in replicate 1 / number of reads in replicate 1
renkonen_distance_quantile = 0.9
read_count_df <- FilterRenkonen(read_count_df, write_csv=T, outdir=outdir, renkonen_distance_quantile=renkonen_distance_quantile, sep=sep)
stat_df <- get_stat(read_count_df, stat_df, stage="FilerRenkonen", params=renkonen_distance_quantile)


###
### PoolReplicates
###
digits = 0
read_count_samples_df <- PoolReplicates(read_count_df, digits=digits, write_csv=T, outdir=outdir, sep=sep)

###
### write ASV table and stat file
###
outfile=paste(outdir, "FROGS_asvtable_optimized.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)
# write sequence and variant counts after each step
write.csv(stat_df, file = paste(outdir, "count_stat.csv", sep=""))

###
### Count FP, FN, TP
###
know_occurrences = paste(outdir, "known_occurrences_after.csv", sep="")
missing_occurrences = paste(outdir, "missing_occurrences_after.csv", sep="")
TP_count <- make_known_occurrences2(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=know_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=0.5)


##############
# Read earlier vtam wide output into and transform it to a long read_count_df format
#############
setwd("/home/meglecz/ITS_pipeline/")
fileinfo<- "FROGS/filter_vtamR/user_input/fileinfo_its21_bqt1.csv"
mock_composition <- "FROGS/filter_vtamR/user_input/mock_composition_its21_fw.csv"

outdir <- "vtam_fw/fw_filter2_stringent/"
outdir <- check_dir(dir=outdir)
vtam_asvtable_file <- "vtam_fw/fw_filter2_stringent/asvtable2_stringent.tsv"


"/home/meglecz/novaseq_metadata/novaseq_2023_03/Epi12_Bqt/ITS/asvtables/its_bqt_asv_taxassigned/asvtable2_relaxed_taxa_UNITE.tsv"
"/home/meglecz/novaseq_metadata/novaseq_2023_03/Epi12_Bqt/ITS/asvtables/its_bqt_asv_taxassigned/asvtable2_stringent_bis_taxa_UNITE.tsv"
"/home/meglecz/novaseq_metadata/novaseq_2023_03/Epi12_Bqt/ITS/asvtables/its_bqt_asv_taxassigned/asvtable2_stringent_taxa_UNITE.tsv"


df <- read.csv(vtam_asvtable_file, sep="\t")

#read_count_samples_df <- df %>%
#  select(-read_count, -"keep_MSA1010_1_Bqt1",-"clusterid",-"clustersize",-"chimera_borderline",-"kingdom",-"kingdom_bootstr",-"phylum",-"phylum_bootstr",-"class",-"class_bootstr",-"order",-"order_bootstr",-"family",-"family_bootstr",-"genus",-"genus_bootstr",-"species",-"species_bootstr",-"variant",-"sequence_length")

read_count_samples_df <- df %>%
  select(-variant, -sequence_length, -read_count, -read_count, -"keep_MSA1010_1_Bqt1",-"clusterid",-"clustersize",-"chimera_borderline")


read_count_samples_df <- pivot_longer(read_count_samples_df, cols=-c("sequence", run, marker), names_to="sample", values_to="read_count")
read_count_samples_df <- read_count_samples_df %>%
  filter(read_count > 0)

read_count_samples_df <- read_count_samples_df %>%
  select("asv" = sequence, "plate"=run,"marker","sample", "mean_read_count"=read_count)

read_count_samples_df_backup <- read_count_samples_df

# nb ASV
length(unique(read_count_samples_df$asv))
sum(read_count_samples_df$mean_read_count)

###
### Count FP, FN, TP
###
know_occurrences = paste(outdir, "known_occurrences_after.csv", sep="")
missing_occurrences = paste(outdir, "missing_occurrences_after.csv", sep="")
TP_count <- make_known_occurrences2(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=know_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=0.5)

###
### write ASV table and stat file
###
outfile=paste(outdir, "vtam_asvtable_stringent_vtamR.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)



##############
# Read earlier vtam wide output into and transform it to a long read_count_df format
#############
outdir <- "vtam_fw/fw_filter2_relaxed/"
outdir <- check_dir(dir=outdir)
vtam_asvtable_file <- "vtam_fw/fw_filter2_relaxed/asvtable2_relaxed.tsv"

df <- read.csv(vtam_asvtable_file, sep="\t")

read_count_samples_df <- df %>%
  select(-variant, -sequence_length, -read_count, -read_count, -"keep_MSA1010_1_Bqt1",-"clusterid",-"clustersize",-"chimera_borderline")


read_count_samples_df <- pivot_longer(read_count_samples_df, cols=-c("sequence", run, marker), names_to="sample", values_to="read_count")
read_count_samples_df <- read_count_samples_df %>%
  filter(read_count > 0)

read_count_samples_df <- read_count_samples_df %>%
  select("asv" = sequence, "plate"=run,"marker","sample", "mean_read_count"=read_count)

read_count_samples_df_backup <- read_count_samples_df

# nb ASV
length(unique(read_count_samples_df$asv))
sum(read_count_samples_df$mean_read_count)

###
### Count FP, FN, TP
###
know_occurrences = paste(outdir, "known_occurrences_after.csv", sep="")
missing_occurrences = paste(outdir, "missing_occurrences_after.csv", sep="")
TP_count <- make_known_occurrences2(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=know_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=0.5)

###
### write ASV table and stat file
###
outfile=paste(outdir, "vtam_asvtable_relaxed_vtamR.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)



##############
# Read earlier vtam wide output into and transform it to a long read_count_df format
#############

outdir <- "vtam_fw/fw_filter2_stringent_bis/"
outdir <- check_dir(dir=outdir)
vtam_asvtable_file <- "vtam_fw/fw_filter2_stringent_bis/asvtable2_stringent_bis.tsv"

df <- read.csv(vtam_asvtable_file, sep="\t")

read_count_samples_df <- df %>%
  select(-variant, -sequence_length, -read_count, -read_count, -"keep_MSA1010_1_Bqt1",-"clusterid",-"clustersize",-"chimera_borderline")


read_count_samples_df <- pivot_longer(read_count_samples_df, cols=-c("sequence", run, marker), names_to="sample", values_to="read_count")
read_count_samples_df <- read_count_samples_df %>%
  filter(read_count > 0)

read_count_samples_df <- read_count_samples_df %>%
  select("asv" = sequence, "plate"=run,"marker","sample", "mean_read_count"=read_count)

read_count_samples_df_backup <- read_count_samples_df

# nb ASV
length(unique(read_count_samples_df$asv))
sum(read_count_samples_df$mean_read_count)

###
### Count FP, FN, TP
###
know_occurrences = paste(outdir, "known_occurrences_after.csv", sep="")
missing_occurrences = paste(outdir, "missing_occurrences_after.csv", sep="")
TP_count <- make_known_occurrences2(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=know_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=0.5)

###
### write ASV table and stat file
###
outfile=paste(outdir, "vtam_asvtable_stringent_bis_vtamR.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)



##############
# Read vtam on rv reads
#############
setwd("/home/meglecz/ITS_pipeline/")
fileinfo<- "FROGS/filter_vtamR/user_input/fileinfo_its21_bqt1.csv"
mock_composition <- "FROGS/filter_vtamR/user_input/mock_composition_its21_rv.csv"

outdir <- "vtam_rv/rv_filter2/"
outdir <- check_dir(dir=outdir)
vtam_asvtable_file <- "vtam_rv/rv_filter2/asvtable2.tsv"

df <- read.csv(vtam_asvtable_file, sep="\t")

read_count_samples_df <- df %>%
  select(-variant, -sequence_length, -read_count, -read_count, -"keep_MSA1010_1_Bqt1",-"clusterid",-"clustersize",-"chimera_borderline")


read_count_samples_df <- pivot_longer(read_count_samples_df, cols=-c("sequence", run, marker), names_to="sample", values_to="read_count")
read_count_samples_df <- read_count_samples_df %>%
  filter(read_count > 0)

read_count_samples_df <- read_count_samples_df %>%
  select("asv" = sequence, "plate"=run,"marker","sample", "mean_read_count"=read_count)

read_count_samples_df_backup <- read_count_samples_df

# nb ASV
length(unique(read_count_samples_df$asv))
sum(read_count_samples_df$mean_read_count)

###
### Count FP, FN, TP
###
know_occurrences = paste(outdir, "known_occurrences_after.csv", sep="")
missing_occurrences = paste(outdir, "missing_occurrences_after.csv", sep="")
TP_count <- make_known_occurrences2(read_count_samples_df, fileinfo=fileinfo, mock_composition=mock_composition, sep=sep, out=know_occurrences, missing_occurrences=missing_occurrences, habitat_proportion=0.5)

###
### write ASV table and stat file
###
outfile=paste(outdir, "vtam_asvtable_rv_optimized_vtamR.csv", sep="")
write_asvtable(read_count_samples_df, outfile=outfile, fileinfo=fileinfo, add_empty_samples=T, add_sums_by_sample=T, add_sums_by_asv=T, add_expected_asv=T, mock_composition=mock_composition, sep=sep)

