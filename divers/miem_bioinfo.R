
library(vtamR)
library(dplyr)

library("devtools")
library("roxygen2")
setwd("/home/meglecz/vtamR")
load_all(".")
roxygenise()
usethis::use_roxygen_md()

log = "~/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/vtamR_log_zfzr.csv"
r_versions = "~/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/R_package_versions.csv"

miem_bioinformatics <- function(log, r_versions, outfile = "", sep = ","){
  
  if(is.character(log)){
    # read known occurrences
    log_df <- read.csv(log, header=T, sep=sep)
  }else{
    log_df <- log
  }
  
  if(is.character(r_versions)){
    # read known occurrences
    r_versions_df <- read.csv(r_versions, header=T, sep=sep)
  }else{
    r_versions_df <- r_versions
  }
  ### package versions ##########################################
  rownames(r_versions_df) <- r_versions_df$Package
  R_version <- r_versions_df["base", "Version"]
  vtamR_version <- r_versions_df["vtamR", "Version"]
  rRDP_version <- r_versions_df["rRDP", "Version"]
  rRDPData_version <- r_versions_df["rRDPData", "Version"]
  
  ### Third Party program versions ##########################################
  tpp <- log_df %>%
    select(argument_name, version) %>%
    filter(!is.na(version)) %>%
    distinct()
  
  if(nrow(tpp) != length(unique(tpp$argument_name))){
    tpp <- tpp %>%
      group_by(argument_name) %>%
      summarize(version = last(version))
    msg <- "WARNING: For some third party programs more then one version has been used. Only the last one will be reported in the MIEM file."
    warning(msg)
  }
  tpp <- as.data.frame(tpp)
  rownames(tpp) <- tpp$argument_name
  
  blast_version <- tpp["blast_path", "version"]
  cutadapt_version <- tpp["cutadapt_path", "version"]
  pigz_version <- tpp["pigz_path", "version"]
  swarm_version <- tpp["swarm_path", "version"]
  vsearch_version <- tpp["vsearch_path", "version"]
  
  # set df ##########################################
  fields = c("Database creation: Source of sequences and steps to identify locus of interest",
             "Database creation: Method for sequence curation",
             "Database creation: Link to database or repository",
             "Primer removal (trimming): program, version, parameters",
             "QC program: program, version, parameters",
             "Read pair merging: program, version, parameters",
             "Chimera removal: program, version, parameters",
             "Clustering: OTUs or ASVs (and thresholds)",
             "Additional filtering: removal of singletons or other methods",
             "Additional filtering: decontamination using sequenced controls",
             "Taxonomic assignment method",
             "Taxonomic assignment parameters (and thresholds)",
             "Read normalization: methods)")
  
  report_req <- c("Report",
                  "Report",
                  "Report",
                  "Report",
                  "Report",
                  "Report",
                  "Report",
                  "Report",
                  "If Applicable",
                  "If Applicable",
                  "Report",
                  "Report",
                  "If Applicable")
  
  miem <- data.frame("Step" = fields,
                     "Information" = NA_character_,
                     "Reporting Requirements" = report_req
  )
  rownames(miem) <- fields
  
  
# database #############################################
  # Database creation: Source of sequences and steps to identify locus of interest
  # Database creation: Method for sequence curation
  # Database creation: Link to database or repository
  # Taxonomic assignment method
  # Taxonomic assignment parameters (and thresholds)
  
# LTG --------------------------------------------------
  
  tmp <- log_df %>%
    filter(function_name == "assign_taxonomy_ltg") %>%
    select(function_name, argument_name, value) %>%
    distinct()
    
  if(nrow(tmp) > 0){
    
    if(nrow(tmp) != length(unique(tmp$argument_name))){
      tmp <- tmp %>%
        group_by(argument_name) %>%
        summarize(function_name = last(function_name), value = last(value))
      msg <- "WARNING: The assign_taxonomy_ltg function was run more than once. Only the last one will be reported in the MIEM file."
      warning(msg)
    }
    
    rownames(tmp) <- tmp$argument_name
    
    db <- basename(tmp["blast_db", "value"])
    tax <- basename(tmp["taxonomy", "value"])
    ltg_par <- tmp["ltg_params", "value"]
    
    # Taxonomic assignment method
    msg <- paste0("The LTG method implemented in vtamR R package (v", vtamR_version, ") was used to assign sequences to taxa [Meglécz, 2023](https://link.springer.com/article/10.1007/s42977-024-00201-x).
  This is a BLAST based Lowest Common Ancestor method using NCBI-BLAST (v", blast_version, ").")
    miem["Taxonomic assignment method", "Information"] <- msg
    
    # Taxonomic assignment parameters (and thresholds)
    if(ltg_par == "NULL"){
      msg <- paste0("Default values of the vtamR package (v", vtamR_version, ") assign_taxonomy_ltg function were used.")
    }else{
      msg <- paste0("TO BE COMPLETED: Report the custom ltg_params used")
    }
    miem["Taxonomic assignment parameters (and thresholds)", "Information"] <- msg
    
    # Database creation: Source of sequences and steps to identify locus of interest
    msg <- paste0("The ", db, " database was used with the ", tax, " taxonomy file.  
    TO BE COMPLETED!!! If you have used the COInr database: 
    COInr includes sequences of COI from ncbi_nt and BOLD databases [Meglécz, 2023](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756)")
    miem["Database creation: Source of sequences and steps to identify locus of interest", "Information"] <- msg
    
    # Database creation: Link to database or repository
    msg <- paste0("TO BE COMPLETED!!! 
    If you have used the COInr database, refer to https://zenodo.org/records/20020232 (UPDATE to the version you have used) 
    downloaded by the download_osf function of vtamR (v", vtamR_version,").")
    miem["Database creation: Link to database or repository", "Information"] <- msg
    
    # Database creation: Method for sequence curation
    msg <- paste0("TO BE COMPLETED!!!! 
    If you have used the COInr database: 
    In COInr, sequence redundancy within taxa is eliminated, but the database is not currecated for potential mis-labeling")
    miem["Database creation: Method for sequence curation", "Information"] <- msg
  }
    
# RDP --------------------------------------------------
  tmp <- log_df %>%
    filter(function_name == "assign_taxonomy_rdp") %>%
    select(function_name, argument_name, value) %>%
    distinct()
  
    if(nrow(tmp) > 0){
      
      if(nrow(tmp) != length(unique(tmp$argument_name))){
        tmp <- tmp %>%
          group_by(argument_name) %>%
          summarize(function_name = last(function_name), value = last(value))
        msg <- "WARNING: The assign_taxonomy_rdp function was run more than once. Only the last one will be reported in the MIEM file."
        warning(msg)
      }
      
      rownames(tmp) <- tmp$argument_name
      
      db <- tmp["dir", "value"]
      confidence <- tmp["confidence", "value"]
      chloroplast <- tmp["rm_chloroplast", "value"]
    
      
      # Taxonomic assignment method
      msg <- paste0("The RDP classifier implemented in rRDP package (v", rRDP_version, ";https://bioconductor.org/packages/3.21/bioc/html/rRDP.html) was used to assign sequences to taxa.")
      miem["Taxonomic assignment method", "Information"] <- msg
      
      # Taxonomic assignment parameters (and thresholds)
      msg <- paste0("Default parameters of the rRDP package were used, accepting only assignement with greater or equal then ", confidence, " bootstrap values.")
      if(chloroplast){
        msg <-paste0(msg, "Taxonomic assignments are set to NA when the class was Chloroplast.")
      }
      miem["Taxonomic assignment parameters (and thresholds)", "Information"] <- msg
      
      # Database creation: Source of sequences and steps to identify locus of interest
      if(db == "NULL"){
        msg <- paste0("The rRDPData (v",  rRDPData_version, ";https://bioconductor.org/packages/3.21/data/experiment/html/rRDPData.html) package was used as a database.")
      }else{
        msg <- paste0("TO BE COMPLETED!!!")
      }
      miem["Database creation: Source of sequences and steps to identify locus of interest", "Information"] <- msg
      
      # Database creation: Link to database or repository
      if(db == "NULL"){
        msg <- paste0("The rRDPData v",  rRDPData_version, "; https://bioconductor.org/packages/3.21/data/experiment/html/rRDPData.html)")
      }else{
        msg <- paste0("TO BE COMPLETED!!!")
      }
      miem["Database creation: Link to database or repository", "Information"] <- msg
    }
  
# Primer removal (trimming): program, version, parameters #############################################
  trim_functions <- c("trim_primers", "demultiplex_and_trim", "demultiplex_fastq_pairs")
  
  params <- c("check_reverse",
              "primer_to_end",
              "cutadapt_error_rate",
              "cutadapt_minimum_length",
              "cutadapt_maximum_length",
              "tag_to_end")
  
  tmp <- log_df %>%
    filter(function_name %in% trim_functions) %>%
    filter(argument_name %in% params) %>%
    select(function_name, argument_name, value) %>%
    distinct()
  
  if(nrow(tmp) >0){
    
    txt <- paste(capture.output(print(tmp, row.names = FALSE)), collapse = "\n")
    
    msg <- paste0("Primer/Tag trimming was performed by Cutadapt (v",cutadapt_version,") integrated to the following vtamR (v",vtamR_version,") function(s) and parameters:", "\n", txt)
    miem["Primer removal (trimming): program, version, parameters", "Information"] <- msg
  }
    
  # QC program: program, version, parameters #############################################
    params <- c("fastq_maxee",
                "fastq_minlen",
                "fastq_maxlen",
                "fastq_minmergelen",
                "fastq_maxmergelen",
                "fastq_maxns",
                "fastq_truncqual",
                "fastq_maxdiffs",
                "fastq_minovlen")
    
    tmp <- log_df %>%
      filter(function_name == "merge_fastq_pairs") %>%
      filter(argument_name %in% params) %>%
      select(function_name, argument_name, value) %>%
      distinct()
  if(nrow(tmp) >0){
    txt <- paste(capture.output(print(tmp, row.names = FALSE)), collapse = "\n")
    
    msg <- paste0("Quality filtering was done by the fastq_mergepairs function of vsearch (v", vsearch_version, ") integrated to the merge_fastq_pairs function of vtamR (v",vtamR_version,") with the following parameters:", "\n", txt)
    miem["QC program: program, version, parameters", "Information"] <- msg
    
  # Read pair merging: program, version, parameters #############################################
    params <- c("fastq_maxdiffs",
                "fastq_minovlen",
                "fastq_minlen",
                "fastq_maxlen",
                "fastq_minmergelen",
                "fastq_maxmergelen",
                "fastq_allowmergestagger")
    
    tmp <- log_df %>%
      filter(function_name == "merge_fastq_pairs") %>%
      filter(argument_name %in% params) %>%
      select(function_name, argument_name, value) %>%
      distinct()
    
    txt <- paste(capture.output(print(tmp, row.names = FALSE)), collapse = "\n")
    
    msg <- paste0("Read pair merging was done by the fastq_mergepairs function of vsearch (v", vsearch_version, ") integrated to the merge_fastq_pairs function of vtamR (v",vtamR_version,") with the following parameters:", "\n", txt)
    miem["Read pair merging: program, version, parameters", "Information"] <- msg
  }
  
# Chimera removal: program, version, parameters #############################################
  
  params <- c("abskew",
              "by_sample",
              "filter_occurrence",
              "sample_prop")
  
  tmp <- log_df %>%
    filter(function_name == "filter_chimera") %>%
    filter(argument_name %in% params) %>%
    select(function_name, argument_name, value) %>%
    distinct()
  
  if(nrow(tmp) >0){
    txt <- paste(capture.output(print(tmp, row.names = FALSE)), collapse = "\n")
    
    msg <- paste0("Chimera removal was done by the uchime3_denovo function of vsearch (v", vsearch_version, ") integrated to the filter_chimera function of vtamR (v",vtamR_version,") with the following parameters:", "\n", txt)
    miem["Chimera removal: program, version, parameters", "Information"] <- msg
  }

  # Clustering: OTUs or ASVs (and thresholds) #############################################
  
  params <- c("method",
              "by_sample",
              "swarm_d",
              "fastidious",
              "identity")
  
  tmp <- log_df %>%
    filter(function_name == "cluster_asv") %>%
    filter(argument_name %in% params) %>%
    select(function_name, argument_name, value) %>%
    distinct()
  if(nrow(tmp) >0){
  
    if(nrow(tmp) != length(unique(tmp$argument_name))){
      tmp <- tmp %>%
        group_by(argument_name) %>%
        summarize(function_name = last(function_name), value = last(value))
  
      msg <- "WARNING: The cluster_asv was run more than once. Only the last one will be reported in the MIEM file."
      warning(msg)
    }
  
    rownames(tmp) <- tmp$argument_name
    cluster_method <- tmp["method", "value"]
    
    if(cluster_method == "vsearch"){
      tmp <- tmp %>%
        filter(argument_name %in% c("by_sample", "identity"))
      
      txt <- paste(capture.output(print(tmp, row.names = FALSE)), collapse = "\n")
      
      msg <- paste0("Clustering ASV to mOTUs was done by the cluster_size function of vsearch (v", vsearch_version, ") integrated to the cluster_asv function of vtamR (v",vtamR_version,") with the following parameters:", "\n", txt)
      miem["Clustering: OTUs or ASVs (and thresholds)", "Information"] <- msg
    }else{
      tmp <- tmp %>%
        filter(argument_name %in% c("by_sample", "swarm_d", "fastidious"))
        
      txt <- paste(capture.output(print(tmp, row.names = FALSE)), collapse = "\n")
      
      msg <- paste0("Clustering ASV to mOTUs was done by swarm (v", swarm_version, ") integrated to the cluster_asv function of vtamR (v",vtamR_version,") with the following parameters:", "\n", txt)
      miem["Clustering: OTUs or ASVs (and thresholds)", "Information"] <- msg
    }
  }
  
# Additional filtering: removal of singletons or other methods #############################################
  
  functions <- c("denoise_by_swarm",
                 "filter_asv_global",
                 "filter_stop_codon",
                 "filter_indel",
                 "filter_min_replicate",
                 "filter_replicate")
  
  delete_params <- c("read_count",
                     "swarm_path",
                     "vsearch_path",
                     "pigz_path",
                     "num_threads",
                     "outfile",
                     "sep",
                     "quiet",
                     "log_file",
                     "sampleinfo",
                     "conta_file")
  
  tmp <- log_df %>%
    filter(function_name %in% functions) %>%
    filter(!(argument_name %in% delete_params)) %>%
    select(function_name, argument_name, value) 
  
  if(nrow(tmp) >0){
    
    txt <- paste(capture.output(print(tmp, row.names = FALSE)), collapse = "\n")
    
    msg <- paste0("The additional filtering were done by the following vtamR (v",vtamR_version,") functions, 
   using with the parameters below. The denoise_by_swarm function uses swarm (v", swarm_version, ").
   For detailed information on the order of filtering, see the log_file provided in the Supplementary Materials.", "\n", txt)
    miem["Additional filtering: removal of singletons or other methods", "Information"] <- msg
  }
    
  
# Additional filtering: decontamination using sequenced controls #############################################
  
  functions <- c("filter_contaminant",
                 "filter_pcr_error",
                 "filter_occurrence_read_count",
                 "filter_occurrence_sample",
                 "filter_occurrence_variant",
                 "SuggestFilterParametersPCRError",
                 "SuggestFilterParametersReadCountVariant",
                 "SuggestFilterParametersSample")
  
  delete_params <- c("read_count",
                     "swarm_path",
                     "vsearch_path",
                     "pigz_path",
                     "num_threads",
                     "outfile",
                     "sep",
                     "quiet",
                     "log_file",
                     "sampleinfo",
                     "conta_file")
  
  tmp <- log_df %>%
    filter(function_name %in% functions) %>%
    filter(!(argument_name %in% delete_params)) %>%
    select(function_name, argument_name, value) 
  
  if(nrow(tmp) >0){
    
    txt <- paste(capture.output(print(tmp, row.names = FALSE)), collapse = "\n")
    
    msg <- paste0("The parameters for below mentionned 'filter_xxx' functions were chosen based on output of 'suggest_xxx' functions. 
    These functions suggest parameter values based on the composition of control samples to minimise fanse positives and false negatives. 
    All of them are implemented in vtamR (v",vtamR_version,"). The filter_pcr_error function uses the usearch_global function of vsearch (v", vsearch_version, ").
    For detailed information on the order of filtering, see the log_file provided in the Supplementary Materials.", "\n", txt)
    miem["Additional filtering: decontamination using sequenced controls", "Information"] <- msg
  }
  
    #### write csv
    if(outfile != ""){
      check_dir(outfile, is_file=TRUE)
      write.table(miem, file = outfile,  row.names = F, sep=sep)
    }
    invisible(miem)
}
  
  log = "~/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/vtamR_log_zfzr.csv"
  r_versions = "~/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/R_package_versions.csv"
  outfile = "~/vtamR/vignettes/vtamR_demo_out/zfzr_plate1/MIEM_bioinformatics.csv"
  miem_bioinformatics(log = log, r_versions = r_versions, outfile = outfile)
  
  log = "/home/meglecz/vtamR_benchmark/pipeline_output_11/shark/vtamR/vtamR_log.csv"
  r_versions = "/home/meglecz/vtamR_benchmark/pipeline_output_11/shark/vtamR/R_package_versions.csv"
  outfile = "/home/meglecz/vtamR_benchmark/pipeline_output_11/shark/vtamR/MIEM_bioinformatics.csv"
  miem_bioinformatics(log = log, r_versions = r_versions, outfile = outfile)
  
  colnames(log_df)
  colnames(r_versions_df)
  #sessionInfo()
  
  