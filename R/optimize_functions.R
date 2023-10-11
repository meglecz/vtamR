#' make_known_occurrences
#' 
#' Prepare a file that list all occurrences that are clearly a TP (expected variants in mock) or FP (unexpected variants in mocks, all variants in negative control samples, variants present in a wrong habitat)
#'  
#' @param read_count_samples_df data frame with the following variables: asv, plate, marker, sample, read_count
#' @param fileinfo csv file with columns: plate, marker, sample, sample_type(mock/negative/real), habitat, replicate, (optional: file)
#' @param mock_composition csv file with columns: plate, marker, sample, action (keep/tolerate), asv
#' @param sep separator used in fileinfo and mock_composition
#' @param out name of the output file containing known occurrences (TP in mock, and FP)
#' @param missing_occurrences name of the output file containing the missing occurrences (FN); file is written only the name has been specified
#' @param habitat_proportion for each asv, if the proportion of reads in a habitat is below this cutoff, it is considered as an artifact in all samples of the habitat
#' @export
#' 

make_known_occurrences <- function(read_count_samples_df, fileinfo="", mock_composition="", sep=",", out="", missing_occurrences="", habitat_proportion=0.5){
  
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
  # flag occurrences as keep with low read count in habitat, compared to the others habitats
  occurrence_df <- flag_from_habitat(occurrence_df, fileinfo_df, habitat_proportion=habitat_proportion) 
  
  # keep only relevant columns and lines, sort data
  occurrence_df <- occurrence_df %>%
    select(plate,marker,sample,action,asv) %>%
    filter(!is.na(action)) %>%
    arrange(plate, marker, sample, action)
  # write to outfile
  write.table(occurrence_df, file=out, row.names = F, sep=sep)
  
  if(missing_occurrences != ""){
    make_missing_occurrences(read_count_samples_df, mock_composition=mock_composition, sep=sep, out=missing_occurrences)
  }
  
}

#' flag_from_mock
#' 
#' Flag all occurrences mock samples. 
#' Expected variants as 'keep', unexpected ASVs as 'delete', tolerate ASVs as NA. 
#' Those are ASVs that can be in the mock, but the filtering should not be optimized to keep them. 
#' (e.g. badly amplified species present in the mock)
#'  
#' @param occurrence_df data frame with the following variables: asv, plate, marker, sample, mean_read_count, sample_type, habitat, action  
#' @param fileinfo_df csv file with columns: plate, marker, sample, sample_type(mock/negative/real), habitat
#' @param mock_composition csv file with columns: plate, marker, sample, action (keep/tolerate), asv
#' @param sep separator used in mock_composition
#' @export
#' 
flag_from_mock <- function(occurrence_df, mock_composition="", fileinfo_df, sep=","){
  # read mock composition
  mock_composition_df <- read.csv(mock_composition, header=T, sep=sep) %>%
    rename(action_mock=action)
  # add action_mock to occurrence_df; use full join, so expected ASV missing from occurrence_df will be added
  occurrence_df <- full_join(occurrence_df, mock_composition_df, by=c("plate", "marker", "sample", "asv"))
  # so expected ASV was missing from occurrence_df, complete the sample_type as mock
  occurrence_df$sample_type[which(is.na(occurrence_df$sample_type) & occurrence_df$action_mock=="keep")] <- "mock"
  # set the action to delete, keep of tolerate in function of the mock composition
  occurrence_df$action[which((is.na(occurrence_df$action)) & occurrence_df$sample_type == "mock")] <- "delete"
  occurrence_df$action[which(occurrence_df$action_mock =="keep")] <- "keep"
  occurrence_df$action[which(occurrence_df$action_mock =="tolerate")] <- NA
  # select original columns
  occurrence_df <- occurrence_df %>%
    select(asv, plate, marker, sample, mean_read_count, sample_type, habitat, action)
  
  return(occurrence_df)
}


#' flag_from_negative_controls
#' 
#' flag all occurrences in all negative controls as "delete" (FP)
#'  
#' @param occurrence_df data frame with the following variables: asv, plate, marker, sample, mean_read_count, sample_type, habitat, action  
#' @param fileinfo_df csv file with columns: plate, marker, sample, sample_type(mock/negative/real), habitat
#' @export
#' 
flag_from_negative_controls <- function(occurrence_df, fileinfo_df){
  
  occurrence_df$action[which(occurrence_df$sample_type=="negative")] <- "delete"
  return(occurrence_df)
}



#' flag_from_habitat
#' 
#' Flag FP occurrences in samples based on the habitat the samples are from.
#' All ASVs present in more than one habitat are checked. 
#' For each of these ASVs, if the proportion of reads in a habitat is below a cutoff (habitat_proportion), 
#' it is considered as an artifact in all samples of the habitat.
#'  
#' @param occurrence_df data frame with the following variables: asv, plate, marker, sample, mean_read_count, sample_type, habitat, action  
#' @param fileinfo_df csv file with columns: plate, marker, sample, sample_type(mock/negative/real), habitat
#' @param habitat_proportion For each of ASVs, if the proportion of reads in a habitat is below this cutoff it is considered as an artifact in all samples of the habitat.
#' @export
#' 
flag_from_habitat <- function(occurrence_df, fileinfo_df, habitat_proportion=0.5){
  
  # group by asv and habitat and count the total number of reads for each habitat-asv combination
  tmp <- occurrence_df %>%
    group_by(habitat, asv) %>%
    summarize(habitat_read_count=sum(mean_read_count), .groups="drop_last") %>%
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
    select(-hab_action)
  
  return(occurrence_df)
}

#' make_missing_occurrences
#' 
#' Prepare a file that list all expected occurrences that are missing (False negatives)
#'  
#' @param read_count_samples_df data frame with the following variables: asv, plate, marker, sample, read_count
#' @param mock_composition csv file with columns: plate, marker, sample, action (keep/tolerate), asv
#' @param sep separator used in fileinfo and mock_composition
#' @param out name of the output file
#' @export
#'
make_missing_occurrences <- function(read_count_samples_df, mock_composition="", sep=",", out=""){
  
  # read mock composition to a df
  mock_comp <- read.csv(mock_composition, header=T, sep=sep) %>%
    filter(action=="keep")
  # add mean_read_count to df from read_count_samples_df, and keep only if value is NA
  df <- left_join(mock_comp, read_count_samples_df,  by=c("plate", "marker", "sample", "asv")) %>%
    filter(is.na(mean_read_count)) %>%
    select(-"mean_read_count")
  
  # write to outfile
  write.table(df, file=out, row.names = F, sep=sep)
}

#' OptimizePCRError
#' 
#' Prepare a file that lists all pairs of expected and unexpected ASVs in mock samples with a single difference between them,
#'  and their read counts. The pcr_error_var_prop parameter should be above the highest pcr_error_var_prop value in the table.
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param mock_composition csv file with columns: plate, marker, sample, action (keep/tolerate), asv
#' @param sep separator used in mock_composition
#' @param outdir name of the output directory
#' @param min_read_count occurrences under this read_count limits are ignored
#' @export
#'

OptimizePCRError <- function(read_count_df, mock_composition="", sep=",", outdir="", min_read_count=2){
  
  # check outdir and make tmp dir
  outdir <- check_dir(outdir)
  out = paste(outdir, "OptimizePCRError.csv", sep="")
  outdir_tmp <- paste(outdir, 'tmp_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  outdir_tmp <- check_dir(outdir_tmp)
  
  # read the mock composition file and keep only lines with keep and tolerate
  mock_composition_df <- read.csv(mock_composition, header=T, sep=sep) %>%
    filter(action=="keep" | action=="tolerate")
  # add pms column
  mock_composition_df$pms <- paste(mock_composition_df$plate, mock_composition_df$marker, mock_composition_df$sample, sep=".")
  # get unique liste of pms
  unique_mock_list <- unique(mock_composition_df$pms)
  
  # sum read_counts over replicates 
  df <- read_count_df %>%
    group_by(plate, marker, sample, asv) %>%
    summarize(read_count_sample=sum(read_count), .groups="drop_last") %>%
    filter(read_count_sample >=min_read_count)
  # add pms to df and keep only mock samples
  df$pms <- paste(df$plate, df$marker, df$sample, sep=".")
  df <- df %>%
    filter(pms %in% unique_mock_list)
  
  # define an empty dataframe for the output
  asv_pairs <- data.frame(
    plate = character(),
    marker= character(),
    sample= character(),
    expected_read_count= numeric(),
    unexpected_read_count= numeric(),
    pcr_error_var_prop= numeric(),
    expected_asv= character(),
    unexpected_asv= character())
  ###
  # loop over all mock samples
  ###
  for(mock in unique_mock_list){
    # get the list of keep ASV in the given mock sample from mock_composition
    tmp_mock <- mock_composition_df %>%
      filter(pms==mock) %>%
      filter(action=="keep")
    asv_list_keep <- unique(tmp_mock$asv)
    # make fasta file with unique mock variants; use sequences as ids
    fas_keep <- paste(outdir_tmp, mock, "_keepASV.fas", sep="")
    write_fasta(asv_list_keep, fas_keep, seq_as_id=T)
    
    # get the list of tolerate ASV in the given mock sample from mock_composition
    tmp_mock <- mock_composition_df %>%
      filter(pms==mock) %>%
      filter(action=="tolerate")
    asv_list_tolerate <- unique(tmp_mock$asv)
    
    # get list of ASVs present in the mock sample in read_count_df that they and neither keep nor tolerate 
    tmp <- df %>%
      filter(pms==mock) %>%
      filter(!(asv %in% asv_list_keep)) %>%
      filter(!(asv %in% asv_list_tolerate))
    asv_list_delete <- unique(tmp$asv)   
    # make fasta file with unique variants that are neither keep nor tolerate in mock; use sequences as ids
    fas_delete <- paste(outdir_tmp, mock, "_deleteASV.fas", sep="")
    write_fasta(asv_list_delete, fas_delete, seq_as_id=T)
    
    
    # vsearch --usearch_global to find highly similar sequence pairs
    vsearch_out <- paste(outdir_tmp, mock, '_vsearch_out.out', sep="")
#    vsearch <- paste(vsearch_path, "vsearch --usearch_global ", fas_delete, " --db ", fas_keep, " --quiet --iddef 1 --self --id 0.90 --maxaccepts 0 --maxrejects 0 --userfields 'query+target+ids+aln' --userout ", vsearch_out, sep="")
    vsearch <- paste(vsearch_path, "vsearch --usearch_global ", fas_delete, " --db ", fas_keep, ' --quiet --iddef 1 --self --id 0.90 --maxaccepts 0 --maxrejects 0 --userfields "query+target+ids+aln" --userout ', vsearch_out, sep="")
    #https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/system
    system(vsearch)
    
    
    if(file.size(vsearch_out) == 0){
      # Delete the temp directory
      unlink(outdir_tmp, recursive = TRUE)
    }else{
      # read vsearch results
      results_vsearch<- read.csv(vsearch_out, header = FALSE, sep="\t")
      colnames(results_vsearch) <- c("query","target","nb_ids","aln")
      # none of the values easily outputted by vsearch take into the external gaps as a diff => correct this, based on the alnlen and the number of identities
      results_vsearch$nb_diff <- nchar(results_vsearch$aln) - results_vsearch$nb_ids
      # keep only pairs with 1 difference 
      results_vsearch <- results_vsearch %>%
        filter(nb_diff == 1)
      # delete unnecessary columns and add plate, marker, sample
      results_vsearch <- select(results_vsearch, -c(nb_ids, aln, nb_diff))
      pms_mock <- strsplit(mock, "[.]")
      results_vsearch$plate <- pms_mock[[1]][1]
      results_vsearch$marker <- pms_mock[[1]][2]
      results_vsearch$sample <- pms_mock[[1]][3]
      
      # add read_counts to results_vsearch
      results_vsearch <- rename(results_vsearch, asv = target)
      results_vsearch <- left_join(results_vsearch, df, by=c("plate", "marker", "sample", "asv"))
      results_vsearch <- rename(results_vsearch, expected_asv = asv)
      results_vsearch <- rename(results_vsearch, expected_read_count = read_count_sample)
      results_vsearch <- rename(results_vsearch, asv = query)
      results_vsearch <- left_join(results_vsearch, df, by=c("plate", "marker", "sample", "asv"))
      results_vsearch <- rename(results_vsearch, unexpected_asv = asv)
      results_vsearch <- rename(results_vsearch, unexpected_read_count = read_count_sample)
      results_vsearch$pcr_error_var_prop <- results_vsearch$unexpected_read_count / results_vsearch$expected_read_count
      results_vsearch <- results_vsearch %>%
        arrange(desc(pcr_error_var_prop)) %>%
        select(plate, marker, sample, expected_read_count, unexpected_read_count, pcr_error_var_prop, expected_asv, unexpected_asv)
      # append results to existing asv_pairs
      asv_pairs <- rbind(asv_pairs, results_vsearch)
    }
  }
  ###
  # end loop 
  ### 
  
  asv_pairs <- asv_pairs %>%
    arrange(desc(pcr_error_var_prop))
  
  # Delete the temp directory
  unlink(outdir_tmp, recursive = TRUE)
  
  write.table(asv_pairs, file=out, sep=sep, row.names = F)
}

#' OptimizeLFNsampleReplicate
#' 
#' Prepare a file that lists all expected occurrences in all mock sample replicates their read_counts and the proportion of 
#' read_counts to the total number of reads in the sample-replicate. 
#' The lfn_variant_replicate parameter should be bellow the smallest proportion in order to keep all expected ASVs in the dataset.
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param mock_composition csv file with columns: plate, marker, sample, action (keep/tolerate), asv
#' @param sep separator used in mock_composition
#' @param outdir name of the output directory
#' @export
#'

OptimizeLFNsampleReplicate <- function(read_count_df, mock_composition="", sep=",", outdir=""){
  
  # check outdir and make tmp dir
  outdir <- check_dir(outdir)
  out = paste(outdir, "OptimizeLFNsampleReplicate.csv", sep="")
  
  # read the mock composition file and keep only lines with keep
  mock_composition_df <- read.csv(mock_composition, header=T, sep=sep) %>%
    filter(action=="keep")
  
  # get a complete and unique list of plate, marker, sample, replicate
  sample_replicate_list <- read_count_df %>%
    select(plate, marker, sample, replicate) %>%
    unique
  
  mock_composition_df <- left_join(mock_composition_df, sample_replicate_list, by=c("plate", "marker", "sample"), relationship = "many-to-many")
  
  unique_asv_keep <-  unique(mock_composition_df$asv)
  
  
  # get the total number of reads for each sample replicate
  sample_replicate_rc <- read_count_df %>%
    group_by(plate, marker, sample, replicate) %>%
    summarize(read_count_sample_replicate= sum(read_count), .groups="drop_last") %>%
    filter(sample %in% mock_composition_df$sample)
  
  # sum read_counts over replicates 
  asv_keep_df <- left_join(mock_composition_df, read_count_df, by=c("plate", "marker", "sample", "replicate", "asv"))
  asv_keep_df <- left_join(asv_keep_df, sample_replicate_rc, by=c("plate", "marker", "sample", "replicate"))
  asv_keep_df$lfn_variant_replicate <- asv_keep_df$read_count/asv_keep_df$read_count_sample_replicate
  asv_keep_df$lfn_variant_replicate <- round(asv_keep_df$lfn_variant_replicate-0.00005, digits=4)
  
  asv_keep_df <- asv_keep_df %>%
    arrange(lfn_variant_replicate) %>%
    select(plate, marker, sample, replicate, everything())
  
  write.table(asv_keep_df, file=out, sep=sep, row.names = F)
  
}


#' OptimizeLFNReaCountAndLFNvariant
#' 
#' Suggest optimal parametres for lfn_read_count_cutoff and lfn_read_count_cutoff. 
#' This script will run LFN_sample_replicate, FilterPCRerror and FilterMinReplicateNumber on the input 
#' read_count_df using parameters set by the user (ideally optimized ones), to eliminate part of the noise. 
#' The the LFN_read_count and LFN_variant is run for a series of parameter value combinations, end after each the number of FN, TP, and FP is reported. 
#' The results are written to the OptimizeLFNReaCountAndLFNvariant.csv.
#' Users should chose the parameter setting that minimizes, FN and FP.
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param known_occurrences file produced by make_known_occurrences function, with known FP and TP
#' @param sep separator used in fileinfo and mock_composition
#' @param outdir name of the output directory
#' @param min_lfn_read_count_cutoff the lowest cutoff value for LFN_read_count function to start from (10 by default). Values from this value to 100 by increments of 5 and tested
#' @param min_lnf_variant_cutoff the lowest cutoff value for LFN_variant function to start from (0.001 by default). Values from this value to 0.05 by increments of 0.001 and tested
#' @param by_replicate T/F (False by default); see LFN_variant function
#' @param lfn_sample_replicate_cutoff cutoff value for LFN_sample_replicate (see LFN_sample_replicate function; default 0.001)
#' @param pcr_error_var_prop cutoff value for FilterPCRerror (see FilterPCRerror function; default 0.1)
#' @param vsearch_path path to vsearch executables
#' @param max_mismatch parameter for FilterPCRerror (see FilterPCRerror function; default 1)
#' @param by_sample,parameter for FilterPCRerror (see FilterPCRerror function; default T)
#' @param sample_prop for FilterPCRerror (see FilterPCRerror function; default 0.8)
#' @param min_replicate_number for FilterMinReplicateNumber (see FilterMinReplicateNumber function; default 2)
#' @export
#'

OptimizeLFNReaCountAndLFNvariant <- function(read_count_df, known_occurrences="", sep=sep, outdir="", min_lfn_read_count_cutoff=10, min_lnf_variant_cutoff=0.001, by_replicate=FALSE, lfn_sample_replicate_cutoff=0.001, pcr_error_var_prop=0.1, vsearch_path="", max_mismatch=1, by_sample=T, sample_prop=0.8, min_replicate_number=2){
#  read_count_df = optimize_read_count_df
#  min_lfn_read_count_cutoff = 10
#  min_lnf_variant_cutoff = 0.001
#  rc_cutoff = 10
#  var_cutoff = 0.001
#  by_replicate = T
  
  
  outdir <- check_dir(outdir)
  out = paste(outdir, "OptimizeLFNReaCountAndLFNvariant.csv", sep="")
  
  # read known occurrences
  known_occurrences_df <- read.csv(known_occurrences, header=T, sep=sep)
  
  # filter by sample-replicate
  df <- LFN_sample_replicate(read_count_df, lfn_sample_replicate_cutoff, write_csv=F, outdir = outdir, sep=sep)
  # FilterPCRerror
  df <- FilterPCRerror(df, write_csv=F, outdir=outdir, vsearch_path=vsearch_path, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch, by_sample=by_sample, sample_prop=sample_prop, sep=sep)
  # FilterMinReplicateNumber
  df <- FilterMinReplicateNumber(df, min_replicate_number, write_csv=F, outdir = outdir, sep=sep)
  
  # make a series of cutoff values for LFN_read_count
  rc_cutoff_list <- seq(from=min_lfn_read_count_cutoff, to=100, by=10)
  # make a series of cutoff values for LFN_read_count
  var_cutoff_list <- seq(from=min_lnf_variant_cutoff, to=0.04, by=0.01)
  
  out_df <- data.frame( 
    lfn_sample_replicate_cutoff =numeric(),
    pcr_error_var_prop =numeric(),
    lfn_read_count_cutoff=numeric(),
    lnf_variant_cutoff=numeric(),
    FN=numeric(),
    TP=numeric(),
    FP=numeric()
  )
  # go through all parameter combination and count the number of TP and FN

  for(rc_cutoff in rc_cutoff_list){
    df_tmp <- df
    #LFN_read_count
    df_tmp <- LFN_read_count(df_tmp, rc_cutoff, write_csv=F, outdir = outdir, sep=sep)
    for(var_cutoff in var_cutoff_list){
      # LFN_variant
      df_tmp <- LFN_variant(df_tmp, var_cutoff, by_replicate=by_replicate, write_csv=F, outdir = outdir, sep=sep)
      # FilterMinReplicateNumber
      df_tmp <- FilterMinReplicateNumber(df_tmp, min_replicate_number, write_csv=F, outdir = outdir, sep=sep)
      # PoolReplicates
      df_tmp_sample <- PoolReplicates(df_tmp, digits=digits, write_csv=F, outdir=outdir, sep=sep)
      # pool readcount info and known occurrences info
      ko <- full_join(df_tmp_sample, known_occurrences_df, by=c("plate", "marker", "sample", "asv")) %>%
        filter(!is.na(action)) %>% # keep only lines mentioned in the known occurrences
        filter(!(is.na(mean_read_count) & action=="delete")) # delete lines if asv is not present (mean_read_count==NA) and the action is delete
      # get the number of FN (misssing expected occurrences) 
      missing <- ko %>%
        filter(is.na(mean_read_count) & action=="keep")
      FN_count <- nrow(missing)
      # get the number of TP and FP
      ko <- ko %>%
        filter(!(is.na(mean_read_count) & action=="keep")) %>%
        group_by(action) %>%
        summarise(TPFP=length(action))
      
      TP_count <- ko %>%
        filter(action=="keep") %>%
        pull(TPFP)
      FP_count <- ko %>%
        filter(action=="delete") %>%
        pull(TPFP)
      
      print(lfn_sample_replicate_cutoff)
      print(pcr_error_var_prop)
      print(rc_cutoff)
      print(var_cutoff)
      print(FN_count)
      print(TP_count)
      print(FP_count)
      
      new_line <- data.frame(lfn_sample_replicate_cutoff=lfn_sample_replicate_cutoff, pcr_error_var_prop=pcr_error_var_prop, lfn_read_count_cutoff=rc_cutoff, lnf_variant_cutoff=var_cutoff ,FN=FN_count, TP=TP_count, FP=FP_count)
      print(rc_cutoff)
      print(var_cutoff)
      print(new_line)
      print(dim(new_line))
      print(out_df)
      print(dim(out_df))
#      out_df <- bind_rows(out_df, new_line )
    }
  }
  
  out_df <- out_df %>%
    arrange(FN, FP, lnf_variant_cutoff, lfn_read_count_cutoff)
  
  write.table(out_df, file=out, sep=sep, row.names = F)
}