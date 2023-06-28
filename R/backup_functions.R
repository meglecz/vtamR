#' FilerPCRerror_NW
#' 
#' Filter out all ASVs if they very similar (max_mismatch) to another more frequent ASV (pcr_error_var_prop)
#' The whole plate can be analysed et once (by_sampe=F)
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param write_csv T/F; write read_counts to csv file; default=FALSE
#' @param outdir name of the output directory
#' @param pcr_error_var_prop if the proportion of read_counts of two similar ASVs is bellow pcr_error_var_prop, the less abundant is flagged as a PCR error
#' @param max_mismatch maximum number of mismatches (gaps included) to consider two ASVs as similar
#' @param by_sample T/F, if T ASVs are flagged as an PCR error separately for each sample
#' @param sample_prop if by_sample=T, the ASV must be flagged as an PCRerror in sample_prop of the cases to be eliminated
#' @export
#' 
FilterPCRerror_NW <- function(read_count_df, by_sample=T, write_csv=F, outdir=outdir, pcr_error_var_prop=0.1, max_mismatch=1, sample_prop=0.8){
  
  # get unique list of ASVs with their total read_count in the run
  unique_asv_df <- read_count_df %>%
    group_by(asv) %>%
    summarize(read_count = sum(read_count)) %>%
    arrange(desc(read_count))
  
  if(by_sample){ # sample by sample
    sample_list <- unique(read_count_df$sample)
    # loop over samples
    for(sample_loc in sample_list){
      # get unique list of ASVs with their total read_count in the sample
      unique_asv_df_sample <- read_count_df %>%
        filter(sample == sample_loc) %>%
        group_by(asv)%>%
        summarize(read_count = sum(read_count))%>%
        arrange(desc(read_count))
      # flag PCR errors
      unique_asv_df_sample <- flagPCRerrors(unique_asv_df_sample, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch)
      # remove read_count column
      unique_asv_df_sample <- unique_asv_df_sample[, -which(names(unique_asv_df_sample) == "read_count")]
      # add a column for for each sample to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
      unique_asv_df <- left_join(unique_asv_df, unique_asv_df_sample, by = "asv")
    }
  }
  else{ # whole dataset
    # add a column for for each sample to unique_asv_df, with 1 if ASV is flagged in the sample, 0 otherwise
    unique_asv_df <- flagPCRerrors(unique_asv_df_sample, pcr_error_var_prop=pcr_error_var_prop, max_mismatch=max_mismatch)
  }
  
  # count the number of times each ASV has been flagged and when it has not. Ignore NA, when the ASV is not present in the sample
  unique_asv_df$yes <- rowSums(unique_asv_df[3:ncol(unique_asv_df)] == 1, na.rm = TRUE)
  unique_asv_df$no <- rowSums(unique_asv_df[3:(ncol(unique_asv_df)-1)] == 0, na.rm = TRUE)
  # keep only ASVs, that are not flagged in sample_prop proportion of the samples where they are present  
  unique_asv_df <- unique_asv_df %>%
    filter(no/(yes+no) >= (1-sample_prop))
  
  # eliminate potential PCRerrors from read_count_df
  read_count_df <- read_count_df %>%
    filter(asv %in% unique_asv_df$asv)
  
  if(write_csv){
    write.csv(read_count_df, file = paste(outdir, "Input.csv", sep=""))
  }
  
  return(read_count_df)
}

#' flagPCRerrors
#' 
#' Flag unique ASVs if they can be a PCR error
#' 
#'  
#' @param read_count_df data frame with the following variables: asv, plate, marker, sample, replicate, read_count
#' @param pcr_error_var_prop if the proportion of read_counts of two similar ASVs is bellow pcr_error_var_prop, the less abundant is flagged as a PCR error
#' @param max_mismatch maximum number of mismatches (gaps included) to consider two ASVs as similar
#' @export
#' 
flagPCRerrors <- function(unique_asv_df, pcr_error_var_prop=0.1, max_mismatch=1){
  
  # define DNAfull substitution matrix
  # https://rosalind.info/glossary/dnafull/
  scores <- matrix(c(5,-4,-4,-4,-4, 1, 1,-4,-4, 1,-4,-1,-1,-1,-2,
                     -4, 5,-4,-4,-4, 1,-4, 1, 1,-4,-1,-4,-1,-1,-2,
                     -4,-4, 5,-4, 1,-4, 1,-4, 1,-4,-1,-1,-4,-1,-2,
                     -4,-4,-4, 5, 1,-4,-4, 1,-4, 1,-1,-1,-1,-4,-2,
                     -4,-4, 1, 1,-1,-4,-2,-2,-2,-2,-1,-1,-3,-3,-1,
                     1, 1,-4,-4,-4,-1,-2,-2,-2,-2,-3,-3,-1,-1,-1,
                     1,-4, 1,-4,-2,-2,-1,-4,-2,-2,-3,-1,-3,-1,-1,
                     -4, 1,-4, 1,-2,-2,-4,-1,-2,-2,-1,-3,-1,-3,-1,
                     -4, 1, 1,-4,-2,-2,-2,-2,-1,-4,-1,-3,-3,-1,-1,
                     1,-4,-4, 1,-2,-2,-2,-2,-4,-1,-3,-1,-1,-3,-1,
                     -4,-1,-1,-1,-1,-3,-3,-1,-1,-3,-1,-2,-2,-2,-1,
                     -1,-4,-1,-1,-1,-3,-1,-3,-3,-1,-2,-1,-2,-2,-1,
                     -1,-1,-4,-1,-3,-1,-3,-1,-3,-1,-2,-2,-1,-2,-1 ,
                     -1,-1,-1,-4,-3,-1,-1,-3,-1,-3,-2,-2,-2,-1,-1,
                     -2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1), nrow=15, ncol=15)
  nt_list <- c("A" ,  "T" ,  "G" ,  "C"  , "S" ,  "W"  , "R"   ,"Y" , "K" ,  "M" ,  "B" ,  "V"  , "H" ,  "D"  , "N")
  DNAfull <-  matrix(scores, nrow = 15, ncol = 15, dimnames = list(nt_list, nt_list))
  
  # add an error column
  unique_asv_df$PCRerror <- rep(0, length(unique_asv_df$asv))
  
  # go through all pairs of ASVs, starting with the 2nd most abundant
  for(fille in 2:length(unique_asv_df$asv)){ # potential daughter asvs
    for(parent in 1:(fille-1)){ # potential parents (higher of equal read_counts; then daughter)
      if((unique_asv_df$read_count[parent] * pcr_error_var_prop) >= unique_asv_df$read_count[fille] ){ # read_count of parent is high compared to daughter
        # align the two sequences
        #methods(class = class(alignment))
        alignment <- pairwiseAlignment(pattern = DNAString(unique_asv_df$asv[parent]), subject = DNAString(unique_asv_df$asv[fille]), substitutionMatrix = DNAfull, type = "global")
        # get he number of mismatches (gaps included)
        nmismatch <- (width(alignedPattern(alignment)) - nmatch(alignment))
        if(nmismatch <= max_mismatch){ # mark the daughter sequence as a probable error
          unique_asv_df$PCRerror[fille] <- 1
          break # we found a parent, stop searching
        }
      }
    }
  }
  return(unique_asv_df)
}