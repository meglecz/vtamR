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