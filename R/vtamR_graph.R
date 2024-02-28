#' graph_read_count_by_sample
#' 
#' Create barplot with the number of reads in each sample or sample-replicate;
#' If information is given on sample types, bars are colored in function of them
#' 
#' @param read_count_df data frame with sample, read_count and replicate (optional) columns
#' @param sample_types file with sample and sample_type (real/mock/negative) columns
#' @param sep separator used in csv files
#' @param sample_replicate [T/F] if true barplot is made by sample-replicates, by sample otherwise
#' @export
#' 
graph_read_count_by_sample <- function(read_count_df, sample_types="", sample_replicate=T, sep=","){
  
  if(sample_types != ""){
    sample_types_df <- read.csv(sample_types, sep=sep)
    # get sample type for each sample
    sample_types_df <- sample_types_df %>%
      select(sample, sample_type) %>%
      unique()
  }else{
    sample_types_df <- data.frame(sample=unique(read_count_df$sample),
                               sample_type=rep("sample_type", length(unique(read_count_df$sample)))
                                 )
  }
  
  if(sample_replicate){ # make a graph for each sample-replicate
    
    df <- read_count_df %>%
      group_by(sample, replicate) %>%
      summarize("Number_of_reads" = sum(read_count), .groups="drop_last") %>%
      arrange(desc(Number_of_reads))
    df$sample_replicate <- paste(df$sample, df$replicate, sep="-")
    # Convert 'sample_replicate' to a factor with the desired order
    df$sample_replicate <- factor(df$sample_replicate, levels = unique(df$sample_replicate))
    # add sample_type
    df <- left_join(df, sample_types_df, by="sample")
    
    p <- ggplot(df, aes(x = sample_replicate, y = Number_of_reads, fill = sample_type)) +
      geom_bar(stat = "identity") +
      labs(title = "Barplot of Read Counts by Sample-Replicate",
           x = "Sample-Replicate",
           y = "Read Count",
           fill = "Sample Type") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels by 45 degrees
            plot.title = element_text(hjust = 0.5))  # Center the title
    
  }else{ # make a graph for each sample
    
    df <- read_count_df %>%
      group_by(sample) %>%
      summarize("Number_of_reads" = sum(read_count)) %>%
      arrange(desc(Number_of_reads))
    df <- left_join(df, sample_types_df, by="sample")
    # Convert 'sample' to a factor with the desired order
    df$sample <- factor(df$sample, levels = unique(df$sample))
    
    p <- ggplot(df, aes(x = sample, y = Number_of_reads, fill = sample_type)) +
      geom_bar(stat = "identity") +
      labs(title = "Barplot of Read Counts by Sample",
           x = "Sample",
           y = "Read Count",
           fill = "Sample Type") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels by 45 degrees
            plot.title = element_text(hjust = 0.5))  # Center the title
  }
  return(p)
}

#' graph_read_count_by_variant
#' 
#' Create histogram with the number of reads in each sample or sample-replicate;
#' If information is given on sample types, bars are colored in function of them
#' 
#' @param read_count_df data frame with asv and read_count columns
#' @param min_read_count filter out ASVs with less than min_read_count reads, before making the graph
#' @param binwidth width of the read count intervals 
#' @export
#' 

graph_read_count_by_variant <- function(read_count_df, min_read_count=0, binwidth=100){
  
  # get read_count for each asv
  df <- read_count_df %>%
    group_by(asv) %>%
    summarize("Number_of_reads"= sum(read_count))
  # filter out low read_count
  df <- subset(df, Number_of_reads > min_read_count)

  
  p <- ggplot(df, aes(x = Number_of_reads)) +
      geom_histogram(binwidth = binwidth, fill = "blue", color = "blue", aes(y = after_stat(count))) +
      labs(title = "Distribution of Read Counts",
           x = "Read Count",
           y = "Frequency") +
      theme_minimal()  
  return(p)
}

#' make_renkonen_df_all
#' 
#' Calculate the Renkonen distance among all pairs of sample-replicates
#' Returns a data frame with the following columns: sample1,sample2,replicate1,replicate2,renkonen_d,sample_comp (within, if sample1 and sample2 are identical, between otherwise)
#' 
#' @param read_count_df data frame with asv, sample, replicate, and read_count columns
#' @export
#' 

make_renkonen_df_all <- function(read_count_df){
  
  df <- read_count_df %>%
    select(asv, sample, replicate, read_count)
  df$sr <- paste(df$sample, df$replicate, sep="-")
  # list of samples
  sample_replicate_df <- df %>%
    select(sample, replicate, sr) %>%
    unique()
  
  # fine empty dataframe
  renkonen_df <- data.frame("sample1" = character(),
                            "sample2" = character(),
                            "sr1"  = character(),
                            "sr2"  = character(),
                            "renkonen_d" = numeric())
  
  # loop over all pairs of sample-replicates within sample
  for(i in 1:(nrow(sample_replicate_df)-1)){
    sri <- sample_replicate_df$sr[i]
    sampi <- sample_replicate_df$sample[i]
    repli <- sample_replicate_df$replicate[i]
    dfi <- filter(df, sr == sri)
    for(j in ((i+1):nrow(sample_replicate_df))){
      srj <- sample_replicate_df$sr[j]
      sampj <- sample_replicate_df$sample[j]
      replj <- sample_replicate_df$replicate[j]
      dfj <- filter(df, sr == srj)
      rdist <- calculate_renkonen_dist(dfi, dfj)
      # add line to renkonen_df
      new_line <- data.frame(sample1 = sampi, sample2 = sampj, replicate1 = sri, replicate2 = srj, renkonen_d = rdist)
      renkonen_df <- rbind(renkonen_df, new_line)
    }
  }
  
  renkonen_df$comparaison <- ifelse(renkonen_df$sample1 == renkonen_df$sample2, "within", "between")
  return(renkonen_df)
}
