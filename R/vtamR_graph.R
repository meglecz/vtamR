#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first if_else
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym :=
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL

#' Barplot_ReadCountBySample
#' 
#' Create barplot with the number of reads in each sample or sample-replicate;
#' If information is given on sample types, bars are colored in function of them
#' 
#' @param read_count_df data frame with sample, read_count and replicate (optional) columns
#' @param sample_types file with sample and sample_type (real/mock/negative) columns
#' @param sep separator used in csv files
#' @param sample_replicate Boolean. If TRUE barplot is made by sample-replicates,
#' by sample otherwise
#' @param x_axis_label_size size of labels in x axis
#' @export
#' 
#' 
Barplot_ReadCountBySample <- function(read_count_df, 
                                      sample_types="", 
                                      sample_replicate=T, 
                                      sep=",", 
                                      x_axis_label_size=6
                                      ){
  
  if(sample_types != ""){
    sample_types_df <- read.csv(sample_types, sep=sep)
    # get sample type for each sample
    sample_types_df <- sample_types_df %>%
      select(sample, sample_type) %>%
      unique()
  }else{
    sample_types_df <- data.frame(sample=unique(read_count_df$sample),
                               sample_type=rep("sample_type", 
                                               length(unique(read_count_df$sample))
                                               )
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
      # Rotate x-axis labels by 45 degrees
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = x_axis_label_size),  
            plot.title = element_text(hjust = 0.5)) + # Center the title
            scale_y_continuous(expand=c(0,0)) # avoid space between labels and x axis
    
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
      # Rotate x-axis labels by 45 degrees
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = x_axis_label_size),  
            plot.title = element_text(hjust = 0.5)) + # Center the title
      scale_y_continuous(expand=c(0,0)) # avoid space between labels and x axis
  }
  return(p)
}

#' Histogram_ReadCountByVariant
#' 
#' Create histogram with the number of reads in each sample or sample-replicate;
#' If information is given on sample types, bars are colored in function of them
#' 
#' @param read_count_df data frame with asv and read_count columns
#' @param min_read_count filter out ASVs with less than min_read_count reads, 
#' before making the graph
#' @param binwidth width of the read count intervals 
#' @export
#' 

Histogram_ReadCountByVariant <- function(read_count_df, min_read_count=0, binwidth=100){
  
  # get read_count for each asv
  df <- read_count_df %>%
    group_by(asv) %>%
    summarize("Number_of_reads"= sum(read_count))
  # filter out low read_count
  df <- subset(df, Number_of_reads > min_read_count)

  
  p <- ggplot(df, aes(x = Number_of_reads)) +
      geom_histogram(binwidth = binwidth, fill = "blue", color = "blue", 
                     aes(y = after_stat(count))) +
      labs(title = "Distribution of Read Counts",
           x = "Read Count",
           y = "Frequency") +
      theme_minimal()  
  return(p)
}

#' Barplot_RenkonenDistance
#' 
#' Create barplot with renkonen distances between pairs of sample-replicates
#' If information is given on sample types, bars are colored in function of them
#' 
#' @param df data frame with the following columns: 
#' sample1,sample2,replicate1,replicate2,renkonen_d 
#' (can be produced by make_renkonen_distance_matrix)
#' @param sample_types Data frame or CSV file with the following columns:
#' sample, sample_type (real/mock/negative)
#' @param sep separator used in csv files
#' @param x_axis_label_size size of labels in x axis
#' @export
#' 
Barplot_RenkonenDistance <- function(df, 
                                     sample_types=NULL, 
                                     sep=",", 
                                     x_axis_label_size=6
                                     ){
  
  if(is.character(sample_types)){ # input file
    sample_types_df <- read.csv(sample_types, sep=sep)
    # get sample type for each sample
    sample_types_df <- sample_types_df %>%
      select(sample, sample_type) %>%
      unique()
  }else if (!is.null(sample_types)){
    sample_types_df <- sample_types %>%
      select(sample, sample_type) %>%
      unique()
  }
  else{
    sample_types_df <- data.frame(sample=unique(read_count_df$sample),
                                  sample_type=rep("sample_type", 
                                                  length(unique(read_count_df$sample))
                                                  )
                                )
  }
  
  df <- left_join(df, sample_types_df, by=c("sample1"="sample")) %>%
    arrange(renkonen_d)
  # make replicate pairs (replicate is a concatenation of sample and replicate)
  df$replicate_pair <- paste(df$sample1, df$replicate1, df$replicate2, sep = ":")
  # Convert 'replicate_pair' to a factor with the desired order
  df$replicate_pair <- factor(df$replicate_pair, levels = df$replicate_pair)
  
  p <- ggplot(df, aes(x = replicate_pair, y = renkonen_d, fill = sample_type)) +
    geom_bar(stat = "identity") +
    labs(title = "Renkonen distances between pairs of replicates of the same sample",
         x = "Replicate pair",
         y = "Renkonen distance",
         fill = "Sample Type") +
    # Rotate x-axis labels by 45 degrees
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = x_axis_label_size), 
          plot.title = element_text(hjust = 0.5)) +  # Center the title
    scale_y_continuous(expand=c(0,0)) # avoid space between labels and x axis
  return(p)
}

#' DensityPlot_RenkonenDistance
#' 
#' Create density plot with Renkonen distances between pairs of sample-replicates
#' 
#' @param df data frame with the following columns: 
#' sample1,sample2,replicate1,replicate2,renkonen_d 
#' (can be produced by MakeRenkonenDistances)
#' @export
#' 
DensityPlot_RenkonenDistance <- function(df){
  
  df$comparison <- ifelse(df$sample1 == df$sample2, "within samples", "between samples")
  
  ggplot(df, aes(x = renkonen_d, fill = comparison)) +
    geom_density(alpha = 0.5) +  # Add transparency to the density plot
    labs(title = "Density of Renkonen distances",
         x = "Distribution of Renkonen Distances between pairs of replicates",
         y = "Density") +
    theme_minimal()
}
