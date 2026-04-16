#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first if_else
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat scale_y_log10
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym :=
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL


#' Plot Read Counts by Sample
#'
#' Creates a bar plot showing the number of reads per sample or sample–replicate.
#' If sample metadata is provided, bars are colored according to sample type.
#'
#' @param read_count_df Data frame with at least `sample` and `read_count`
#'   columns, and optionally a `replicate` column.
#' @param sampleinfo Data frame or CSV file with `sample` and `sample_type`
#'   columns (e.g. `real`, `mock`, `negative`).
#' @param sep Field separator used in CSV files.
#' @param sample_replicate Logical; if `TRUE`, the plot is generated per
#'   sample–replicate combination, otherwise per sample.
#' @param x_axis_label_size Numeric; size of x-axis labels.
#' @param plotfile Character string: name of the output PNG file. If empty,
#'   no file is written.
#'
#' @return A bar plot object.
#'
#' @export
#' 
plot_read_count_by_sample <- function(read_count_df, 
                                      sampleinfo="", 
                                      sample_replicate=T, 
                                      sep=",", 
                                      x_axis_label_size=6,
                                      plotfile=""
){
  
  if(sampleinfo != ""){
    sampleinfo_df <- read.csv(sampleinfo, sep=sep)
    # get sample type for each sample
    sampleinfo_df <- sampleinfo_df %>%
      select(sample, sample_type) %>%
      unique()
  }else{
    sampleinfo_df <- data.frame(sample=unique(read_count_df$sample),
                                sampleinfo=rep("sample_type", 
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
    df <- left_join(df, sampleinfo_df, by="sample")
    
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
    df <- left_join(df, sampleinfo_df, by="sample")
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
  
  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile, width = 2000, height = 1500, res = 300) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  return(p)
}


#' Plot Read Count Histogram
#'
#' Creates a histogram of read counts per ASV.
#'
#' This function filters low-abundance ASVs (if requested) and then plots the
#' distribution of read counts across ASVs.
#'
#' @param read_count_df Data frame with `asv` and `read_count` columns.
#' @param min_read_count Numeric; minimum read count threshold. ASVs with fewer
#'   reads are filtered out before plotting.
#' @param binwidth Numeric; width of histogram bins for read count intervals.
#' @param plotfile Character string: name of the output PNG file. If empty,
#'   no file is written.
#'
#' @return A histogram plot object.
#'
#' @export
#' 
plot_read_count_histogram <- function(read_count_df, 
                                      min_read_count=0, 
                                      binwidth=100,
                                      plotfile= ""
){
  
  # get read_count for each asv
  df <- read_count_df %>%
    group_by(asv) %>%
    summarize("Number_of_reads"= sum(read_count))
  # filter out low read_count
  df <- subset(df, Number_of_reads > min_read_count)
  
  
  p <- ggplot(df, aes(x = Number_of_reads)) +
    geom_histogram(binwidth = binwidth, fill = "blue", color = "blue", 
                   aes(y = after_stat(count))) +
    scale_y_log10() +
    labs(title = "Distribution of Read Counts",
         x = "Read Count",
         y = "Frequency") +
    theme_minimal()  
  
  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile, width = 2000, height = 1500, res = 300) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  return(p)
}

#' Plot Renkonen Distance Barplot
#'
#' Creates a bar plot of Renkonen distances between pairs of samples or
#' sample–replicates.
#'
#' If sample metadata is provided, bars are colored according to sample type.
#'
#' @param df Data frame with the following columns:
#'   `sample1`, `sample2`, `replicate1`, `replicate2`, `renkonen_d`.
#'   This output can be generated using `make_renkonen_distance_matrix()`.
#' @param sampleinfo Data frame or CSV file with columns `sample` and
#'   `sample_type` (e.g. `real`, `mock`, `negative`).
#' @param sep Field separator used in CSV files.
#' @param x_axis_label_size Numeric; size of x-axis labels.
#' @param plotfile Character string: name of the output PNG file. If empty,
#'   no file is written.
#'
#' @return A bar plot object.
#'
#' @export
#' 
plot_renkonen_distance_barplot <- function(df, 
                                           sampleinfo=NULL, 
                                           sep=",", 
                                           x_axis_label_size=6,
                                           plotfile=""
){
  
  if(is.character(sampleinfo)){ # input file
    sampleinfo_df <- read.csv(sampleinfo, sep=sep)
    # get sample type for each sample
    sampleinfo_df <- sampleinfo_df %>%
      select(sample, sample_type) %>%
      unique()
  }else if (!is.null(sampleinfo)){
    sampleinfo_df <- sampleinfo %>%
      select(sample, sample_type) %>%
      unique()
  }
  else{
    sampleinfo_df <- data.frame(sample=unique(read_count_df$sample),
                                sample_type=rep("sample_type", 
                                                length(unique(read_count_df$sample))
                                )
    )
  }
  
  df <- left_join(df, sampleinfo_df, by=c("sample1"="sample")) %>%
    arrange(renkonen_d)
  # make replicate pairs (replicate is a concatenation of sample and replicate)
  df$replicate_pair <- paste(df$sample1, df$replicate1, df$replicate2, sep = ":")
  # Convert 'replicate_pair' to a factor with the desired order
  df$replicate_pair <- factor(df$replicate_pair, levels = df$replicate_pair)
  
  p <- ggplot(df, aes(x = replicate_pair, y = renkonen_d, fill = sample_type)) +
    geom_bar(stat = "identity") +
    labs(title = "Renkonen distances \n between pairs of replicates of the same sample",
         x = "Replicate pair",
         y = "Renkonen distance",
         fill = "Sample Type") +
    # Rotate x-axis labels by 45 degrees
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = x_axis_label_size), 
          plot.title = element_text(hjust = 0.5)) +  # Center the title
    scale_y_continuous(expand=c(0,0)) # avoid space between labels and x axis
  
  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile, width = 2000, height = 1500, res = 300) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  return(p)
}


#' Plot Renkonen Distance Density
#'
#' Creates a density plot of Renkonen distances between pairs of samples or
#' sample–replicates.
#'
#' @param df Data frame with the following columns:
#'   `sample1`, `sample2`, `replicate1`, `replicate2`, `renkonen_d`.
#'   This output can be produced using `compute_renkonen_distances()`.
#' @param plotfile Character string: name of the output PNG file. If empty,
#'   no file is written.
#'
#' @return A density plot object.
#'
#' @export
#' 
plot_renkonen_distance_density <- function(df, plotfile=""){
  
  df$comparison <- ifelse(df$sample1 == df$sample2, "within samples", "between samples")
  
  p <- ggplot(df, aes(x = renkonen_d, fill = comparison)) +
    geom_density(alpha = 0.5) +  # Add transparency to the density plot
    labs(title = "Density of Renkonen distances",
         x = "Distribution of Renkonen Distances between pairs of replicates",
         y = "Density") +
    theme_minimal()
  
  if(plotfile != ""){
    check_dir(plotfile, is_file=TRUE)
    png(filename=plotfile, width = 2000, height = 1500, res = 300) # one png file per plot
    print(p) # print plot to file
    dev.off()
  }
  
  return(p)
}