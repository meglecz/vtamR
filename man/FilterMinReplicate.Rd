% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{FilterMinReplicate}
\alias{FilterMinReplicate}
\title{FilterMinReplicate}
\usage{
FilterMinReplicate(read_count, cutoff = 2, outfile = "", sep = ",")
}
\arguments{
\item{read_count}{Data frame or csv file with the following variables:
asv_id, sample, replicate, read_count, asv.}

\item{cutoff}{Positive integer: minimum number of replicates.}

\item{outfile}{Character string: csv file name to print the output data
frame if necessary. If empty, no file is written.}

\item{sep}{Field separator character in input and output csv files.}
}
\value{
Filtered read_count_df data frame.
}
\description{
Filter out all occurrences where the asv in not present in at least
\code{cutoff} number of replicates of the sample.
}
\examples{
\dontrun{
filtered_read_count_df <- FilterMinReplicate(read_count_df, cutoff=3)
}
}
