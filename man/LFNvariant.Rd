% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{LFNvariant}
\alias{LFNvariant}
\title{LFNvariant}
\usage{
LFNvariant(
  read_count,
  cutoff = 0.001,
  by_replicate = FALSE,
  outfile = "",
  sep = ",",
  min_read_count_prop = 0.7
)
}
\arguments{
\item{read_count}{Data frame or csv file with the following variables:
asv_id, sample, replicate, read_count, asv.}

\item{cutoff}{Real (0-1): minimum proportion of the read count of
an occurrence within all reads of the asv or asv-replicate. Bellow this cutoff
the occurrence is deleted.}

\item{by_replicate}{logical: Compare read count of the occurrence to the
read counts of the ASV-replicate.}

\item{outfile}{Character string: csv file name to print the output data
frame if necessary. If empty, no file is written.}

\item{sep}{Field separator character in input and output csv files.}

\item{min_read_count_prop}{Real (0-1): If the proportion of the read count
of a variant in the output compared to the input is less then
min_read_count_prop, prints out a warning, since it suggest a
to high cutoff value}
}
\value{
Filtered read_count_df data frame.
}
\description{
If by_replicate is FALSE: Eliminate occurrences where the
(read_count/read_count of the asv in the data set) is less than cutoff.
If by_replicate is TRUE: Eliminate occurrences where the
(read_count/read_count of the asv in its replicate) is less than cutoff.
}
\details{
Issues a warning if the total read count of an ASV has been reduced
bellow min_read_count_prop, since it can indicate a to high cutoff value.
}
\examples{
\dontrun{
filtered_read_count_df <- LFNvariant(read_count_df, cutoff=0.005, min_read_count_prop=0.8)
}
}
