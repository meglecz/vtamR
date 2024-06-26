% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{FilterChimera}
\alias{FilterChimera}
\title{FilterChimera}
\usage{
FilterChimera(
  read_count,
  outfile = "",
  vsearch_path = "",
  by_sample = T,
  sample_prop = 0.8,
  abskew = 2,
  sep = ","
)
}
\arguments{
\item{read_count}{Data frame or csv file with the following variables:
asv_id, sample, replicate, read_count, asv.}

\item{outfile}{Character string: csv file name to print the output data
frame if necessary. If empty, no file is written.}

\item{vsearch_path}{Character string: path to vsearch executables.
Can be empty if vsearch in the the PATH.}

\item{by_sample}{Boolean: ASVs are flagged as chimera separately
for each sample.}

\item{sample_prop}{Boolean: if by_sample=TRUE, the ASV deleted if
they are flagged as chimera in at least \code{sample_prop} of the
samples among the sample they are present.}

\item{abskew}{Positive integer: a chimera must be at least \code{abskew}
times less frequent that the parental ASVs.}
}
\value{
Filtered read_count_df data frame.
}
\description{
Filter out Chimeras.
}
\examples{
filtered_read_count_df <- FilterChimera(read_count_df, vsearch_path=vsearch_path, by_sample=T, sample_prop=0.7, abskew=4)
filtered_read_count_df <- FilterChimera(read_count_df, vsearch_path=vsearch_path, by_sample=F, abskew=4)
}
