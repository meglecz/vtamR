% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{flagChimera}
\alias{flagChimera}
\title{Flag chimeric sequences}
\usage{
flagChimera(unique_asv_df, vsearch_path = "vsearch", abskew = 2)
}
\arguments{
\item{unique_asv_df}{Data frame with the following variables:
asv, read_count; ASVs must be unique.}

\item{vsearch_path}{Character string: path to vsearch executables.}

\item{abskew}{Positive integer: a chimera must be at least \code{abskew}
times less frequent that the parental ASVs.}
}
\value{
The input data frame completed by chimera \code{column}.
1 if potential chimera, 0 otherwise.
}
\description{
Select chimeras in a data frame of unique ASVs.
Add a column to the input data frame with 1 if the ASV is a probable
chimera and 0 otherwise.
}
\examples{
\dontrun{
unique_asv_df <- read_count_df \%>\%
  group_by(asv) \%>\%
  summarize(read_count = sum(read_count))
flagChimera(unique_asv_df, vsearch_path=vsearch_path, abskew=2)
}
}
