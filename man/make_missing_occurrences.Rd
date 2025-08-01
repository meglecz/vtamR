% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{make_missing_occurrences}
\alias{make_missing_occurrences}
\title{Make Missing Occurrences}
\usage{
make_missing_occurrences(
  read_count_samples,
  mock_composition,
  sep = ",",
  out = ""
)
}
\arguments{
\item{read_count_samples}{Data frame or csv file with the following variables:
asv, plate, marker, sample, read_count.}

\item{mock_composition}{Data frame or csv file with columns:
sample, action (keep/tolerate), asv.}

\item{sep}{Field separator character in input and output csv files.}

\item{out}{Character string: output file. If empty, no file is written.}
}
\value{
Data frame with the following columns:
sample, action, asv, asv_id
}
\description{
Prepare a data fram that list all expected occurrences that are missing
(False negatives).
}
\examples{
\dontrun{
make_missing_occurrences(read_count_samples=read_count_samples_df, 
    mock_composition="data/mock_composition.csv"
    )
}
}
