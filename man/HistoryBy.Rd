% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{HistoryBy}
\alias{HistoryBy}
\title{History By}
\usage{
HistoryBy(dir, feature, value, sep = ",")
}
\arguments{
\item{dir}{Character string: directory containing the output files
of intermediate filtering steps. File names should start with a number
followed by an underscore (e.g. 5_LFNsampleReplicate.csv).}

\item{feature}{Character string with the following possible values:
"asv_id", "asv", "sample", "replicate", "read_count".}

\item{value}{Character string: values of feature that should selected.
The output data frame contains all lines where value is present
in feature in the input files.}

\item{sep}{Field separator character in input and output csv files.}
}
\value{
Data frame with the selected lines of the input files.
}
\description{
Filters a feature (asv_id/asv/sample/replicate),
in all output files of intermediate filtering steps to retain only
lines corresponding a value.
}
\examples{
\dontrun{
HistoryBy(dir="out", feature="asv_id", value=1)
HistoryBy(dir="out", feature="sample", value="tpos1")
}
}
