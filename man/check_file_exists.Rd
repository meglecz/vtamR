% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{check_file_exists}
\alias{check_file_exists}
\title{Check existence of files}
\usage{
check_file_exists(file_list, dir = "")
}
\arguments{
\item{file_list}{Vector of file names.}

\item{dir}{Character string: name of the directory containing the input files}
}
\value{
Error message if some files do not exist and stop program.
}
\description{
Check if all files in the input vector exist.
}
\examples{
\dontrun{
file_list <- c("14ben01-1.fasta", "14ben01-2.fasta")
check_file_exists(file_list=file_list, dir="vtamR_test/out_mfzr/sorted")
}
}
