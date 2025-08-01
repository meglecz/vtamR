% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{CountReadsDir}
\alias{CountReadsDir}
\title{Count reads in directory}
\usage{
CountReadsDir(
  dir,
  pattern = ".",
  file_type = "fasta",
  outfile = "",
  sep = ",",
  quiet = T
)
}
\arguments{
\item{dir}{Character string: name of the input directory.}

\item{pattern}{Regular expression: pattern in the name of the files in the
input directory. Read only files with pattern in their names.}

\item{file_type}{Character string with the following values: "fasta", "fastq".
For all other values, the number of lines is returned.}

\item{outfile}{Character string: csv file name to print the output data
frame if necessary. If empty, no file is written.}

\item{sep}{Field separator character in input and output csv files.}

\item{quiet}{logical: If TRUE, suppress informational messages and only
show warnings or errors.}
}
\value{
Data frame with 2 columns: filename, read_count.
}
\description{
Counts the number of sequences in fasta or fastq files, or the number of
lines in other files. Reads all files in the input directory with
\code{pattren} in their names.
}
\details{
Efficient in linux-like systems, but slow on Windows for large files.
Can handle gz compressed and uncompressed files, but not zip files.
}
\examples{
\dontrun{
CountReadsDir(dir="out", pattern=".fastq", file_type="fastq")
CountReadsDir(dir="out", pattern="^mfzr", file_type="fasta")
}
}
