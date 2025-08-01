% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{RandomSeq}
\alias{RandomSeq}
\title{Select random sequences}
\usage{
RandomSeq(
  fastainfo,
  n,
  fasta_dir = "",
  outdir = "",
  vsearch_path = "vsearch",
  randseed = 0,
  compress = F,
  sep = ",",
  quiet = T
)
}
\arguments{
\item{fastainfo}{Data frame or csv file with a fasta column containing input
fasta file names. Files can be gzip compressed.}

\item{n}{Positive integer: the number of randomly selected sequences from each input file.}

\item{fasta_dir}{Character string: directory that contains the
input fasta files}

\item{outdir}{Character string: output directory.}

\item{vsearch_path}{Character string: path to vsearch executables.}

\item{randseed}{Positive integer: seed for random sampling.
0 (default value) means to use a pseudo-random seed.
A given non-zero seed produces always the same result.}

\item{compress}{logical: Compress output files to gzip format.}

\item{sep}{Field separator character in input and output csv files.}

\item{quiet}{logical: If TRUE, suppress informational messages and only
show warnings or errors.}
}
\value{
The input data frame with updated file names and read counts.
}
\description{
Random select n sequences from each input fasta file.
}
\details{
Do not work on Windows! If using Windows, please, use RandomSeqWindows
}
\examples{
\dontrun{
fastainfo_df <- RandomSeq(fastainfo, 
   fasta_dir="data/fasta", 
   outdir="data/randomseq"
   )
}
}
