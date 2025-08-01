% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{read_fasta_seq}
\alias{read_fasta_seq}
\title{Read sequences from fasta}
\usage{
read_fasta_seq(filename = filename, dereplicate = F)
}
\arguments{
\item{filename}{Character string: input fasta file name including full path.}

\item{dereplicate}{logical: If TRUE, return ASVs with read counts.}
}
\value{
Data frame with one read in each line (dereplicate==F), or a data frame with
asv and read_count columns (dereplicate==T).
}
\description{
Read sequences from a fasta file.
}
\details{
Input fasta can be gzip compressed or uncompressed.
}
\examples{
\dontrun{
read_df <- read_fasta_seq(filename="data/test.fasta", dereplicate=F)
asv_df <- read_fasta_seq(filename="data/test.fasta", dereplicate=T)
}
}
