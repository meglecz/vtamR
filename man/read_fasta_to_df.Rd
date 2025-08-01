% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{read_fasta_to_df}
\alias{read_fasta_to_df}
\title{Read fasta file to data frame}
\usage{
read_fasta_to_df(file, dereplicate = F)
}
\arguments{
\item{file}{Character string: name of the input fasta file.
Can be  uncompressed or compressed in gz format (zip files are not supported).}

\item{dereplicate}{logical: If TRUE returns a data frame with asv and
read_count columns. If FALSE data frame with header and sequence columns.}
}
\value{
Data frame with two columns:
\itemize{
\item headers and sequences if dereplicate==FALSE
\item asv and read_count if dereplicate==TRUE
}
}
\description{
Read a fasta file to a data frame
}
\examples{
\dontrun{
read_fasta_to_df(file="data/test.fasta", dereplicate=F)
read_fasta_to_df(file="data/test.fasta", dereplicate=T)
}
}
