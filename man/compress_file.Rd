% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{compress_file}
\alias{compress_file}
\title{Compress file}
\usage{
compress_file(filename = "", remove_input = F)
}
\arguments{
\item{filename}{Character string: uncompressed input file.}

\item{remove_input}{logical: Remove the input uncompressed file.}
}
\value{
Character string: output gz compressed file.
}
\description{
Compress input file to gzip format.
}
\details{
This function work in all operating systems, but might not work with very large files,
since it reads the file to memory.
}
\examples{
\dontrun{
outfile <- compress_file(filename="data/test.fasta", remove_input=T)
}
}
