% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{read_blast_res}
\alias{read_blast_res}
\title{Read BLAST results to data frame}
\usage{
read_blast_res(file)
}
\arguments{
\item{file}{Character string naming the output of BLAST.
Tab separated colums: qseqid,pident,qcovhsp,staxids}
}
\value{
Data frame with the following columns: qseqid,pident,qcovhsp,staxids
}
\description{
Read BLAST result to a data frame.
If more than one taxid for a hit, make a separate line for each taxid.
}
\examples{
\dontrun{
read_blast_res(file="blastout.txt")
}
}
