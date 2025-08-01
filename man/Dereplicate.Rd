% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{Dereplicate}
\alias{Dereplicate}
\title{Read all fasta files to a data frame and dereplicate}
\usage{
Dereplicate(
  sortedinfo,
  dir = "",
  outfile = "",
  sep = ",",
  asv_list = "",
  updated_asv_list = "",
  quiet = T
)
}
\arguments{
\item{sortedinfo}{Data frame or csv file with columns:
sample,replicate,fasta, optional: sample_type,habitat.
The fasta column contains the names of the fasta file to be dereplicated.}

\item{dir}{Character string naming of the directory with input fasta files.}

\item{outfile}{Character string: csv for the output data frame
(asv_id, sample, replicate, read_count). If empty, no file is written.}

\item{sep}{Field separator character in input and output csv files.}

\item{asv_list}{Character string: file with asvs and asv_ids
from earlier analyses. Optional. It is used to homogenize asv_ids
between different data sets.}

\item{updated_asv_list}{Character string naming of the output file
with updated asv_id - asv pairs. Optional.}

\item{quiet}{logical: If TRUE, suppress informational messages and only
show warnings or errors.}
}
\value{
Data frame with the following columns:
asv_id,sample,replicate,read_count,asv
}
\description{
Read all fasta files in the fasta columns of sortedinfo data frame (or csv file).
Dereplicate reads to ASVs.
Count the number of reads of each ASV in each input file.
Add a unique asv_id to each asv.
}
\details{
If asv_list is given (containing earlier asv and asv_id pairs),
uses already existing asv_id when possible for the present data
and adds new unique asv_id, to new ASVs.
If updated_asv_list is given, writes an updated file containing
all asv - asv_id pairs.
}
\examples{
\dontrun{
Dereplicate(sortedinfo=sortedinfo, dir="data/sorted", asv_list="data/asv_list.csv")
}
}
