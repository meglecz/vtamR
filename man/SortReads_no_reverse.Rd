% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{SortReads_no_reverse}
\alias{SortReads_no_reverse}
\title{Demultiplex and trim off tags and primers without checking reverse strand}
\usage{
SortReads_no_reverse(
  fastainfo,
  fasta_dir,
  outdir = "",
  cutadapt_path = "cutadapt",
  tag_to_end = T,
  primer_to_end = T,
  cutadapt_error_rate = 0.1,
  cutadapt_minimum_length = 50,
  cutadapt_maximum_length = 500,
  sep = ",",
  compress = F,
  quiet = T
)
}
\arguments{
\item{fastainfo}{Data frame or csv file with columns:
tag_fw,primer_fw,tag_rv,primer_rv,sample,sample_type(mock/negative/real),habitat(optional),replicate,fasta}

\item{fasta_dir}{Character string: directory with input fasta files
(listed in the fasta columns of fastainfo).}

\item{outdir}{Character string: output directory.}

\item{cutadapt_path}{Character string: path to cutadapt executables.}

\item{tag_to_end}{logical: tags are at the extremity of the reads
(starting at the first base).}

\item{primer_to_end}{logical: primers follow directly the tags
(no heterogeneity spacer).}

\item{cutadapt_error_rate}{Real (0-1): maximum proportion of errors
between primers and reads (for tags, exact match is required).}

\item{cutadapt_minimum_length}{Positive integer: minimum length of the
trimmed sequence.}

\item{cutadapt_maximum_length}{Positive integer: maximum length of the
trimmed sequence.}

\item{sep}{Field separator character in input and output csv files.}

\item{compress}{logical: Compress output files to gzip format.}

\item{quiet}{logical: If TRUE, suppress informational messages and only
show warnings or errors.}
}
\value{
Data frame corresponding to the output sortedinfo.csv file
and one fasta file for each tag combination for each input fasta file,
containing trimmed reads.
}
\description{
Same as SortReads, but do not check the reverse complement of the sequences.
Demultiplex each input fasta file using the tag combinations
at the extremities of the merged reads.
Trim primers from demultiplexed reads.
}
\details{
The output sortedinfo.csv file is similar to the fastainfo file,
but the but do not have tag and primer columns.
}
\examples{
\dontrun{
fastainfo_df <- SortReads_no_reverse(fastainfo=fastainfo_df, 
    fasta_dir="data/fasta", 
    outdir="data/sorted", 
    tag_to_end=T, 
    primer_to_end=T, 
    cutadapt_minimum_length=300,
    cutadapt_maximum_length=350, 
    sep=","
    )
}
}
