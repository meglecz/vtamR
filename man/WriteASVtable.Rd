% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{WriteASVtable}
\alias{WriteASVtable}
\title{Write ASV table}
\usage{
WriteASVtable(
  read_count_samples_df,
  outfile = "",
  asv_tax = NULL,
  sortedinfo = "",
  add_empty_samples = F,
  add_sums_by_sample = F,
  add_sums_by_asv = F,
  add_expected_asv = F,
  mock_composition = "",
  sep = ","
)
}
\arguments{
\item{read_count_samples_df}{Data frame with the following variables:
asv_id, sample, read_count, asv.}

\item{outfile}{Character string: csv file name to print the output data
frame if necessary. If empty, no file is written.}

\item{asv_tax}{Data frame or CSV file with taxonomic assignments with the following columns:
asv_id,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,
superkingdom_taxid,superkingdom,kingdom_taxid,kingdom,
phylum_taxid,phylum,class_taxid,class,order_taxid,order,
family_taxid,family,genus_taxid,genus,species_taxid,species.
pid,pcov,phit,taxn,seqn,refres,ltgres,asv.
If given, the output is completed with taxonomic assignment of each ASV.}

\item{sortedinfo}{Data frame or csv file with columns: sample, sample_type.
Only necessary if add_empty_samples==T or add_expected_asv==T.}

\item{add_empty_samples}{logical: add a column for each samples
in the original data set, even if they do not have reads after filtering.}

\item{add_sums_by_sample}{logical: add a line with the total number of reads
in each sample, and another with the number of ASVs in each sample.}

\item{add_sums_by_asv}{logical: add a column with the total number of reads
for each ASV, and another with the number of samples, where the ASV is present.}

\item{add_expected_asv}{logical: add a column for each mock sample in which
keep and tolerate ASVs are flagged.}

\item{mock_composition}{Data frame or CSV file with the following columns:
sample,action,asv. Action can take the following values: keep/tolerate.
Only necessary if add_expected_asv==T.}

\item{sep}{Field separator character in input and output csv files.}
}
\value{
Data frame corresponding to the output file.
Samples in columns, ASVs in lines, read_counts in cells, plus
other information according to the input parameters.
}
\description{
Write csv file with samples in columns, ASVs in lines, read_counts in cells.
}
\examples{
\dontrun{
WriteASVtable(read_count_samples_df, 
    outfile="out/asv_table.csv", 
    asv_tax=asv_tax, 
    sortedinfo=sortedinfo_df,
    add_empty_samples=T, 
    add_sums_by_sample=T,
    add_sums_by_asv=T, 
    add_expected_asv=T,
    mock_composition="data/mock_compostion.csv"
    )
}
}
