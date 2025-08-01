% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{run_swarm}
\alias{run_swarm}
\title{Run Swarm on one file}
\usage{
run_swarm(
  read_count_df,
  swarm_path = "swarm",
  num_threads = 1,
  swarm_d = 1,
  fastidious = T,
  quiet = T
)
}
\arguments{
\item{read_count_df}{Data frame with the following variables:
asv, plate, marker, sample, replicate, read_count.}

\item{swarm_path}{Character string: path to swarm executables.}

\item{num_threads}{Positive integer: Number of CPUs.}

\item{swarm_d}{Positive integer: d parameter for swarm.
Maximum number of differences allowed between two ASVs,
meaning that two ASVs will be grouped if they have d (or less) differences.}

\item{fastidious}{logical: when working with d = 1, perform a second
clustering pass to reduce the number of small clusters.}

\item{quiet}{logical: If TRUE, suppress informational messages and only
show warnings or errors.}
}
\value{
Data frame with the same structure as the input, but ASVs of
the same cluster pooled to one row.
}
\description{
Runs swarm \url{https://github.com/torognes/swarm}
on input read_count data frame or csv file.
Pools variants of the same cluster and sums read counts of the
underlying ASVs.
}
\details{
\code{run_swarm} runs swarm on the whole data set in one go, while
\code{Swarm} can run swarm sample by sample or in one go.
}
\examples{
\dontrun{
run_swarm(read_count=read_count_df, swarm_path=swarm_path, num_threads=4)
}
}
