% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtamR_fun.R
\name{check_dir_old}
\alias{check_dir_old}
\title{Check directory}
\usage{
check_dir_old(dir)
}
\arguments{
\item{dir}{Character string naming a directory.}
}
\value{
Character string naming a directory appended by slash if necessary.
}
\description{
Create directory if does not exists.
Add slash to the end of the directory name.
}
\examples{
\dontrun{
check_dir(dir="data")
}
}
