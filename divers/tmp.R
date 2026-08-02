
library(vtamR)
library(dplyr)

library("devtools")
library("roxygen2")
setwd("/home/meglecz/vtamR")
load_all(".")
roxygenise()
usethis::use_roxygen_md()



sessionInfo()

