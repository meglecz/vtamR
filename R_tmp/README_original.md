# vtamR

vtamR is functional and ready to be tested, but the R package has not yet been released.

At the moment, you have to clone de [github repository](https://github.com/meglecz/vtamR) and install third party programs. Please, report bugs and add comments in the [issues](https://github.com/meglecz/vtamR/issues). The release of a package its on it way.


## Install third party programs 

### Linux

* `vsearch` [https://github.com/torognes/vsearch](https://github.com/torognes/vsearch)
* `BLAST` [https://www.ncbi.nlm.nih.gov/books/NBK52640/](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
* `cutadapt` [https://cutadapt.readthedocs.io/en/stable/installation.html](https://cutadapt.readthedocs.io/en/stable/installation.html)
* `swarm` [https://github.com/torognes/swarm](https://github.com/torognes/swarm)
* `git` [https://github.com/git-guides/install-git](https://github.com/git-guides/install-git)

### Windows

Download binaries and save it a convenient place on your computer (path without space, e.g. `C:/Users/Public/`)

**vsearch**

* Download binaries from https://github.com/torognes/vsearch/releases/tag/v2.23.0
* Decompress the zip file
* The executable are found in `vsearch-x.xx.x-win-x86_64/bin/`

**cutadapt**

* Download `cutadapt.exe` from [https://github.com/marcelm/cutadapt/releases](https://github.com/marcelm/cutadapt/releases)

**swarm**

* Download binaries from [https://github.com/torognes/swarm/releases](https://github.com/torognes/swarm/releases)
* Decompress the zip file
* The executable is in the `swarm-x.x.x-win-x86_64/bin` directory

**BLAST**

* Detailed instructions: [https://www.ncbi.nlm.nih.gov/books/NBK52637/](https://www.ncbi.nlm.nih.gov/books/NBK52637/)
* Download executable (ncbi-blast-x.xx.x+-win64.exe) from [https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)
* Double click on the .exe file, accept the license agreement and specify the install location when prompted. Attention! Do not accept the standard location (C:/Program Files/...) since it contains a space. Chose a directory with a path without space (e.g. `C:/Users/Public/`).
* The executables are found in `blast-x.xx.x+/bin/`

**git**

* `git` [https://github.com/git-guides/install-git](https://github.com/git-guides/install-git)



## Clone github repository

```{bash clone}
git clone https://github.com/meglecz/vtamR.git
cd vtamR
```


## Install R packages

```{r install_pck}
install.packages("devtools")
install.packages("roxygen2")
install.packages("seqinr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
```

## Test your installation

**Set path to third party programs and specify the number of CPUs**

Adapt the paths according to your installation.

```{r set_path, eval=TRUE}
cutadapt_path <- "C:/Users/Public"
vsearch_path <- "C:/Users/Public/vsearch-2.23.0-win-x86_64/bin"
blast_path <- "C:/Users/Public/blast-2.14.1+/bin/"
swarm_path <- "C:/Users/Public/swarm-3.1.4-win-x86_64/bin"
num_threads <- 4
vtam_dir <- "C:/Users/emese/vtamR"
setwd(vtam_dir)
```


**Load libraries**
```{r load_libraries}
library("devtools")
library("roxygen2")
library("seqinr")
library("dplyr")
library("tidyr")
library("ggplot2") 
load_all(".")
roxygenise()
usethis::use_roxygen_md()
```


**Run test functions on test data**
```{r test_vtamR}
Test_MergeSortReads(vsearch_path=vsearch_path, cutadapt_path=cutadapt_path)
Test_Filters(vsearch_path=vsearch_path, swarm_path=swarm_path)
Test_MakeKnownOccurrences()
Test_Optimize(vsearch_path=vsearch_path)
Test_TaxAssign(blast_path=blast_path, num_threads=num_threads)
```
You should see "PASS" after the tested function names.

## Troubleshooting

### install devtools

`install.packages("devtools")` can lead to many error messages. Try to follow these.
Most of the time you need to install separately some of the dependencies using `install.packages("pkg_name")` commands.
On linux systems, sometimes it is necessary to install applications in a 
terminal and not in R.

