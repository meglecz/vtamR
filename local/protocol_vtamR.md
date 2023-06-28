# vtamR

## get test data

~~~
mkdir vtamR
cd vtamR

wget -nc https://github.com/aitgon/vtam/releases/latest/download/fastq.tar.gz -O fastq.tar.gz
tar zxvf fastq.tar.gz
rm fastq.tar.gz
~~~


## merge

~~~
conda activate vtam
vtam merge --fastqinfo user_input/fastqinfo_mfzr.tsv --fastqdir fastq --fastainfo mfzr/fastainfo.tsv --fastadir mfzr/merged -v

vtam sortreads --fastainfo mfzr/fastainfo.tsv --fastadir mfzr/merged --sorteddir mfzr/sorted -v
~~~



## make a package

https://ourcodingclub.github.io/tutorials/writing-r-package/

~~~
install.packages("devtools")
install.packages("roxygen2")
~~~

make a directory

~~~
mkdir -p ~/vtamR/R
~~~
Make an R file with the funtions and put it to /home/meglecz/vtamR/R/vtamR_functions.R

Create a new file called `DESCRIPTION` in the `vtamR` directory 
~~~
touch  ~/vtamR/DESCRIPTION
~~~
Write minimal descriton to this file
~~~
Package: vtamR
Type: Package
Title: Validation and Taxonomic Assigment of Metabarcoding data with R
Version: 0.0.1.0
~~~

make a new project from ~/vtamR. this creates .RprojUser dir .Rbuildignore, vtamR.Rproj.

To ignore the local dir that contains my prinate file, complete   .Rbuildignore

~~~
^.*\.Rproj$
^\.Rproj\.user$
local*
~~~



## load the package and make help

~~~
setwd("~/vtamR")
library(devtools)
load_all(".")

library(roxygen2); # Read in the roxygen2 R package
roxygenise();      # Builds the help files
~~~

## documentation

https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html


## changes between vtam and vtamR

run => plate
vtamR can deal with just one marker and plate in one go



## versioning with git

https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github

~~~
cd ~/vtamR

git init -b main
git add .

git config --global user.email emese.meglecz@imbe.fr
git config --global user.name meglecz

git commit -m "First commit"

~~~

create a repository on github
https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository






