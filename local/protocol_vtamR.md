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


~~~
cd ~/vtamR
git remote add origin https://github.com/meglecz/vtamR.git
git remote -v
~~~





## Windows

### vsearch

Download binaries from

https://github.com/torognes/vsearch/releases/download/v2.22.1/vsearch-2.22.1-win-x86_64.zip

Dezipped to C:

It is possible to ron from a terminal using vsearch.exe at the beginning of the command

 





# make vtam test data

Use unzipped files



~~~
conda activate vtam_2
cd /home/meglecz/vtamR/

vtam merge --fastqinfo vtamR_test/vtam/user_input/fastqinfo_mfzr.tsv --fastqdir vtamR_test/data --fastainfo vtamR_test/vtam/merged/fastainfo.tsv --fastadir  vtamR_test/vtam/merged

vtam sortreads --fastainfo vtamR_test/vtam/merged/fastainfo.tsv --fastadir vtamR_test/vtam/merged --sorteddir vtamR_test/vtam/sorted
~~~







!!!!!! not use yet

~~~
vtam filter --db vtamR_test/vtam/db.sqlite --sortedinfo vtamR_test/vtam/sorted/sortedinfo.tsv --sorteddir vtamR_test/vtam/sorted --asvtable vtamR_test/vtam/filter/asvtable_default.tsv
~~~



~~~
!!!!!!!!!!!!!!!!! TODO  vtam make_known_occurrences --asvtable asvtable filter/asvtable_default.tsv --sample_types asper1/user_input/sample_types.tsv --mock_composition asper1/user_input/mock_composition_mfzr.tsv --known_occurrences asper1/run1_mfzr/known_occurrences_mfzr.tsv --missing_occurrences asper1/run1_mfzr/missing_occurrences_mfzr.tsv -v
~~~





## Windows

## install miniconda => very much unix centered => give up

https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html#installing-on-windows

Points 1-5

### install dependencies to a conda env

~~~~
conda create --name vtamR -y
conda activate vtamR

conda install -c bioconda blast -y
conda install -c bioconda vsearch -y
conda install cutadapt -y
~~~~

No blast available for conda in windows 

No vsearch available for conda in windows

No cutadapt available for vsearch in windows

**Forget conda on windows**



## pip

~~~
pip install blast 
blastn
Nucleotide-Nucleotide BLAST 2.13.0+

pip install cutadapt
Cutadapt --version 
4.5

pip install vsearch => only version 1.1 and it does not work
~~~
OK and easy for blast and cutadapt, but not for vsearch + python and pip should be installed

## istall dependencies individually


### vsearch
Download binaries from https://github.com/torognes/vsearch/releases/tag/v2.23.0
Extraire zip file => use 
C:/Users/Public/vsearch-2.23.0-win-x86_64/bin/vsearch to start a command

### BLAST
Detailed instructions: https://www.ncbi.nlm.nih.gov/books/NBK52637/
Download executable (ncbi-blast-x.xx.x+-win64.exe) from https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

Double click on the exe file, accept the license agreement and specify the install location in a new prompt. Attention! Do not accept the standard location (C:/Program Files/...) since it contains a space. Chose a directory with a path without space.

C:/Users/Public/blast-2.14.1+/bin/blastn

### cutadapt

download cutadapt.exe and save it a convenient place on your computer

C:/Users/Public/cutadapt to start a command



## swarm

download binaries from https://github.com/torognes/swarm/releases

unzip

The executable is in the C:\swarm-3.1.4-win-x86_64\bin directory





# File compression

## windows

### gzip files 

- can be read and written in windows
- vsearch cannot deal with them directly in windows, 
- the output of vsearch is uncompressed
- cutadapt can deal with gzip files directly and produces output in unzipped or gzipped format according to the file extension.

### zip files

-  can be read and created by R but it is a bit difficult to handle file path
- if path is given to zip => the zip archive contains reproduce the directory structure, even it is a single file. This is annoying when dezipping
- it is possible to overcome of this problem by changing the wd, but it is tiresome
- wen zipping and dezipping files, the name is not always coherent. The orignal file extention is usually replaced instead of appended by .zip. this makes it more difficult to guess the name of the dezipped files
- HTS sequencing does not come in zip format usually

## linux

### gzip files 

- easy to handle
- vsearch and cutadapt can use them directly
- output is uncompressed for vsearch
- output in unzipped or gzipped format according to the file extension
- Check if bz2 files are handled by vsearch, cutadapt and R!!!! 

### zip files

- same problems as in windows

## Conclusion
 - vtamR does not support zipped files
 - On windows an input gz files are dezipped and passed to consecutive steps. It is tested automatically if it is windows and if input is gz
 - on linux, gzip files are handled directly by vsearch and cutadapt
 - the compress option (T/F), handles wether output files should be compressed or not, for cutadapt output on linux, check if compression corresponds to the output, and change it if necessary





# Rmarkdown

https://ourcodingclub.github.io/tutorials/rmarkdown/











