
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vtamR

<!-- badges: start -->
<!-- badges: end -->

**vtamR** is a revised, completed version of
[VTAM](https://www.csbj.org/article/S2001-0370(23)00034-X/fulltext)
(Validation and Taxonomic Assignation of Metabarcoding Data) rewritten
in R. It is a complete metabarcoding pipeline:

- Sequence analyses **from raw fastq files** of amplicon sequences till
  Amplicon Sequence Variant
  ([ASV](https://people.imbe.fr/~emeglecz/vtamR/tutorial-vtamr-pipeline.html#glossary))
  table of **validated ASVs assigned to taxonomic groups**.
- Handles technical or biological **replicates** of the same sample.
- Uses positive and negative **control samples to fine tune the
  filtering** and reduce [false
  positive](https://people.imbe.fr/~emeglecz/vtamR/tutorial-vtamr-pipeline.html#glossary)
  and [false
  negative](https://people.imbe.fr/~emeglecz/vtamR/tutorial-vtamr-pipeline.html#glossary)
  occurrences.
- Can pool multiple data sets (results of earlier analyses)
- Can pool results from overlapping markers

**Novelties compared to VTAM:**

- As it is a series of R functions, `vtamR` is highly adaptable to
  include/exclude and order different steps of the analyses
- Includes SWARM for denoising
- Graphic options
- Include functions to get statistics of each filtering steps (read and
  variant count etc.)
- ASV clustering to mOTU
- The notion of marker and run has been dropped to simplify the analyses

# Tutorial

A detailed **[Tutorial](https://people.imbe.fr/~emeglecz/vtamR)** on how
to construct a full metabarcoding pipeline is avalibale online.

# Installation

## Implementation

The heavy lifting of sequence analyses are done by the following third
party programs. They should be installed on your computer:

- BLAST([Altschul et al.,
  1990](https://pubmed.ncbi.nlm.nih.gov/2231712/)) (v2.11.0+)
- vsearch ([Rognes et al., 2016](https://peerj.com/articles/2584/))
  (v2.7.0)
- cutadapt([Martin,
  2011](https://journal.embnet.org/index.php/embnetjournal/article/view/200/479))
  (v4.0)
- swarm ([Mahé et al., 2015](https://peerj.com/articles/1420/))
  (v2.1.12)
- pigz (optional)

`vtamR` has been tested using the above mentioned versions, but it
should work with later versions. `vtamR` was tested on **Windows** and
**Linux**, and should work in all operating systems.

## Install vtamR

Install `vtamR` and its CRAN dependencies by `pak`.

``` r
if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install("meglecz/vtamR@develop", dependencies=TRUE)
```

Alternatively, you can use `devtools`.

``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("meglecz/vtamR@develop", dependencies = TRUE)
```

You also need to install Bioconductor packages manually if you intend to
use `TaxAssignRDP` function.

``` r
# Install Bioconductor suggested packages manually if needed (only for TaxAssignRDP)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("rRDP", quietly = TRUE)) BiocManager::install("rRDP")
if (!requireNamespace("rRDPData", quietly = TRUE)) BiocManager::install("rRDPData")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
```

Read the **[Manual](https://people.imbe.fr/~emeglecz/vtamR)** on how to
construct a full metabarcoding pipeline.

## Install third party programs

### Linux

- `vsearch` <https://github.com/torognes/vsearch>
- `BLAST` <https://www.ncbi.nlm.nih.gov/books/NBK52640/>
- `cutadapt`
  <https://cutadapt.readthedocs.io/en/stable/installation.html>
- `swarm` <https://github.com/torognes/swarm>
- `pigz` (optional) <https://github.com/madler/pigz>

**Add the executables to your PATH** so you don’t need to type their
full location each time. If you prefer **not** to add them, make sure to
provide the full path to each program through the appropriate `vtamR`
function argument (e.g., `vsearch_path`).

You can add a directory to your PATH so that R and RStudio can find the
executables you want to use.

1.  Edit (or create) your `~/.Renviron` file

``` bash
nano ~/.Renviron
```

2.  Add the directory to the PATH

``` bash
PATH="/home/user/tools/myprogram/bin:${PATH}"
```

**Note:** The PATH must contain a **directory**, not an executable file.
For example, use: `/home/username/miniconda3/envs/vtamr/bin` and
**not**: `/home/username/miniconda3/envs/vtamr/bin/cutadapt`

3.  Restart R or RStudio

The new PATH is loaded at startup.

4.  Verify that R sees the correct executables

``` r
system("cutadapt --version")
system("swarm --version")
system("vsearch --version")
system("blastn -version")
system("pigz --version")
```

### Windows

Download binaries and save the to a convenient place on your computer
(path without space, e.g. `C:/Users/Public/`)

**vsearch**

- Download binaries from
  <https://github.com/torognes/vsearch/releases/tag/v2.23.0>
- Decompress the zip file
- The executable are found in `vsearch-x.xx.x-win-x86_64/bin/`

**cutadapt**

- Download `cutadapt.exe` from
  <https://github.com/marcelm/cutadapt/releases>

**swarm**

- Download binaries from <https://github.com/torognes/swarm/releases>
- Decompress the zip file
- The executable is in the `swarm-x.x.x-win-x86_64/bin` directory

**BLAST**

- Download the installer (ncbi-blast-x.xx.x+-win64.exe) from:
  <https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/>
- Double-click the .exe file, accept the license agreement, and keep the
  suggested installation path. This will automatically add BLAST to your
  system PATH, so you won’t need to specify its location when using it.
- The executables will be installed in:
  `C:\Program Files\NCBI\blast-x.xx.x+\bin`

**pigz** (Optional)

- Download
  [pigz-win32.zip](https://sourceforge.net/projects/pigz-for-windows/files/)
- Decompress the zip file
- The executable (pigz.exe) is in the pigz-win32 directory

**Add executables to your PATH**

This step is optional. If third-party executables are included in your
system PATH, you won’t need to manually specify their locations when
running functions that use them.

If you prefer **not** to add them, make sure to provide the full path to
each program through the appropriate `vtamR` function argument (e.g.,
`vsearch_path`).

You can use either **GUI** or **PowerShell** to add executables to your
PATH.

1.  Using the Windows GUI (System Settings)

- Press **Win + R**, type `sysdm.cpl`, and press **Enter**. *(Opens
  System Properties.)*
- Go to the **Advanced** tab.
- Click **Environment Variables…**
- In the **System variables** section, find and select **Path**.
- Click **Edit…**
- Click **New**, then paste the directory you want to add (for example:
  `C:\Users\YourName\Tools\bin`)
- Click **OK** on all windows to save.

2.  Using PowerShell

- Press **Win + R**, type `powershell`, and press **Enter**. *(Opens
  Powershell)*

<!-- -->

    setx PATH "$($env:Path);C:\My\New\Folder\vsearch:C:\My\New\Folder\cutadapt;C:\My\New\Folder\pigz;C:\My\New\Folder\swarm"

- Close de Powershell, since the change applies only to **new**
  terminals, not the current one.

**Verify Program Versions**

Open a Command Prompt on Windows (**Win + R**, type `cmd`)

``` bash
cutadapt --version
swarm --version
vsearch --version
pigz --version
blastn -version
```

**Note** *These commands will display the version of each program **if
the executables are in your PATH**. Otherwise, don’t forget to provide
the full path to the executables.*

Example:

    C:/Users/Public/vsearch-x.xx.x-win-x86_64/bin/vsearch --version

## TaxAssign reference data base

### TaxAssignLTG

For `TaxAssignLTG` (BLAST based Lowest Common Ancestor method), a
ready-to-use COI database is available from OSF
[OSF](https://osf.io/vrfwz/), ([Meglécz,
2023](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756)).

Download and extract using:

``` r
library(vtamR)

download_osf(filename = "COInr_for_vtam_2025_05_23_dbV5.tar.gz",
    url = "https://osf.io/download/jyhz6/",
    dest_dir = "COInr_db",
    untar = TRUE,
    quiet = FALSE
    )
```

To check for newer versions:

- Go to: <https://osf.io/vrfwz/files/osfstorage>
- In **Files** open the the **dbV5** folder.
- Choose your version → Click on the tree-dot menu → Right-click
  Download → Copy Link Address.

Alternatively, you can also create a custom version from the **TSV
format** database available at
[Zenodo](https://zenodo.org/records/15515860) and make a custom version
using [mkCOInr](https://github.com/meglecz/mkCOInr) ([Meglécz,
2023](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756)).

**For other markers** you will need a [database formatted to
BLAST](tutorial-vtamr-pipeline.html#reference-database-for-taxonomic-assignments),
containing taxIDs and a [taxonomy
file](tutorial-vtamr-pipeline.html#reference-database-for-taxonomic-assignments).
If all sequences are extracted from NCBI-nt. The taxonomy file provided
in [OSF](https://osf.io/vrfwz/) can be used without modification.

### TaxAssignRDP

`TaxAssignRDP` function uses de RDP algorithm to assign bacterial 16S
sequences. The database is provided in the `rRDPData` R data package.

It is also possible to provide to provide a trained classifier object
created with `trainRDP()`.
