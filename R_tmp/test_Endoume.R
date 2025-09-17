if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install("meglecz/vtamR@develop")

environment_vtamR_yml_path <- system.file("environment_vtamR.yml", package = "vtamR")
cat(environment_vtamR_yml_path)

# In a terminal
### uriginal yml try to install many R packages with incompatibilities
# environment_vtamR_modif.yml contains only TPP s adependencies and making the env goes OK
conda env create -f /home/emese/Bureau/environment_vtamR_modif.yml   

conda activate vtamRenv

which vsearch
which blastn
which cutadapt
which swarm

/home/emese/anaconda3/envs/vtamRenv/bin/vsearch
/home/emese/anaconda3/envs/vtamRenv/bin/blastn
/home/emese/anaconda3/envs/vtamRenv/bin/cutadapt
/home/emese/anaconda3/envs/vtamRenv/bin/swarm

### in R
library(vtamR)
download_osf(filename = "COInr_for_vtam_2025_05_23_dbV5.tar.gz",
            url = "https://osf.io/download/jyhz6/",
            dest_dir = "COInr_db",
            untar = TRUE,
            quiet = FALSE
)



################################################"

library(dplyr)
library(ggplot2)

### set up
cutadapt_path <- "/home/emese/anaconda3/envs/vtamRenv/bin/cutadapt"
vsearch_path <- "/home/emese/anaconda3/envs/vtamRenv/bin/vsearch"
blast_path <- "/home/emese/anaconda3/envs/vtamRenv/bin/blastn"
swarm_path <- "/home/emese/anaconda3/envs/vtamRenv/bin/swarm"
num_threads <- 8
sep <- ","
outdir <- "~/vtamR_demo_out"

fastq_dir <- system.file("extdata/demo/fastq", package = "vtamR")
fastqinfo <-  system.file("extdata/demo/fastqinfo.csv", package = "vtamR")
mock_composition <-  system.file("extdata/demo/mock_composition.csv", package = "vtamR")
asv_list <-  system.file("extdata/demo/asv_list.csv", package = "vtamR")
taxonomy <- system.file("extdata/db_test/taxonomy_reduced.tsv", package = "vtamR")
blast_db <- system.file("extdata/db_test", package = "vtamR")
blast_db <- file.path(blast_db, "COInr_reduced")


### Merge
merged_dir <- file.path(outdir, "merged")
fastainfo_df <- Merge(fastqinfo, 
                      fastq_dir=fastq_dir, 
                      vsearch_path=vsearch_path, 
                      outdir=merged_dir,
                      fastq_maxee=1,
                      fastq_maxns=0,
                      fastq_allowmergestagger=F,
                      num_threads=num_threads
                      
)

### demultiplex

sorted_dir <- file.path(outdir, "sorted")
sortedinfo_df <- SortReads(fastainfo_df, 
                           fasta_dir=merged_dir, 
                           outdir=sorted_dir, 
                           check_reverse=TRUE, 
                           cutadapt_path=cutadapt_path, 
                           vsearch_path=vsearch_path,
                           num_threads=num_threads
)

### Error
"/home/emese/anaconda3/envs/vtamRenv/bin/cutadapt --cores=8 --quiet -e 0.1 --no-indels --trimmed-only --minimum-length 50 --maximum-length 500 -g ^TCCACTAATCACAARGATATTGGTAC...GGAGGATTTGGWAATTGATTAGTW$ --output ~/vtamR_demo_out/sorted/14ben01-2.fasta /tmp/Rtmp0F5NsC/mfzr_2_fw.fasta_175811470397/tagtrimmed-GTCGATCATGTCA-ACATCGACGTACG.fasta"

### PB with cutadapt version, which is 2.6 in vtalRenv => see how to chaang this => add a teste of version to install vignette.

### OK
rc_dir <- paste('rc_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
rc_dir <- file.path(tempdir(), rc_dir)
check_dir(rc_dir)
