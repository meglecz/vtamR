
fastq_dir <- system.file("extdata/demo/fastq", package = "vtamR")
fastqinfo <-  system.file("extdata/demo/fastqinfo_mfzr_plate1.csv", package = "vtamR")
outdir <- "/home/meglecz/vtamR/tmp"
cutadapt_path <- "~/miniconda3/envs/vtam/bin/cutadapt"
vsearch_path <- "~/miniconda3/envs/vtam/bin/vsearch"
blast_path <- "~/miniconda3/envs/vtam/bin/blastn"
swarm_path <- "swarm" # swarm is in the PATH
pigz_path <- "pigz"   # optional; pigz is in the PATH
check_reverse=TRUE
compress_method="R"
num_threads=0
tag_to_end=TRUE
primer_to_end=TRUE
cutadapt_error_rate=0.1
sep=","
compress=FALSE
quiet=TRUE

fastqinfo_delultiplexed <- demultiplex_fastq_pairs(fastqinfo, 
                                    fastq_dir, 
                                    outdir, 
                                    cutadapt_path="cutadapt",
                                    check_reverse=TRUE, 
                                    num_threads=0,
                                    tag_to_end=FALSE, 
                                    primer_to_end=FALSE, 
                                    cutadapt_error_rate=0.1,
                                    sep=",",
                                    compress=FALSE,
                                    quiet=T)

fastqinfo <- file.path(outdir, "fastqinfo.csv")
merged_dir <- "/home/meglecz/vtamR/merged"
sample_info <- merge_fastq_pairs(
  fastqinfo = fastqinfo,
  fastq_dir = outdir,
  outdir = merged_dir,
  vsearch_path = "vsearch",
  compress_method = "R",
  pigz_path = "pigz",
  num_threads = 0,
  fastq_ascii = 33,
  fastq_maxdiffs = 10,
  fastq_maxee = 1,
  fastq_minlen = 50,
  fastq_maxlen = 500,
  fastq_minmergelen = 50,
  fastq_maxmergelen = 1000,
  fastq_maxns = 0,
  fastq_truncqual = 10,
  fastq_minovlen = 50,
  fastq_allowmergestagger = TRUE,
  sep = ",",
  compress = FALSE,
  quiet = T
)




headers <- data.frame(
  header = as.character(),
  file = as.character()
)

for(i in 1:nrow(fastqinfo_delultiplexed)){
  
  filename <- fastqinfo_delultiplexed$fastq_fw[i]
  filename <- file.path(outdir, filename)
  file_connection <- file(filename, "r")
  # read file to a vector. Each element is a line
  file_contents <- readLines(file_connection, warn = FALSE)
  close(file_connection)
  
  # Identify lines starting with '>'
  header_indices <- grepl("^@", file_contents)
  
  fastq_headers <- as.data.frame(file_contents[header_indices])
  colnames(fastq_headers) <- c("header")
  
  fastq_headers <- fastq_headers %>%
    mutate(header = sub(" .*$", "", header)) %>%
    mutate(file = filename)
  
  headers <- rbind(headers, fastq_headers)
}

nrow(headers)
headers <- headers %>%
  group_by(header) %>%
  summarize(count = n())



read_count_input <- count_reads_in_dir(
  dir=fastq_dir, 
  pattern="_fw.fastq", 
  file_type="fastq"
)

read_count_input <- read_count_input %>%
  filter(startsWith(filename, "mfzr"))

sum(read_count_input$read_count)
sum(fastqinfo_delultiplexed$read_count)


read_count_fw <- count_reads_in_dir(
  dir="/tmp/RtmpZdJIT1/fw_178306060884", 
  pattern="_fw.fastq", 
  file_type="fastq"
)

read_count_rv <- count_reads_in_dir(
  dir="/tmp/RtmpZdJIT1/rv_178306060859", 
  pattern="_fw.fastq", 
  file_type="fastq"
)

tmp <- left_join(fastqinfo_delultiplexed, read_count_fw, by=c("fastq_fw" = "filename"))
tmp <- left_join(tmp, read_count_rv, by=c("fastq_fw" = "filename"))

tmp <- tmp %>%
  mutate(diff = read_count.x - read_count.y - read_count)

fasta_df <- read_fasta_to_df("/home/meglecz/vtamR/tmp/14ben01-1_fw.fastq")



case_3a <- read.csv("/home/meglecz/vtamR/vignettes/vtamR_demo_case3a/1_before_filter.csv")

case_3b <- read.csv("/home/meglecz/vtamR/vignettes/vtamR_demo_case3b/1_before_filter.csv")

nrow(case_3a)
nrow(case_3b)
