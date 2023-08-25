#' run_blast
#' 
#' Run BLAST
#'  
#' @param seqs vector containing the quary sequences
#' @param blast_db BLAST DB incuding path
#' @param blast_path path to BLAST executables
#' @param qcov_hsp_perc minimum query coverage
#' @param perc_identity minimum percentage of identity
#' @param num_threads number of threads
#' @export
#'

run_blast <- function(seqs, blast_db="", blast_path="", out="", qcov_hsp_perc=70, perc_identity=70, num_threads=1){
  
  # make fasta file with unique reads; use sequences as ids
  outdir <- outdir_tmp
  fas <- paste(outdir, 'unique.fas', sep="")
  write_fasta(seqs, fas, seq_as_id=F)
  
  task = "megablast"
  e = 1e-20
  dust = "yes"
  max_target_seqs=500
  
  
  blast <- paste(blast_path, "blastn -task ", task, " -db ",blast_db ," -query ",fas," -evalue ",e," -out ",out," -outfmt '6 qseqid sseqid pident length qcovhsp staxids evalue' -dust ",dust," -qcov_hsp_perc ",qcov_hsp_perc," -perc_identity ",perc_identity," -num_threads ",num_threads," -max_target_seqs ",max_target_seqs, sep="")
  system(blast)
}