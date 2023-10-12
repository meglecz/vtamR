#' run_blast
#' 
#' Run BLAST
#'  
#' @param seqs vector containing the query sequences
#' @param blast_db BLAST DB incuding path
#' @param blast_path path to BLAST executables
#' @param outdir output directory
#' @param qcov_hsp_perc minimum query coverage
#' @param perc_identity minimum percentage of identity
#' @param num_threads number of threads
#' @export
#'
run_blast <- function(df, blast_db="", blast_path="", outdir="", qcov_hsp_perc=70, perc_identity=70, num_threads=8){
 
   # make fasta file with unique reads; use numbers as ids  
  seqs <- unique(df$asv)
  fas <- paste(outdir, 'unique.fas', sep="")
  write_fasta(seqs, fas, seq_as_id=F)
  # define the name of the output file
  blast_out <- paste(outdir, 'blast.out', sep="")
  
  task = "megablast"
  e = 1e-20
  dust = "yes"
  max_target_seqs=500
  
  
#  blast <- paste(blast_path, "blastn -task ", task, " -db ",blast_db ," -query ",fas," -evalue ",e," -out ",blast_out," -outfmt '6 qseqid pident qcovhsp staxids' -dust ",dust," -qcov_hsp_perc ",qcov_hsp_perc," -perc_identity ",perc_identity," -num_threads ",num_threads," -max_target_seqs ",max_target_seqs, sep="")
  blast <- paste(blast_path, "blastn -task ", task, " -db ",blast_db ," -query ",fas," -evalue ",e," -out ",blast_out,' -outfmt "6 qseqid pident qcovhsp staxids" -dust ',dust," -qcov_hsp_perc ",qcov_hsp_perc," -perc_identity ",perc_identity," -num_threads ",num_threads," -max_target_seqs ",max_target_seqs, sep="")
    system(blast)
  
  # read BLAST results 
  blast_res <- read.delim(blast_out, header=F, sep="\t", fill=T, quote="")
  colnames(blast_res) <- c("qseqid","pident","qcovhsp","staxids") 
  return(blast_res)
}

#' delete_1_by_row
#' 
#' Delete all 1 from the beginning of a raw , shift the other values and replace de missing ones at the end by NA
#'  
#' @param row vector
#'
delete_1_by_row <- function(row) {
  
  n <- length(row) 
  # Remove all occurrences of 1
  row <- row[row != 1]
  
  # Create a new row with NA at the end
  new_row <- c(row, rep(NA, n - length(row)))
  
  return(new_row)
}
#' get_lineage_ids
#' 
#' Get the complete taxonomic lineage of each taxID in the input vector; return taxIDs of each taxa in the lineage
#'  
#' @param taxids vector taxIDs (taxonomic IDs)
#' @param tax_df data frame with the following columns: tax_id, parent_tax_id, rank, name_txt, taxlevel (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 2: kingdom, 1: superkingdom, 0: root)
#' @export
#'
get_lineage_ids <- function(taxids, tax_df){
  
  # taxids is a vector of taxids; there can be duplicated values
  lineages <- as.data.frame(taxids)
  colnames(lineages) <- c("tax_id")
  
  # if input is dataframe with blast output
  #lineages <- df_intern %>%
  #  select(staxids) %>%
  #  rename(tax_id=staxids)
  
  i <- 1
  while(i < 100){
    # use i as name instead of tax_id
    new_colname <- as.character(i)
    # add parent_tax_id and rename columns
    lineages <- left_join(lineages, tax_df, by="tax_id")%>%
      select(-rank, -name_txt, -taxlevel) %>%
      # !! = interprent the variable
      rename(!!new_colname :=tax_id, "tax_id"=parent_tax_id)
    
    i <- i+1
    # stop if all lines has the same value (usually 1)
    tid_list <- unique(lineages$tax_id)
    if(length(tid_list) == 1 && tid_list[1] ==1){
      break
    }
  }
  # delete the last column, where all values are 1
  lineages <- lineages %>%
    select(-tax_id)
  # reverse order of columns
  lineages <- lineages[, ncol(lineages):1]
  # Apply the function to each row of the lineages data frame: 
  # delete all 1, shift the remaining elements of each row to the beginning, 
  # and replace missing values at the end of the row bu NA
  lineages <- as.data.frame(t(apply(lineages, 1, delete_1_by_row)))
  # add as a first column the taxid, so thay can easily me accessed
  lineages <- cbind(taxids, lineages)
  
  return(lineages)
}
#' make_ltg
#' 
#' Determine the Lowest Taxonomic Group (LTG) that contains phit percentage of the input taxids
#'  
#' @param taxids vector taxIDs (taxonomic IDs); there can be duplicated values
#' @param lineages data frame: taxids in the first column followed by their lineages represented taxids (starting at the lowest resolution)
#' @export
#'
make_ltg <- function(taxids, lineages, phit=70){
  # taxids is a vector of taxids; there can be duplicated values
  # mak a data frame from the vector
  lin <- as.data.frame(taxids)
  colnames(lin) <- c("staxid")
  
  # add lineage to each taxid
  lin <- left_join(lin, lineages, by=c("staxid" = "taxids")) %>%
    select(-where(~all(is.na(.)))) # delete columns if all elements are NA
  
  ltg <- NA
  if(length(unique(lin$staxid)) == 1){ # only one taxid among hits; avoid loops
    ltg <-lin$staxid[1]
  }else{
    for(i in 2:ncol(lineages)){ # start from low resolution
      tmp <- as.data.frame(lineages[,i])
      colnames(tmp) <- c("tax_id")
      # get unique taxids, and their numbers in the i-th column
      tmp <- tmp %>%
        group_by(tax_id) %>%
        summarize(nhit=length(tax_id)) %>%
        arrange(desc(nhit))
      
      # if the taxid with the highest number of sequences does not contain at least phit percent of the hits, stop
      max_hitn <- as.numeric(tmp[1,"nhit"])
      if(max_hitn/sum(tmp[,"nhit"]) < phit/100){
        break
      }
      ltg <- as.numeric(tmp[1,"tax_id"])
    }
  }
  return(ltg)
}

#' update_taxids
#' 
#' Replaces old taxids by valid ones
#'  
#' @param df data frame with the following columns: qseqid,pident,qcovhsp,staxids
#' @param old_taxid data frame with the following columns:  tax_id,old_tax_id
#' @export
#'

update_taxids <- function(df, old_taxid){
  # df is a data frame with the following columns: qseqid,pident,qcovhsp,staxids
  # old_taxid is a data frame with the following columns:  tax_id,old_tax_id
  
  # replace old taxids (if any) in df by up to date ones 
  df <- left_join(df, old_taxid, by=c("staxids" = "old_tax_id"))
  df$staxids[which(!is.na(df$tax_id))] <- df$tax_id[which(!is.na(df$tax_id))]
  # delete tax_id column since the values (if non NA were used to replace staxids)
  df <- df %>%
    select(-tax_id)
  return(df)
}

#' get_ranked_lineages
#' 
#' Returns  a data frame with the ranked lineages of taxids (columns: ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,
#' superkingdom_taxid,superkingdom,kingdom_taxid,kingdom,phylum_taxid,phylum,class_taxid,class,order_taxid,order,
#' family_taxid,family,genus_taxid,genus,species_taxid,species,)
#'  
#' @param taxids vector of taxIDs
#' @param tax_df data frame with the following columns: tax_id, parent_tax_id, rank, name_txt, taxlevel (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 2: kingdom, 1: superkingdom, 0: root)
#' @export
#'
get_ranked_lineages <- function(taxids, tax_df){
  
  # taxids is a vector of taxids; there can be duplicated values
  ranked_lineages <- as.data.frame(taxids)%>%
    filter(!is.na(taxids)) %>%
    rename(tax_id=taxids)
  # make tmp data frame to keep a list of taxids
  tmp <- ranked_lineages
  # make tmp_lin data frame to keep a list of taxids and the lineage of each taxid (including, names, taxid, taxlevel)
  tmp_lin <- ranked_lineages
  
  # define first colums, with taxid, name, taxlevel
  ranked_lineages <- left_join(ranked_lineages, tax_df, by="tax_id")%>%
    select(-parent_tax_id)%>%
    rename(ltg_taxid=tax_id, ltg_name=name_txt, ltg_rank=rank, ltg_rank_index=taxlevel) %>%
    select(ltg_taxid, ltg_name, ltg_rank, ltg_rank_index)
  # add columns for each major taxlevel (taxid and name)
  now_cols <- c(
    "superkingdom_taxid", "superkingdom",
    "kingdom_taxid", "kingdom",
    "phylum_taxid", "phylum",
    "class_taxid", "class",
    "order_taxid", "order",
    "family_taxid", "family",
    "genus_taxid", "genus",
    "species_taxid", "species"
  )
  ranked_lineages[now_cols] <- NA
  
  # get linegaes of each taxid
  i <- 1
  while(i < 100){
    # get the tax name, and tax rank for each taxid
    tmp <- left_join(tmp, tax_df, by="tax_id")
    # tock info in tmp_lin
    tmp_lin <- cbind(tmp_lin, tmp$tax_id, tmp$name_txt, tmp$rank )
    # re-initilize tmp
    tmp <- tmp %>%
      select(parent_tax_id)%>%
      rename(tax_id=parent_tax_id)
    # stop if all linage ends with root
    if(all(tmp$tax_id ==1)){
      break
    }
    i<- i+1
  }
  
  # select only major taxonomic ranks from each line of tmp_lin; keep the results in ranked_lineages
  for (c in seq(from=6, to=ncol(ranked_lineages), by=2)) {# go though all major taxlevel
    taxrank <- colnames(ranked_lineages[c])
    for (i in 1:nrow(tmp_lin)) {
      row <- tmp_lin[i, ]  # Extract the current row
      col_index <- which(row == taxrank)  # Find the column index containing "species"
      
      if (length(col_index) > 0) {
        # Add taxon name and taxid to ranked_lineages
        ranked_lineages[i,c-1] <- tmp_lin[i,col_index-2]
        ranked_lineages[i,c] <- tmp_lin[i,col_index-1]
      }
    }
  }
  return(ranked_lineages)
}

#' adjust_ltgres
#' 
#' If the ltg has a higher resolution than the ltgres parameter, adjust the ltg and stop lineage at ltgres level; Returns the adjusted data frame
#'  
#' @param taxres_df data frame with the following columns: asv,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,superkingdom_taxid,
#' superkingdom,kingdom_taxid,kingdom,phylum_taxid,phylum,class_taxid,class,order_taxid,order,family_taxid,family,genus_taxid,genus,species_taxid,species,pid,
#' pcov,phit,taxn,seqn,refres,ltgres
#' @param tax_df data frame with the following columns: tax_id, parent_tax_id, rank, name_txt, taxlevel (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 2: kingdom, 1: superkingdom, 0: root)
#' @export
#'
adjust_ltgres <- function(taxres_df, tax_df){
  
  # link taxlevel index and tax rank
  taxlevel_index = data.frame(taxlevel_index=c(1,2,3,4,5,6,7,8),
                              taxrank=c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  )
  
  # add the name of the tax rank equivalent to the index in ltgref
  taxres_df <- left_join(taxres_df, taxlevel_index, by=c("ltgres" = "taxlevel_index"))
  
  for(i in 1:nrow(taxres_df)){ # all rows
    
    if(!is.na(taxres_df[i,"ltg_taxid"]) & taxres_df[i,"ltg_rank_index"] > taxres_df[i,"ltgres"]){ # if resolution of ltg is higher then ltgres
      # get the taxrank that corresponds to ltgres 
      tl <- taxres_df[i,"taxrank"]
      # get the index of the column that corresponds to the ltgres
      col_index <- which(colnames(taxres_df) == tl)
      # make a data frame with taxid, and get taxinfo from tax_df
      new_taxid <- as.data.frame(taxres_df[i, col_index-1]) 
      colnames(new_taxid) <- c("tax_id")
      new_taxid <- left_join(new_taxid, tax_df, by="tax_id") %>%
        select(tax_id, name_txt, rank, taxlevel)
      
      # replace ltg taxid and associated info
      taxres_df[i, 2:5] <- new_taxid[1,]
      # replace tax lineage over the ltgref by NA
      taxres_df[i, (col_index+1):(ncol(taxres_df)-9)] <- NA
    }# end if
  }# end for i
  
  taxres_df <- taxres_df %>%
    select(-taxrank)
  
  return(taxres_df)
}

#' TaxAssign
#' 
#' Find LTG for each asv in the input dataframe
#'  
#' @param taxres_df data frame with the following columns: asv,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,superkingdom_taxid,
#' superkingdom,kingdom_taxid,kingdom,phylum_taxid,phylum,class_taxid,class,order_taxid,order,family_taxid,family,genus_taxid,genus,species_taxid,species,pid,
#' pcov,phit,taxn,seqn,refres,ltgres
#' 
#' 
#' @param df a data frame contining and asv column
#' @param ltg_params_df data frame with a list of ercentage of identity values (pid) and associated parameters (pcov,phit,taxn,seqn,refres,ltgres)
#' @param taxonomy file containing the following columns: tax_id,parent_tax_id,rank,name_txt,old_tax_id(has been mered to tax_id),taxlevel (8: species, 7: genus, 6: family, 5: order, 4: class, 3: phylum, 2: kingdom, 1: superkingdom, 0: root)
#' @param blast_db BLAST database
#' @param blast_path path to BBAST executables
#' @param outdir name of the output directory
#' @export
#'
TaxAssign <- function(df, ltg_params_df="", taxonomy="", blast_db="", blast_path="", outdir="", num_threads=8){
  
  # default value for ltg_params_df
  if(nrow(ltg_params_df)==0){
    ltg_params_df = data.frame( pid=c(100,97,95,90,85,80),
                                pcov=c(70,70,70,70,70,70),
                                phit=c(70,70,70,70,70,70),
                                taxn=c(1,1,2,3,4,4),
                                seqn=c(1,1,2,3,4,4),
                                refres=c(8,8,8,7,6,6),
                                ltgres=c(8,8,8,8,7,7)
    )
  }
  
  #### Read taxonomy info 
  # read taxonomy file; quote="" is important, since some of the taxon names have quotes and this should be ignored
  tax_df <- read.delim(taxonomy, header=T, sep="\t", fill=T, quote="")
  # make data frame with old taxids as line numbers and taxids in a columns
  old_taxid <- tax_df %>%
    filter(!is.na(old_tax_id)) %>%
    select(tax_id, old_tax_id)
  # delete old_tax_ids from tax_df and make taxids unique
  tax_df <- tax_df %>%
    select(-old_tax_id)
  tax_df <- unique(tax_df)
  
  ####
  # create a tmp directory for temporary files using time and a random number
  outdir_tmp <- paste(outdir, 'tmp_', trunc(as.numeric(Sys.time())), sample(1:100, 1), sep='')
  outdir_tmp <- check_dir(outdir_tmp)
  
  ### run blast and clean/complete results
  # run blast and read results to data frame (blast_res columns: "qseqid","pident","qcovhsp","staxids")
  blast_res <- run_blast(df, blast_db=blast_db, blast_path=blast_path, outdir=outdir_tmp, qcov_hsp_perc=min(ltg_params_df$pcov), perc_identity=min(ltg_params_df$pid), num_threads=num_threads)
  # add update old taxids to valid ones
  blast_res <- update_taxids(blast_res, old_taxid)
  # add taxlevel
  blast_res <- left_join(blast_res, tax_df, by=c("staxids" = "tax_id")) %>%
    select(-parent_tax_id, -rank, -name_txt)
  
  ### make a lineage for each taxid in blastres
  lineages <- get_lineage_ids(unique(blast_res$staxids), tax_df)
  # initialize data frame with asv and NA for all other cells
  taxres_df <- data.frame(asv = unique(df$asv), ltg_taxid = NA, pid=NA, pcov=NA, phit=NA, taxn=NA, seqn=NA, refres=NA, ltgres=NA)
  for(i in 1:nrow(taxres_df)){ # go through all sequences 
    for(p in 1:nrow(ltg_params_df)){ # for each pid
      pidl <- ltg_params_df[p,"pid"]
      pcovl <- ltg_params_df[p,"pcov"]
      phitl <- ltg_params_df[p,"phit"]
      taxnl <- ltg_params_df[p,"taxn"]
      seqnl <- ltg_params_df[p,"seqn"]
      refresl <- ltg_params_df[p,"refres"]
      ltgresl <- ltg_params_df[p,"ltgres"]
      
      # filter the blastres according to  pid, pcov, refres
      df_intern <- blast_res %>%
        filter(qseqid==i & pident>=pidl & qcovhsp>=pcovl & taxlevel>=refresl)
      
      # check if enough taxa and seq among validated hits
      tn <- length(unique(df_intern$staxids))
      if(tn >= taxnl & nrow(df_intern) >= seqnl ){
        # make ltg if all conditions are met
        ltg <- make_ltg(df_intern$staxids, lineages, phit = phitl)
        # fill out line with the ltg and the parmeters that were used to get it
        taxres_df[i,2:ncol(taxres_df)] <- c(ltg, pidl, pcovl, phitl, taxnl, seqnl, refresl, ltgresl)
        break
      } # end if
    } # end p (pids)
  } # end i (asvs)
  
  # get the ranked lineage for each taxid in taxres_df
  ranked_lineages <- get_ranked_lineages(unique(taxres_df$ltg_taxid), tax_df)
  # add lineage to taxres_df
  taxres_df <- left_join(taxres_df, ranked_lineages, by="ltg_taxid") %>%
    select(asv,ltg_taxid,ltg_name,ltg_rank,ltg_rank_index,superkingdom_taxid,superkingdom,kingdom_taxid,kingdom,phylum_taxid,phylum,class_taxid,class,order_taxid,order,family_taxid,family,genus_taxid,genus,species_taxid,species,pid,pcov,phit,taxn,seqn,refres,ltgres)
  # adjust resolution if it is higher than ltgres
  taxres_df <- adjust_ltgres(taxres_df, tax_df)
  
  # delete temporary  dir
  unlink(outdir_tmp, recursive = TRUE)
  
  return(taxres_df)
}
