

library("dada2")

# Use getAnywhere() to access the code for isBimeraDenovo()
code <- getAnywhere(isBimeraDenovo)

# Print the code
print(code)

?isBimeraDenovo

###########################################
function (unqs, minFoldParentOverAbundance = 2, minParentAbundance = 8, 
          allowOneOff = FALSE, minOneOffParentDistance = 4, maxShift = 16, 
          multithread = FALSE, verbose = FALSE) 
{
  if (any(duplicated(getSequences(unqs)))) 
    message("Duplicate sequences detected.")
  unqs.int <- getUniques(unqs, silence = TRUE)
  abunds <- unname(unqs.int)
  seqs <- names(unqs.int)
  seqs.input <- getSequences(unqs)
  rm(unqs)
  gc(verbose = FALSE)
  if (is.logical(multithread)) {
    if (multithread == TRUE) {
      mc.cores <- getOption("mc.cores", detectCores())
    }
  }
  else if (is.numeric(multithread)) {
    mc.cores <- multithread
    multithread <- TRUE
  }
  else {
    warning("Invalid multithread parameter. Running as a single thread.")
    multithread <- FALSE
  }
  loopFun <- function(i, unqs.loop, minFoldParentOverAbundance, 
                      minParentAbundance, allowOneOff, minOneOffParentDistance, 
                      maxShift) {
    sq <- names(unqs.loop)[[i]]
    abund <- unqs.loop[[i]]
    pars <- names(unqs.loop)[(unqs.loop > (minFoldParentOverAbundance * 
                                             abund) & unqs.loop > minParentAbundance)]
    if (length(pars) < 2) {
      return(FALSE)
    }
    else {
      isBimera(sq, pars, allowOneOff = allowOneOff, minOneOffParentDistance = minOneOffParentDistance, 
               maxShift = maxShift)
    }
  }
  if (multithread) {
    mc.indices <- sample(seq_along(unqs.int), length(unqs.int))
    bims <- mclapply(mc.indices, loopFun, unqs.loop = unqs.int, 
                     allowOneOff = allowOneOff, minFoldParentOverAbundance = minFoldParentOverAbundance, 
                     minParentAbundance = minParentAbundance, minOneOffParentDistance = minOneOffParentDistance, 
                     maxShift = maxShift, mc.cores = mc.cores)
    bims <- bims[order(mc.indices)]
  }
  else {
    bims <- lapply(seq_along(unqs.int), loopFun, unqs.loop = unqs.int, 
                   allowOneOff = allowOneOff, minFoldParentOverAbundance = minFoldParentOverAbundance, 
                   minParentAbundance = minParentAbundance, minOneOffParentDistance = minOneOffParentDistance, 
                   maxShift = maxShift)
  }
  bims <- unlist(bims)
  bims.out <- seqs.input %in% seqs[bims]
  names(bims.out) <- seqs.input
  if (verbose) 
    message("Identified ", sum(bims.out), " bimeras out of ", 
            length(bims.out), " input sequences.")
  return(bims.out)
}

###########################################

# Use getAnywhere() to access the code for isBimeraDenovo()
code <- getAnywhere(isBimera)
# Print the code
print(code)
?isBimera

###########################################
function (sq, parents, allowOneOff = FALSE, minOneOffParentDistance = 4, 
          maxShift = 16) 
{
  rval <- C_is_bimera(sq, parents, allowOneOff, minOneOffParentDistance, 
                      getDadaOpt("MATCH"), getDadaOpt("MISMATCH"), getDadaOpt("GAP_PENALTY"), 
                      maxShift)
  return(rval)
}
###########################################



# Use getAnywhere() to access the code for isBimeraDenovo()
code <- getAnywhere(C_is_bimera)
# Print the code
print(code)
###########################################
function (sq, pars, allow_one_off, min_one_off_par_dist, match, 
          mismatch, gap_p, max_shift) 
{
  .Call("_dada2_C_is_bimera", PACKAGE = "dada2", sq, pars, 
        allow_one_off, min_one_off_par_dist, match, mismatch, 
        gap_p, max_shift)
}
###########################################

