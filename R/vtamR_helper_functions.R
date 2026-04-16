

#' @importFrom dplyr filter mutate group_by select summarize summarise arrange 
#' @importFrom dplyr desc left_join full_join inner_join %>% n_distinct distinct 
#' @importFrom dplyr bind_rows ungroup rename rename_with rowwise n do first if_else
#' @importFrom ggplot2 ggplot geom_bar labs theme element_text scale_y_continuous 
#' @importFrom ggplot2 aes geom_density theme_minimal geom_histogram after_stat
#' @importFrom utils read.csv write.table read.table read.delim count.fields
#' @importFrom tidyr everything pivot_wider gather separate 
#' @importFrom tidyselect where
#' @importFrom rlang sym :=
#' @importFrom magrittr %>%
#' @importFrom seqinr splitseq
NULL

#' Download from Zenodo
#'
#' Downloads a file from Zenodo and optionally extracts it if it is a
#' `.tar.gz` archive.
#'
#' @param filename Character string specifying the name of the file to download.
#' @param url Character string; URL of the Zenodo record (copy from browser address bar).
#' @param dest_dir Character string; path to the directory where the file should be downloaded.
#' @param untar Logical; if `TRUE`, the archive is extracted after download.
#' @param quiet Logical; if `TRUE`, suppress informational messages and show only warnings or errors.
#'
#' @return Invisible full path to the directory where the file was downloaded.
#'
#' @examples
#' \dontrun{
#' download_zenodo(
#'   filename = "COInr_2025_05_23.tar.gz",
#'   url = "https://zenodo.org/records/15515860",
#'   dest_dir = "~/vtamR",
#'   untar = TRUE,
#'   quiet = FALSE
#' )
#' }
#'
#' @export
#' 
download_zenodo <- function(
    filename,
    url,
    dest_dir = ".",
    untar = TRUE,
    quiet = FALSE
) {
  
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest_dir <- normalizePath(dest_dir, mustWork = FALSE)
  
  url <- paste0(url, "/files/", filename, "?download=1")
  destfile <- file.path(dest_dir, filename)
  
  if (!quiet) message("Downloading ", filename, " to ", dest_dir)
  utils::download.file(url, destfile)
  
  if (!quiet) message("Extracting archive...")
  utils::untar(destfile, exdir = dest_dir)
  
  unlink(destfile)
  
  if (!quiet) message("Extraction complete: ", dest_dir)
  return(invisible(dest_dir))
}

#' Download from OSF
#'
#' Downloads a file from the Open Science Framework (OSF).
#' If the file is a `.tar.gz` archive, it can optionally be extracted after download.
#'
#' @param filename Character string specifying the name of the file to download.
#' @param url Character string; direct download URL of the OSF file.
#' @param dest_dir Character string; path to the directory where the file should be downloaded.
#' @param untar Logical; if `TRUE`, the archive is extracted after download.
#' @param quiet Logical; if `TRUE`, suppress informational messages and show only warnings or errors.
#'
#' @return Invisible full path to the directory where the file was downloaded.
#'
#' @examples
#' \dontrun{
#' download_osf(
#'   filename = "COInr_for_vtam_2025_05_23_dbV5.tar.gz",
#'   url = "https://osf.io/download/jyhz6/",
#'   dest_dir = "~/vtamR/OSF",
#'   untar = TRUE,
#'   quiet = FALSE
#' )
#' }
#'
#' @export
#' 
download_osf <- function(
    filename,
    url,
    dest_dir = ".",
    untar = TRUE,
    quiet = FALSE
) {
  
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest_dir <- normalizePath(dest_dir, mustWork = FALSE)
  
  destfile <- file.path(dest_dir, filename)
  
  if (!quiet) message("Downloading ", filename, " to ", dest_dir)
  utils::download.file(url, destfile)
  
  if (!quiet) message("Extracting archive...")
  utils::untar(destfile, exdir = dest_dir)
  
  unlink(destfile)
  
  if (!quiet) message("Extraction complete: ", dest_dir)
  return(invisible(dest_dir))
}