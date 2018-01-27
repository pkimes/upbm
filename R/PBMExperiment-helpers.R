#' @author Patrick Kimes
makePBMExperiment <- function(tab) {
    ## currently only support all scans with same design
    if (length(unique(tab$vers)) > 1) {
        stop("All samples/scans must have the same assay design version")
    }

    ## define arrays to be unique combination of date/id, scan-type, laser power
    tab <- tidyr::unite(tab, "assay", c("date", "id", "scan", "lp"),
                        remove = FALSE)

    ## assays, colData, rowData #################################################
    ## start by reading in each assay
    assay_metadata <- dplyr::select(tab, assay, date, id, scan, lp)
    assay_metadata <- dplyr::distinct(assay_metadata)
    assay_list <- lapply(assay_metadata$assay, readAssay, stab = tab)

    ## merge SummarizedExperiments?
    ## - or I should have been using data.frames up until now...

    ## assayData ################################################################
    ## then, need to add assay data
}


#' Read GPR Files as Assay
#'
#' Helper function for reading in each assay, i.e. collection of
#' GPR files with same laser intensity and scan (Alex488 or Cy3)
#' across different conditions.
#'
#' @param x name of assay
#' @param stab table with `assay`, and `gpr` columns providing the
#' paths to all sample GPR files in the assay.
#'
#' @return
#' SummarizedExperiment object with one assay with multiple columns
#' corresponding to different conditions.
#'
#' @author Patrick Kimes
readAssay <- function(x, stab) {
    stab <- filter(stab, assay == x)
    gpr_list <- lapply(stab$gpr, readGPR)
    gpr_list <- mapply(function(g, n) { colnames(g) <- n; g },
                       gpr_list, stab$condition)
    gprs <- do.call(cbind, gpr_list)
    names(assays(gprs)) <- x
    gprs
}


#' Read GPR File as Assay
#' 
#' Helper function for reading in each GPR file.
#'
#' @param x path to GPR fle 
#'
#' @return
#' SummarizedExperiment object with one assay with one column,
#' and corresponding rowData
#' 
#' @details
#' Since the name of the value column can vary across samples,
#' we pull it out using its index (5th column)
#' 
#' @author Patrick Kimes
readGPR <- function(x) {
    colt <- rep("-", 45)
    colt[c(2, 3)] <- 'i'
    colt[c(4, 5, 38, 40)] <- 'c'
    colt[c(34)] <- 'd'
    vals <- read_tsv(x, skip = 35, col_types = paste(colt, collapse = ""))
    SummarizedExperiment(assays = list(a = DataFrame(val = dplyr::pull(vals, 5))),
                         rowData = DataFrame(select(vals, Column, Row, Name, ID)))
}

