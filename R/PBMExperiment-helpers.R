#' Read PBMExperiment Data
#'
#' @import SummarizedExperiment
#' @importFrom purrr reduce
#' @importFrom dplyr select
#' @export
#' @author Patrick Kimes
makePBMExperiment <- function(tab) {
    ## currently only support all scans with same design
    if (length(unique(tab$vers)) > 1) {
        stop("All samples/scans must have the same assay design version")
    }

    ## read in all gpr scans
    assay_table <- lapply(tab$gpr, readGPR)
    assay_table <- purrr::reduce(assay_table, left_join,
                                 by = c("Column", "Row", "Name", "ID"))

    ## assay data
    assaydat <- DataFrame(dplyr::select(assay_table, -Column, -Row, -Name, -ID))
    names(assaydat) <- paste0("s", 1:ncol(assaydat))

    ## row/probe-level metadata
    rowdat <- DataFrame(dplyr::select(assay_table, Column, Row, Name, ID))
    
    ## column/condition-level metadata
    coldat <- DataFrame(dplyr::select(tab, -gpr))
    rownames(coldat) <- names(assaydat)

    ## SummarizedExperiment
    SummarizedExperiment(assays = list(gpr = assaydat),
                         rowData = rowdat, colData = coldat)
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
#' @importFrom readr read_tsv
#' @importFrom dplyr select
#' @author Patrick Kimes
readGPR <- function(x) {
    colt <- rep("-", 45)
    colt[c(2, 3)] <- 'i'
    colt[c(4, 5, 38, 40)] <- 'c'
    colt[c(34)] <- 'd'
    vals <- readr::read_tsv(x, skip = 35, col_types = paste(colt, collapse = ""))
    names(vals)[5] <- 'intensity'
    vals %>% dplyr::select(Column, Row, Name, ID, intensity)
}

