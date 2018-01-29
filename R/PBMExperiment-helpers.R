#' @author Patrick Kimes
makePBMExperiment <- function(tab) {
    ## currently only support all scans with same design
    if (length(unique(tab$vers)) > 1) {
        stop("All samples/scans must have the same assay design version")
    }

    ## define arrays to be unique combination of date/id, scan-type, laser power
    tab <- tidyr::unite(tab, "assay", c("date", "id", "scan", "lp"),
                        remove = FALSE)

    ## read in all assays
    assay_metadata <- dplyr::select(tab, assay, date, id, scan, lp)
    assay_metadata <- dplyr::distinct(assay_metadata)
    assay_tables <- lapply(assay_metadata$assay, readAssay, stab = tab)
    assay_tables <- purrr::reduce(assay_tables, left_join,
                                  by = c("Column", "Row", "Name", "ID"))

    ## make list of asssays
    assay_list <- assay_tables %>%
        mutate_if(is.list, funs(ifelse(sapply(., length) > 0, ., NA_real_))) %>%
        gather(array, vals, -Column, -Row, -Name, -ID) %>%
        select(., -Column, -Row, -Name, -ID) %>%
        group_by(array) %>%
        do(l = DataFrame(unnest(select(., -array))))

    ## add names to assays
    assay_names <- assay_list$array
    assay_list <- assay_list$l
    names(assay_list) <- assay_names

    ## row/probe-level metadata
    rowdat <- DataFrame(select(assay_tables, Column, Row, Name, ID))
    
    ## column/condition-level metadata
    coldat <- DataFrame(select(assay_metadata, -assay))
    rownames(coldat) <- assay_metadata$assay

    ## SummarizedExperiment
    SummarizedExperiment(assays = assay_list, rowData = rowdat,
                         metadata = list(assayinfo = coldat))
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
    gpr_list <- mapply(function(g, n) { names(g)[5] <- n; g },
                       gpr_list, stab$condition, SIMPLIFY = FALSE)
    gprs <- purrr::reduce(gpr_list, left_join,
                          by = c("Column", "Row", "Name", "ID"))
    tidyr::nest(gprs, -Column, -Row, -Name, -ID, .key = UQ(x))
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
    names(vals)[5] <- 'intensity'
    vals %>% select(Column, Row, Name, ID, intensity)
}

