#' @title Tidy SummarizedExperiment object
#'
#' @description
#' Simple helper function to convert a single assay data in
#' a SummarizedExperiment object into a tidy tibble. If \code{long = FALSE},
#' a tibble with the same number of rows as the input SummarizedExperiment
#' is constructed with columns corresponding to the columns of the
#' SummarizedExperiment and any rowData. If \code{long = TRUE}, the assay
#' data is "tidied" such that each row of the tibble corresponds to a single
#' value of the assay. When \code{long = TRUE}, in addition to rowData, any
#' colData is also added as additional columns. 
#' 
#' @param x a SummarizedExperiment object.
#' @param assay a numeric index or string name specifying the assay to tidy.
#'        If multiple assays are specified, assays will be combined
#'        as separate columns.
#'        (default = \code{SummarizedExperiment::assayNames(x)[1]})
#' @param long a logical whether to transform data to long format and
#'        include colData in output rather than default wide format with
#'        dimension similar to original PBMExperiment object.
#'        Ignored and set to TRUE if \code{assay} specifies more than a
#'        single assay. (default = FALSE)
#' @param ... other parameters for the \code{tidy} generic function. 
#'
#' @return
#' tibble containing values from a single SummarizedExperiment assay
#' along with rowData and optionally colData.
#'
#' @name tidy-SummarizedExperiment
#' @importFrom dplyr as_tibble bind_cols left_join n 
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom broom tidy
#' @export 
#' @author Patrick Kimes
tidy.SummarizedExperiment <- function(x, assay = SummarizedExperiment::assayNames(x)[1],
                                      long = FALSE, ...) {

    assay <- match.arg(assay, SummarizedExperiment::assayNames(x),
                       several.ok = TRUE)
    
    ## extract row data
    rowdat <- as.data.frame(rowData(x), optional = TRUE)
    rowdat <- dplyr::as_tibble(rowdat)
    
    ## extract intensities
    pdat <- lapply(assay, SummarizedExperiment::assay, x = x)
    names(pdat) <- assay
    pdat <- lapply(pdat, as.data.frame, optional = TRUE)
    pdat <- lapply(pdat, dplyr::as_tibble)
    pdat[[1]] <- dplyr::bind_cols(pdat[[1]], rowdat)
    
    if (long || length(assay) > 1) {
        ## combine assays
        pdat <- lapply(pdat, dplyr::mutate, `__rowid` = 1:dplyr::n())
        pdat <- lapply(pdat, tidyr::pivot_longer, names_to = "cname",
                       values_to = "value", cols = colnames(x))
        pdat <- mapply(function(tib, nam) {
            dplyr::rename(tib, !! nam := value)
        }, tib = pdat, nam = names(pdat), SIMPLIFY = FALSE) 
        pdat <- purrr::reduce(pdat, dplyr::left_join,
                              by = c("cname", "__rowid"))
        pdat <- dplyr::select(pdat, -`__rowid`)

        ## extract column data
        coldat <- as.data.frame(colData(x), optional = TRUE)
        coldat <- tibble::rownames_to_column(coldat, "cname")
        pdat <- dplyr::left_join(pdat, coldat, by = "cname")
    } else {
        pdat <- pdat[[1]]
    }
    
    return(pdat)
}


#' @title Tidy PBMExperiment object
#'
#' @description
#' Simple helper function to convert assay format data in PBMExperiment
#' object into a tidy tibble. If 
#' If \code{process = TRUE}, probe filtering and trimming are performed
#' according to the \code{probeFilter} and \code{probeTrim} slots of the PBMExperiment
#' object before passing the object to the SummarizedExperiment
#' \code{\link[=tidy-SummarizedExperiment]{tidy}} method.
#' If \code{process = FALSE}, the PBMExperiment is passed directly to the
#' \code{\link[=tidy-SummarizedExperiment]{tidy}} method without any
#' processing.
#'
#' @param x a PBMExperiment object.
#' @param assay a numeric index or string name specifying the assay to tidy.
#'        (default = \code{SummarizedExperiment::assayNames(x)[1]})
#' @param long a logical value whether to transform data to long format and
#'        include colData in output rather than default wide format with
#'        dimension similar to original PBMExperiment object.
#'        (default = FALSE)
#' @param process a logical value whether to filter probes and trim probe sequences
#'        according to PBMExperiment rules in \code{probeFilter} and \code{probeTrim}
#'        slots. (default = TRUE)
#' @param ... other parameters for the \code{tidy} generic function. 
#'
#' @return
#' tibble containing values from a single PBMExperiment assay
#' along with rowData and optionally colData.
#'
#' @seealso \code{\link{tidy-SummarizedExperiment}}
#' @name tidy-PBMExperiment
#' @importFrom dplyr as_tibble bind_cols left_join
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom broom tidy
#' @export 
#' @author Patrick Kimes
tidy.PBMExperiment <- function(x, assay = SummarizedExperiment::assayNames(x)[1],
                               long = FALSE, process = TRUE, ...) {
    ## filter probes
    if (process) {
        x <- pbmFilterProbes(x)
        x <- pbmTrimProbes(x)
    }

    ## treat as regular SummarizedExperiment
    tidy.SummarizedExperiment(x, assay, long, ...)
}

