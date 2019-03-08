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
#'        (default = \code{SummarizedExperiment::assayNames(x)[1]})
#' @param long a logical whether to transform data to long format and
#'        include colData in output rather than default wide format with
#'        dimension similar to original PBMExperiment object.
#'        (default = FALSE)
#' @param ... other parameters for the \code{tidy} generic function. 
#'
#' @return
#' tibble containing values from a single SummarizedExperiment assay
#' along with rowData and optionally colData.
#'
#' @name tidy-SummarizedExperiment
#' @importFrom dplyr as_tibble bind_cols left_join
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom broom tidy
#' @export 
#' @author Patrick Kimes
tidy.SummarizedExperiment <- function(x, assay = SummarizedExperiment::assayNames(x)[1],
                                      long = FALSE, ...) {
    ## extract row data
    rowdat <- as.data.frame(rowData(x), optional = TRUE)
    rowdat <- dplyr::as_tibble(rowdat)
    
    ## extract intensities
    pdat <- SummarizedExperiment::assay(x, assay)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- dplyr::as_tibble(pdat)
    pdat <- dplyr::bind_cols(pdat, rowdat)

    if (long) {
        ## extract column data
        coldat <- as.data.frame(colData(x), optional = TRUE)
        coldat <- tibble::rownames_to_column(coldat, "cname")

        pdat <- tidyr::gather(pdat, cname, value, colnames(x))
        pdat <- dplyr::left_join(pdat, coldat, by = "cname")
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
#' @importFrom tidyr gather
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
