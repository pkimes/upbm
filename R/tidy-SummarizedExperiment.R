#' Tidy SummarizedExperiment Object
#'
#' @description
#' Simple helper function to convert assay format data in SummarizedExperiment
#' object into a tidy tibble for interactive analysis.
#' 
#' @param x SummarizedExperiment object.
#' @param assay string name of the assay to use.
#'        (default = \code{assayNames(x)[1]})
#' @param long logical whether to transform data to long format and
#'        include colData in output rather than default wide format with
#'        dimension similar to original SummarizedExperiment object.
#'        (default = FALSE)
#' @param .filter integer specifying level of probe filtering. (default = 1L)
#'
#' @return
#' tibble with assay data.
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
                                      long = FALSE, .filter = 1L, ...) {
    ## filter probes
    x <- pbmFilterProbes(x, .filter) 

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
