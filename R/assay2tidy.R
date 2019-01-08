#' Tidy SummarizedExperiment Object
#'
#' Simple helper function to convert assay format data in SummarizedExperiment
#' object into a tidy tibble for interactive analysis.
#' 
#' @param se SummarizedExperiment object.
#' @param assay string name of the assay to use.
#'        (default = \code{assayNames(se)[1]})
#' @param long logical whether to transform data to long format and
#'        include colData in output rather than default wide format with
#'        dimension similar to original SummarizedExperiment object.
#'        (default = FALSE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to plotting. (default = 1)
#'
#' @return
#' tibble with assay data.
#'
#' @importFrom dplyr as_tibble bind_cols left_join
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom SummarizedExperiment assayNames assay
#' @export
#' @author Patrick Kimes
assay2tidy <- function(se, assay = SummarizedExperiment::assayNames(se)[1],
                       long = FALSE, .filter = 1L) {

    ## filter probes
    se <- pbmFilterProbes(se, .filter) 

    ## extract row data
    rowdat <- as.data.frame(rowData(se), optional = TRUE)
    rowdat <- dplyr::as_tibble(rowdat)
    
    ## extract intensities
    pdat <- SummarizedExperiment::assay(se, assay)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- dplyr::as_tibble(pdat)
    pdat <- dplyr::bind_cols(pdat, rowdat)

    if (long) {
        ## extract column data
        coldat <- as.data.frame(colData(se), optional = TRUE)
        coldat <- tibble::rownames_to_column(coldat, "cname")

        pdat <- tidyr::gather(pdat, cname, value, colnames(se))
        pdat <- dplyr::left_join(pdat, coldat, by = "cname")
    }
    
    return(pdat)
}
