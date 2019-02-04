#' @title Create a new PBMExperiment object
#'
#' @description
#' Create a PBMExperiment object by directly specifying optional probe annotation
#' information (\code{probeFilter}, \code{probeTrim}, \code{probeCols}), along with
#' a based SummarizedExperiment object containing probe intensity data. 
#' 
#' @param probeFilter an optional named list of probe filters to be used to subset
#'        probes during data analysis steps. List names must correspond to columns in
#'        rowData. List entries must be single-parameter functions to be called on the
#'        corresponding column to return a logical vector of probes to keep (TRUE) or
#'        drop (FALSE) during analysis. (default = \code{list()})
#' @param probeTrim an optional integer vector of length 2 specifying start and end
#'        positions in probe `Sequence' to use in analysis steps. (default = \code{numeric()})
#' @param probeCols an optional character vector of rowData column names corresponding
#'        to probe design information. (default = \code{character()})
#' @param ... an optional SummarizedExperiment containing PBM probe intensity data and
#'        corresponding annotations or parameters to be passed to the SummarizedExperiment
#'        constructor.
#'
#' @return
#' \code{\link[=PBMExperiment-class]{PBMExperiment}} object.
#' 
#' @seealso \code{\link{PBMExperiment-class}}, \code{\link[SummarizedExperiment:RangedSummarizedExperiment-class]{SummarizedExperiment-class}}, \code{\link{gpr2PBMExperiment}}
#' @export
#' @author Patrick Kimes
PBMExperiment <- function(probeFilter = list(),
                          probeTrim = numeric(),
                          probeCols = character(),
                          ...) {
    se <- SummarizedExperiment(...)
    .PBMExperiment(se, probeFilter = probeFilter,
                   probeTrim = probeTrim,
                   probeCols = probeCols)
}
