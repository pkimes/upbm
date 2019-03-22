#' @title Create a new PBMExperiment object
#'
#' @description
#' Create a PBMExperiment object by directly specifying optional probe annotation
#' information (\code{probeFilter}, \code{probeTrim}, \code{probeCols}), along with
#' a based SummarizedExperiment object containing probe intensity data. 
#' 
#' @param ... an optional SummarizedExperiment containing PBM probe intensity data and
#'        corresponding annotations or parameters to be passed to the SummarizedExperiment
#'        constructor.
#' @param pbmDesign an optional PBMDesign containing design information. If specified,
#'        \code{probeFilter}, \code{probeTrim}, and \code{probeCols} parameters will be 
#'        ignored. Ignored if NULL. (default = NULL)
#' @param probeFilter an optional named list of probe filters to be used to subset
#'        probes during data analysis steps. List names must correspond to columns in
#'        rowData. List entries must be single-parameter functions to be called on the
#'        corresponding column to return a logical vector of probes to keep (TRUE) or
#'        drop (FALSE) during analysis. (default = \code{list()})
#' @param probeTrim an optional integer vector of length 2 specifying start and end
#'        positions in probe `Sequence' to use in analysis steps. (default = \code{numeric()})
#' @param probeCols an optional character vector of rowData column names corresponding
#'        to probe design information. (default = \code{c("Sequence", "probeID")})
#'
#' @return
#' \code{\link[=PBMExperiment-class]{PBMExperiment}} object.
#' 
#' @seealso \code{\link{PBMExperiment-class}}, \code{\link[SummarizedExperiment:RangedSummarizedExperiment-class]{SummarizedExperiment-class}}, \code{\link{gpr2PBMExperiment}}
#' @export
#' @author Patrick Kimes
PBMExperiment <- function(...,
                          pbmDesign = NULL,
                          probeFilter = list(),
                          probeTrim = numeric(),
                          probeCols = c("Sequence", "probeID")) {

    args <- list(...)
    if (length(args) == 1L && is(args[[1L]], "SummarizedExperiment")) {
        se <- args[[1L]]
    } else {
        se <- do.call(SummarizedExperiment, args)
    }
    rd <- rowData(se)
    if (is.null(pbmDesign)) {
        if (! "Sequence" %in% colnames(rd)) {
            rowData(se)$Sequence <-rep(NA_character_, nrow(se))
        }
        if (! "probeID" %in% colnames(rd)) {
            rowData(se)$probeID <- rep(NA_character_, nrow(se))
        }
    } else {
        if (!(is(pbmDesign, "PBMDesign"))) {
            stop("If specified, 'pbmDesign' must be a PBMDesign object.")
        }
        if (nrow(design(pbmDesign)) != nrow(se)) {
            stop("Dimensions of SummarizedExperiment and PBMDesign do not match.")
        }
        if ("Sequence" %in% colnames(rd)) {
            warning("Existing 'Sequence' column of SummarizedExperiment will be overwritten.")
            rowData(se)$Sequence <- NULL
        }
        if ("probeID" %in% colnames(rd)) {
            warning("Existing 'probeID' column of SummarizedExperiment will be overwritten.")
            rowData(se)$probeID <- NULL
        }
        rowData(se) <- rowData(se)[, setdiff(colnames(rowData(se)), colnames(design(pbmDesign)))]
        rowData(se) <- cbind(rowData(se), design(pbmDesign))
        probeCols <- probeCols(pbmDesign)
        probeTrim <- probeTrim(pbmDesign)
        probeFilter <- probeFilter(pbmDesign)
    }
    
    .PBMExperiment(se,
                   probeFilter = probeFilter,
                   probeTrim = probeTrim,
                   probeCols = probeCols)
}
