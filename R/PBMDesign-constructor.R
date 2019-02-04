#' Create a new PBMDesign object
#'
#' @description
#' Create a new PBMDesign object of protein binding microarray probe design information.
#' Alternatively, the function can be called on a PBMExperiment to extract the
#' probe design information associated with experimental data.
#'
#' @param x data.frame with each row corresponding to a probe on the array.
#'        Must include `Sequence' and (unique) `probeID' columns, along with any
#'        other metadata for probes, e.g. array `Row' or `Column' spatial coordinates.
#'        Alternatively, a \code{\link[=PBMExperiment-class]{PBMExperiment}} object to
#'        return the associated \code{PBMDesign} object.
#' @param ... optional probe design parameters to be defined as part of the \code{PBMDesign}
#'        object. See the \code{\link[=PBMDesign-class]{PBMDesign}} class definition for
#'        a list of probe design parameters. Important parameters are as described in
#'        the Details section below.
#'
#' @details
#' Probe design parameters can be specified by name. The following are a couple important
#' parameters which are defined by default for universal PBM designs in the \pkg{upbmAux}
#' package.
#' \enumerate{
#' \item \code{probeFilter}: optional named list of probe filters to be used to subset
#' probes during data analysis steps. List names must correspond to columns in `design'.
#' List entries must be single-parameter functions to be called on the corresponding column
#' to return a logical vector of probes to keep (TRUE) and drop (FALSE) during analysis.
#' \item \code{probeTrim}: optional integer vector of length 2 specifying start and end
#' positions in probe `Sequence' to use in analysis steps.
#' }
#'
#' @return
#' \code{PBMDesign} object.
#' 
#' @seealso \code{\link{PBMDesign-class}}
#' @name PBMDesign
#' @export
#' @author Patrick Kimes
NULL

.PBMDesign.table <- function(x, ...) {
    new("PBMDesign", design = x, ...)
}

.PBMDesign.PBMExperiment <- function(x) {
    new("PBMDesign", design = rowData(x)[, x@probeCols],
        probeFilter = x@probeFilter, probeTrim = x@probeTrim)
}

#' @rdname PBMDesign
#' @exportMethod "PBMDesign"
setMethod("PBMDesign", signature(x = "data.frame"), .PBMDesign.table)
setMethod("PBMDesign", signature(x = "DataFrame"), .PBMDesign.table)

#' @rdname PBMDesign
#' @exportMethod "PBMDesign"
setMethod("PBMDesign", signature(x = "PBMExperiment"), .PBMDesign.PBMExperiment)
