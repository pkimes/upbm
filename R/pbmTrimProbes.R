#' @title Trim probe sequences 
#'
#' @description
#' Not all nucleotides in the probe sequences should be used for analysis.
#' Given a PBMExperiment or a PBMDesign object, this function returns a
#' modified object with probe sequences trimmed according to the
#' associated \code{probeTrim} slot value.
#' 
#' @param pe a PBMExperiment or PBMDesign object.
#'
#' @return
#' Original PBMExperiment or PBMDesign object with trimmed sequences.
#'
#' @importFrom stringr str_sub
#' @export
#' @author Patrick Kimes
pbmTrimProbes <- function(pe) {
    stopifnot(is(pe, "PBMExperiment") || is(pe, "PBMDesign"))

    if (length(pe@probeTrim) == 0L) {
        return(pe)
    }

    if (is(pe, "PBMExperiment")) {
        rowData(pe)$Sequence <- stringr::str_sub(rowData(pe)$Sequence,
                                                 pe@probeTrim[1], pe@probeTrim[2])
    } else {
        pe@design$Sequence <- stringr::str_sub(pe@design$Sequence,
                                               pe@probeTrim[1], pe@probeTrim[2])
    }

    return(pe)
}
