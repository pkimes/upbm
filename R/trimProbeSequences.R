#' @title Trim probe sequences 
#'
#' @description
#' Simple helper function to trim probe sequences in the rowData of
#' a PBMExperiment or PBMDesign object.
#' 
#' @param pe a PBMExperiment object or PBMDesign object.
#'
#' @return
#' Original PBMExperiment or PBMDesign object with trimmed sequences.
#'
#' @importFrom stringr str_sub
#' @export
#' @author Patrick Kimes
trimProbeSequences <- function(pe) {
    stopifnot(is(pe, "PBMExperiment") || is(pe, "PBMDesign"))
    
    if (length(pe@probeTrim) == 0L) {
        return(pe)
    }

    if (is(pe, "PBMExperiment")) {
        pd <- rowData(pe)
    } else {
        pd <- pe@design
    }
    pd$Sequence <- stringr::str_sub(pd$Sequence, pe@probeTrim[1], pe@probeTrim[2])

    return(pe)
}
