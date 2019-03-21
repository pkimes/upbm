.show.PBMExperiment <- function(object) {
    methods::callNextMethod(object)
    
    cat(paste0(stringr::str_trunc(paste0("probeCols(", length(object@probeCols), "): ",
                                         paste0(object@probeCols, collapse = " ")),
                                  width = 60),
               "\n"))
    cat(paste0(stringr::str_trunc(paste0("probeFilter names(", length(object@probeFilter), "): ",
                                         paste0(names(object@probeFilter), collapse = " ")),
                                  width = 60),
               "\n"))
    cat(paste0("probeTrim: ", paste0(object@probeTrim, collapse = " "), "\n"))
}

#' @importFrom methods callNextMethod
#' @rdname show-methods
#' @export
setMethod("show", signature(object = "PBMExperiment"), .show.PBMExperiment)
