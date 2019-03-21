.show.PBMDesign <- function(object) {
    cat("class: PBMDesign\n")
    cat(paste0("design dim: ", nrow(object@design), " ", ncol(object@design), "\n"))
    cat(paste0(stringr::str_trunc(paste0("design colnames(", ncol(object@design), "): ",
                                         paste0(colnames(object@design), collapse = " ")),
                                  width = 60),
               "\n"))
    cat(paste0(stringr::str_trunc(paste0("probeFilter names(", length(object@probeFilter), "): ",
                                         paste0(names(object@probeFilter), collapse = " ")),
                                  width = 60),
               "\n"))
    cat(paste0("probeTrim: ", paste0(object@probeTrim, collapse = " "), "\n"))
}

#' @title Show upbm objects
#' 
#' @param object a PBMDesign or PBMExperiment object to show.
#' 
#' @importFrom stringr str_trunc
#' @name show-methods
#' @aliases show,PBMDesign-method
#' @export
setMethod("show", signature(object = "PBMDesign"), .show.PBMDesign)
