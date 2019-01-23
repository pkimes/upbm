#' PBMExperiment class
#'
#' @description
#' Extension of the \code{\link{SummarizedExperiment}} class to store
#' universal protein binding microarray (uPBM) data for use with functions
#' in the \pkg{upbm} package. 
#'
#' @aliases PBMExperiment-class
#' @import SummarizedExperiment
#' @export
#' @exportClass PBMExperiment
#' @name PBMExperiment-class
#' @author Patrick Kimes
setClass("PBMExperiment", contains = "SummarizedExperiment",
         slots = c(probeFilter = "list", probeTrim = "numeric"))

setValidity("PBMExperiment",
            function(object) {
                rd <- rowData(object)
                if (!all(c("Sequence", "ID") %in% colnames(rd))) {
                    stop("PBMExperiment must include probe sequence information ",
                         "as rowData columns: 'Sequence' and 'ID'.")
                }
                if (!is.vector(rd$ID, mode = "character")) {
                    stop("PBMExperiment probe IDs must be character strings.\n",
                         "Please check 'ID' rowData values.")
                }
                if (!is.vector(rd$Sequence, mode = "character")) {
                    stop("PBMExperiment probe sequences must be character strings.\n",
                         "Please check 'Sequence' rowData values.")
                }
                ## if (any(nchar(rd$Sequence) != nchar(rd$Sequence[1]))) {
                ##     stop("PBMExperiment probe sequences must all be same length.\n",
                ##          "Please check 'Sequence' rowData values.")
                ## }
                if (any(duplicated(rd$ID))) {
                    stop("PBMExperiment probe IDs must be unique.\n",
                         "Please check 'ID' rowData values.")
                }
                ## check probeFilter specification
                if (length(object@probeFilter) > 1) {
                    if (is.null(names(object@probeFilter)) | any(names(object@probeFilter) == "")) {
                        stop("PBMExperiment 'probeFilter' must be a named list.")
                    }
                    if (!all(names(object@probeFilter) %in% colnames(rd))) {
                        stop("PBMExperiment 'probeFilter' must be a named list matching ",
                             "rowData columns.")
                    }
                }
                ## check probeTrim specification
                if (length(probeTrim) != 0L || length(probeTrim) != 2L) {
                    stop("PBMExperiment 'probeTrim' is too long. \n",
                         "If specified, 'probeTrim' must be a vector of length 2 corresponding ",
                         "to the start and end of the probe sequence to keep for analysis.")
                }
                TRUE
            })


