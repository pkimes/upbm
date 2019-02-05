#' Filter PBMExperiment and PBMDesign objects
#'
#' Not all probes in a PBM experiment should be used for analysis.
#' Given a PBMExperiment or a PBMDesign object with defined probe
#' filters, this function returns the subset of probes passing all
#' filtering criteria.
#' 
#' @param x a PBMExperiment or PBMDesign object
#'
#' @return
#' Filtered object of same class as input \code{se}.
#'
#' @export
#' @author Patrick Kimes
pbmFilterProbes <- function(x) {
    if (!is(x, "PBMExperiment") & !is(x, "PBMDesign"))
        stop("Probe filtering can only be performed with 'PBMExperiment' and 'PBMDesign' objects.")
    
    filters <- x@probeFilter
    if (is(x, "PBMExperiment")) {
        pddata <- rowData(x)[, x@probeCols, drop = FALSE]
    } else {
        pddata <- x@design
    }

    fres <- mapply(function(fn, cl) { fn(pddata[[cl]]) },
                   fn = filters, cl = names(filters),
                   SIMPLIFY = FALSE)

    if (any(vapply(fres, length, numeric(1)) != nrow(pddata))) {
        stop("probeFilter functions should return logical vectors of length equal to number of probes.\n",
             "number of probes = ", nrow(pddata), "\n",
             "length of probeFilter output. \n",
             paste0(paste0("    ", names(filters), " = ", vapply(fres, length, numeric(1))),
                    collapse = "\n"))
    }
    if (!all(vapply(fres, is, logical(1), class2 = "logical"))) {
        stop("probeFilter functions should return logical vectors of length equal to number of probes.\n",
             "number of probes = ", nrow(pddata), "\n",
             "length of probeFilter output. \n",
             paste0(paste0("    ", names(filters), " = ", vapply(fres, class, character(1))),
                    collapse = "\n"))
    }

    ## combine filtering results
    fres <- simplify2array(fres)
    fres <- apply(fres, 1, all)
    
    ## filter according to probeFilter rules
    if (is(x, "PBMExperiment")) {
        x <- x[fres, , drop = FALSE]
    } else {
        x@design <- x@design[fres, , drop = FALSE]
    }
    return(x)
}
