#' fit probe distribution
#' 
#' This is a collection of several models to fit the probe intensities. It applys
#' normal+exponential convolutional model from \code{limma::normexp.fit} and 
#' normal+gamma convolutional model from \code{NormalGamma::normgam.fit}. This function returns the
#' parameter estimate of corresponding mode.
#' 
#' @param se SummarizedExperiment object containing GPR intensity information. Notice usually the intensity should
#'        should be raw intensity in non-log scale.
#' @param assay_name string name of the assay to normalize. (default = "fore")
#' @param model character string specifying model for fitting intensityies. 
#'        Must be one of "NormGam", "NormExp". (default = "NormGam")
#' @param method method used to estimate the parameters for normal+exponential model. See more details from 
#'        \code{limma::normexp.fit}. (default = "mle")
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to normalization. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#'        
#' @return A dataframe with rows as condition and cols as parameter
#' 
#' @import dplyr tidyr
#' @importFrom limma normexp.fit
#' @importFrom NormalGamma normgam.fit
#' @export
#' @author Dongyuan Song, Patrick Kimes
fitProbeDistribution <- function(se,
                                 assay_name = "fore",
                                 model = c("NormGam", "NormExp"),
                                 method = "mle",
                                 .filter = 1L) {
  stopifnot(assay_name %in% assayNames(se))
  
  match.arg(model)
  
  se <- pbmFilterProbes(se, .filter)
  
  # extract intensity matrix
  new_assay <- as.matrix(assay(se, assay_name))
  assay_vector <- apply(new_assay, 2, na.omit)
  
  if (model == "NormGam") {
    par_e <- sapply(assay_vector, function(x) {NormalGamma::normgam.fit(x)$par}, simplify = TRUE)
    par_e <- as.data.frame(t(par_e))
    par_e <- cbind(colData(se)$condition, par_e)
    colnames(par_e) <- c("condition", "norm_mean", "norm_sigma", "gam_shape", "gam_scale")
    
  } else {
    
    par_e <- sapply(assay_vector, function(x) {limma::normexp.fit(x, method = method)$par}, simplify = TRUE)
    par_e[2, ] <- exp(par_e[2, ])
    par_e[3, ] <- exp(par_e[3, ])
    par_e <- as.data.frame(t(par_e))
    par_e <- cbind(colData(se)$condition, par_e)
    colnames(par_e) <- c("condition", "norm_mean", "norm_sigma", "exp_mean")
    }
  
  return(par_e)
  
}