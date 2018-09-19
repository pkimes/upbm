#' Fit Probe-Level Convolution Distributions
#' 
#' This is a collection of several models to fit a normal+exponential or
#' normal+gamma convolution model to the probe-level intensities of a PBM experiment.
#' The normal+exponential convolutional model is fit using \code{limma::normexp.fit} and 
#' the normal+gamma convolutional model is fit using \code{NormalGamma::normgam.fit}.
#' This function returns the model parameter estimates for each sample in the
#' SummarizedExperiment object.
#' 
#' @param se SummarizedExperiment object containing GPR intensity information.
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
  
  model <- match.arg(model)
  
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
