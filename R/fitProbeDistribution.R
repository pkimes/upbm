#' @title Fit Probe-Level Convolution Distributions
#' 
#' @description
#' This function fits either a normal+exponential or normal+gamma convolution model to
#' the probe-level intensities of a PBMExperiment.
#' The normal+exponential convolution model is fit using \code{limma::normexp.fit} and 
#' the normal+gamma convolution model is fit using \code{NormalGamma::normgam.fit}.
#' This function returns the model parameter estimates for each sample in a data.frame.
#'
#' For more details on the estimated parameters, see the underlying functions.
#' 
#' @param pe a PBMExperiment object containing GPR intensity information.
#' @param assay string name of the assay. (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param model character string specifying model for fitting intensityies. 
#'        Must be one of "NormGam", "NormExp". (default = "NormGam")
#' @param method method used to estimate the parameters for normal+exponential model. See more details from 
#'        \code{limma::normexp.fit}. (default = "mle")
#'
#' @return
#' A data.frame with scans as rows and parameters in columns.
#'
#' @references
#' \itemize{
#' \item Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47.
#' \item Plancade, S., Rozenholc, Y., & Lund, E. (2012). Generalization of the normal-exponential model: exploration of a more accurate parametrisation for the signal distribution on Illumina BeadArrays. BMC Bioinformatics, 13(1), 329.
#' }
#' 
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom limma normexp.fit
#' @importFrom NormalGamma normgam.fit
#' @importFrom stats na.omit
#' @export
#' @author Dongyuan Song
fitProbeDistribution <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                                 model = c("NormGam", "NormExp"), method = "mle") {
    stopifnot(is(pe, "PBMExperiment"))
    stopifnot(assay %in% SummarizedExperiment::assayNames(pe))
    
    model <- match.arg(model)
    
    pe <- pbmFilterProbes(pe)
    
    ## extract intensity matrix
    new_assay <- as.matrix(SummarizedExperiment::assay(pe, assay))
    assay_vector <- apply(new_assay, 2, stats::na.omit)
    
    if (model == "NormGam") {
        par_e <- sapply(assay_vector, function(x) {NormalGamma::normgam.fit(x)$par}, simplify = TRUE)
        par_e <- as.data.frame(t(par_e))
        par_e <- cbind(colData(pe)$condition, par_e)
        colnames(par_e) <- c("condition", "norm_mean", "norm_sigma", "gam_shape", "gam_scale")
        
    } else {
        par_e <- sapply(assay_vector, function(x) {limma::normexp.fit(x, method = method)$par}, simplify = TRUE)
        par_e[2, ] <- exp(par_e[2, ])
        par_e[3, ] <- exp(par_e[3, ])
        par_e <- as.data.frame(t(par_e))
        par_e <- cbind(colData(pe)$condition, par_e)
        colnames(par_e) <- c("condition", "norm_mean", "norm_sigma", "exp_mean")
    }
    
    return(par_e)
}
