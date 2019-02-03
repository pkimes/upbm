## source: `metafor` CRAN package, R/misc.func.hidden.r (git hash: f037e1b)
.invcalc <- function(X, W, k) {
   sWX <- sqrt(W) %*% X
   res.qrs <- qr.solve(sWX, diag(k))
   return(tcrossprod(res.qrs))
}

#' DerSimonian and Laird Estimator
#'
#' @description
#' This is an implementation of the DerSimonian and Laird one step estimator of
#' cross-study variance, originally proposed in the context of meta analysis,
#' adapted from the \code{metafor} package. The function is used to estimate
#' the cross-probe variance for each k-mer probe set.  
#' 
#' @param Y probe effect sizes
#' @param vi probe variances
#' @param k number of probes
#'
#' @return
#' list of estimates:
#' \itemize{
#' \item betaFE: effect size with no cross-study variance
#' \item varFE: total variance with no cross-study variance
#' \item betaME: effect size with cross-study variance
#' \item varME: total variance with cross-study variance
#' \item tau2: cross-study variance 
#' }
#'
#' @references
#' \itemize{
#' \item DerSimonian, R., & Laird, N. (1986). Meta-analysis in clinical trials. Controlled Clinical Trials, 7(3), 177-188.
#' \item Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. URL: http://www.jstatsoft.org/v36/i03/
#' }
#'
#' @keywords internal
#' @author Patrick Kimes
dl_estimator <- function(Y, vi, k) {
    X     <- rep(1, k)
    p     <- 1
    
    wi    <- 1/vi
    W     <- diag(wi, nrow = k, ncol = k)
    stXWX <- .invcalc(X = X, W = W, k = k)
    stXWX_tXW <- stXWX %*% crossprod(X, W)
    P     <- W - W %*% X %*% stXWX_tXW
    RSS   <- crossprod(Y, P) %*% Y
    trP   <- sum(diag(P))

    tau2  <- max((RSS - k + p) / trP, 0)
    
    betaFE  <- as.numeric(stXWX_tXW %*% Y)
    varFE <- 1 / sum(wi)

    if (tau2 > 0) {
        W_     <- diag(1 / (vi + tau2), nrow = k, ncol = k)
        M_     <- diag(vi + tau2, nrow = k, ncol = k)
        stXWX_ <- .invcalc(X = X, W = W_, k = k)
        betaME  <- as.numeric(stXWX_ %*% crossprod(X, W_) %*% Y)
        varME <- diag(stXWX_)
    } else {
        betaME <- betaFE
        varME <- varFE
    }
    list(betaFE = betaFE, varFE = varFE,
         betaME = betaME, varME = varME,
         tau2 = tau2)
}

#' Two-Step DerSimonian aand Kacker Estimator
#'
#' @description
#' This is an implementation of the DerSimonian and Kacker two-step estimator of
#' cross-study variance, originally proposed as an improvement over the one-step
#' DerSimonian and Laird estimator in the context of meta analysis.
#' 
#' @param Y probe effect sizes
#' @param vi probe variances
#' @param k number of probes
#'
#' @return
#' list of estimates:
#' \itemize{
#' \item betaFE: effect size with no cross-study variance
#' \item varFE: total variance with no cross-study variance
#' \item betaME: effect size with cross-study variance
#' \item varME: total variance with cross-study variance
#' \item tau2: cross-study variance 
#' }
#' 
#' @references
#' \itemize{
#' \item DerSimonian, R., & Kacker, R. (2007). Random-effects model for meta-analysis of clinical trials: an update. Contemporary Clinical Trials, 28(2), 105-114.
#' }
#'
#' @keywords internal
#' @author Patrick Kimes
dl2_estimator <- function(Y, vi, k) {
    X   <- rep(1, k)
    p   <- 1
    res <- dl_estimator(Y, vi, k)
    
    wi    <- 1 / (vi + res$tau2)
    W     <- diag(wi, nrow = k, ncol = k)
    stXWX <- .invcalc(X = X, W = W, k = k)
    stXWX_tXW <- stXWX %*% crossprod(X, W)
    P     <- W - W %*% X %*% stXWX_tXW
    RSS   <- crossprod(Y, P) %*% Y
    trP   <- sum(diag(P))

    SSE <- sum(wi * vi) - sum(wi^2 * vi) / sum(wi)
    tau2  <- max((RSS - SSE) / trP, 0)

    betaFE  <- as.numeric(stXWX_tXW %*% Y)
    varFE <- 1 / sum(wi)
    
    if (tau2 > 0) {
        W_     <- diag(1 / (vi + tau2), nrow = k, ncol = k)
        M_     <- diag(vi + tau2, nrow = k, ncol = k)
        stXWX_ <- .invcalc(X = X, W = W_, k = k)
        betaME  <- as.numeric(stXWX_ %*% crossprod(X, W_) %*% Y)
        varME <- diag(stXWX_)
    } else {
        betaME <- betaFE
        varME <- varFE
    }
    list(betaFE = betaFE, varFE = varFE,
         betaME = betaME, varME = varME,
         tau2 = tau2)
}


