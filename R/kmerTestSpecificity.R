#' @title Test for k-mer differential specificities
#'
#' @description
#' Using estimated k-mer level affinities returned by \code{kmerFit}, this
#' function tests for differential specificities across conditions for each k-mer
#' separately. The call to \code{kmerFit} must have been made with
#' \code{contrasts = TRUE}. Here, shared specificity is defined as a common ordering of
#' k-mer affinities across conditions. Therefore, differential specificity is tested by
#' assessing the statistical significance of deviations from a common trend
#' of k-mer affinities observed between conditions.
#'
#' For each condition, the common trend is estimated by fitting either a LOESS curve
#' or B-spline smoother to the corresponding columns of the \code{"contrastAverage"} (x)
#' and \code{"contrastDifference"} (y) assays of the input SummarizedExperiment object.
#' If \code{useref = TRUE}, instead of the \code{"contrastAverage"}, the
#' \code{"affinityEstimate"} values of the reference condition are used as the x-axis.
#' The difference between the observed differential affinity and the
#' estimated trend is taken as the estimate of differential specificity.
#' 
#' When \code{method = "subset"} (default), two trends are fit using weighted
#' orinary linear regression and a B-spline basis. The two trends are meant to capture
#' any differential specificities that may be smoothed out by fitting a single trend
#' to the data. The two trends are identified by first splitting the data along
#' the x-axis into 20 equal-width bins and fitting a
#' two-component normal mixture to the \code{"contrastDifference"} values of k-mers
#' in each bin. The upper and lower components in each cluster are grouped across
#' the bins. The two trends are fit separately by weighing each k-mer based
#' on the likelihood of belonging to either the upper or lower clusters. Finally,
#' a single trend is selected to be used for testing by taking the trend which
#' minimizes the mean residual across the top 50\% of k-mers.
#' 
#' @description
#' After estimating k-mer level affinities using probe-set aggregate, this
#' function tests for differential specificity between conditions. 
#'
#' @param se SummarizedExperiment of k-mer results from \code{kmerFit}.
#' @param method string value specifying method to use for fitting
#'        specificity trend. Must be one of \code{"subset"}, \code{"bs"},
#'        or \code{"loess"}. (default = "subset")
#' @param span span parameter to be passed to \code{limma::loessFit}. Only
#'        used when \code{method = "loess"}. (default = 0.05)
#' @param useref logical value whether to fit specificity trend as function of
#'        reference intensities rather than the average of reference and
#'        variant intensities. Must be TRUE if \code{method = "subset"}.
#'        (defualt = TRUE)
#' @param ... other parameters to pass to \code{limma::loessFit}.
#' 
#' @return
#' SummarizedExperiment of k-mer differential specificity results with the following
#' assays:
#' 
#' \itemize{
#' \item \code{"contrastAverage"}: input k-mer average affinities.
#' \item \code{"contrastDifference"}: input k-mer differential affinities.
#' \item \code{"contrastVariance"}: input k-mer differential affinity variances.
#' \item \code{"contrastFit"}: estimated specificity trend.
#' \item \code{"contrastResidual"}: residual from estimated specificity trend (\code{contrastDifference - contrastFit}).
#' \item \code{"specificityZ"}: studentized differences (\code{contrastResidual / sqrt(contrastVariance)}).
#' \item \code{"specificityP"}: two-sided tail p-values for studentized differences.
#' \item \code{"specificityQ"}: FDR-controlling Benjamini-Hochberg adjusted p-values.
#' }
#'
#' @import mclust
#' @importFrom splines bs
#' @importFrom broom tidy
#' @importFrom limma loessFit
#' @importFrom stats p.adjust pnorm
#' @importFrom dplyr select group_by left_join ungroup mutate as_tibble
#' @importFrom tidyr nest_legacy unnest_legacy
#' @export
#' @author Patrick Kimes
kmerTestSpecificity <- function(se, method = c("subset", "bs", "loess"), span = 0.05, useref = TRUE, ...) {

    stopifnot(is(se, "SummarizedExperiment"))
    if (!all(c("contrastAverage", "contrastDifference", "contrastVariance") %in%
             assayNames(se))) {
        stop("Input SummarizedExperiment is missing k-mer contrast estimates.\n",
             "SummarizedExperiment should be created by calling kmerFit(..) with 'contrasts=TRUE'.")
    }
    method <- match.arg(method)
    
    kmers <- rowData(se)$seq

    ## gather data
    cdat <- broom::tidy(se, c("contrastAverage", "contrastDifference", "contrastVariance"))
    cdat <- dplyr::rename(cdat, condition = cname)

    ## use either average between ref/var or just ref as x-axis of trend
    if (useref) {
        bl <- metadata(se)$baseline
        if (is.null(bl)) {
            bl <- grep("-REF$", colnames(se), value = TRUE)
        }
        
        caxis <- broom::tidy(se, "affinityEstimate", long = TRUE)
        caxis <- dplyr::filter(caxis, cname == bl)
        caxis <- dplyr::select(caxis, seq, specificityAxis = affinityEstimate)
        cdat <- dplyr::left_join(cdat, caxis, by = "seq")
    } else {
        cdat <- dplyr::mutate(cdat, specificityAxis = contrastAverage)
    }
    
    cdat <- tidyr::nest_legacy(cdat, -condition)

    if (method == "loess") {
        cdat <- dplyr::mutate(cdat,
                              data = lapply(data, function(x) {
                                  fit <- limma::loessFit(x$contrastDifference, x$specificityAxis,
                                                         span = span, ...)
                                  x$contrastFit <- fit$fitted
                                  x$contrastResidual <- fit$residuals
                                  x
                              }))
    } else if (method == "bs") {
        cdat <- dplyr::mutate(cdat,
                              data = lapply(data, function(x) {
                                  if (all(is.na(x$contrastDifference))) {
                                      x$contrastFit <- NA
                                      x$contrastResidual <- NA
                                  } else {
                                      xknot <- min(x$specificityAxis, na.rm = TRUE) +
                                          diff(range(x$specificityAxis, na.rm = TRUE)) / 2
                                      xns <- splines::bs(x$specificityAxis, knots = xknot,
                                                         degree = 3, intercept = FALSE)
                                      fit <- lm.fit(x = xns, y = x$contrastDifference)
                                      x$contrastFit <- fit$fitted.values
                                      x$contrastResidual <- fit$residuals
                                  }
                                  x
                              }))
    } else if (method == "subset") {
        cdat <- dplyr::mutate(cdat,
                              data = lapply(data, function(x) {
                                  if (all(is.na(x$contrastDifference))) {
                                      x$contrastFit <- NA
                                      x$contrastResidual <- NA
                                  } else {
                                      ## cut into nbin bins for clustering
                                      x <- dplyr::mutate(x, affbin = cut_interval(x$specificityAxis, 20))
                                      x <- dplyr::mutate(x, affbin = LETTERS[as.numeric(affbin)])
                                      x <- tidyr::nest_legacy(x, -affbin)
                                      x <- dplyr::mutate(x, n = vapply(data, nrow, numeric(1L)))
                                      x <- dplyr::arrange(x, affbin)
                                      ## merge small bins at beginning end (assume middle bins are fine)
                                      x <- dplyr::mutate(x,
                                                         nsum_rev = rev(cumsum(rev(n))),
                                                         nsum_fwd = cumsum(n))
                                      x <- dplyr::mutate(x,
                                                         affbin = ifelse(nsum_rev < 50, tail(affbin[nsum_rev >= 50], 1), affbin),
                                                         affbin = ifelse(nsum_fwd < 50, head(affbin[nsum_fwd >= 50], 1), affbin))
                                      x <- dplyr::select(x, -n, -nsum_rev, -nsum_fwd)
                                      x <- tidyr::unnest_legacy(x)
                                      x <- tidyr::nest_legacy(x, -affbin)
                                      ## cluster within bins
                                      x <- dplyr::mutate(x,
                                                         mcl = lapply(data, function(z) {
                                                             if (nrow(z) > 1L) {
                                                                 mclust::Mclust(z$contrastDifference, G = 2,
                                                                                modelNames = "E", verbose = FALSE)
                                                             } else {
                                                                 NA
                                                             }
                                                         }))
                                      ncomps <- vapply(x$mcl, function(z) { z$parameters$variance$G },
                                                       numeric(1L))
                                      if (any(ncomps != 2L)) {
                                          stop("Some bins split into fewer than 2 clusters.")
                                      }
                                      x <- dplyr::mutate(x, 
                                                         data = mapply(function(d, m) {
                                                             if (is(m, "Mclust")) {
                                                                 d$cidx <- m$classification
                                                                 d$var <- m$parameters$variance$sigmasq[1]
                                                                 d$p1 <- dnorm(d$contrastDifference, m$parameters$mean[1],
                                                                               m$parameters$variance$sigmasq[1]^.5)
                                                                 d$p2 <- dnorm(d$contrastDifference, m$parameters$mean[2],
                                                                               m$parameters$variance$sigmasq[1]^.5)
                                                                 d$p1 <- d$p1 / (d$p1 + d$p2)
                                                                 d$p2 <- 1 - d$p1
                                                             } else {
                                                                 d$cidx <- 1L
                                                                 d$p1 <- d$p2 <- d$var <- NA
                                                             }
                                                             d
                                                         }, d = data, m = mcl, SIMPLIFY = FALSE))
                                      x <- dplyr::select(x, -mcl)
                                      x <- tidyr::unnest(x, data)
                                      ## don't use multiple clusters in lower 50% of 8-mers
                                      x <- dplyr::mutate(x,
                                                         p1 = ifelse(specificityAxis <
                                                                     quantile(specificityAxis, .5),
                                                                     1, p1),
                                                         p2 = ifelse(specificityAxis <
                                                                     quantile(specificityAxis, .5),
                                                                     1, p2))
                                      ##
                                      xknot <- min(x$specificityAxis, na.rm = TRUE) + 
                                          diff(range(x$specificityAxis, na.rm = TRUE)) / 2
                                      xns <- splines::bs(x$specificityAxis, knots = xknot, 
                                                         degree = 3, intercept = FALSE)
                                      fit1 <- lm.wfit(x = xns, y = x$contrastDifference,
                                                      w = x$p1)
                                      x$contrastFit1 <- fit1$fitted.values
                                      x$contrastResidual1 <- fit1$residuals
                                      fit2 <- lm.wfit(x = xns, y = x$contrastDifference,
                                                      w = x$p2)
                                      x$contrastFit2 <- fit2$fitted.values
                                      x$contrastResidual2 <- fit2$residuals
                                      ## determine better fit
                                      if (mean(abs(x$contrastResidual1[x$specificityAxis >
                                                                       quantile(x$specificityAxis, .5)]), na.rm = TRUE) <
                                          mean(abs(x$contrastResidual2[x$specificityAxis >
                                                                       quantile(x$specificityAxis, .5)]), na.rm = TRUE)) {
                                          x <- dplyr::mutate(x,
                                                             contrastFit = contrastFit1,
                                                             contrastResidual = contrastResidual1)
                                      } else {
                                          x <- dplyr::mutate(x,
                                                             contrastFit = contrastFit2,
                                                             contrastResidual = contrastResidual2)
                                      }
                                  }
                                  x
                              }))
    }
    
    cdat <- tidyr::unnest_legacy(cdat)

    ## compute z-scores, p-values, and adjusted p-values
    cdat <- dplyr::mutate(cdat, specificityZ = contrastResidual / sqrt(contrastVariance))
    cdat <- dplyr::mutate(cdat, specificityP = 2 * stats::pnorm(-abs(specificityZ)))
    cdat <- dplyr::group_by(cdat, condition)
    cdat <- dplyr::mutate(cdat, specificityQ = stats::p.adjust(specificityP, method = "BH"))
    cdat <- dplyr::ungroup(cdat)

    ## tidy results to assays
    assaylist <- list(contrastAverage = .tidycol2mat(cdat, "contrastAverage", kmers, colnames(se)),
                      contrastDifference = .tidycol2mat(cdat, "contrastDifference", kmers, colnames(se)),
                      contrastVariance = .tidycol2mat(cdat, "contrastVariance", kmers, colnames(se)),
                      contrastFit = .tidycol2mat(cdat, "contrastFit", kmers, colnames(se)),
                      contrastResidual = .tidycol2mat(cdat, "contrastResidual", kmers, colnames(se)),
                      specificityZ = .tidycol2mat(cdat, "specificityZ", kmers, colnames(se)),
                      specificityP = .tidycol2mat(cdat, "specificityP", kmers, colnames(se)),
                      specificityQ = .tidycol2mat(cdat, "specificityQ", kmers, colnames(se)),
                      specificityAxis = .tidycol2mat(cdat, "specificityAxis", kmers, colnames(se)))

    rdat <- dplyr::select(cdat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]
    
    SummarizedExperiment(assays = assaylist, rowData = rdat)
}
