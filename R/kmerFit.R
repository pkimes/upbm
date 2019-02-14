#' @title Fit k-mer probe set models
#'
#' @description
#' Given a PBMExperiment of probe-level summaries returned by \code{probeFit} and a list of k-mers,
#' this function applies probe set aggregation to obtain k-mer level estimates of affinity and
#' variance on the scale of log2 signal intensities. Additionally, if \code{contrasts = TRUE},
#' effect size and variance estimates are also returned for differential k-mer affinities against
#' a baseline condition specified with \code{baseline=}.
#'
#' The output can be passed to \code{\link{kmerTestContrast}}, \code{\link{kmerTestAffinity}},
#' \code{\link{kmerTestSpecificity}} to perform various statistical tests with the estimated
#' k-mer statistics.
#'
#' @param pe a PBMExperiment object containing probe-level summarized intensity data
#'        returned by \code{\link{probeFit}}.
#' @param kmers a character vector of k-mers. (default = \code{uniqueKmers(8L)})
#' @param positionbias a logical value whether to correct for bias due to position
#'        of k-mer along probe sequence. (default = TRUE)
#' @param method a character name specifying the method to use for estimating cross-probe
#'        variance in each k-mer probe set. Currently, the two-step DerSimonian-Kacker ("dl2") and
#'        non-iterative DerSimonian-Laird ("dl") methods are supported. (default = "dl2")
#' @param contrasts a logical value whether to compute contrasts for all columns against a
#'        specified \code{baseline} column. (default = TRUE)
#' @param baseline a character string specifying the baseline condition across \code{pe} columns to
#'        use when calculating contrasts. If not specified and set to NULL, the baseline
#'        value is guessed by looking for ``ref" in the column names of \code{pe}. If a unique
#'        matching value is not found, an error is thrown.
#'        This parameter is ignored when \code{contrasts = FALSE}.
#'        (default = NULL)
#' @param outlier_cutoff a numeric threshold used for filtering probes from k-mer
#'        probe sets before fitting each k-mer level model. The threshold is
#'        applied to the absolute value of an approximate robust studentized residual
#'        computed for each probe in each probe set and can be turned off by
#'        setting the value to NULL. By default, approximate 0.5% tails are trimmed.
#'        (default = \code{stats::qnorm(0.995)})
#' @param outlier_maxp a numeric threshold on the maximum proportion of probes to filter
#'        for each k-mer probe set according to \code{outlier_cutoff}. This should be
#'        set to a reasonably small value to avoid over-filtering based on the approximate
#'        residual threshold. (default = 0.2)
#' @param verbose a logical value whether to print verbose output during analysis. (default = FALSE)
#' 
#' @return
#' SummarizedExperiment of estimated k-mer affinities and differences with some or all
#' of the following assays:
#'
#' \itemize{
#' \item \code{"affinityEstimate"}: k-mer affinities.
#' \item \code{"affinityVariance"}: k-mer affinity variances.
#' \item \code{"contrastDifference"}: (optional) k-mer differential affinities with \code{baseline} condition.
#' \item \code{"contrastAverage"}: (optional) k-mer average affinities with \code{baseline} condition.
#' \item \code{"contrastVariance"}: (optional) k-mer differential affinity variances.
#' }
#'
#' If computed, the values of the \code{"contrast"} assays will be NA for the specified
#' baseline condition.
#' 
#' @details
#' By default, probe intensities are corrected within each k-mer probe set
#' to account for biases introduced by where the k-mer lies along the probe sequence. Bias
#' correction is performed such that the mean cross-probe intensity for each k-mer is
#' (typically) unchanged. This bias correction step only serves to reduce the cross-probe variance 
#' and improve downstream inference for each k-mer.
#'
#' For many low affinity k-mers, probe sets may include several probes with high intensity due to
#' the k-mer sharing a probe with a separate high affinity k-mer. These probes do not provide
#' meaningful affinity information for the lower affinity k-mer. To adjust for this possibility,
#' outlier probes are filtered from each k-mer probe set prior after position bias correction, but
#' before aggregation. Probes with large approximate studentized residuals are
#' filtered from each probe set according to a user-specified threshold (\code{outlier_cutoff}).
#' However, to prevent overfiltering, a maximum proportion of probes to filter from any probe set
#' should also be specified (\code{outlier_maxp}).
#' 
#' After bias correction and probe filtering, a meta analysis model is fit to each probe set.
#' Under this model, cross-probe variances are estimated using either the DerSimonian and Kacker (2007)
#' or DerSimonian and Laird (1986) estimator. The estimated k-mer affinities and variances are
#' included in the returned SummarizedExperiment as two assays, \code{"affinityEstimate"} and
#' \code{"affinityVariance"}.
#'
#' If \code{contrast = TRUE}, k-mer differential affinities, the corresponding variances, and
#' average affinities are also returned as three assays, \code{"contrastDifference"},
#' \code{"contrastVariance"}, and \code{"contrastAverage"}. Positive differential affinities indicate
#' higher affinity relative to the baseline condition. 
#' 
#' @references
#' If using \code{method = "dl2"} cross-probe variance estimator:
#' \itemize{
#' \item DerSimonian, R., & Kacker, R. (2007). Random-effects model for meta-analysis of clinical trials: an update. Contemporary Clinical Trials, 28(2), 105-114.
#' }
#' If using \code{method = "dl"} cross-probe variance estimator:
#' \itemize{
#' \item DerSimonian, R., & Laird, N. (1986). Meta-analysis in clinical trials. Controlled Clinical Trials, 7(3), 177-188.
#' }
#' Cross-probe variance estimation code adapted from:
#' \itemize{
#' \item Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. URL: http://www.jstatsoft.org/v36/i03/
#' }
#' 
#' @seealso \code{\link{probeFit}}, \code{\link{uniqueKmers}}
#' @importFrom dplyr select_ group_by left_join ungroup do mutate arrange one_of
#' @importFrom tidyr unnest spread
#' @importFrom stats qnorm
#' @export
#' @author Patrick Kimes
kmerFit <- function(pe, kmers = uniqueKmers(8L), positionbias = TRUE,
                    method = c("dl2", "dl"), contrasts = TRUE, baseline = NULL,
                    outlier_cutoff = stats::qnorm(0.995), outlier_maxp = 0.2,
                    verbose = FALSE) {
    
    stopifnot(is(pe, "PBMExperiment"))
    stopifnot(c("beta", "sd") %in% assayNames(pe))
    method <- match.arg(method)
    
    stopifnot(is.null(outlier_cutoff) ||
              (is.numeric(outlier_cutoff) & outlier_cutoff > 0))
    stopifnot(is.numeric(outlier_maxp) && outlier_maxp >= 0 && outlier_maxp <= 1)

    if (verbose) {
        cat("|| upbm::kmerFit \n")
        cat("|| - Starting k-mer probe set model fitting for", ncol(pe), "conditions.\n")
    }
    
    ## check if specificed 'baseline' is valid
    if (contrasts) {
        if (is.null(baseline)) {
            baseline <- grep("ref$", colnames(pe), value = TRUE, ignore.case = TRUE)
            if (length(baseline) > 1) {
                stop("Too many candidate baseline states in column names: ",
                     paste0(baseline, collapse = ", "), ".\n",
                     "Specify correct baseline column name w/ 'baseline ='.")
            }
        } else {
            if (! baseline %in% colnames(pe)) {
                stop(baseline, " is not a column name of the SummarizedExperiment.\n",
                     "Specify correct baseline column name w/ 'baseline ='.")
            }
        }
    }
    
    ## check kmers specified
    if (!is.vector(kmers, mode = "character")) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences as character strings.")
    }
    if (length(unique(nchar(kmers))) != 1L) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences of equal length.")
    }

    ## filter probes
    if (verbose) {
        cat("|| - Filtering probes according to", length(pe@probeFilter),
            "probeFilter rule(s).\n")
        ntotal <- nrow(pe)
    }

    ## filter using rules
    pe <- pbmFilterProbes(pe)
    
    if (verbose) {
        cat("|| - Data filtered from", ntotal, "probes to", nrow(pe), "probes.\n")
    }

    if (verbose && length(pe@probeTrim > 0L)) {
        cat("|| - Trimming probes according to probeTrim settings:",
            paste0("[", paste0(pe@probeTrim, collapse = ", "), "]"), "\n")
    }

    ## trim probe sequences
    pe <- pbmTrimProbes(pe)

    ## find mapping between kmers and probes - fine to assume cols exist
    kmermap <- mapKmers(PBMDesign(pe), kmers)
    kmermap <- dplyr::select_at(kmermap,
                                c(setdiff(names(kmermap), pe@probeCols),
                                  "probeID"))
    
    ## create table of probe-level betas
    bdat <- broom::tidy(pe, "beta", long = FALSE)
    bdat <- dplyr::left_join(kmermap, bdat, by = "probeID")
    bdat <- dplyr::select(bdat, seq, pos, probeID, dplyr::one_of(colnames(pe)))
    
    if (positionbias) {
        if (verbose) {
            cat("|| - Estimating sequence position biases for probe sets (positionbias = TRUE).\n")
        }
        ## reshape table and compute probe set bias per k-mer
        bdat <- tidyr::gather(bdat, sample, value, -pos, -probeID, -seq)
        
        bdat <- dplyr::group_by(bdat, seq, sample)
        bdat <- dplyr::mutate(bdat, pmean = mean(value, na.rm = TRUE),
                              pbias = value - pmean)
        bdat <- dplyr::ungroup(bdat)

        ## compute average bias over 2% bins
        bdat <- dplyr::group_by(bdat, sample)
        bdat <- dplyr::mutate(bdat, qbin = as.numeric(ggplot2::cut_number(pmean, 1 / .02)))
        bdat <- dplyr::group_by(bdat, sample, qbin, pos)
        bdat <- dplyr::mutate(bdat, pbias = mean(pbias, na.rm = TRUE))
        bdat <- dplyr::ungroup(bdat)
        bdat <- dplyr::mutate(bdat, value = value - pbias)
        bdat <- dplyr::select(bdat, -pmean, -qbin)

        ## create table of adjusted beta estimates, samples as cols (so slow..)
        bdat_beta <- dplyr::select(bdat, seq, pos, probeID, sample, value)
        bdat_beta <- tidyr::spread(bdat_beta, sample, value)

        ## ## create table of adjusted beta estimates, samples as cols (so slow..)
        ## bdat_pbias <- dplyr::select(bdat, seq, pos, ID, sample, value)
        ## bdat_pbias <- tidyr::spread(bdat_pbias, sample, value)
        ## bdat_pbias <- dplyr::arrange(bdat_pbias, seq, ID, pos)
    } else {
        if (verbose) {
            cat("|| - Skipping sequence position bias estimation (positionbias = FALSE).\n")
        }
        bdat_beta <- bdat
    }

    bdat_beta <- dplyr::arrange(bdat_beta, seq, probeID, pos)

    ## extract consistent table columns
    rowdat <- dplyr::select(bdat_beta, -dplyr::one_of(colnames(pe)))

    ## create table of sd estimates, samples as cols
    bdat_sd <- broom::tidy(pe, "sd", long = FALSE)
    bdat_sd <- dplyr::select(bdat_sd, -dplyr::one_of(setdiff(names(rowData(pe)), "probeID")))
    bdat_sd <- dplyr::left_join(rowdat, bdat_sd, by = "probeID")
    bdat_sd <- dplyr::arrange(bdat_sd, seq, probeID, pos)

    ## turn beta, sd values into tall table and join
    bdat_beta <- tidyr::gather(bdat_beta, condition, beta, -seq, -pos, -probeID)
    bdat_sd <- tidyr::gather(bdat_sd, condition, sd, -seq, -pos, -probeID)

    ## tidyr join operations takes time, so check order of rows and throw error if not same 
    stopifnot(bdat_beta$condition == bdat_sd$condition)
    stopifnot(bdat_beta$seq == bdat_sd$seq)
    stopifnot(bdat_beta$pos == bdat_sd$pos)
    stopifnot(bdat_beta$probeID == bdat_sd$probeID)

    ## simple merge if rows are same order
    adat <- bdat_beta
    adat$sd <- bdat_sd$sd

    ## only keep necessary columns
    adat <- dplyr::select(adat, condition, seq, beta, sd)
    
    if (!is.null(outlier_cutoff)) {
        if (verbose) {
            cat("|| - Filtering outlier probes from probe sets (outlier_cutoff =", paste0(round(outlier_cutoff, 3), ")."), "\n")
        }
        ## compute quick studentized residuals
        ## -- cross-probe var estimate as MAD
        ## -- cross-probe point estimate as median
        adat <- dplyr::group_by(adat, condition, seq)
        adat <- dplyr::mutate(adat, probeZ = (beta - median(beta, na.rm = TRUE)) /
                                        sqrt(mad(beta, na.rm = TRUE)^2 + sd^2))
        if (outlier_maxp < 1) {
            if (verbose) {
                cat("|| - Maximum proportion of filtered probes per probe set was specified (outlier_maxp =", paste0(round(outlier_maxp, 3), ")."), "\n")
            }
            ## compute quantiles of residuals
            adat <- dplyr::mutate(adat,
                                  probeZq = rank(-abs(probeZ), na.last = TRUE, ties.method = "first"),
                                  probeZq = probeZq / sum(!is.na(probeZ)))
            ## prevent more than maxp to be rejected at any cutoff
            adat <- dplyr::mutate(adat, probeZ = ifelse(probeZq > outlier_maxp, 0, probeZ))
            adat <- dplyr::select(adat, -probeZq)
        }
        adat <- dplyr::ungroup(adat)

        ## filter out probes
        if (verbose) {
            cat("|| - Filtering",
                paste0(round(100 * sum(abs(adat$probeZ) > outlier_cutoff, na.rm = TRUE) / nrow(adat), 3), "%"),
                "of probes from k-mer probe sets.\n")
        }
        adat <- dplyr::mutate(adat, beta = ifelse(abs(probeZ) > outlier_cutoff, NA, beta))
        adat <- dplyr::select(adat, -probeZ)
    }

    ## compute probe set mixed effects model for each k-mer and condition
    adat <- tidyr::nest(adat, -condition, -seq)
    if (method == "dl") {
        if (verbose) {
            cat("|| - Estimating cross-probe variance using DerSimonian-Laird estimator.\n")
        }
        adat <- dplyr::mutate(adat,
                              res = lapply(data, function(x) {
                                  ## filter out NA probe results
                                  x <- dplyr::filter(x, !is.na(beta), !is.na(sd))
                                  dl_estimator(x$beta, x$sd^2, nrow(x))
                              }))
    } else if (method == "dl2") {
        if (verbose) {
            cat("|| - Estimating cross-probe variance using DerSimonian-Kacker estimator.\n")
        }
        adat <- dplyr::mutate(adat,
                              res = lapply(data, function(x) {
                                  ## filter out NA probe results
                                  x <- dplyr::filter(x, !is.na(beta), !is.na(sd))
                                  dl2_estimator(x$beta, x$sd^2, nrow(x))
                              }))
    } else {
        stop("specified method is invalid")
    }

    adat <- dplyr::mutate(adat, beta = vapply(res, `[[`, numeric(1), "betaME"))
    adat <- dplyr::mutate(adat, betaVar = vapply(res, `[[`, numeric(1), "varME"))

    ## tidy results to assays
    assaylist <- list(affinityEstimate = .tidycol2mat(adat, "beta", kmers, colnames(pe)),
                      affinityVariance = .tidycol2mat(adat, "betaVar", kmers, colnames(pe)))

    ## rearrange data to compute contrast covariances if needed
    if (contrasts) {
        if (verbose) {
            cat("|| - Estimating cross-condition covariances for contrasts (contrasts = TRUE).\n")
        }
        ares <- dplyr::mutate(adat, data = mapply(function(d, r) {
            d$resid <- d$beta - r$betaME
            d$tau2 <- r$tau2
            d$totvar <- d$sd^2 + d$tau2
            d
        }, d = data, r = res, SIMPLIFY = FALSE))
        ares <- dplyr::select(ares, condition, seq, data)
        ares <- tidyr::spread(ares, condition, data)
        ares <- dplyr::rename(ares, "bldata" = !!baseline)
        ares <- tidyr::gather(ares, condition, data, -seq, -bldata)

        ## compute empirical covariance and upperbound (assuming indep sampling error)
        ares <- dplyr::mutate(ares,
                              ecov = mapply(function(x, y) {
                                  sum(x$resid * y$resid, na.rm = TRUE) /
                                      (sum(!is.na(x$resid) & !is.na(y$resid)) - 1)
                              }, x = data, y = bldata))
        ares <- dplyr::mutate(ares, 
                              covmax = mapply(function(x, y) {
                                  sqrt(x$tau2[1] * y$tau2[1])
                              }, x = data, y = bldata))
        ares <- dplyr::mutate(ares, ecovT = pmin(ecov, covmax, na.rm = TRUE))

        ## compute total variance of difference using error covariance estimate
        ares <- dplyr::mutate(ares, 
                              dvar = mapply(function(d1, d2, e) {
                                  v1 <- 1 / sum(1/d1$totvar, na.rm = TRUE)
                                  v2 <- 1 / sum(1/d2$totvar, na.rm = TRUE)
                                  ccv <- sum(1 / (d1$totvar * d2$totvar), na.rm = TRUE) * e * v1 * v2
                                  v1 + v2 - 2 * ccv
                              }, d1 = data, d2 = bldata, e = ecovT))
        ares <- dplyr::select(ares, seq, condition, dvar)

        ## clean up and tidy contrasts results
        cdat <- dplyr::select(adat, condition, seq, beta)
        cdat <- tidyr::spread(cdat, condition, beta)
        cdat <- dplyr::rename(cdat, "blbeta" = !!baseline)
        cdat <- tidyr::gather(cdat, condition, beta, -seq, -blbeta)
        cdat <- dplyr::mutate(cdat, M = beta - blbeta)
        cdat <- dplyr::mutate(cdat, A = (beta + blbeta) / 2)
        assaylist <- c(assaylist,
                       list(contrastDifference = .tidycol2mat(cdat, "M", kmers, colnames(pe)),
                            contrastAverage = .tidycol2mat(cdat, "A", kmers, colnames(pe)),
                            contrastVariance = .tidycol2mat(ares, "dvar", kmers, colnames(pe))))
    }
    
    rdat <- dplyr::select(adat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]

    se <- SummarizedExperiment(assays = assaylist, rowData = rdat)
    if (verbose) {
        cat("|| - Finished k-mer probe set model fitting.\n")
        cat("|| - Returning SummarizedExperiment with", nrow(se), "rows and", ncol(se), "columns.\n")
    }
    return(se)
}


## helper to turn tidy table column into matrix for SE
.tidycol2mat <- function(x, cn, km, s) {
    x <- dplyr::select(x, seq, condition, !!cn)
    x <- tidyr::spread(x, condition, !!cn)
    x <- x[match(km, x$seq), sort(names(x)), ]
    x <- as.matrix(dplyr::select(x, -seq))
    if (!all(colnames(x) %in% s))
        stop("column names do not match specified samples")
    if (ncol(x) < length(s)) {
        x <- cbind(matrix(NA, nrow(x), length(s) - ncol(x)), x)
        colnames(x)[seq_len(length(setdiff(s, colnames(x))))] <- setdiff(s, colnames(x))
    }
    x[, s, drop = FALSE]
}
