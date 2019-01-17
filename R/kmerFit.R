#' Fit k-mer Probe Set Models
#'
#' @description
#' After performing probe-level aggregation across samples, this function applies probe
#' set aggregation to obtain k-mer level estimates of affinity. 
#'
#' @param se SummarizedExperiment of k-mer, probe pairs by \code{kmerAggregate}.
#' @param method character name of method to use for estimating cross-probe variance
#'        in each k-mer probe set. Currently, the non-iterative DerSimonian-Laird ("dl")
#'        and two-step Dersimonian-Laird ("dl2") methods are supported. (default = "dl")
#' @param baseline character name of baseline condition to use for calculating
#'        contrasts. (default = NULL)
#' 
#' @return
#' SummarizedExperiment of estimated k-mer affinities and differences.
#'
#' @importFrom dplyr select_ group_by left_join ungroup do mutate
#' @importFrom tidyr unnest spread
#' @export
#' @author Patrick Kimes
kmerFit <- function(se, method = c("dl", "dl2"), baseline = NULL) {
    stopifnot(is(se, "SummarizedExperiment"))
    method <- match.arg(method)
    
    ## kmers should be in rowData as "seq"
    if (! "seq" %in% names(rowData(se))) {
        stop("k-mers must be in rowData as column 'seq'")
    }
    kmers <- levels(rowData(se)$seq)
    
    ## turn assays into single tibble
    rd <- rowData(se)
    rd <- as.tibble(as.data.frame(rd, optional = TRUE))
    adat <- replicate(ncol(se), rd, simplify = FALSE)
    names(adat) <- colnames(se)
    adat <- bind_rows(adat, .id = "condition")
    adat$beta <- as.numeric(as.matrix(assay(se, "beta")))
    adat$sd <- as.numeric(as.matrix(assay(se, "sd")))

    ## only keep necessary columns
    adat <- dplyr::select(adat, condition, seq, beta, sd)

    ## filter out NA probe results
    adat <- dplyr::filter(adat, !is.na(beta), !is.na(sd))
    
    ## compute probe set mixed effects model for each k-mer and condition
    adat <- tidyr::nest(adat, -condition, -seq)
    
    if (method == "dl") {
        adat <- dplyr::mutate(adat,
                              res = lapply(data, function(x) {
                                  dl_estimator(x$beta, x$sd^2, nrow(x))
                              }))
    } else if (method == "dl2") {
        adat <- dplyr::mutate(adat,
                              res = lapply(data, function(x) {
                                  dl2_estimator(x$beta, x$sd^2, nrow(x))
                              }))
    } else {
        stop("specified method is invalid")
    }

    adat <- dplyr::mutate(adat, beta = vapply(res, `[[`, numeric(1), "betaME"))
    adat <- dplyr::mutate(adat, betaVar = vapply(res, `[[`, numeric(1), "varME"))
    
    ## rearrange data to compute difference covariances
    ares <- dplyr::mutate(adat, data = mapply(function(d, r) {
        d$resid <- d$beta - r$betaME
        d$tau2 <- r$tau2
        d$totvar <- d$sd^2 + d$tau2
        d
    }, d = data, r = res, SIMPLIFY = FALSE))
    ares <- dplyr::select(ares, condition, seq, data)
    ares <- tidyr::spread(ares, condition, data)
    ares <- dplyr::rename(ares, "bldata" =  !!baseline)
    ares <- tidyr::gather(ares, condition, data, -seq, -bldata)

    ## compute empirical covariance and upperbound (assuming indep sampling error)
    ares <- dplyr::mutate(ares,
                          ecov = mapply(function(x, y) { sum(x$resid * y$resid, na.rm = TRUE) / (nrow(x) - 1) },
                                        x = data, y = bldata))
    ares <- dplyr::mutate(ares, 
                          covmax = mapply(function(x, y) { sqrt(x$tau2[1] * y$tau2[1]) },
                                          x = data, y = bldata))
    ares <- dplyr::mutate(ares, ecovT = pmin(ecov, covmax, na.rm = TRUE))

    ## compute total variance of estimator using error covariance estimate
    ares <- dplyr::mutate(ares, 
                          dvar = mapply(function(d1, d2, e) {
                              v1 <- 1 / sum(1/d1$totvar, na.rm = TRUE)
                              v2 <- 1 / sum(1/d2$totvar, na.rm = TRUE)
                              ccv <- sum(1 / (d1$totvar * d2$totvar), na.rm = TRUE) * e * v1 * v2
                              v1 + v2 - 2 * ccv
                          }, d1 = data, d2 = bldata, e = ecovT))
    ares <- dplyr::select(ares, seq, condition, dvar)

    ## clean up all results
    bdat <- dplyr::select(adat, condition, seq, beta)
    bdat <- tidyr::spread(bdat, condition, beta)
    bdat <- dplyr::rename(bdat, "blbeta" =  !!baseline)
    bdat <- tidyr::gather(bdat, condition, beta, -seq, -blbeta)
    bdat <- dplyr::mutate(bdat, M = beta - blbeta)
    bdat <- dplyr::mutate(bdat, A = (beta + blbeta) / 2)

    ## tidy results to assays
    assaylist <- list(affinityEstimate = .tidycol2mat(adat, "beta", kmers, colnames(se)),
                      affinityVariance = .tidycol2mat(adat, "betaVar", kmers, colnames(se)),
                      contrastDifference = .tidycol2mat(bdat, "M", kmers, colnames(se)),
                      contrastAverage = .tidycol2mat(bdat, "A", kmers, colnames(se)),
                      contrastVariance = .tidycol2mat(ares, "dvar", kmers, colnames(se)))
    
    rdat <- dplyr::select(adat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]

    SummarizedExperiment(assays = assaylist, rowData = rdat)
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
