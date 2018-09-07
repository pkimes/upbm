#' qsmooth Normalization (TESTING)
#' 
#' Given a list of SummarizedExperiment objects containing GPR intensity 
#' information, this function performs normalization across samples using qsmooth's
#' \code{qsmooth} function, and returns the same list of
#' SummarizedExperiment objects with the normalized intensities.
#' Currently, since qsmooth does not allow inputs with missing values, this function just removes all
#' probes with NAs and shows a warning message. 
#' \emph{This code is still under testing}.
#'
#' @param se_list a list object of SummarizedExperiment objects containing GPR intensity 
#'        information.
#' @param ref_rep integer specifying which replicate will be used for scaling before normalization. (default = 1)
#' @param assay_name string name of the assay to normalize. (default = "fore")
#' @param scaling_between_rep logical whether to scale two reps into same median. (default = TRUE)
#' @param scaling_within_rep logical whether to perform scaling within arrays in advance. (default = TRUE)
#' @param q percentile between 0 and 1 specifying the quantile to align different tfs within a rep.
#'        See \code{lowertailNormalization}. (default = 0.5)
#' @param .force logical whether to run normalization even if data
#'        has already been normalized within arrays. (default = FALSE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to normalization. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @param plot_weight logical whether to plot \code{qsmooth::qsmoothPlotWeights}. (default = TRUE)
#' @param ... parameters to pass to \code{qsmooth::qsmooth}.
#'        See details below for more information on main parameters.
#'
#' @return
#' list of SummarizedExperiment object with normalized intensities in
#' new assay
#' 
#' @import SummarizedExperiment
#' @importFrom qsmooth qsmooth
#' @export
#' @author Dongyuan Song, Patrick Kimes 
qsmoothNormalization <- function(se_list, 
                                 ref_rep = 1,
                                 assay_name = "fore",
                                 scaling_between_rep = TRUE,
                                 scaling_within_rep = TRUE,
                                 q = 0.5,
                                 .force = FALSE, .filter = 1L,
                                 plot_weight = TRUE, 
                                 ...) {
  ## Check num of replicates
  rep_num <- length(se_list)
  stopifnot(rep_num >= 2)
  
  ## check if already normalized
  if (!.force) {
    sapply(se_list, 
           function(x) {stopifnot(is.null(metadata(x)$betweenArrayNormalization))}, simplify = TRUE)
  }
  
  
  ## perform median scaling within reps
  if (scaling_within_rep) {
    new_assay_list <- lapply(se_list, lowertailNormalization, 
                             assay_name = assay_name,  q = q,  shift = FALSE, method = "quantile", .filter = .filter)
    new_assay_list <- lapply(new_assay_list, function(x) {as.matrix(assay(x, "scaled"))})
  }
  else {
    new_assay_list <- lapply(se_list, function(x) {as.matrix(assay(x, "scaled"))})
  }
  
  ## perform median scaling between reps (with geometric mean of ratios)
  if (scaling_between_rep) {
    median_ref <- matrixStats::colMedians(new_assay_list[[ref_rep]], na.rm = TRUE)
    
    new_assay_list <- lapply(new_assay_list, function(x) {
      x * exp(mean(log(median_ref/matrixStats::colMedians(x, na.rm = TRUE))))})
    
  }
  
  ## bind list of assays
  new_assay_list <- lapply(new_assay_list, as.data.frame)
  new_assay <- dplyr::bind_cols(new_assay_list)
  
  
  row.names(new_assay) <- rowData(se_list[[ref_rep]])$Sequence
  new_assay_cmpl <- new_assay[complete.cases(new_assay), ]
  
  ## deal with missing value
  
  warning(paste0("Notice: because of missing value ",  round((1 - dim(new_assay_cmpl)[1]/dim(new_assay)[1])*100, 2), 
                 "% rows will be removed."))
  
  ## variant numbers
  var_num <- dim(new_assay)[2]/length(se_list)
  
  ## perform quantile normalization
  q_assay <- qsmooth::qsmooth(object = new_assay[complete.cases(new_assay), ], 
                              groupFactor = rep(1:var_num, rep_num),
                              ...)
  if(plot_weight) qsmooth::qsmoothPlotWeights(q_assay)
  
  q_assay <- as.matrix(q_assay@qsmoothData)
  
  row.names(q_assay) <- row.names(new_assay_cmpl)
  
  ## match back to the original new_assay
  new_assay[row.names(new_assay) %in% row.names(q_assay), ] <- q_assay[ ,]  
  new_assay[!complete.cases(new_assay), ] <- NA
  
  
  new_assay <- DataFrame(new_assay)
  names(new_assay) <- rep(rownames(colData(se_list[[ref_rep]])), length(se_list))
  row.names(new_assay) <- NULL
  
  ## modify input SummarizedExperiment (need modified)
  
  for(i in 1:length(se_list)) {
    temp_df <- new_assay[, 1:var_num + (i-1)*var_num]
    dimnames(temp_df) <- dimnames(assay(se_list[[i]], assay_name))
    assay(se_list[[i]], assay_name) <- temp_df
  }
  
  
  ## add step to metadata
  method_str <- paste("qsmooth::normalizeBetweenArrays ->", assay_name)
  se_list <- lapply(se_list, function(x) {metadata(x)$steps <- c(metadata(x)$steps, method_str)
  return(x)})
  
  if (.force) {
    se_list <- lapply(se_list, function(x) {metadata(x)$betweenArrayNormalization <-
      c(metadata(x)$betweenArrayNormalization, method_str) 
    return(x)})
  } else {
    se_list <- lapply(se_list, function(x) {metadata(x)$betweenArrayNormalization <-
      method_str 
    return(x)})
  }
  
  return(se_list)
}