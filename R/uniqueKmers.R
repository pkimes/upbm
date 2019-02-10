#' @title Generate Unique K-mer Sequences
#'
#' @description
#' This simple function returns a vector of all unique k-mers of
#' a specified length, filtered to remove reverse complement redundancy
#' such that, e.g. only one of "AA" and "TT" is kept. When deciding
#' between a k-mer and the reverse complement, the alphabetically earlier
#' sequence is kept. In the example above, "AA" will be kept because "AA"
#' occurs before "TT".
#' 
#' Note, this function is not meant to be used for \code{k > 10} as the
#' computational cost of the approach can quickly grow to be too much.
#' By default, an error willbe thrown if k > 10 is specified. This 
#' can be overriden by setting \code{largek = TRUE}.
#' 
#' @param k an integer length of oligonucleotide sequences to be
#'        returned. (default = 8L)
#' @param .largek a logical value whether to allow specifying a large \code{k},
#'        i.e. \code{k > 10L}. (default = FALSE)
#' 
#' @return
#' Vector of unique k-mer strings.
#' 
#' @export
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet reverseComplement
#' @importFrom dplyr tibble mutate filter
#' @author Patrick Kimes
uniqueKmers <- function(k = 8L, .largek = FALSE) {
    if (!.largek & k > 10L) 
        stop("Specified 'k' too large, specify a 'k' no greater than 10.")
    kmers <- Biostrings::oligonucleotideFrequency(Biostrings::DNAStringSet(""), k)
    kmers <- dplyr::tibble(fwd = colnames(kmers))
    kmers <- dplyr::mutate(kmers, rev = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(fwd))))
    kmers <- dplyr::filter(kmers, fwd <= rev)
    kmers$fwd
}
