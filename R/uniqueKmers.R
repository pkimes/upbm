#' Generate Unique K-mer Sequences
#'
#' This simple function returns a vector of all unique k-mers of
#' a specified length, filtered to remove reverse complement redundancy
#' such that, e.g. only one of "AA" and "TT" is kept. When deciding
#' between a k-mer and the reverse complement, the "smaller" sequence
#' is kept, where simple alphabetical ordering is used. In the example
#' above, "AA" will be kept because "AA" <= "TT". For \code{k = 8}, 
#' this will return the same set of 8-mers as is used in the Universal
#' PBM Analysis Suite software.
#'
#' Note, this function is not meant to be used for k > 10, and will return
#' an error if k > 10 is specified unless, \code{largek = TRUE}.
#' 
#' @param k integer length of oligonucleotide sequences to be
#'        returned. (default = 8)
#' @param .largek logical whether to allow specifying a large k,
#'        i.e. k > 10L. (default = FALSE)
#' 
#' @return
#' vector of unique k-mer strings.
#' 
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
