#' Calculate consensus sequence from multiple sequences
#'
#' Returns the most frequent character at each position. Sequences are converted to
#' uppercase before processing. Positions without a character meeting the frequency
#' threshold are marked with "?".
#'
#' @param sequences Character vector of aligned sequences
#' @param excluded_characters Characters to ignore when calculating frequencies
#'   (e.g., "-", "X", "N")
#' @param min_freq Minimum frequency (0-1) required for a character to be used.
#'   Positions below this threshold become "?" (default: 0.5)
#'
#' @return Single consensus sequence as a character string
#'
#' @examples
#' # Basic consensus
#' seqs <- c("ACGT", "ACGT", "ACAT")
#' get_consensus(seqs)  # Returns "ACGT" (G appears in 2/3 sequences)
#'
#' # With frequency threshold
#' seqs <- c("ACGT", "ACAT", "ATTT")
#' get_consensus(seqs, min_freq = 0.7)  # Returns "A??T"
#'
#' # Exclude gaps from calculation
#' seqs <- c("AC-T", "ACGT", "AC-T")
#' get_consensus(seqs, excluded_characters = "-")
#'
#' @importFrom Biostrings consensusMatrix
#' @export
get_consensus = function(sequences, excluded_characters = c(), min_freq = 0.5) {
  sequences = toupper(sequences)
  cm = Biostrings::consensusMatrix(sequences)
  cm = cm[!rownames(cm) %in% excluded_characters, ]

  get_aa = function(x) {
    names(x)[which.max(x)]
  }

  get_freq = function(x) {
    max(x) / sum(x)
  }

  highest_frequencies = setNames(
    apply(cm, 2, get_freq),
    apply(cm, 2, get_aa)
  )

  aas = names(highest_frequencies)
  aas[!is.finite(highest_frequencies) | highest_frequencies < min_freq] = "?"
  paste0(aas, collapse = "")
}
