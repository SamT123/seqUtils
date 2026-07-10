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
#'   Positions below this threshold become `fill` (default: 0.5)
#' @param fill Character written at positions where no character clears
#'   `min_freq`, or where every character is excluded (default: "?")
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
#' # Plurality consensus of a nucleotide alignment, gaps where nothing is called
#' seqs <- c("AC-T", "AC-T")
#' get_consensus(seqs, excluded_characters = "-", min_freq = 0, fill = "-")
#'
#' @importFrom Biostrings consensusMatrix
#' @export
get_consensus <- function(
  sequences,
  excluded_characters = c(),
  min_freq = 0.5,
  fill = "?"
) {
  sequences <- toupper(sequences)
  cm <- Biostrings::consensusMatrix(sequences)
  cm <- cm[!rownames(cm) %in% excluded_characters, , drop = FALSE]

  # ties.method = "first" matches which.max, so a tied column resolves the same
  # way on every run.
  winners <- max.col(t(cm), ties.method = "first")
  counts <- cm[cbind(winners, seq_len(ncol(cm)))]
  highest_frequencies <- counts / colSums(cm)

  aas <- rownames(cm)[winners]
  aas[!is.finite(highest_frequencies) | highest_frequencies < min_freq] <- fill
  paste0(aas, collapse = "")
}
