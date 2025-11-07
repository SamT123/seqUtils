#' Get the consensus sequence from a vector of sequences
#'
#' P
#'
#' @param sequences a vector of sequences
#' @param excluded_characters characters to disregard for calculating frequencies; perhaps "-", "X", "N"
#' @param min_freq positions where there is no aa/nt at frequency >= 'min_freq' will be "?"
#'
#' @export
getConsensus = function(sequences, excluded_characters = c(), min_freq = 0.5) {
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
  aas[highest_frequencies < min_freq] = "?"
  paste0(aas, collapse = "")
}
