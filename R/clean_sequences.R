#' Remove non-standard characters from sequences
#'
#' @param sequences A character vector of sequences
#' @param type "aa" or "nt"
#' @param alphabet (optional) Provide a custom alphabet (defaults for "aa" and "nt" are the alphabet or "-")
#' @param replacement_character (optional) Provide a custom replacement character (defaults for "aa" and "nt" are "X" and "N"; you may want "-")
#'
#' @return a character vector of cleaned sequences
#' @export
#'
#' @importFrom stringr str_replace_all
#' @importFrom Biostrings AA_STANDARD
clean_sequences = function(
  sequences,
  type,
  alphabet = NULL,
  replacement_character = NULL
) {
  sequences = toupper(sequences)

  if (is.null(replacement_character)) {
    replacement_character = c(aa = "X", nt = "N")[[type]]
  }

  if (is.null(alphabet)) {
    alphabet = list(
      aa = c(Biostrings::AA_STANDARD, "-"),
      nt = c("A", "T", "C", "G", "-")
    )[[type]]
  }

  regex = paste0(
    "[^",
    paste(alphabet, collapse = ""),
    "]"
  )

  sequences = setNames(
    str_replace_all(
      sequences,
      regex,
      replacement_character
    ),
    names(sequences)
  )

  sequences
}
