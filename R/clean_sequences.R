#' Remove non-standard characters from sequences
#'
#' Replaces any characters not in the specified alphabet with a replacement character.
#' All sequences are converted to uppercase before processing.
#'
#' @param sequences A character vector of sequences (names are preserved if present)
#' @param type Either "aa" (amino acid) or "nt" (nucleotide)
#' @param alphabet (optional) Character vector defining valid characters. Defaults:
#'   \itemize{
#'     \item For type "aa": The 20 standard amino acids (ARNDCQEGHILKMFPSTWYV) plus "-"
#'     \item For type "nt": A, T, C, G, and "-"
#'   }
#' @param replacement_character (optional) Character to replace non-standard characters with.
#'   Defaults to "X" for amino acids and "N" for nucleotides. Use "-" to replace with gaps.
#'
#' @return A character vector of cleaned sequences with names preserved
#'
#' @examples
#' # Clean nucleotide sequences
#' clean_sequences(c("ATCGXYZ", "GGGNNN"), type = "nt")  # X,Y,Z become N
#'
#' # Clean amino acid sequences
#' clean_sequences(c("ACDEF123"), type = "aa")  # 1,2,3 become X
#'
#' # Use custom replacement
#' clean_sequences("ATC-G", type = "nt", replacement_character = "-")
#'
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
  # Input validation
  if (!type %in% c("aa", "nt")) {
    stop("type must be either 'aa' or 'nt', got '", type, "'")
  }

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
