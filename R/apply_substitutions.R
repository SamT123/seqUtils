#' Apply substitutions to a sequence
#'
#' @param sequence Character string of amino acids or nucleotides
#' @param substitutions Character vector of substitutions in format "K2N" (change K at position 2 to N)
#' @param type Sequence type: "aa" for amino acids or "nt" for nucleotides
#'
#' @return Modified sequence with substitutions applied
#'
#' @export
apply_substitutions = function(sequence, substitutions, type = "aa") {
  alphabet = switch(
    type,
    "aa" = c(Biostrings::AA_STANDARD, "X", "-"),
    "nt" = c(Biostrings::DNA_BASES, "N", "-"),
    stop("type must be 'aa' or 'nt'")
  )

  new_sequence = sequence
  for (s in substitutions) {
    fat = convert_substitution(s, alphabet)
    old_char = stringr::str_sub(new_sequence, fat$at, fat$at)
    if (!old_char == fat$from) {
      # fmt: skip
      stop("sub = ", s, "; from_char mismatch, expected: ", fat$from, "; got ", old_char)
    }

    stringr::str_sub(new_sequence, fat$at, fat$at) = fat$to
  }

  new_sequence
}

convert_substitution = function(s, alphabet) {
  invalid = function(s) {
    stop("Invalid substitution? ", s)
  }
  from = stringr::str_sub(s, 1, 1)
  to = stringr::str_sub(s, -1, -1)
  at = tryCatch(
    as.integer(stringr::str_sub(s, 2, -2)),
    warning = function(msg) invalid(s)
  )

  if (!all(c(from, to) %in% alphabet) | is.na(at)) {
    invalid(s)
  }

  list(from = from, at = at, to = to)
}
