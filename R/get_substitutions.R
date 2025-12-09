#' Get substitutions between reference and query sequences
#'
#' Compares sequences and returns substitutions in the format "A123T" (reference amino acid,
#' position, query amino acid). Sequences are trimmed to the shortest length if they differ.
#'
#' @param sequence_1 Reference sequence (single sequence)
#' @param sequence_2 Query sequence(s) to compare (character vector, can be multiple). If named,
#'   names will be preserved in the output for multi-sequence cases.
#' @param position_map Custom position numbering (default: 1:nchar(sequence_1)). Useful when
#'   sequences represent a subset of positions from a larger sequence.
#' @param exclude Characters to ignore when comparing (e.g., "X" for unknown amino acids,
#'   "N" for ambiguous nucleotides)
#' @param simplify If TRUE and only one query sequence provided, returns a character vector
#'   instead of a list (default: TRUE)
#' @param allow_empty_sequence_2 If TRUE, allows sequence_2 to be empty (length 0),
#'   returning an empty list. If FALSE (default), throws an error for empty sequence_2.
#'
#' @return List of character vectors (one per query sequence) containing substitutions with
#'   names from sequence_2 (if present), or a single character vector if simplify = TRUE and
#'   only one query sequence provided. Single-sequence output is always unnamed.
#'
#' @examples
#' # Single comparison
#' get_substitutions("ACDEF", "ACDEG")  # Returns "F5G"
#'
#' # Multiple sequences
#' get_substitutions("ACDEF", c("ACDEG", "TCDEG"))
#'
#' # Custom position numbering
#' get_substitutions("ACE", "TCE", position_map = c(145, 155, 156))
#'
#' # Exclude ambiguous characters
#' get_substitutions("ACDEF", "XCDXF", exclude = "X")
#'
#' @importFrom Biostrings AAStringSet subseq
#' @importFrom purrr map2 pmap
#' @importFrom stringr str_split
#' @export
get_substitutions = function(
  sequence_1,
  sequence_2,
  position_map = 1:nchar(sequence_1),
  exclude = c(),
  simplify = TRUE,
  allow_empty_sequence_2 = FALSE
) {
  # Input validation
  if (length(sequence_1) != 1) {
    stop("sequence_1 must be a single sequence, got ", length(sequence_1))
  }
  if (length(sequence_2) < 1) {
    if (!allow_empty_sequence_2) {
      stop(
        "sequence_2 must contain at least one sequence unless `allow_empty_sequence_2` is set to TRUE"
      )
    } else {
      return(list())
    }
  }

  lengths = nchar(c(sequence_1, sequence_2))

  if (!all(lengths == min(lengths))) {
    message('Trimming to shortest length: ', min(lengths))
    sequence_1 = substr(sequence_1, 1, min(lengths))
    sequence_2 = substr(sequence_2, 1, min(lengths))
  }

  if (length(sequence_1) == 1 & length(sequence_2) == 1 & simplify) {
    return(get_substitutions_CASE_SINGLE(
      sequence_1,
      sequence_2,
      position_map,
      exclude
    ))
  }

  sequence_2_uniques = unique(sequence_2)
  idxs = match(sequence_2, sequence_2_uniques)

  # Store original names from sequence_2
  sequence_2_names <- names(sequence_2)
  has_names <- !is.null(sequence_2_names) && !all(is.na(sequence_2_names))

  unique_substitutions = get_substitutions_CASE_MULTIPLE(
    sequence_1,
    sequence_2_uniques,
    position_map,
    exclude
  )

  all_substitutions = unique_substitutions[idxs]

  # Preserve names from sequence_2 in output (only for multi-sequence case)
  if (has_names) {
    names(all_substitutions) <- sequence_2_names
  }

  all_substitutions
}


get_substitutions_CASE_SINGLE = function(
  sequence_1,
  sequence_2,
  position_map = 1:max(nchar(sequence_1), nchar(sequence_2)),
  exclude = c()
) {
  min_len = min(
    stringr::str_length(sequence_1),
    stringr::str_length(sequence_2)
  )
  if (stringr::str_length(sequence_1) != stringr::str_length(sequence_2)) {
    stop('Sequence lengths differ')
  }

  # Handle empty sequences
  if (min_len == 0) {
    return(character(0))
  }

  sequence_1_split <- stringr::str_split(sequence_1, '')[[1]][1:min_len]
  sequence_2_split <- stringr::str_split(sequence_2, '')[[1]][1:min_len]

  diffs = unlist(
    purrr::pmap(
      list(
        s1 = sequence_1_split,
        s2 = sequence_2_split,
        n = 1:length(sequence_1_split)
      ),
      function(s1, s2, n) {
        if (s1 != s2 & !any(c(s1, s2) %in% exclude)) {
          return(paste0(s1, position_map[n], s2))
        }
      }
    )
  )

  # Convert NULL to character(0) for consistent return type
  if (is.null(diffs)) {
    diffs = character(0)
  }

  diffs
}


get_substitutions_CASE_MULTIPLE = function(
  sequence_1,
  sequence_2,
  position_map = 1:nchar(sequence_1),
  exclude = c()
) {
  # Handle empty sequences
  if (nchar(sequence_1) == 0) {
    return(rep(list(character(0)), length(sequence_2)))
  }

  diffs = rep(list(NULL), length(sequence_2))

  sequence_1_split = stringr::str_split(sequence_1, pattern = "")[[1]]
  sequence_2_split = Biostrings::AAStringSet(sequence_2)

  for (i in seq_along(sequence_1_split)) {
    sequence_2_pos_i = as.character(Biostrings::subseq(sequence_2_split, i, i))
    incl = sequence_1_split[[i]] != sequence_2_pos_i &
      (!(sequence_2_pos_i %in% exclude | sequence_1_split[[i]] %in% exclude))

    diffs[incl] = purrr::map2(
      diffs[incl],
      paste0(
        sequence_1_split[[i]],
        position_map[[i]],
        sequence_2_pos_i[incl]
      ),
      append
    )
  }

  # Convert NULL to character(0) for consistent return type
  diffs = lapply(diffs, function(x) {
    if (is.null(x)) character(0) else x
  })

  diffs
}
