#' Get substitutions which differ between two sequences
#'
#' @param sequence_1 "from" sequnce
#' @param sequence_2 "to" sequence(s); can be multiple sequences
#' @param position_map "rename" positions, if you want n'th position of sequences to be referred to with a number other than "n". e.g. if sequence_1 & sequence_2 are just the Koel-7 positions, use `position_map = c(145, 155, 156, 158, 159, 189, 193)`
#' @param exclude characters to exclude from substitution list (perhaps "X" for aas, "N" for nts)
#' @param simplify should output of length 1 (i.e. `length(sequence_2)` = 1) be unlisted?
#'
#' @return list of character vectors
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
    simplify = T
) {

  lengths = nchar(c(sequence_1, sequence_2))

  if (!all(lengths == min(lengths))){
    message('Trimming to shortest length: ', min(lengths))
    sequence_1 = substr(sequence_1, 1, min(lengths))
    sequence_2 = substr(sequence_2, 1, min(lengths))
  }

  if (length(sequence_1) == 1 & length(sequence_2) == 1){
    return(get_substitutions_CASE_SINGLE(sequence_1, sequence_2, position_map, exclude))
  }

  sequence_2_uniques = unique(sequence_2)
  idxs = match(sequence_2, sequence_2_uniques)

  unique_substitutions = get_substitutions_CASE_MULTIPLE(
    sequence_1,
    sequence_2_uniques,
    position_map,
    exclude
  )

  all_substitutions = unique_substitutions[idxs]

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
  if (stringr::str_length(sequence_1) != stringr::str_length(sequence_2) ){
    stop('Sequence lengths differ')
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
      function(s1, s2, n){
        if (s1!=s2 & !any(c(s1, s2) %in% exclude)){
          return(paste0(s1,position_map[n],s2))}
        }
    )
  )

  diffs
}


get_substitutions_CASE_MULTIPLE = function(
    sequence_1,
    sequence_2,
    position_map = 1:nchar(sequence_1),
    exclude = c()
) {
  diffs = rep(list(NULL), length(sequence_2))

  sequence_1_split = stringr::str_split(sequence_1, pattern = "")[[1]]
  sequence_2_split = Biostrings::AAStringSet(sequence_2)

  for (i in seq_along(sequence_1_split)){

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
  diffs
}
