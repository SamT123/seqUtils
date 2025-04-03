
#' Read a fasta file quite quickly using `readr::read_lines()`
#'
#' @param path path of .fasta to read
#'
#' @return Named character vector of sequences
#'
#' @importFrom readr read_lines
#' @export
fast_fasta = function(path){
  lines = readr::read_lines(path, skip_empty_rows = F)
  name_lines = which(substr(lines, 1, 1) == '>')
  message("Loaded ", length(name_lines), " sequences")

  if (length(name_lines) == length(lines)/2 && all(name_lines == seq(1, length(lines), 2))){
    # case: single-line sequences
    sequences = lines[-name_lines]
    names(sequences) = lines[name_lines]
  } else {
    # case: multi-line sequences
    sequences = rep(NA, length(name_lines))
    names(sequences) = rep(NA, length(name_lines))

    name_lines = c(name_lines, length(lines)+1)

    for (i in seq_along(name_lines[-1])){
      names(sequences)[[i]] = lines[[name_lines[[i]]]]
      sequences[[i]] = paste(
        lines[(name_lines[[i]]+1):(name_lines[[i+1]]-1)],
        collapse = ""
      )
    }
  }

  names(sequences) = substring(names(sequences), 2, 1e6)
  sequences
}


#' Write a fasta file quite quickly using `readr::write_lines()`
#'
#' @param seqs a named or unnamed (in which case `names` must be provided) character vector of sequences
#' @param names (optional) names of sequences, required if `sequences` is unnamed.
#' @param path path to write .fasta
#'
#' @return NULL
#'
#' @importFrom readr write_lines
#' @export
write_fast_fasta = function(seqs, names = NULL, path){
  if (is.null(names)) names = names(seqs)

  stopifnot(length(seqs) == length(names))

  x <- vector(class(seqs), 2*length(seqs))
  x[c(TRUE, FALSE)] <- paste0('>',names)
  x[c(FALSE, TRUE)] <- seqs

  readr::write_lines(x, file = path)

  invisible(NULL)
}
