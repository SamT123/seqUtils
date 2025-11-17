# Helper function to add a path to the system PATH
add_to_PATH = function(path) {
  current_path = Sys.getenv("PATH")
  new_path = paste(path, current_path, sep = .Platform$path.sep)
  Sys.setenv(PATH = new_path)
}

#' Build a phylogenetic tree using CMAPLE
#'
#' @description
#' Constructs a maximum likelihood phylogenetic tree from DNA/RNA sequences using
#' CMAPLE (Comprehensive Maximum Likelihood Phylogeny Estimation). The tree is
#' optionally ladderized and rooted before being returned.
#'
#' @details
#' CMAPLE is a fast maximum likelihood phylogenetic inference tool. This function
#' provides a convenient R wrapper that handles input/output files, parameter
#' configuration, and tree post-processing (ladderizing and rooting).
#'
#' If \code{freeze_starting_tree = TRUE}, CMAPLE uses NORMAL search mode which
#' preserves relationships in the starting tree. If \code{FALSE}, EXHAUSTIVE
#' search mode is used which may rearrange the starting tree topology.
#'
#' @param sequences Named character vector of DNA or RNA sequences. All sequences
#'   must have unique names.
#' @param tree_path Character string specifying the file path where the tree will
#'   be saved. A log file will also be created at \code{paste0(tree_path, ".log")}.
#' @param multifurcating Logical. If \code{TRUE} (default), the output tree may
#'   contain multifurcations (polytomies). If \code{FALSE}, the tree will be
#'   strictly bifurcating.
#' @param starting_tree_path Optional character string specifying a path to a
#'   starting tree file (Newick format) containing a subset of the sequences.
#'   CMAPLE will add the remaining sequences to this tree. Default is \code{NULL}.
#' @param freeze_starting_tree Logical. If \code{TRUE}, relationships in the
#'   starting tree are preserved (uses NORMAL search). If \code{FALSE} (default),
#'   CMAPLE may rearrange the topology (uses EXHAUSTIVE search). Only relevant
#'   when \code{starting_tree_path} is provided.
#' @param cmaple_path Optional character string specifying the directory containing
#'   the CMAPLE executable. If provided, this path will be added to the system PATH.
#'   If \code{NULL} (default), assumes CMAPLE is already in PATH.
#' @param out_sequence Optional character string specifying the name of a sequence
#'   to use as the outgroup for rooting the tree. Must match one of the sequence
#'   names. If \code{NULL} (default), the tree is ladderized but not explicitly rooted.
#' @param num_threads Character string or numeric specifying the number of threads
#'   for CMAPLE to use. Default is \code{"AUTO"} which lets CMAPLE choose automatically.
#' @param model Character string specifying the substitution model. Default is
#'   \code{"GTR"} (General Time Reversible).
#' @param keep_files Character vector specifying which files to keep after CMAPLE
#'   execution. Valid options are \code{"nwk"} (tree file), \code{"fasta"}
#'   (alignment file), and \code{"log"} (log file). Default is \code{c("nwk", "log")}.
#'   Note: \code{"nwk"} can only be omitted when \code{return_tree = TRUE}.
#' @param return_tree Logical. If \code{TRUE} (default), returns the tree object
#'   (phylo class). If \code{FALSE}, returns the file path as a character string.
#'
#' @return
#' If \code{return_tree = TRUE}, returns a phylogenetic tree object of class
#' \code{phylo} (from the \code{ape} package). If \code{return_tree = FALSE},
#' returns a character string containing the path to the saved tree file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' sequences <- c(seq1 = "ATCG", seq2 = "ATCC", seq3 = "TTCG")
#' tree <- make_cmaple_tree(sequences, tree_path = "my_tree.nwk")
#'
#' # With a starting tree and custom parameters
#' tree <- make_cmaple_tree(
#'   sequences,
#'   tree_path = "my_tree.nwk",
#'   starting_tree_path = "starting.nwk",
#'   freeze_starting_tree = TRUE,
#'   out_sequence = "seq3",
#'   num_threads = 4,
#'   model = "GTR"
#' )
#'
#' # Keep only the tree file (delete fasta and log)
#' tree <- make_cmaple_tree(
#'   sequences,
#'   tree_path = "my_tree.nwk",
#'   keep_files = "nwk"
#' )
#'
#' # Delete all files, only return tree object
#' tree <- make_cmaple_tree(
#'   sequences,
#'   tree_path = "my_tree.nwk",
#'   keep_files = character(0),
#'   return_tree = TRUE
#' )
#' }
make_cmaple_tree = function(
  sequences,
  tree_path,
  multifurcating = TRUE,
  starting_tree_path = NULL,
  freeze_starting_tree = FALSE,
  cmaple_path = NULL,
  out_sequence = NULL,
  num_threads = "AUTO",
  model = "GTR",
  keep_files = c("nwk", "log"),
  return_tree = TRUE
) {
  # Validate inputs
  if (!is.character(sequences)) {
    stop("sequences must be a character vector")
  }
  if (is.null(names(sequences)) || any(names(sequences) == "")) {
    stop("sequences must be a named vector with all sequences having names")
  }

  # Check that tree_path directory exists, create if needed
  tree_dir = fs::path_dir(tree_path)
  if (!fs::dir_exists(tree_dir)) {
    message("Creating directory: ", tree_dir)
    fs::dir_create(tree_dir, recurse = TRUE)
  }

  # Validate keep_files
  if (!is.null(keep_files)) {
    valid_files = c("nwk", "fasta", "log")
    invalid_files = setdiff(keep_files, valid_files)
    if (length(invalid_files) > 0) {
      stop(
        "Invalid keep_files values: ",
        paste(invalid_files, collapse = ", "),
        ". Valid options are: ",
        paste(valid_files, collapse = ", ")
      )
    }
  }

  # Check that "nwk" is in keep_files if return_tree = FALSE
  if (!return_tree && !("nwk" %in% keep_files)) {
    stop("keep_files must include 'nwk' when return_tree = FALSE")
  }

  if (!is.null(cmaple_path)) {
    add_to_PATH(cmaple_path)
  }

  # Check if CMAPLE is available
  if (system("command cmaple", ignore.stdout = T) != 0) {
    stop(
      "CMAPLE is not found in PATH. Please install CMAPLE or provide cmaple_path parameter."
    )
  }

  if (!is.null(starting_tree_path)) {
    if (!file.exists(starting_tree_path)) {
      stop("starting_tree_path does not exist: ", starting_tree_path)
    }
    starting_tree = castor::read_tree(file = starting_tree_path)
    message(
      sum(names(sequences) %in% starting_tree$tip.label),
      " / ",
      length(sequences),
      " sequences already in starting tree"
    )

    message(
      sum(starting_tree$tip.label %in% names(sequences)),
      " / ",
      length(starting_tree$tip.label),
      " starting tree tips in sequence list"
    )
  }

  fasta_path = fs::path_ext_set(tree_path, ".fasta")

  seqUtils::write_fast_fasta(
    sequences,
    names(sequences),
    path = fasta_path
  )

  # fmt: skip
  cmaple_call = c(
    "cmaple",
    "-aln", shQuote(fasta_path)
  )

  if (multifurcating) {
    # fmt: skip
    cmaple_call = c(
      cmaple_call,
      "--out-mul-tree"
    )
  }

  if (!is.null(starting_tree_path)) {
    # fmt: skip
    cmaple_call = c(
      cmaple_call,
      "-t", shQuote(starting_tree_path)
    )
  }

  if (!freeze_starting_tree) {
    # fmt: skip
    cmaple_call = c(
      cmaple_call,
      "-search", "EXHAUSTIVE"
    )
  } else {
    # fmt: skip
    cmaple_call = c(
      cmaple_call,
      "-search", "NORMAL"
    )
  }

  # fmt: skip
  cmaple_call = c(
      cmaple_call,
      "-nt", num_threads,
      "-m", model
    )

  exit_code = system(paste(cmaple_call, collapse = " "))
  if (exit_code != 0) {
    stop("CMAPLE execution failed with exit code ", exit_code)
  }

  # Check if output files were created
  treefile_path = paste0(fasta_path, ".treefile")
  if (!file.exists(treefile_path)) {
    stop("CMAPLE did not create expected output file: ", treefile_path)
  }

  # Move output files to desired locations
  exit_code = system(paste("mv", shQuote(treefile_path), shQuote(tree_path)))
  if (exit_code != 0) {
    stop("Failed to move tree file to ", tree_path)
  }

  logfile_path = paste0(fasta_path, ".log")
  if (file.exists(logfile_path)) {
    system(paste(
      "mv",
      shQuote(logfile_path),
      shQuote(paste0(tree_path, ".log"))
    ))
  }

  tree = castor::read_tree(file = tree_path)
  tree = ladderizeAndMaybeRoot(tree, out_sequence)
  castor::write_tree(tree, tree_path)

  # Clean up files based on keep_files
  if (!("fasta" %in% keep_files)) {
    if (fs::file_exists(fasta_path)) {
      fs::file_delete(fasta_path)
    }
  }

  if (!("log" %in% keep_files)) {
    log_path = paste0(tree_path, ".log")
    if (fs::file_exists(log_path)) {
      fs::file_delete(log_path)
    }
  }

  if (!("nwk" %in% keep_files) && return_tree) {
    if (fs::file_exists(tree_path)) {
      fs::file_delete(tree_path)
    }
  }

  if (return_tree) {
    return(tree)
  } else {
    return(tree_path)
  }
}


ladderizeAndMaybeRoot = function(tree, out_sequence = NULL) {
  if (!is.null(out_sequence)) {
    stopifnot(out_sequence %in% tree$tip.label)

    tree = ape::root(ape::unroot(tree), outgroup = out_sequence)
  }

  tree = ape::ladderize(tree)

  f = fs::file_temp()
  castor::write_tree(tree, file = f)
  tree_read = castor::read_tree(file = f)

  tree_read
}
