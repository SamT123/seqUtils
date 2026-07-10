#' Build a phylogenetic tree using IQ-TREE
#'
#' @description
#' Constructs a maximum likelihood phylogenetic tree from DNA/RNA sequences
#' using IQ-TREE 2. The tree is ladderized before being returned.
#'
#' @param sequences Named character vector of DNA or RNA sequences. All
#'   sequences must have unique names.
#' @param tree_path Character string specifying the file path where the tree
#'   will be saved. A log file is also created at
#'   `paste0(tree_path, ".log")`.
#' @param iqtree_path Optional character string specifying the directory
#'   containing the IQ-TREE executable. If provided, this path is prepended to
#'   `PATH`. If `NULL` (default), assumes `iqtree2` (or `iqtree`) is already
#'   in `PATH`.
#' @param num_threads Character string or numeric specifying the number of
#'   threads. Passed to IQ-TREE via `-T`. Default `"AUTO"`.
#' @param model Character string specifying the substitution model. Default
#'   `"GTR"`.
#' @param seed Integer seed for the IQ-TREE RNG. `NULL` lets IQ-TREE choose.
#' @param overwrite Logical. If `TRUE`, passes `-redo` to IQ-TREE so it
#'   overwrites existing output. Default `FALSE`.
#' @param keep_files Character vector specifying which files to keep after
#'   IQ-TREE execution. Valid options are `"nwk"` (tree file), `"fasta"`
#'   (alignment file), and `"log"` (log file). Default `c("nwk", "log")`.
#' @param return_tree Logical. If `TRUE` (default), returns the tree as a
#'   `phylo`. If `FALSE`, returns the file path to the saved tree.
#'
#' @return A `phylo` object (default) or a character path.
#' @export
#'
#' @examples
#' \dontrun{
#' seqs <- c(seq1 = "ATCG", seq2 = "ATCC", seq3 = "TTCG")
#' tree <- make_iqtree_tree(seqs, tree_path = "iq_tree.nwk", seed = 1)
#' }
make_iqtree_tree <- function(
  sequences,
  tree_path,
  iqtree_path = NULL,
  num_threads = "AUTO",
  model = "GTR",
  seed = NULL,
  overwrite = FALSE,
  keep_files = c("nwk", "log"),
  return_tree = TRUE
) {
  if (!is.character(sequences)) {
    stop("sequences must be a character vector")
  }
  if (is.null(names(sequences)) || any(names(sequences) == "")) {
    stop("sequences must be a named vector with all sequences having names")
  }

  tree_dir <- fs::path_dir(tree_path)
  if (!fs::dir_exists(tree_dir)) {
    message("Creating directory: ", tree_dir)
    fs::dir_create(tree_dir, recurse = TRUE)
  }

  valid_files <- c("nwk", "fasta", "log")
  invalid_files <- setdiff(keep_files, valid_files)
  if (length(invalid_files) > 0) {
    stop(
      "Invalid keep_files values: ",
      paste(invalid_files, collapse = ", "),
      ". Valid options are: ",
      paste(valid_files, collapse = ", ")
    )
  }
  if (!return_tree && !("nwk" %in% keep_files)) {
    stop("keep_files must include 'nwk' when return_tree = FALSE")
  }

  if (!is.null(iqtree_path)) {
    add_to_PATH(iqtree_path)
  }

  iqtree_bin <- Sys.which("iqtree2")
  if (iqtree_bin == "") {
    iqtree_bin <- Sys.which("iqtree")
  }
  if (iqtree_bin == "") {
    stop(
      "IQ-TREE is not found in PATH. Install iqtree2 or pass iqtree_path."
    )
  }

  fasta_path <- fs::path_ext_set(tree_path, ".fasta")
  prefix <- fs::path_ext_remove(fasta_path)

  seqUtils::write_fast_fasta(
    sequences,
    names(sequences),
    path = fasta_path
  )

  # fmt: skip
  iqtree_call = c(
    shQuote(iqtree_bin),
    "-s", shQuote(fasta_path),
    "-m", model,
    "-T", num_threads,
    "--prefix", shQuote(prefix)
  )

  if (!is.null(seed)) {
    iqtree_call <- c(iqtree_call, "--seed", shQuote(seed))
  }

  if (overwrite) {
    iqtree_call <- c(iqtree_call, "-redo")
  }

  exit_code <- system(paste(iqtree_call, collapse = " "))
  if (exit_code != 0) {
    stop("IQ-TREE execution failed with exit code ", exit_code)
  }

  treefile_path <- paste0(prefix, ".treefile")
  if (!file.exists(treefile_path)) {
    stop("IQ-TREE did not create expected output file: ", treefile_path)
  }

  fs::file_move(treefile_path, tree_path)

  logfile_path <- paste0(prefix, ".log")
  if (file.exists(logfile_path)) {
    fs::file_move(logfile_path, paste0(tree_path, ".log"))
  }

  # IQ-TREE also produces .iqtree, .bionj, .mldist, .ckp.gz, .splits.nex,
  # .contree. None are downstream inputs here; remove them so the run leaves
  # only nwk / fasta / log according to keep_files.
  aux_exts <- c(
    ".iqtree",
    ".bionj",
    ".mldist",
    ".ckp.gz",
    ".splits.nex",
    ".contree",
    ".parstree",
    ".uniqueseq.phy"
  )
  unlink(paste0(prefix, aux_exts))

  tree <- castor::read_tree(file = tree_path)
  tree <- ladderizeAndMaybeRoot(tree, out_sequence = NULL)
  castor::write_tree(tree, tree_path)

  if (!("fasta" %in% keep_files)) {
    unlink(fasta_path)
  }
  if (!("log" %in% keep_files)) {
    unlink(paste0(tree_path, ".log"))
  }
  if (!("nwk" %in% keep_files)) {
    unlink(tree_path)
  }

  if (return_tree) {
    tree
  } else {
    tree_path
  }
}
