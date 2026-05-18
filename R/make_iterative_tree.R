#' Build a phylogenetic tree iteratively with IQ-TREE bootstrap + CMAPLE growth
#'
#' @description
#' Builds a high-quality phylogeny by bootstrapping a small backbone with
#' IQ-TREE, then growing it through a cascade of CMAPLE runs at increasing
#' per-year sample sizes, each freezing the previous tree. The final CMAPLE
#' run uses every row of `alignment`.
#'
#' The cascade goes:
#' \enumerate{
#'   \item IQ-TREE on `initial_iqtree_size` sequences/year.
#'   \item For each `n` in `cascade_sizes`: CMAPLE on `n`/year, with
#'         `starting_tree_path` = previous step's tree and
#'         `freeze_starting_tree = TRUE`. Every tip from the previous tree is
#'         forced into the new subsample.
#'   \item A final CMAPLE on the entirety of `alignment`, again freezing the
#'         previous tree.
#' }
#' Intermediate trees are written under `<work_dir>/cascade/`.
#
#' @param alignment A `data.frame`/tibble with columns
#'   `Isolate_unique_identifier` (character), `dna_sequence` (character), and
#'   `Collection_date` (Date or coercible). The cascade samples from this and
#'   the final step uses it whole.
#' @param tree_path Character path for the final tree (Newick).
#' @param initial_iqtree_size Integer. Sequences per year for the IQ-TREE
#'   bootstrap step.
#' @param cascade_sizes Integer vector of sequences-per-year for the CMAPLE
#'   cascade steps. Must be strictly increasing.
#' @param work_dir Optional character. Directory for intermediate trees. If
#'   `NULL` (default), uses `fs::path_dir(tree_path)`.
#' @param unique_only Logical. If `TRUE` (default), duplicate `dna_sequence`
#'   entries are deduplicated before subsampling (forced-through tips
#'   excepted).
#' @param seed Integer. RNG seed for subsampling and for IQ-TREE / CMAPLE.
#'   Default 100.
#' @param num_threads Threads to pass to IQ-TREE and CMAPLE. Default
#'   `"AUTO"`.
#' @param model Substitution model for both tools. Default `"GTR"`.
#' @param return_tree Logical. If `TRUE` (default), returns the final phylo;
#'   else returns `tree_path`.
#'
#' @return A `phylo` object (default) or a character path.
#' @export
#'
#' @seealso [make_iqtree_tree()], [make_cmaple_tree()]
#'
#' @examples
#' \dontrun{
#' aln <- readRDS("aln.RDS") # tibble with the required columns
#' tree <- make_iterative_tree(
#'   alignment           = aln,
#'   tree_path           = "tree.nwk",
#'   initial_iqtree_size = 10,
#'   cascade_sizes       = c(20, 40, 80, 160, 320, 640, 1280, 2560),
#'   work_dir            = "tree"
#' )
#' }
make_iterative_tree = function(
  alignment,
  tree_path,
  initial_iqtree_size,
  cascade_sizes,
  work_dir = NULL,
  unique_only = TRUE,
  seed = 100,
  num_threads = "AUTO",
  model = "GTR",
  return_tree = TRUE
) {
  required_cols = c(
    "Isolate_unique_identifier",
    "dna_sequence",
    "Collection_date"
  )
  missing_cols = setdiff(required_cols, colnames(alignment))
  if (length(missing_cols) > 0) {
    stop(
      "alignment is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!is.numeric(initial_iqtree_size) || initial_iqtree_size < 1) {
    stop("initial_iqtree_size must be a positive integer")
  }
  if (length(cascade_sizes) == 0) {
    stop("cascade_sizes must contain at least one value")
  }
  if (is.unsorted(cascade_sizes, strictly = TRUE)) {
    stop("cascade_sizes must be strictly increasing")
  }
  if (cascade_sizes[1] <= initial_iqtree_size) {
    stop(
      "cascade_sizes must all be greater than initial_iqtree_size (",
      initial_iqtree_size,
      ")"
    )
  }

  if (is.null(work_dir)) {
    work_dir = fs::path_dir(tree_path)
  }
  cascade_dir = fs::path(work_dir, "cascade")
  fs::dir_create(cascade_dir, recurse = TRUE)

  alignment$year = lubridate::year(alignment$Collection_date)

  to_named_sequences = function(aln_subset) {
    stats::setNames(
      aln_subset$dna_sequence,
      aln_subset$Isolate_unique_identifier
    )
  }

  message(
    "[make_iterative_tree] IQ-TREE bootstrap on ",
    initial_iqtree_size,
    " sequences/year"
  )
  init_aln = subsamplePerYear(
    alignment,
    n_per_year = initial_iqtree_size,
    forced_identifiers = character(0),
    unique_only = unique_only,
    seed = seed
  )
  init_tree_path = fs::path(cascade_dir, "iqtree_initial", ext = "nwk")
  make_iqtree_tree(
    sequences = to_named_sequences(init_aln),
    tree_path = init_tree_path,
    num_threads = num_threads,
    model = model,
    seed = seed,
    overwrite = TRUE,
    keep_files = c("nwk", "log"),
    return_tree = FALSE
  )

  prev_tree_path = init_tree_path
  prev_tips = init_aln$Isolate_unique_identifier

  for (n_per_year in cascade_sizes) {
    message(
      "[make_iterative_tree] CMAPLE cascade step: ",
      n_per_year,
      " sequences/year"
    )
    step_aln = subsamplePerYear(
      alignment,
      n_per_year = n_per_year,
      forced_identifiers = prev_tips,
      unique_only = unique_only,
      seed = seed
    )
    step_tree_path = fs::path(
      cascade_dir,
      paste0("size_", n_per_year),
      ext = "nwk"
    )
    make_cmaple_tree(
      sequences = to_named_sequences(step_aln),
      tree_path = step_tree_path,
      multifurcating = TRUE,
      starting_tree_path = prev_tree_path,
      freeze_starting_tree = TRUE,
      num_threads = num_threads,
      model = model,
      seed = seed,
      overwrite = TRUE,
      keep_files = c("nwk", "log"),
      return_tree = FALSE
    )
    prev_tree_path = step_tree_path
    prev_tips = step_aln$Isolate_unique_identifier
  }

  message(
    "[make_iterative_tree] Final CMAPLE on full alignment (",
    nrow(alignment),
    " sequences)"
  )
  final_tree = make_cmaple_tree(
    sequences = to_named_sequences(alignment),
    tree_path = tree_path,
    multifurcating = TRUE,
    starting_tree_path = prev_tree_path,
    freeze_starting_tree = TRUE,
    num_threads = num_threads,
    model = model,
    seed = seed,
    overwrite = TRUE,
    keep_files = c("nwk", "log", "fasta"),
    return_tree = return_tree
  )

  final_tree
}


# Internal: stratified subsample by year, force-include specific isolates.
# Mirrors the inspiration code at predicting_influenza_using_convergence_ms_public
# /make_trees/R/make_CMAPLE_tree.R::ss(). Caller must add the `year` column.
subsamplePerYear = function(
  alignment,
  n_per_year,
  forced_identifiers = character(0),
  unique_only = TRUE,
  seed = NULL
) {
  stopifnot("year" %in% colnames(alignment))

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  # Put forced rows first so de-duplication keeps them.
  alignment = alignment[
    order(!alignment$Isolate_unique_identifier %in% forced_identifiers),
  ]

  if (unique_only) {
    keep = !duplicated(alignment$dna_sequence) |
      alignment$Isolate_unique_identifier %in% forced_identifiers
    alignment = alignment[keep, ]
  }

  picked = unlist(
    lapply(split(seq_len(nrow(alignment)), alignment$year), function(idx) {
      idx[sample.int(length(idx), size = min(length(idx), n_per_year))]
    }),
    use.names = FALSE
  )

  is_forced = alignment$Isolate_unique_identifier %in% forced_identifiers
  keep_rows = sort(union(picked, which(is_forced)))
  alignment[keep_rows, , drop = FALSE]
}
