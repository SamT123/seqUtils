# Helper function to check if CMAPLE is available
cmaple_available = function() {
  result = suppressWarnings(system(
    "cmaple --version",
    ignore.stdout = TRUE,
    ignore.stderr = TRUE
  ))
  return(result == 0)
}

# Input validation tests ----

test_that("make_cmaple_tree errors with unnamed sequences", {
  sequences = c("ATCGATCG", "GCTAGCTA", "TTAACCGG")
  tree_path = tempfile(fileext = ".nwk")

  expect_error(
    make_cmaple_tree(sequences, tree_path),
    "sequences must be a named vector with all sequences having names"
  )
})

test_that("make_cmaple_tree errors with partially named sequences", {
  sequences = c(seq1 = "ATCGATCG", "GCTAGCTA", seq3 = "TTAACCGG")
  tree_path = tempfile(fileext = ".nwk")

  expect_error(
    make_cmaple_tree(sequences, tree_path),
    "sequences must be a named vector with all sequences having names"
  )
})

test_that("make_cmaple_tree errors with non-character sequences", {
  sequences = c(seq1 = 1, seq2 = 2, seq3 = 3)
  tree_path = tempfile(fileext = ".nwk")

  expect_error(
    make_cmaple_tree(sequences, tree_path),
    "sequences must be a character vector"
  )
})

test_that("make_cmaple_tree errors with non-existent starting tree", {
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  sequences = c(seq1 = "ATCGATCG", seq2 = "GCTAGCTA", seq3 = "TTAACCGG")
  tree_path = tempfile(fileext = ".nwk")

  expect_error(
    make_cmaple_tree(
      sequences,
      tree_path,
      starting_tree_path = "/nonexistent/path/tree.nwk"
    ),
    "starting_tree_path does not exist"
  )
})

test_that("make_cmaple_tree errors when CMAPLE is not installed", {
  # Temporarily clear PATH to simulate CMAPLE not being installed
  old_path = Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path))
  Sys.setenv(PATH = "")

  sequences = c(seq1 = "ATCGATCG", seq2 = "GCTAGCTA", seq3 = "TTAACCGG")
  tree_path = tempfile(fileext = ".nwk")

  expect_error(
    make_cmaple_tree(sequences, tree_path),
    "CMAPLE is not found in PATH"
  )
})

# Basic functionality tests ----

test_that("make_cmaple_tree builds a tree with minimal input", {
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  sequences = c(
    seq1 = "ATCGATCGATCGATCG",
    seq2 = "ATCGATTGATCGATCG",
    seq3 = "ATCGATCGAACGATCG",
    seq4 = "TTCGATCGATCGATCG"
  )
  tree_path = tempfile(fileext = ".nwk")
  on.exit(unlink(c(tree_path, paste0(tree_path, ".log"))))

  result = make_cmaple_tree(sequences, tree_path)

  expect_s3_class(result, "phylo")
  expect_equal(length(result$tip.label), 4)
  expect_true(all(names(sequences) %in% result$tip.label))
  expect_true(file.exists(tree_path))
  expect_true(file.exists(paste0(tree_path, ".log")))
})

test_that("make_cmaple_tree returns tree object by default", {
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  sequences = c(
    seq1 = "ATCGATCGATCGATCG",
    seq2 = "ATCGATTGATCGATCG",
    seq3 = "ATCGATCGAACGATCG"
  )
  tree_path = tempfile(fileext = ".nwk")
  on.exit(unlink(c(tree_path, paste0(tree_path, ".log"))))

  result = make_cmaple_tree(sequences, tree_path)

  expect_s3_class(result, "phylo")
})

test_that("make_cmaple_tree returns path when return_tree = FALSE", {
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  sequences = c(
    seq1 = "ATCGATCGATCGATCG",
    seq2 = "ATCGATTGATCGATCG",
    seq3 = "ATCGATCGAACGATCG"
  )
  tree_path = tempfile(fileext = ".nwk")
  on.exit(unlink(c(tree_path, paste0(tree_path, ".log"))))

  result = make_cmaple_tree(sequences, tree_path, return_tree = FALSE)

  expect_type(result, "character")
  expect_equal(result, tree_path)
})

# Parameter tests ----

test_that("make_cmaple_tree respects multifurcating parameter", {
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  sequences = c(
    seq0 = "ATCGATCGATCGATCG",
    seq1 = "ATCGATCGATCGATCG",
    seq2 = "ATCGATCGATCGATCG",
    seq3 = "ATCGATCGATCGATCG",
    seq4 = "TTCGATCGATCGATCG",
    seq5 = "ATCGATCGATCGAACG"
  )

  # Test with multifurcating = FALSE - should be strictly bifurcating
  tree_path_bifurcating = tempfile(fileext = ".nwk")
  on.exit(unlink(c(
    tree_path_bifurcating,
    paste0(tree_path_bifurcating, ".log")
  )))

  result_bifurcating = make_cmaple_tree(
    sequences,
    tree_path_bifurcating,
    multifurcating = FALSE
  )
  expect_s3_class(result_bifurcating, "phylo")
  # Check that the tree is strictly bifurcating (binary)
  expect_true(ape::is.binary(result_bifurcating))

  tree_path_multi = tempfile(fileext = ".nwk")
  on.exit(
    unlink(c(tree_path_multi, paste0(tree_path_multi, ".log"))),
    add = TRUE
  )

  result_multi = make_cmaple_tree(
    sequences,
    tree_path_multi,
    multifurcating = TRUE
  )
  expect_s3_class(result_multi, "phylo")
  expect_false(ape::is.binary(result_multi))
})

test_that("make_cmaple_tree respects keep_files parameter", {
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  sequences = c(
    seq1 = "ATCGATCGATCGATCG",
    seq2 = "ATCGATTGATCGATCG",
    seq3 = "ATCGATCGAACGATCG"
  )
  tree_path = tempfile(fileext = ".nwk")
  fasta_path = fs::path_ext_set(tree_path, ".fasta")
  log_path = paste0(tree_path, ".log")
  on.exit(unlink(c(tree_path, fasta_path, log_path)))

  result = make_cmaple_tree(
    sequences,
    tree_path,
    keep_files = c("nwk", "log")
  )

  expect_true(fs::file_exists(tree_path))
  expect_true(fs::file_exists(log_path))
  expect_false(fs::file_exists(fasta_path))
})
