# Helper function to check if CMAPLE is available
cmaple_available <- function() {
  result <- suppressWarnings(system(
    "cmaple --help",
    ignore.stdout = TRUE,
    ignore.stderr = TRUE
  ))
  return(result == 0)
}

test_that("root_tree_using_outsequence roots a tree correctly and removes outgroup", {
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  # Create a simple unrooted tree
  tree <- ape::read.tree(text = "((seq1:1,seq2:1):1,(seq3:1,seq4:1):1);")

  sequences <- c(
    seq1 = "ATCGATCGATCGATCG",
    seq2 = "ATCGATTGATCGATCG",
    seq3 = "ATCGATCGGGGGGGGG",
    seq4 = "TTCGATCGATCGATCG"
  )

  outsequence <- "GGGGGGGGGGGGGGGG"

  result <- root_tree_using_outsequence(tree, sequences, outsequence)

  expect_s3_class(result, "phylo")
  expect_equal(length(result$tip.label), 4)
  expect_true(all(c("seq1", "seq2", "seq3", "seq4") %in% result$tip.label))
  expect_false("outsequence" %in% result$tip.label)
  expect_true(!is.null(result$root.edge) || result$Nnode > 0)
  root <- castor::find_root(result)
  expect_true(
    "seq3" %in% result$tip.label[result$edge[, 2][result$edge[, 1] == root]]
  )
})

test_that("root_tree_using_outsequence keeps outgroup when requested", {
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  tree <- ape::read.tree(text = "((seq1:1,seq2:1):1,(seq3:1,seq4:1):1);")

  sequences <- c(
    seq1 = "ATCGATCGATCGATCG",
    seq2 = "ATCGATTGATCGATCG",
    seq3 = "ATCGATCGAACGATCG",
    seq4 = "TTCGATCGATCGATCG"
  )

  outsequence <- "GGGGGGGGGGGGGGGG"

  result <- root_tree_using_outsequence(
    tree,
    sequences,
    outsequence,
    remove_outsequence = FALSE
  )

  expect_s3_class(result, "phylo")
  expect_equal(length(result$tip.label), 5)
  expect_true("outsequence" %in% result$tip.label)
})

test_that("root_tree_using_outsequence forwards CMAPLE settings", {
  captured <- NULL
  local_mocked_bindings(
    make_cmaple_tree = function(sequences, tree_path, ...) {
      captured <<- list(...)
      castor::write_tree(
        ape::read.tree(
          text = "(((seq1:1,seq2:1):1,(seq3:1,seq4:1):1):1,outsequence:1);"
        ),
        file = tree_path
      )
      invisible(NULL)
    }
  )

  tree <- ape::read.tree(text = "((seq1:1,seq2:1):1,(seq3:1,seq4:1):1);")
  sequences <- c(
    seq1 = "ATCGATCGATCGATCG",
    seq2 = "ATCGATTGATCGATCG",
    seq3 = "ATCGATCGAACGATCG",
    seq4 = "TTCGATCGATCGATCG"
  )

  root_tree_using_outsequence(
    tree,
    sequences,
    outsequence = "GGGGGGGGGGGGGGGG",
    model = "UNREST",
    site_rate = "SCALAR",
    seed = 99,
    num_threads = 4
  )

  expect_equal(captured$model, "UNREST")
  expect_equal(captured$site_rate, "SCALAR")
  expect_equal(captured$seed, 99)
  expect_equal(captured$num_threads, 4)
})

test_that("root_tree_using_outsequence defaults reproduce CMAPLE defaults", {
  captured <- NULL
  local_mocked_bindings(
    make_cmaple_tree = function(sequences, tree_path, ...) {
      captured <<- list(...)
      castor::write_tree(
        ape::read.tree(text = "((seq1:1,seq2:1):1,outsequence:1);"),
        file = tree_path
      )
      invisible(NULL)
    }
  )

  root_tree_using_outsequence(
    ape::read.tree(text = "(seq1:1,seq2:1);"),
    c(seq1 = "ATCGATCG", seq2 = "ATCGATTG"),
    outsequence = "GGGGGGGG"
  )

  expect_equal(captured$model, "GTR")
  expect_equal(captured$site_rate, "NONE")
  expect_null(captured$seed)
  expect_equal(captured$num_threads, "AUTO")
})
