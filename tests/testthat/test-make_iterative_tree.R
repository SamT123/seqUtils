iqtree_available = function() {
  ok = function(bin) {
    res = suppressWarnings(
      system(paste(bin, "--help"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    )
    res == 0
  }
  ok("iqtree2") || ok("iqtree")
}

cmaple_available = function() {
  res = suppressWarnings(
    system("cmaple --help", ignore.stdout = TRUE, ignore.stderr = TRUE)
  )
  res == 0
}

# Build a tiny synthetic alignment: 3 years x 8 seqs, all of length 60, with a
# few SNPs so that tree inference is well-defined.
make_synthetic_alignment = function() {
  base = strsplit(
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    ""
  )[[
    1
  ]]
  vary = function(positions, alts) {
    s = base
    s[positions] = alts
    paste(s, collapse = "")
  }
  rows = list()
  i = 1
  for (yr in c(2018L, 2019L, 2020L)) {
    for (k in 1:8) {
      positions = ((yr - 2017L) * 5L + k) %% 60L + 1L
      alt = c("A", "T", "C", "G")[((yr * k) %% 4L) + 1L]
      rows[[i]] = list(
        Isolate_unique_identifier = paste0("s", yr, "_", k),
        dna_sequence = vary(positions, alt),
        Collection_date = as.Date(paste0(yr, "-06-15"))
      )
      i = i + 1
    }
  }
  do.call(rbind, lapply(rows, as.data.frame))
}


test_that("subsamplePerYear caps per-year count and honours forced ids", {
  aln = make_synthetic_alignment()
  aln$year = lubridate::year(aln$Collection_date)
  out = subsamplePerYear(
    aln,
    n_per_year = 3,
    forced_identifiers = "s2018_1",
    seed = 1
  )
  # 3 years * 3/year = 9 picks; forced id is one of those (already in year
  # 2018), so total should be <= 9
  expect_lte(nrow(out), 9)
  expect_true("s2018_1" %in% out$Isolate_unique_identifier)
})


test_that("make_iterative_tree errors on bad cascade_sizes", {
  aln = make_synthetic_alignment()
  expect_error(
    make_iterative_tree(
      aln,
      tree_path = tempfile(fileext = ".nwk"),
      initial_iqtree_size = 2,
      cascade_sizes = integer(0)
    ),
    "cascade_sizes"
  )
  expect_error(
    make_iterative_tree(
      aln,
      tree_path = tempfile(fileext = ".nwk"),
      initial_iqtree_size = 4,
      cascade_sizes = c(3, 5)
    ),
    "cascade_sizes must all be greater"
  )
  expect_error(
    make_iterative_tree(
      aln,
      tree_path = tempfile(fileext = ".nwk"),
      initial_iqtree_size = 2,
      cascade_sizes = c(5, 3)
    ),
    "strictly increasing"
  )
})


test_that("make_iterative_tree errors when required columns missing", {
  bad = data.frame(a = 1, b = 2)
  expect_error(
    make_iterative_tree(
      bad,
      tree_path = tempfile(fileext = ".nwk"),
      initial_iqtree_size = 2,
      cascade_sizes = c(4)
    ),
    "missing required columns"
  )
})


test_that("make_iterative_tree builds a tree end-to-end on a tiny alignment", {
  skip_if_not(iqtree_available(), "IQ-TREE not installed")
  skip_if_not(cmaple_available(), "CMAPLE not installed")

  aln = make_synthetic_alignment()
  out_dir = tempfile()
  fs::dir_create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  tree_path = fs::path(out_dir, "final", ext = "nwk")
  tree = make_iterative_tree(
    alignment = aln,
    tree_path = tree_path,
    initial_iqtree_size = 2,
    cascade_sizes = c(3, 5),
    work_dir = out_dir,
    seed = 1
  )

  expect_s3_class(tree, "phylo")
  # Final step uses the entire alignment
  expect_setequal(tree$tip.label, aln$Isolate_unique_identifier)
  # Cascade intermediate files exist
  expect_true(fs::file_exists(fs::path(
    out_dir,
    "cascade",
    "iqtree_initial",
    ext = "nwk"
  )))
  expect_true(fs::file_exists(fs::path(
    out_dir,
    "cascade",
    "size_3",
    ext = "nwk"
  )))
  expect_true(fs::file_exists(fs::path(
    out_dir,
    "cascade",
    "size_5",
    ext = "nwk"
  )))
  expect_true(fs::file_exists(tree_path))
})
