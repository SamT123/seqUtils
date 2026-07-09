# Helper function to check if mafft is available
mafft_available <- function() {
  result <- suppressWarnings(system(
    "mafft --version",
    ignore.stdout = TRUE,
    ignore.stderr = TRUE
  ))
  return(result == 0)
}

test_that("translate works for basic DNA sequence without deletions", {
  seq <- c("ATGGCC")
  result <- translate(seq)

  expect_equal(as.character(result[1]), "MA")
  expect_equal(length(result), 1)
})

test_that("translate preserves sequence names", {
  seqs <- c("ATGGCC", "ATGTTT")
  names(seqs) <- c("seq1", "seq2")

  result <- translate(seqs)

  expect_equal(names(result), c("seq1", "seq2"))
  expect_equal(as.character(result["seq1"]), "MA")
  expect_equal(as.character(result["seq2"]), "MF")
})

test_that("translate errors when deletions present without reference", {
  seq <- c("ATG---GCC")

  expect_error(
    translate(seq),
    "Deletions present but no reference for alignment provided"
  )
})

test_that("translate with reference_aas handles deletions correctly", {
  skip_if_not(mafft_available(), "mafft not installed")

  seqs <- c("ATGGCCAAA", "ATG---AAA")
  names(seqs) <- c("ref_seq", "del_seq")

  ref_aa <- translate(seqs[1])

  result <- suppressMessages(translate(seqs, reference_aas = ref_aa))

  expect_equal(names(result), c("ref_seq", "del_seq"))
  # Both should align to the reference
  expect_equal(as.character(result["ref_seq"]), "MAK")
})

test_that("translate handles leading deletions with reference", {
  skip_if_not(mafft_available(), "mafft not installed")

  seq <- c("-ATGGCCAA")
  ref <- c("GQ")

  result <- suppressMessages(translate(seq, reference_aas = ref))

  expect_equal(as.character(result[1]), "GQ")
})

test_that("translate handles NA sequences with warning", {
  seqs <- c("ATGGCC", NA, "ATGTTT")
  names(seqs) <- c("seq1", "seq2", "seq3")

  expect_warning(
    result <- translate(seqs),
    "Some sequences are NA"
  )

  expect_true(is.na(result["seq2"]))
  expect_equal(as.character(result["seq1"]), "MA")
  expect_equal(as.character(result["seq3"]), "MF")
})

test_that("translate handles sequences not divisible by 3", {
  seq <- c("ATGGCCT")

  result <- translate(seq)

  # Should only translate first 6 bases: ATGGCC -> MA
  expect_equal(as.character(result[1]), "MA")
})

test_that("translate handles fuzzy codons", {
  seq <- c("ATGNNN")

  result <- translate(seq)

  # Should translate ATG to M, NNN to X (unknown)
  expect_equal(as.character(result[1]), "MX")
})

test_that("translate handles stop codons", {
  seq <- c("ATGTAA")

  result <- translate(seq)

  # Should translate ATG to M, TAA to * (stop)
  expect_equal(as.character(result[1]), "M*")
})

test_that("translate handles empty sequence vector", {
  seqs <- character(0)

  result <- translate(seqs)

  expect_equal(length(result), 0)
})


# aligned = TRUE -----------------------------------------------------------

sliceParts <- function(seq, parts) {
  paste0(
    purrr::map_chr(parts, \(p) substr(seq, p$start, p$end)),
    collapse = ""
  )
}

test_that("aligned translation round-trips every CDS in the bundled references", {
  paths <- c(
    setNames(list(references$H1N1$HA), "H1N1/HA"),
    setNames(as.list(references$H3N2), paste0("H3N2/", names(references$H3N2)))
  )

  for (label in names(paths)) {
    gb <- read_genbank(paths[[label]])
    orfs <- extract_orfs(gb, gb$sequence)

    for (feature in purrr::keep(gb$features, \(f) f$type == "CDS")) {
      gene <- as.character(feature$qualifiers$gene)[[1]]
      nt <- sliceParts(gb$sequence, orfs[[gene]]$parts)

      expect_equal(
        unname(translate(nt, aligned = TRUE)),
        paste0(as.character(feature$qualifiers$translation)[[1]], "*"),
        info = paste(label, gene)
      )
    }
  }
})

test_that("aligned translation handles gaps, ambiguity and stops", {
  seqs <- c(
    full = "ATGGCCTAA",
    all_gap_codon = "ATG---TAA",
    partial_gap = "ATGG-CTAA",
    ambiguous = "ATGGNCTAA",
    premature_stop = "ATGTAATAA"
  )

  expect_equal(
    translate(seqs, aligned = TRUE),
    c(
      full = "MA*",
      all_gap_codon = "M-*",
      partial_gap = "MX*",
      ambiguous = "MX*",
      premature_stop = "M**"
    )
  )
})

test_that("aligned translation is coordinate-stable", {
  seqs <- c("ATGGCCTAA", "ATG---TAA", "---GCCTAA")
  result <- translate(seqs, aligned = TRUE)

  expect_equal(unique(nchar(result)), 3L)
})

test_that("aligned translation preserves names and NA sequences", {
  seqs <- c(a = "ATGGCC", b = NA, c = "ATGTTT")

  expect_warning(result <- translate(seqs, aligned = TRUE), "NA")

  expect_equal(names(result), c("a", "b", "c"))
  expect_equal(unname(result), c("MA", NA, "MF"))
})

test_that("aligned translation rejects input that is not an in-frame alignment", {
  expect_error(
    translate(c("ATGGC", "ATGGC"), aligned = TRUE),
    class = "seqUtils_error_not_an_alignment"
  )
  expect_error(
    translate(c("ATGGCC", "ATG"), aligned = TRUE),
    class = "seqUtils_error_not_an_alignment"
  )
})

test_that("aligned translation rejects a reference, which it cannot use", {
  expect_error(
    translate("ATGGCC", reference_aas = "MA", aligned = TRUE),
    class = "seqUtils_error_aligned_with_reference"
  )
})
