gb_cache <- new.env(parent = emptyenv())

cachedGenbank <- function(segment) {
  if (is.null(gb_cache[[segment]])) {
    gb_cache[[segment]] <- read_genbank(references$H3N2[[segment]])
  }
  gb_cache[[segment]]
}

# Trimmed flanks force a non-zero gb_offset, exercising the pwalign remap.
makeForeign <- function(gb_seq) substr(gb_seq, 26, nchar(gb_seq) - 10)

segment_names <- c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2")


# Backward compatibility ---------------------------------------------------

test_that("extract_orfs default output is unchanged from master", {
  skip_if_not_installed("pwalign")
  expected <- readRDS(test_path("fixtures", "extract_orfs_master.rds"))

  for (segment in segment_names) {
    gb <- cachedGenbank(segment)

    expect_equal(
      extract_orfs(gb, gb$sequence),
      expected[[segment]]$identity,
      info = segment
    )

    expect_equal(
      extract_orfs(gb, makeForeign(gb$sequence)),
      expected[[segment]]$foreign,
      info = segment
    )
  }
})

test_that("the identity fast path agrees with the pwalign remap", {
  skip_if_not_installed("pwalign")
  gb <- cachedGenbank("MP")

  expect_equal(
    extract_orfs(gb, gb$sequence),
    extract_orfs(gb, paste0(gb$sequence, "")),
    info = "spliced ORFs, offset 0"
  )
})

test_that("the identity path returns GenBank-native coordinates", {
  gb <- cachedGenbank("HA")
  cds <- purrr::keep(gb$features, \(f) f$type == "CDS")[[1]]

  orfs <- extract_orfs(gb, gb$sequence)

  # gb_offset == 0, so 1-based part starts equal the 0-based feature starts + 1
  expect_equal(orfs$HA$parts[[1]]$start, cds$parts[[1]]$start + 1L)
  expect_equal(orfs$HA$parts[[1]]$end, cds$parts[[1]]$end)
})

test_that("extract_orfs errors when no requested feature is present", {
  gb <- cachedGenbank("NP")
  expect_error(extract_orfs(gb, gb$sequence, types = "mat_peptide"))
})


# The types argument -------------------------------------------------------

test_that("types = CDS is the default and names ORFs by /gene", {
  expect_named(
    extract_orfs(cachedGenbank("PA"), cachedGenbank("PA")$sequence),
    c("PA", "PA-X"),
    ignore.order = TRUE
  )
  expect_named(
    extract_orfs(cachedGenbank("PB1"), cachedGenbank("PB1")$sequence),
    c("PB1", "PB1-F2"),
    ignore.order = TRUE
  )
  expect_named(
    extract_orfs(cachedGenbank("NS"), cachedGenbank("NS")$sequence),
    c("NS1", "NS2"),
    ignore.order = TRUE
  )
})

test_that("mat_peptide and sig_peptide features are named and extracted", {
  gb <- cachedGenbank("HA")
  orfs <- extract_orfs(
    gb,
    gb$sequence,
    types = c("CDS", "mat_peptide", "sig_peptide")
  )

  expect_named(orfs, c("HA", "SigPep", "HA1", "HA2"), ignore.order = TRUE)

  widthOf <- function(name) {
    sum(purrr::map_int(orfs[[name]]$parts, \(p) p$end - p$start + 1L))
  }

  expect_equal(widthOf("HA"), 1701)
  expect_equal(widthOf("SigPep"), 48)
  expect_equal(widthOf("HA1"), 987)
  expect_equal(widthOf("HA2"), 663)

  # sig_peptide + HA1 + HA2 = the HA CDS minus its stop codon
  expect_equal(
    widthOf("SigPep") + widthOf("HA1") + widthOf("HA2"),
    widthOf("HA") - 3L
  )
})

test_that("mat_peptide names come from /product, not /gene", {
  gb <- cachedGenbank("HA")
  mat <- purrr::keep(gb$features, \(f) f$type == "mat_peptide")

  # Both mat_peptides carry /gene="HA", so /gene would collide.
  expect_equal(unique(purrr::map_chr(mat, \(f) f$qualifiers$gene)), "HA")
  expect_setequal(
    purrr::map_chr(mat, \(f) f$qualifiers$product),
    c("HA1", "HA2")
  )
})

test_that("spliced features have exons that must be concatenated", {
  widths <- function(orf) purrr::map_int(orf$parts, \(p) p$end - p$start + 1L)

  m2 <- extract_orfs(cachedGenbank("MP"), cachedGenbank("MP")$sequence)$M2
  ns2 <- extract_orfs(cachedGenbank("NS"), cachedGenbank("NS")$sequence)$NS2
  pax <- extract_orfs(cachedGenbank("PA"), cachedGenbank("PA")$sequence)$`PA-X`

  expect_equal(widths(m2), c(26L, 268L))
  expect_equal(widths(ns2), c(30L, 336L))
  expect_equal(widths(pax), c(573L, 186L))

  # M2's first exon is not codon-aligned, so exons must be joined before the
  # sequence is split into codons. NS2's and PA-X's happen to be.
  expect_true(widths(m2)[[1]] %% 3 != 0)

  for (orf in list(m2, ns2, pax)) {
    expect_equal(sum(widths(orf)) %% 3, 0)
  }
})
