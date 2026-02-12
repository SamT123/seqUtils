test_that("references list is populated with H3N2 segments", {
  expect_type(references, "list")
  expect_true("H3N2" %in% names(references))

  expected_segments = c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2")
  expect_equal(sort(names(references$H3N2)), sort(expected_segments))
})

test_that("reference paths point to existing .gb files", {
  for (segment in names(references$H3N2)) {
    path = references$H3N2[[segment]]
    expect_true(file.exists(path), info = segment)
    expect_match(path, "\\.gb$")
  }
})

test_that("read_genbank parses a GenBank file", {
  gb = read_genbank(references$H3N2$HA)

  expect_type(gb, "list")
  expect_named(
    gb,
    c("id", "name", "description", "sequence", "annotations", "features")
  )

  expect_equal(gb$id, "CY002064.1")
  expect_match(gb$description, "H3N2")
  expect_type(gb$sequence, "character")
  expect_true(nchar(gb$sequence) > 0)
})

test_that("read_genbank extracts CDS features", {
  gb = read_genbank(references$H3N2$HA)

  types = vapply(gb$features, \(f) f$type, character(1))
  expect_true("CDS" %in% types)

  cds = gb$features[types == "CDS"][[1]]
  expect_true(all(
    c("type", "start", "end", "strand", "qualifiers") %in% names(cds)
  ))
  expect_type(cds$start, "integer")
  expect_type(cds$end, "integer")
})
