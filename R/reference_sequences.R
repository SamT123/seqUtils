#' Reference GenBank files
#'
#' A nested list of file paths to reference GenBank files bundled with the
#' package. Access files like `references$H3N2$HA`.
#'
#' @export
references <- NULL

.onLoad <- function(libname, pkgname) {
  ref_dir <- system.file("extdata", "references", package = pkgname)
  subtypes <- list.dirs(ref_dir, recursive = FALSE, full.names = TRUE)

  refs <- lapply(subtypes, function(subtype_dir) {
    files <- list.files(subtype_dir, pattern = "\\.gb$", full.names = TRUE)
    setNames(as.list(files), tools::file_path_sans_ext(basename(files)))
  })
  names(refs) <- basename(subtypes)

  assign("references", refs, envir = parent.env(environment()))
}

#' Read a GenBank file
#'
#' Parses a GenBank file using BioPython's SeqIO, returning a structured list
#' with sequence, annotations, and features (including CDS).
#'
#' Requires a Python installation with BioPython. BioPython will be installed
#' automatically on first use if not already available.
#'
#' @param path Path to a GenBank (.gb) file.
#' @return A list with elements: `id`, `name`, `description`, `sequence`,
#'   `annotations`, and `features`. Each feature includes a `parts` field: a
#'   list of `{start, end}` intervals (0-based, half-open). For simple features
#'   this has length 1; for spliced features (e.g. M2, NEP) it contains one
#'   element per exon. `start`/`end` on the feature itself remain the outer
#'   bounds for backward compatibility.
#' @export
read_genbank <- function(path) {
  reticulate::py_require("biopython")
  SeqIO <- reticulate::import("Bio.SeqIO", as = "SeqIO")
  record <- SeqIO$read(path, "genbank")

  list(
    id = record$id,
    name = record$name,
    description = record$description,
    sequence = record$seq$`__str__`(),
    annotations = as.list(record$annotations),
    features = purrr::map(
      record$features,
      function(f) {
        loc <- f$location
        parts <- if (
          reticulate::py_has_attr(loc, "parts") && length(loc$parts) > 1
        ) {
          purrr::map(
            loc$parts,
            ~ list(start = as.integer(.x$start), end = as.integer(.x$end))
          )
        } else {
          list(list(start = as.integer(loc$start), end = as.integer(loc$end)))
        }
        list(
          type = f$type,
          start = as.integer(loc$start),
          end = as.integer(loc$end),
          strand = as.integer(loc$strand),
          parts = parts,
          qualifiers = as.list(f$qualifiers)
        )
      }
    )
  )
}

#' Extract ORF coordinates from a GenBank entry
#'
#' Remaps annotated feature coordinates from a GenBank reference onto another
#' sequence via pairwise alignment. When \code{seq} is the GenBank sequence
#' itself the alignment is skipped, so \code{pwalign} is not needed.
#'
#' @param gb GenBank object from \code{seqUtils::read_genbank}, with
#'   \code{$features} and \code{$sequence}.
#' @param seq Nucleotide sequence string to remap coordinates onto (e.g. root
#'   sequence from a tree tibble).
#' @param types Character vector of GenBank feature types to extract. Defaults
#'   to \code{"CDS"}. \code{"mat_peptide"} and \code{"sig_peptide"} give the
#'   post-cleavage products, e.g. HA1 and HA2 on influenza HA.
#'
#' @return A named list of ORF objects ordered longest-first. Each ORF has a
#'   \code{$parts} list of \code{list(start, end, cum_nt)} in \code{seq}
#'   coordinates. Spliced ORFs (e.g. M2, NS2) have \code{length(parts) > 1}.
#'   \code{CDS} features are named by their \code{/gene} qualifier,
#'   \code{mat_peptide} by \code{/product} (both HA mat_peptides share
#'   \code{/gene="HA"}), and \code{sig_peptide} by the literal \code{"SigPep"}.
#' @export
extract_orfs <- function(gb, seq, types = "CDS") {
  features <- purrr::keep(gb$features, \(f) f$type %in% types)
  if (length(features) == 0) {
    stop(
      "No ",
      paste(types, collapse = "/"),
      " features found in genbank entry"
    )
  }

  gb_seq <- toupper(gb$sequence)
  gb_offset <- if (identical(toupper(seq), gb_seq)) {
    0L
  } else {
    rlang::check_installed(
      "pwalign",
      reason = "to remap GenBank coordinates onto a different sequence."
    )
    aln <- pwalign::pairwiseAlignment(
      Biostrings::DNAString(toupper(seq)),
      Biostrings::DNAString(gb_seq),
      type = "global-local"
    )
    pwalign::start(pwalign::subject(aln)) - 1L
  }

  seq_len <- nchar(seq)

  orfs_names <- purrr::map_chr(features, featureName)

  orfs <- purrr::map(
    features,
    function(f) {
      codon_start <- as.integer(f$qualifiers$codon_start %||% 1L)
      cum_nt <- 0L
      parts <- purrr::imap(f$parts, function(p, i) {
        gb_s <- p$start + if (i == 1L) (codon_start - 1L) else 0L
        gb_e <- p$end - 1L
        part <- list(
          start = gb_s - gb_offset + 1L,
          end = gb_e - gb_offset + 1L,
          cum_nt = cum_nt
        )
        cum_nt <<- cum_nt + (gb_e - gb_s + 1L)
        part
      })
      list(parts = parts)
    }
  )

  orfs <- setNames(orfs, orfs_names)

  orf_span <- function(o) {
    c(
      min(purrr::map_int(o$parts, ~ .x$start)),
      max(purrr::map_int(o$parts, ~ .x$end))
    )
  }
  orfs <- purrr::keep(
    orfs,
    function(o) {
      sp <- orf_span(o)
      sp[2] >= 1L && sp[1] <= seq_len
    }
  )

  total_len <- purrr::map_int(
    orfs,
    ~ sum(purrr::map_int(.x$parts, ~ .x$end - .x$start + 1L))
  )
  orfs[order(total_len, decreasing = TRUE)]
}

featureName <- function(f) {
  name <- switch(
    f$type,
    mat_peptide = f$qualifiers$product,
    sig_peptide = "SigPep",
    f$qualifiers$gene
  )
  as.character(name)[[1]]
}


# Kept for backwards compatibility

#' 550 aa HA sequence for Alaska/232/2015 (useful for aligning)
#' @export
alaska_232_2015_aas <- "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSSSSFFSRLNWLTHLNYTYPALNVTMPNKEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKHSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGRGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRVQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNETYDHNVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"

#' 1650 nt HA sequence for Alaska/232/2015 (useful for aligning)
#' @export
alaska_232_2015_nts <- "CAAAAAATTCCTGGAAATGACAATAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACAAATGACCGAATTGAAGTTACTAATGCTACTGAGTTGGTTCAGAATTCCTCAATAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAGAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTTCAAAATAAGAAATGGGACCTTTTTGTTGAACGAAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGACTGGAGTCACTCAAAACGGAACAAGTTCTGCTTGCATAAGGAGATCTAGTAGTAGTTTCTTTAGTAGATTAAATTGGTTGACCCACTTAAACTACACATATCCAGCATTGAACGTGACTATGCCAAACAAGGAACAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAAGGACCAAATCTTCCTGTATGCTCAATCATCAGGAAGAATCACAGTATCTACCAAAAGAAGCCAACAAGCTGTAATCCCAAATATCGGATCTAGACCCAGAATAAGGGATATCCCTAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAAGTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCATAGCACTCTGAAATTGGCAACAGGAATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATAGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAAGAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATCGGTTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTAGAAGGAAGAGTTCAAGACCTTGAGAAATATGTTGAGGACACTAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGAAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATAAGAAATGAAACTTATGACCACAATGTGTACAGGGATGAAGCATTAAACAACCGGTTCCAGATCAAGGGAGTTGAGCTGAAGTCAGGGTACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAGGGCAACATTAGATGCAACATTTGCATT"

#' aa HA sequence for Alaska/232/2015 (useful for aligning)
#' @export
california_04_2009_h1_aas <-
  "DTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSIYQNADTYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI"

#' nt HA sequence for Alaska/232/2015 (useful for aligning)
#' @export
california_04_2009_h1_nts <-
  "GACACATTATGTATAGGTTATCATGCGAACAATTCAACAGACACTGTAGACACAGTACTAGAAAAGAATGTAACAGTAACACACTCTGTTAACCTTCTAGAAGACAAGCATAACGGGAAACTATGCAAACTAAGAGGGGTAGCCCCATTGCATTTGGGTAAATGTAACATTGCTGGCTGGATCCTGGGAAATCCAGAGTGTGAATCACTCTCCACAGCAAGCTCATGGTCCTACATTGTGGAAACACCTAGTTCAGACAATGGAACGTGTTACCCAGGAGATTTCATCGATTATGAGGAGCTAAGAGAGCAATTGAGCTCAGTGTCATCATTTGAAAGGTTTGAGATATTCCCCAAGACAAGTTCATGGCCCAATCATGACTCGAACAAAGGTGTAACGGCAGCATGTCCTCATGCTGGAGCAAAAAGCTTCTACAAAAATTTAATATGGCTAGTTAAAAAAGGAAATTCATACCCAAAGCTCAGCAAATCCTACATTAATGATAAAGGGAAAGAAGTCCTCGTGCTATGGGGCATTCACCATCCATCTACTAGTGCTGACCAACAAAGTATCTATCAGAATGCAGATACATATGTTTTTGTGGGGTCATCAAGATACAGCAAGAAGTTCAAGCCGGAAATAGCAATAAGACCCAAAGTGAGGGATCAAGAAGGGAGAATGAACTATTACTGGACACTAGTAGAGCCGGGAGACAAAATAACATTCGAAGCAACTGGAAATCTAGTGGTACCGAGATATGCATTCGCAATGGAAAGAAATGCTGGATCTGGTATTATCATTTCAGATACACCAGTCCACGATTGCAATACAACTTGTCAAACACCCAAGGGTGCTATAAACACCAGCCTCCCATTTCAGAATATACATCCGATCACAATTGGAAAATGTCCAAAATATGTGAAAAGCACAAAATTGAGACTGGCCACAGGATTGAGGAATATCCCGTCTATTCAATCTAGAGGCCTATTTGGGGCCATTGCCGGTTTCATTGAAGGGGGGTGGACAGGGATGGTAGATGGATGGTACGGTTATCACCATCAAAATGAGCAGGGGTCAGGATATGCAGCCGACCTGAAGAGCACACAGAATGCCATTGACGAGATTACTAACAAAGTAAATTCTGTTATTGAAAAGATGAATACACAGTTCACAGCAGTAGGTAAAGAGTTCAACCACCTGGAAAAAAGAATAGAGAATTTAAATAAAAAAGTTGATGATGGTTTCCTGGACATTTGGACTTACAATGCCGAACTGTTGGTTCTATTGGAAAATGAAAGAACTTTGGACTACCACGATTCAAATGTGAAGAACTTATATGAAAAGGTAAGAAGCCAGCTAAAAAACAATGCCAAGGAAATTGGAAACGGCTGCTTTGAATTTTACCACAAATGCGATAACACGTGCATGGAAAGTGTCAAAAATGGGACTTATGACTACCCAAAATACTCAGAGGAAGCAAAATTAAACAGAGAAGAAATAGATGGGGTAAAGCTGGAATCAACAAGGATTTACCAGATTTTGGCGATCTATTCAACTGTCGCCAGTTCATTGGTACTGGTAGTCTCCCTGGGGGCAATCAGTTTCTGGATGTGCTCTAATGGGTCTCTACAGTGTAGAATATGTATT"
