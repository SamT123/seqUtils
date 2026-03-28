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

#' Extract CDS ORF coordinates from a GenBank entry
#'
#' Remaps CDS feature coordinates from a GenBank reference onto the tree root
#' sequence via pairwise alignment.
#'
#' @param gb GenBank object from \code{seqUtils::read_genbank}, with
#'   \code{$features} and \code{$sequence}.
#' @param seq Nucleotide sequence string to remap coordinates onto (e.g. root
#'   sequence from a tree tibble).
#'
#' @return A named list of ORF objects ordered longest-first. Each ORF has a
#'   \code{$parts} list of \code{list(start, end, cum_nt)} in \code{seq}
#'   coordinates. Spliced ORFs (e.g. M2, NEP) have \code{length(parts) > 1}.
#' @export
extract_orfs <- function(gb, seq) {
  cds <- keep(gb$features, ~ .x$type == "CDS")
  if (length(cds) == 0) {
    stop("No CDS features found in genbank entry")
  }

  gb_seq <- toupper(gb$sequence)
  aln <- pwalign::pairwiseAlignment(
    Biostrings::DNAString(toupper(seq)),
    Biostrings::DNAString(gb_seq),
    type = "global-local"
  )
  gb_offset <- pwalign::start(pwalign::subject(aln)) - 1L

  seq_len <- nchar(seq)

  orfs_names = map_chr(cds, \(f) f$qualifiers$gene)

  orfs <- map(
    cds,
    function(f) {
      codon_start <- as.integer(f$qualifiers$codon_start %||% 1L)
      cum_nt <- 0L
      parts <- imap(f$parts, function(p, i) {
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

  orfs = setNames(orfs, orfs_names)

  orf_span <- function(o) {
    c(min(map_int(o$parts, ~ .x$start)), max(map_int(o$parts, ~ .x$end)))
  }
  orfs <- keep(
    orfs,
    function(o) {
      sp <- orf_span(o)
      sp[2] >= 1L && sp[1] <= seq_len
    }
  )

  total_len <- map_int(orfs, ~ sum(map_int(.x$parts, ~ .x$end - .x$start + 1L)))
  orfs[order(total_len, decreasing = TRUE)]
}


# Kept for backwards compatibility

#' 550 aa HA sequence for Alaska/232/2015 (useful for aligning)
#' @export
alaska_232_2015_aas = "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSSSSFFSRLNWLTHLNYTYPALNVTMPNKEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKHSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGRGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRVQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNETYDHNVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"

#' 1650 nt HA sequence for Alaska/232/2015 (useful for aligning)
#' @export
alaska_232_2015_nts = "CAAAAAATTCCTGGAAATGACAATAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACAAATGACCGAATTGAAGTTACTAATGCTACTGAGTTGGTTCAGAATTCCTCAATAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAGAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTTCAAAATAAGAAATGGGACCTTTTTGTTGAACGAAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGACTGGAGTCACTCAAAACGGAACAAGTTCTGCTTGCATAAGGAGATCTAGTAGTAGTTTCTTTAGTAGATTAAATTGGTTGACCCACTTAAACTACACATATCCAGCATTGAACGTGACTATGCCAAACAAGGAACAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAAGGACCAAATCTTCCTGTATGCTCAATCATCAGGAAGAATCACAGTATCTACCAAAAGAAGCCAACAAGCTGTAATCCCAAATATCGGATCTAGACCCAGAATAAGGGATATCCCTAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAAGTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCATAGCACTCTGAAATTGGCAACAGGAATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATAGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAAGAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATCGGTTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTAGAAGGAAGAGTTCAAGACCTTGAGAAATATGTTGAGGACACTAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGAAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATAAGAAATGAAACTTATGACCACAATGTGTACAGGGATGAAGCATTAAACAACCGGTTCCAGATCAAGGGAGTTGAGCTGAAGTCAGGGTACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAGGGCAACATTAGATGCAACATTTGCATT"

#' aa HA sequence for Alaska/232/2015 (useful for aligning)
#' @export
california_04_2009_h1_aas =
  "DTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSIYQNADTYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI"

#' nt HA sequence for Alaska/232/2015 (useful for aligning)
#' @export
california_04_2009_h1_nts =
  "GACACATTATGTATAGGTTATCATGCGAACAATTCAACAGACACTGTAGACACAGTACTAGAAAAGAATGTAACAGTAACACACTCTGTTAACCTTCTAGAAGACAAGCATAACGGGAAACTATGCAAACTAAGAGGGGTAGCCCCATTGCATTTGGGTAAATGTAACATTGCTGGCTGGATCCTGGGAAATCCAGAGTGTGAATCACTCTCCACAGCAAGCTCATGGTCCTACATTGTGGAAACACCTAGTTCAGACAATGGAACGTGTTACCCAGGAGATTTCATCGATTATGAGGAGCTAAGAGAGCAATTGAGCTCAGTGTCATCATTTGAAAGGTTTGAGATATTCCCCAAGACAAGTTCATGGCCCAATCATGACTCGAACAAAGGTGTAACGGCAGCATGTCCTCATGCTGGAGCAAAAAGCTTCTACAAAAATTTAATATGGCTAGTTAAAAAAGGAAATTCATACCCAAAGCTCAGCAAATCCTACATTAATGATAAAGGGAAAGAAGTCCTCGTGCTATGGGGCATTCACCATCCATCTACTAGTGCTGACCAACAAAGTATCTATCAGAATGCAGATACATATGTTTTTGTGGGGTCATCAAGATACAGCAAGAAGTTCAAGCCGGAAATAGCAATAAGACCCAAAGTGAGGGATCAAGAAGGGAGAATGAACTATTACTGGACACTAGTAGAGCCGGGAGACAAAATAACATTCGAAGCAACTGGAAATCTAGTGGTACCGAGATATGCATTCGCAATGGAAAGAAATGCTGGATCTGGTATTATCATTTCAGATACACCAGTCCACGATTGCAATACAACTTGTCAAACACCCAAGGGTGCTATAAACACCAGCCTCCCATTTCAGAATATACATCCGATCACAATTGGAAAATGTCCAAAATATGTGAAAAGCACAAAATTGAGACTGGCCACAGGATTGAGGAATATCCCGTCTATTCAATCTAGAGGCCTATTTGGGGCCATTGCCGGTTTCATTGAAGGGGGGTGGACAGGGATGGTAGATGGATGGTACGGTTATCACCATCAAAATGAGCAGGGGTCAGGATATGCAGCCGACCTGAAGAGCACACAGAATGCCATTGACGAGATTACTAACAAAGTAAATTCTGTTATTGAAAAGATGAATACACAGTTCACAGCAGTAGGTAAAGAGTTCAACCACCTGGAAAAAAGAATAGAGAATTTAAATAAAAAAGTTGATGATGGTTTCCTGGACATTTGGACTTACAATGCCGAACTGTTGGTTCTATTGGAAAATGAAAGAACTTTGGACTACCACGATTCAAATGTGAAGAACTTATATGAAAAGGTAAGAAGCCAGCTAAAAAACAATGCCAAGGAAATTGGAAACGGCTGCTTTGAATTTTACCACAAATGCGATAACACGTGCATGGAAAGTGTCAAAAATGGGACTTATGACTACCCAAAATACTCAGAGGAAGCAAAATTAAACAGAGAAGAAATAGATGGGGTAAAGCTGGAATCAACAAGGATTTACCAGATTTTGGCGATCTATTCAACTGTCGCCAGTTCATTGGTACTGGTAGTCTCCCTGGGGGCAATCAGTTTCTGGATGTGCTCTAATGGGTCTCTACAGTGTAGAATATGTATT"
