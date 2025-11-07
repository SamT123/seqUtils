#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom purrr map_chr
#' @importFrom stats setNames
## usethis namespace: end
NULL

.onAttach <- function(libname, pkgname) {
  if (Sys.which("mafft") == "") {
    packageStartupMessage(
      "Note: MAFFT not found in PATH. Install MAFFT for alignment functions."
    )
  }
}
