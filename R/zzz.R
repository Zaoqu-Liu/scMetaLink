#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "scMetaLink v", utils::packageVersion("scMetaLink"), "\n",
    "Metabolite-mediated cell communication analysis\n",
    "Documentation: https://Zaoqu-Liu.github.io/scMetaLink/"
  )
}

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.scMetaLink <- list(
    scMetaLink.verbose = TRUE,
    scMetaLink.n_cores = 1L,
    scMetaLink.seed = 42L
  )
  toset <- !(names(op.scMetaLink) %in% names(op))
  if (any(toset)) options(op.scMetaLink[toset])
  invisible()
}

#' Show Package Citation
#' @description Display the recommended citation for scMetaLink
#' @return Invisibly returns NULL, prints citation information
#' @export
#' @examples
#' citationScMetaLink()
citationScMetaLink <- function() {
  cat("\n")
  cat("To cite scMetaLink in publications, please use:\n\n")
  cat("  Liu Z (2026). scMetaLink: Inferring metabolite-mediated cell-cell\n")
  cat("  communication from single-cell transcriptomic data.\n")
  cat("  R package version ", as.character(utils::packageVersion("scMetaLink")), ".\n", sep = "")
  cat("  https://github.com/Zaoqu-Liu/scMetaLink\n\n")
  cat("A BibTeX entry for LaTeX users is:\n\n")
  cat("  @Manual{,\n")
  cat("    title = {scMetaLink: Single-Cell Metabolite-Mediated Cell Communication Analysis},\n")
  cat("    author = {Zaoqu Liu},\n")
  cat("    year = {2026},\n")
  cat("    note = {R package version ", as.character(utils::packageVersion("scMetaLink")), "},\n", sep = "")
  cat("    url = {https://github.com/Zaoqu-Liu/scMetaLink},\n")
  cat("  }\n\n")
  invisible(NULL)
}
