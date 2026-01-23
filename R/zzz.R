.onAttach <- function(libname, pkgname) {
  # Check if C++ implementation is available
  cpp_available <- exists("compute_roc_pr_metrics", mode = "function") &&
    exists("predict_ligand_activities_cpp", mode = "function")

  if (cpp_available) {
    packageStartupMessage(
      "NicheNet loaded with C++ implementation enabled\n",
      "To disable: options(nichenetr.disable_cpp = TRUE)"
    )
  } else {
    packageStartupMessage(
      "NicheNet loaded with R implementation\n",
      "For better performance, reinstall: devtools::install_local(force = TRUE)"
    )
  }
}

.onLoad <- function(libname, pkgname) {
  # Set default options
  op <- options()
  op.nichenetr <- list(
    nichenetr.disable_cpp = FALSE,
    nichenetr.verbose = TRUE
  )
  toset <- !(names(op.nichenetr) %in% names(op))
  if (any(toset)) options(op.nichenetr[toset])

  invisible()
}
