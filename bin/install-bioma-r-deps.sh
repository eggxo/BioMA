#!/usr/bin/env sh
set -eu

RSCRIPT=${RSCRIPT:-Rscript}

"$RSCRIPT" - <<'RS'
if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("The remotes R package is required; create the BioMA environment first.")
}

install_exact <- function(package, version, installer) {
  installed <- requireNamespace(package, quietly = TRUE)
  if (installed && identical(as.character(packageVersion(package)), version)) {
    message(package, " ", version, " is already installed")
    return(invisible(TRUE))
  }

  installer()
  installed <- requireNamespace(package, quietly = TRUE)
  actual <- if (installed) as.character(packageVersion(package)) else "missing"
  if (!installed || !identical(actual, version)) {
    stop(package, " installation check failed: expected ", version, ", found ", actual)
  }
  message(package, " ", version, " installed")
  invisible(TRUE)
}

install_exact("extendedForest", "1.6.2", function() {
  remotes::install_url(
    "https://download.r-forge.r-project.org/src/contrib/extendedForest_1.6.2.tar.gz",
    dependencies = FALSE,
    upgrade = "never"
  )
})
install_exact("gradientForest", "0.1.37", function() {
  remotes::install_url(
    "https://download.r-forge.r-project.org/src/contrib/gradientForest_0.1-37.tar.gz",
    dependencies = FALSE,
    upgrade = "never"
  )
})
install_exact("sars", "2.0.0", function() {
  remotes::install_url(
    "https://cran.r-project.org/src/contrib/Archive/sars/sars_2.0.0.tar.gz",
    dependencies = FALSE,
    upgrade = "never"
  )
})
install_exact("mar", "0.2.0", function() {
  remotes::install_github(
    "meixilin/mar@f2a60772504a52e518d6827a8d22de1ef4dd11d4",
    dependencies = FALSE,
    upgrade = "never"
  )
})
RS
