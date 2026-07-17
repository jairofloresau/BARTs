############################################################################
## Instalacion reproducible de dependencias para el flujo QBART.
############################################################################

cran_packages <- c(
  "dplyr", "tidyr", "ggplot2", "quantreg", "quantregForest",
  "digest", "Rcpp", "RcppArmadillo"
)
missing <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0L) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BayesQArt", quietly = TRUE)) {
  commit <- "3f83a7fed1afc28f22dadc07aa4c4ed9f611cc2a"
  archive <- tempfile(fileext = ".zip")
  extract_dir <- tempfile("bayesqart-")
  dir.create(extract_dir)
  url <- paste0(
    "https://github.com/bpkindo/bayesqart/archive/",
    commit,
    ".zip"
  )
  download.file(url, archive, mode = "wb")
  unzip(archive, exdir = extract_dir)
  package_dir <- file.path(
    extract_dir,
    paste0("bayesqart-", commit),
    "BayesQArt"
  )
  status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", shQuote(package_dir))
  )
  if (!identical(status, 0L)) {
    stop(
      "No se pudo compilar BayesQArt. En Windows verifique que Rtools ",
      "corresponda a su version de R."
    )
  }
}

cat("Dependencias QBART instaladas correctamente.\n")
