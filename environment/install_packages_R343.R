legacy_lib <- file.path(getwd(), "environment", "legacy_R343_library")
source_dir <- file.path(getwd(), "environment", "cran-source")

.libPaths(c(legacy_lib, .Library))

archives <- c(
  "backports_1.1.1.tar.gz",
  "checkmate_1.8.5.tar.gz",
  "colorspace_1.3-2.tar.gz",
  "digest_0.6.12.tar.gz",
  "iterators_1.0.8.tar.gz",
  "foreach_1.4.3.tar.gz",
  "Formula_1.2-2.tar.gz",
  "gtable_0.2.0.tar.gz",
  "labeling_0.3.tar.gz",
  "lazyeval_0.2.1.tar.gz",
  "magrittr_1.5.tar.gz",
  "R6_2.2.2.tar.gz",
  "RColorBrewer_1.1-2.tar.gz",
  "Rcpp_0.12.14.tar.gz",
  "rlang_0.1.4.tar.gz",
  "munsell_0.4.3.tar.gz",
  "plyr_1.8.4.tar.gz",
  "stringi_1.1.6.tar.gz",
  "stringr_1.2.0.tar.gz",
  "reshape2_1.4.2.tar.gz",
  "tibble_1.3.4.tar.gz",
  "dichromat_2.0-0.tar.gz",
  "viridisLite_0.2.0.tar.gz",
  "scales_0.5.0.tar.gz",
  "ggplot2_2.2.1.tar.gz",
  "glmnet_2.0-13.tar.gz",
  "MatrixModels_0.4-1.tar.gz",
  "SparseM_1.77.tar.gz",
  "snow_0.4-2.tar.gz",
  "quantreg_5.34.tar.gz",
  "hqreg_1.4.tar.gz",
  "mvtnorm_1.0-6.tar.gz",
  "doSNOW_1.0.16.tar.gz",
  "hdm_0.2.0.tar.gz"
)

for (archive in archives) {
  path <- file.path(source_dir, archive)
  if (!file.exists(path)) stop("Missing archive: ", path)
  package <- sub("_.*$", "", archive)
  expected <- sub("^[^_]+_", "", sub("[.]tar[.]gz$", "", archive))
  current <- tryCatch(
    packageDescription(package, lib.loc = legacy_lib, fields = "Version"),
    error = function(e) NA_character_
  )
  if (!is.na(current) && current == expected) {
    message("Already installed: ", package, " ", expected)
    next
  }
  message("Installing ", package, " ", expected, " from ", path)
  install.packages(path, lib = legacy_lib, repos = NULL, type = "source")
  actual <- tryCatch(
    packageDescription(package, lib.loc = legacy_lib, fields = "Version"),
    error = function(e) NA_character_
  )
  if (is.na(actual) || actual != expected) {
    stop(
      "Installation verification failed for ", package,
      "; requested ", expected,
      "; installed ", ifelse(is.na(actual), "<not installed>", actual)
    )
  }
}

message("All pinned packages installed successfully in ", legacy_lib)
