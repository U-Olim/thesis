legacy_lib <- file.path(getwd(), "environment", "legacy_R343_library")
.libPaths(c(legacy_lib, .Library))
output_file <- file.path(getwd(), "environment", "sessionInfo_R343.txt")
sink(output_file, split = TRUE)
on.exit(sink(), add = TRUE)

expected <- c(
  quantreg = "5.34",
  hdm = "0.2.0",
  hqreg = "1.4",
  mvtnorm = "1.0-6",
  doSNOW = "1.0.16"
)

print(R.version.string)
print(packageVersion("quantreg", lib.loc = legacy_lib))
print(packageVersion("hdm", lib.loc = legacy_lib))
print(packageVersion("hqreg", lib.loc = legacy_lib))
print(packageVersion("mvtnorm", lib.loc = legacy_lib))
print(packageVersion("doSNOW", lib.loc = legacy_lib))

cat("\nLiteral CRAN DESCRIPTION versions:\n")
for (package in names(expected)) {
  actual <- packageDescription(package, lib.loc = legacy_lib, fields = "Version")
  cat(sprintf("%-10s %s\n", package, actual))
  if (actual != expected[[package]]) {
    stop(package, ": expected ", expected[[package]], ", found ", actual)
  }
}

cat("\n")
for (package in names(expected)) {
  library(package, character.only = TRUE, lib.loc = legacy_lib)
}
print(sessionInfo())
