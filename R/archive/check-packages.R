############################################################################
# Check that all required packages are installed to run the repository.
#
#
#
#
#
############################################################################

# Required CRAN packages

cran_pkgs <- c(
  "DHARMa",
  "cowplot",
  "renv",
  "dplyr",
  "ggplot2",
  "ggtext",
  "glmmTMB",
  "gt",
  "knitr",
  "lme4",
  "lmerTest",
  "parallel",
  "purrr",
  "rmarkdown",
  "stats",
  "stringr",
  "tibble",
  "tidyverse"
)

cran_missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(cran_missing) > 0) {
  install.packages(cran_missing, dependencies = TRUE)
}


# Required Bioconductor packages

bioc_pkgs <- c("DESeq2", "edgeR", "clusterProfiler", "org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_missing <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(bioc_missing) > 0) {
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
}


# GitHub packages

github_pkgs <- c(
  "seqwrap",
  "lmerSeq",
  "glmmSeq",
  "ComplexUpset"
)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("seqwrap", quietly = TRUE)) {
  if (as.character(packageVersion("seqwrap")) != "0.6.1") {
    remotes::install_github("trainome/seqwrap", ref = "seqwrap-paper")
  }
}

if (!requireNamespace("lmerSeq", quietly = TRUE)) {
  remotes::install_github("stop-pre16/lmerSeq")
}

if (!requireNamespace("glmmSeq", quietly = TRUE)) {
  remotes::install_github("myles-lewis/glmmSeq")
}

if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
  remotes::install_github("krassowski/complex-upset")
}


# Final check

all_pkgs <- c(cran_pkgs, bioc_pkgs, github_pkgs)
missing <- setdiff(all_pkgs, rownames(installed.packages()))

if (length(missing) > 0) {
  stop(
    "The following packages are still missing: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
} else {
  message("✅ All required packages are installed.")
}
