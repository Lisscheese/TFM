# *********************************************************
# Dependency Management
# This is an optional script that helps to manage all the
# dependencies use in R environment (CRAN and BiocManager)
# *********************************************************
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.es.r-project.org")
  }
}
install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.eu.r-project.org")
}

cran_packages <- c("ggplot2", "plotly", "shiny", "fs", "shinythemes", "shinyFiles")
bioc_packages <- c("ShortRead", "Rsubread", "Rsamtools", "GenomicAlignments",
                   "VariantAnnotation", "vcfR", "biomaRt", "AnnotationDbi",
                   "Gviz", "org.Hs.eg.db", "clusterProfiler", "enrichplot", "GenomicRanges", "S4Vectors", "IRanges")

sapply(cran_packages, install_if_missing)
sapply(bioc_packages, install_bioc_if_missing)
