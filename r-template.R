###############################################################################
## Copyright (c) ZHANG Yang.
## Author: ZHANG Yang <zy0781@connect.hku.hk>
###############################################################################

###############################################################################
## Load Essential Packages
###############################################################################

update.packages(ask = FALSE)

essential_pkgs <- c(
  "devtools",
  "dplyr",
  "ggplot2",
  "openxlsx",
  "org.Mm.eg.db",
  "pathview",
  "pheatmap"
)

pending_pkgs <- unique(essential_pkgs[!essential_pkgs %in% installed.packages()])

if (length(pending_pkgs)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(pending_pkgs)
}

for(pkg in essential_pkgs) library(pkg, character.only = TRUE)

###############################################################################
## Actual Code Start Here
###############################################################################



devtools::session_info()
