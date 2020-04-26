###############################################################################
## Copyright (c) ZHANG Yang.
## Author: ZHANG Yang <zy0781@connect.hku.hk>
###############################################################################

###############################################################################
## Load Essential Packages
###############################################################################


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
  update.packages(ask = FALSE)
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", 
                     lib = normalizePath(paste(find.package("utils"), "/..", sep = "")))
  }
  BiocManager::install(pending_pkgs,
                       lib = normalizePath(paste(find.package("utils"), "/..", sep = "")))
}

for(pkg in essential_pkgs) library(pkg, character.only = TRUE)

###############################################################################
## Actual Code Start Here
###############################################################################



sessionInfo()
