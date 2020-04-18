###############################################################################
## Copyright (c) ZHANG Yang.
## Author: ZHANG Yang <zy0781@connect.hku.hk>
###############################################################################

###############################################################################
## Load Essential Packages
###############################################################################

essential_pkgs <- c(
  "Biobase",
  "clusterProfiler",
  "devtools",
  "dplyr",
  "GEOquery",
  "ggplot2",
  "limma",
  "openxlsx",
  "org.Mm.eg.db",
  "pathview",
  "pheatmap",
  "stringr"
)

pending_pkgs <- unique(essential_pkgs[!essential_pkgs %in% installed.packages()])

if (length(pending_pkgs)) {
  update.packages(ask = FALSE)
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(pending_pkgs)
}

for(pkg in essential_pkgs) {
  library(pkg, character.only = TRUE)
}

###############################################################################
## Actual Code Start Here
###############################################################################


#' load series and platform data from GEO
#' getGEO(): An object of the appropriate class (GDS, GPL, GSM, or GSE) is returned. 
#'  If the GSEMatrix option is used, then a list of *ExpressionSet* objects is returned, 
#'  one for each SeriesMatrix file associated with the GSE accesion.
gset <- getGEO("GSE63514", GSEMatrix =TRUE, AnnotGPL=TRUE)


if (length(gset) > 1) {
  idx <- grep("GPL570", attr(gset, "names"))
} else {
  idx <- 1 # index
}

gset <- gset[[idx]]


#' make proper column names to match toptable 
#' fvarLabels():
#'  Retrieve information on features recorded in eSet-derived classes.
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("000000000000000000000000XXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "1111111111111111111111111111")
sml <- c()
for (i in 1:nchar(gsms)) {
  sml[i] <- substr(gsms,i,i)
}


#' eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

#' log2 transform
#'  Limma expects data values to be in log space
#' exprs(): Retrieve expression data from eSets.
#' quantile(): produces sample quantiles corresponding to the given probabilities.
ex <- exprs(gset)

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

is.log.untransformed <- 
  (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (is.log.untransformed) { 
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) 
}


#' set up the data and proceed with analysis
#' model.matrix(object, ...): creates a design (or model) matrix
#'  ~ operator:
#'    y ~ model: y is modelled by a linear predictor specified symbolically by model
#'  + operators:
#'  : operators:
#'    The terms themselves consist of variable and factor names separated by : operators
#'  * operator:
#'    a*b interpreted as a+b+a:b
#'  ^ operator: indicates crossing to the specified degree
#'    (a+b+c)^2 is identical to (a+b+c)*(a+b+c)
#'  %in% operator: the terms on its left are nested within those on the right
#'    a + b %in% a interpreted as a + a:b
sml <- paste("G", sml, sep="")    # set group names

fl <- as.factor(sml)

gset$description <- fl

design <- model.matrix(~ description + 0, gset)

colnames(design) <- levels(fl)

#' lmFit(): Fit linear model for each gene given a series of arrays
#'  returns An MArrayLM object containing the result of the fits.
fit <- lmFit(gset, design)

#' makeContrasts(): Construct Matrix of Custom Contrasts
#' model.matrix(object, data = environment(object),
#'              contrasts.arg = NULL, xlev = NULL, ...)
#'  returns Matrix which columns corresponding to contrasts.
cont.matrix <- makeContrasts(G1-G0, levels=design)

#' contrasts.fit(): Compute Contrasts from Linear Model Fit
#' contrasts.fit(fit, contrasts=NULL, coefficients=NULL)
fit2 <- contrasts.fit(fit, cont.matrix)

#' eBayes(): Empirical Bayes Statistics for Differential Expression
#' eBayes(fit, proportion = 0.01, stdev.coef.lim = c(0.1,4),
#'        trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05,0.1))
fit2 <- eBayes(fit2, 0.01)

#' topTable(): Table of Top Genes from Linear Model Fit
#' topTable(fit, coef=NULL, number=10, genelist=fit$genes, adjust.method="BH",
#'          sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

write.csv(tT, file="./expression-table.csv", row.names=F, sep="\t")
getwd()

################################################################

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("Normal","Cancer")

# set parameters and draw the plot

palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))

dev.new(width=4+dim(gset)[[2]]/5, height=6)

par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))

title <- paste ("GSE63514", '/', annotation(gset), " selected samples", sep ='')

boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)

legend("topleft", labels, fill=palette(), bty="n")

sessionInfo()
