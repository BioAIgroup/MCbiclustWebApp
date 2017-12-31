library(MCbiclust)
library(ggplot2)
library(gplots)
library(dplyr)
library(gProfileR)
library(MASS)
library(devtools)

data(CCLE_small)
data(Mitochondrial_genes)
mito.loc <- which(row.names(CCLE_small) %in% Mitochondrial_genes)
CCLE.mito <- CCLE_small[mito.loc,]


set.seed(102)
CCLE.seed <- FindSeed(gem = CCLE.mito,
                      seed.size = 10,
                      iterations = 1000,
                      messages = 1000)
