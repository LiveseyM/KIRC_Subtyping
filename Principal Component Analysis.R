# Packages
library(factoextra)
library(FactoMineR)

# Input expression data with Cluster number in last column - match to sample
# Samples as rows and genes as columns
exp <- read.csv2("Top48_genes.csv", header = TRUE, row.names = 1, sep = ";", dec = ".")
exp[1:3,1:3]
exp$Cluster <- as.factor(exp$Cluster)

L <- scale(exp[,-49])
dim(L)

# PCA original data
res.pca <- PCA(L, graph = TRUE)
fviz_pca_ind(res.pca, label = "none", habillage = exp$Cluster, palette = c("#00AFBB", "#E7B800", "#FC4E07"), addEllipses = TRUE, ellipse.level = 0.80) +
  labs(title ="PCA")
