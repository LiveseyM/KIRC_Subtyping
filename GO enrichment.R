# Packages
library("biomartr")
library("clusterProfiler")
library("enrichplot")
library("tidyverse")
library("enrichplot")
library("biomaRt")  

# Input file
df <- read.csv("Gene_list.tsv", sep = "\t", header = TRUE)

# gene ontology = "all" includes CC & MF & BP
ego <- enrichGO(df$gene_list, pvalueCutoff = 0.05,OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE)

# Table: GO results
tab.ego <- as.data.frame(ego)

# Draw dotplot
dotplot(ego, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

