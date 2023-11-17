# Packages
library("biomartr")
library("clusterProfiler")
library("tidyverse")
library("enrichplot")
library("biomaRt")

# Input gene list
df <- read.csv("Gene_list.tsv", sep = "\t", header = TRUE)

# Perform enrichment
ora_analysis_kegg <- enrichKEGG(gene = df$gene_list,
                                organism = "hsa",
                                keyType = "ncbi-geneid",
                                minGSSize = 5,
                                maxGSSize = 500,
                                pAdjustMethod = "fdr",
                                pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                use_internal_data = FALSE)

# Table: KEGG results
tab.kegg <- as.data.frame(ora_analysis_kegg)

# Draw dotplot
dotplot(ora_analysis_kegg, showCategory=5, color = "pvalue", title = "Enrichment of genes")

