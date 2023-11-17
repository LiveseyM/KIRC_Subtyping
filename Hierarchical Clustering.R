library(dendextend)
library(flashClust)
library(dplyr)
library(data.table)
library(magrittr)
library(RColorBrewer)
library(wordspace)

Tumor <- read.csv("L.tsv", sep = "\t", header = TRUE, row.names = 1, dec = ".")
Tumor.t <- t(Tumor)

sampleDist <- Tumor.t %>% dist.matrix(method = "cosine") %>% as.dist

heatDendro <-
    sampleDist %>%
    flashClust(method = "ward") %>%
    as.dendrogram

k <- 
  heatDendro %>%
  find_k

colorpalDend <-
  brewer.pal(k$k,"Set2")

heatDendro %<>% 
  color_branches(k = k$k,
                 col = colorpalDend,
                 groupLabels = TRUE)
groups <-
  heatDendro %>%
  cutree(k = k$k)

plot(heatDendro)

