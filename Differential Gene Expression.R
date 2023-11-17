library(limma)

# Input df
df <- read.csv2("L_ordered.tsv", header = TRUE, row.names = 1, sep = "\t", dec = ".")

# Experimental groups
groups <-  c(rep("C1",42), rep("C2",24), rep("C3",16))
cluster <- factor(groups, levels = c("C1","C2","C3"))

# Design object 
design <- model.matrix(~0+cluster)
colnames(design) <- c("C1","C2","C3")

# Contrast matrix
constrasts <- makeContrasts(C2-C1, C3-C1, C3-C2,levels = design)

# Fit with no intercept 
fit <- eBayes(contrasts.fit(lmFit(df, design = design), constrasts))

#topTable
tabl <- topTable(fit, number = Inf)
head(tabl)

# decideTest
decide<- decideTests(fit,method="separate",adjust.method="BH",p.value=0.05,lfc=0.5)
head(decide)
