# Package
library("ggpubr")

# Input file
bp <- read.table("Gene_expression.csv", sep=";", header = TRUE, dec = ".")
head(bp)

plot_bp <- ggboxplot(bp, x = "Cluster", y = "Gene_name", title = "Gene_name",
              ylab = "Normalized Gene Expression", xlab = "Risk subcategory", add = "jitter", legend = "none",
              fill = "Cluster", palette = c("#00AFBB","#E7B800", "#FC4E07" )) 

# Draw boxplot
plot_bp

# Statistical analysis
compare_means(Gene_name ~ Cluster,  data = bp)

method <- "anova" # one of "anova" or "kruskal.test"
my_comparisons <- list( c("SS", "IS"), c("SS", "LS"), c("IS", "LS") )

print(
  bp + stat_compare_means(comparisons = my_comparisons, method = method, label = "p.format")) 



