The project used a normalization method to track cancer progression throughout the multi-
stage of cancer development from early to advanced-stage using transcriptomic profiles.

The project focused on investigating cancer progression based on RNA-Seq gene expression.
A computational normalization method was employed that normalized advanced-stage
expression samples with early-stage RNA-Seq data. The normalization method corrects for
genes that display less expression variability in advanced-stage cancer samples but display a
high variability in early-stage cancer samples.

Additional bioinformatics analyses performed included; hierarchical clustering, feature
analysis (i.e. differential gene expression analysis and marker gene selection), predictive and
validation of marker genes, survival analysis, statistical analysis, and enrichment analyses.

To start a project create two gene-by-sample matrices (early-stage and advanced-stage).
Where, the rows and columns are composed of the gene expression profiles and samples,
respectively. The matrices are used as an input to the normalization method. Thereafter,
downstream analysis can be performed with normalized gene expression.
