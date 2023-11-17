### X and Y matrices
# Packages
library(magrittr)

setwd(path)

pa <- readr::read_tsv(file.path("data", "pa.tsv")) %>%
  tibble::column_to_rownames("site")
pe <- readr::read_tsv(file.path("data", "pe.tsv")) %>%
  tibble::column_to_rownames("site")
a <- readr::read_tsv(file.path("data", "A.tsv")) %>%
  tibble::column_to_rownames("gene")
e <- readr::read_tsv(file.path("data", "E.tsv")) %>%
  tibble::column_to_rownames("gene")


# Create 3D matrices
df_columns <- function(primary_site_matrix, 
                       gene_expr_matrix) {
  cols <- c()
  for (rn in rownames(primary_site_matrix)) {
    cols <- c(cols, paste0(names(gene_expr_matrix), paste0("_", rn)))
  }
  return(cols)
}

# multiply a primary site vector with a gene 
# expression matrix
mult_ps_ge <- function(primary_site_vector,
                       gene_expr_matrix,
                       suffix = "primary_site") {
  df <- sweep(gene_expr_matrix,
              MARGIN = 2, 
              STATS = primary_site_vector,
              FUN = '*')
  cols <- names(gene_expr_matrix) %>%
    paste0("_", suffix)
  
  return(df %>% 
           setNames(cols))
}

# Function to construct the 3D matrix X and Y
three_dim_matrix <- function(primary_site_matrix,
                             gene_expression_matrix) {
  X <- NULL
  for (i in 1:nrow(primary_site_matrix)) {
    primary_site_vector <- primary_site_matrix[i,] %>%
      as.numeric()
    suffix <- rownames(primary_site_matrix[i,])
    df <- mult_ps_ge(primary_site_vector,
                     gene_expression_matrix,
                     suffix)
    if(i == 1)
      X <- df
    else 
      X <- X %>%
        merge(df, by=0, all=TRUE) %>%
        tibble::column_to_rownames("Row.names")
  }
  return(X)
}

# Construct X and Y matrices
X <- three_dim_matrix(pa, a)
Y <- three_dim_matrix(pe, e)

# save matrices to files
readr::write_tsv(X %>% 
                   tibble::rownames_to_column("gene"),
                 file.path("results", "X_matrix.tsv"))
readr::write_tsv(Y %>% 
                   tibble::rownames_to_column("gene"),
                 file.path("results", "Y_matrix.tsv"))

### Normalization
X <- readr::read_tsv(file.path("results", "X_matrix.tsv")) %>%
  tibble::column_to_rownames("gene")
Y <- readr::read_tsv(file.path("results", "Y_matrix.tsv")) %>%
  tibble::column_to_rownames("gene")
pe <- readr::read_tsv(file.path("data", "pe.tsv")) %>%
  tibble::column_to_rownames("site")
e <- readr::read_tsv(file.path("data", "E.tsv")) %>%
  tibble::column_to_rownames("gene")

mean_primary_site <- function(primary_site_vector) {
  return(primary_site_vector %>%
           sum())
}

# Compute G_tissue matrix
G_tissue <- matrix(0, nrow = nrow(e), ncol = nrow(pe)) %>%
  as.data.frame() %>%
  set_rownames(rownames(e)) %>%
  set_names(rownames(pe))

for (i in rownames(G_tissue)) {
  for (j in names(G_tissue)) {
    m_l <- mean_primary_site(pe[j,])
    G_tissue[i,j] <- 1/m_l * Y[i,] %>%
      dplyr::select(ends_with(j)) %>%
      sum()
  }
}

# Drop columns with all 0 from X
X <- X %>%
  dplyr::select_if(colSums(.) != 0)

L <- matrix(0, nrow = nrow(X), ncol = ncol(X)) %>%
  as.data.frame() %>%
  set_rownames(rownames(X)) %>%
  set_names(colnames(X))

# normalization
for (i in rownames(X)) {
  for (j in colnames(X)) {
    p <- sub(".*_", "", j)
    L[i,j] <- log(X[i,j]/G_tissue[i,p]) 
  }
}

L

readr::write_tsv(L %>% 
                   tibble::rownames_to_column("gene"),
                 file.path("results", "Tissue_cor_matrix.tsv"))


