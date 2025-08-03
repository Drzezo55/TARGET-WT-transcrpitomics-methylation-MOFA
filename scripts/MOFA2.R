# Install MOFA2
# Step 1: Install dependencies
install.packages("devtools")
devtools::install_github("bioFAM/MOFA2")

# Step 2: Load MOFA2
library(MOFA2)
library(matrixStats)
# load exp
expr_var <- rowVars(as.matrix(tpm))
top_expr <- tpm[order(expr_var, decreasing = TRUE)[1:2000], ]
# load meth

mval_mat <- readRDS("TCGA-WT_mvalue_filtered.rds")
meth_var <- rowVars(as.matrix(mval_mat))
top_meth <- mval_mat[order(meth_var, decreasing = TRUE)[1:10000], ]
# Function to remove the last component (e.g., "-01R" or "-01D")
get_core_id <- function(x) sub("-01[A-Z]$", "", x)
# Extract core IDs
expr_ids <- get_core_id(colnames(top_expr))
meth_ids <- get_core_id(colnames(top_meth))
metadata_ids <- get_core_id(row.names(metadata))

# common ids
common_ids <- intersect(expr_ids, meth_ids)
length(common_ids)  # Check how many match
# Subset expression and rename columns to core ID
expr_common <- top_expr[, expr_ids %in% common_ids]
colnames(expr_common) <- expr_ids[expr_ids %in% common_ids]
# Subset methylation and rename columns to core ID
meth_common <- top_meth[, meth_ids %in% common_ids]
colnames(meth_common) <- meth_ids[meth_ids %in% common_ids]
# subset meta data 
metadata_common <- metadata[metadata_ids %in% common_ids, ]
rownames(metadata_common) <- metadata_ids[metadata_ids %in% common_ids]

# Reorder both matrices to the same sample order
expr_common <- expr_common[, common_ids]
meth_common <- meth_common[, common_ids]
metadata_common <- metadata_common[common_ids,]
dim(metadata_common)
# build the list for MOFA 
data_list <- list(
  Expression = as.matrix(expr_common),
  Methylation = as.matrix(meth_common)
)
#
mofa_object <- create_mofa(data_list)
#
# Create a MOFA object from the aligned data matrices
# The list should contain the data views (transcriptomics and methylation)
# in the correct format (e.g., matrices or data frames).
data_list <- list(
  transcriptomics = as.matrix(expr_common),
  methylation = as.matrix(meth_common)
)

MOFAobject <- create_mofa(data_list)

# Define the model options
# The number of factors (k) is a key parameter; start with a reasonable number.
model_options <- get_default_model_options(MOFAobject)
model_options$num_factors <- 15 # Example: set number of factors to 15

# Define the training options
# 'convergence_mode' can be "fast" or "medium" for faster runs.
train_options <- get_default_training_options(MOFAobject)
train_options$convergence_mode <- "fast" 
train_options$maxiter <- 1000

# Prepare the MOFA2 object for training
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  model_options = model_options,
  training_options = train_options
)

# Run the MOFA2 model training
# This is the main computational step.
MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE, outfile = "mofa2_model.hdf5")
# After you have run run_mofa() and have the trained MOFAobject:

# Ensure the rownames of the metadata match the sample names in the MOFA object.
# This should already be the case if you've followed the previous data alignment steps.
all.equal(samples_names(MOFAobject), rownames(metadata_common))
BiocManager::install("MOFA2")
library(MOFA2)
# Add the metadata to the MOFA object
MOFA2::samples_metadata(MOFAobject) <- metadata_common
remove.packages("MOFA2")
devtools::install_github("bioFAM/MOFA2", build_vignettes = FALSE, force = TRUE)
library(MOFA2)
plot_factor(MOFAobject, 
                    factors = "Factor1",  # Specify the factor you want to plot
                    color_by = "ten_year_survival",  # Replace with a column name from your metadata
                    dot_size = 2.5)
plot_weights(MOFAobject, 
             view = "transcriptomics", 
             factor = 1, 
             nfeatures = 20)
plot_variance_explained(MOFAobject)
plot_data_overview(MOFAobject)
factors <- get_factors(MOFAobject)$group1
cor(factors[,1], as.numeric(metadata_common$ten_year_survival))
top_genes <- get_weights(MOFAobject, 
                         view = "transcriptomics", 
                         factors = 3, 
                         as.data.frame = TRUE)

# Sort by absolute weight (most influential)
top_genes <- top_genes[order(abs(top_genes$value), decreasing = TRUE), ]

# Select top 100 genes (or any number you want)
top_gene_list <- head(top_genes$feature, 100)
list <- as.data.frame(top_gene_list)
list <- list %>%
  left_join(gene_mapping, by =c ("top_gene_list" = "ensembl_gene_id"))
go_results <- enrichDO(gene = list$entrezgene_id,
                       ont = "HDO",    # biological process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
barplot(go_results)
# View results
head(go_results)
#
top_dmrs <- get_weights(MOFAobject, 
                         view = "methylation", 
                         factors = 1, 
                         as.data.frame = TRUE)
plot_weights(MOFAobject, 
             view = "methylation", 
             factor = 1, 
             nfeatures = 20)
#
plot_factor_cor(MOFAobject)
#

clusters <- kmeans(factors, centers = 6)
metadata_common$cluster <- factor(clusters$cluster)
samples_metadata(MOFAobject)$cluster <- metadata_common$cluster
plot_factor(MOFAobject, factors = 1, color_by = "cluster")
#
library(uwot)
library(ggplot2)

factors <- get_factors(MOFAobject)$group1
umap_result <- umap(factors, n_neighbors = 15, min_dist = 0.1, metric = "cosine")
#
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$cluster <- metadata_common$cluster
umap_df$survival <- metadata_common$ten_year_survival
#
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster, shape = survival)) +
  geom_point(size = 3) +
  labs(title = "UMAP of MOFA Factors", color = "Cluster", shape = "Survival") +
  theme_minimal()

#
cor(factors[,3], as.numeric(metadata_common$ten_year_survival), method = "spearman")

#
plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "ten_year_survival")  # or any metadata of interest

#
top_genes <- get_weights(MOFAobject, 
                         view = "transcriptomics", 
                         factors = 3, 
                         as.data.frame = TRUE)

top_genes <- top_genes[order(abs(top_genes$value), decreasing = TRUE), ]
top_gene_list <- head(top_genes$feature, 100)
#
top_dmrs <- get_weights(MOFAobject, 
                         view = "methylation", 
                         factors = 3, 
                         as.data.frame = TRUE)

top_dmrs <- top_dmrs[order(abs(top_dmrs$value), decreasing = TRUE), ]
top_dmr_list <- head(top_dmrs$feature, 100)
#
plot_weights(MOFAobject, 
             view = "transcriptomics", 
             factor = 3, 
             nfeatures = 20)

plot_weights(MOFAobject, 
             view = "methylation", 
             factor = 3, 
             nfeatures = 20)
#
genes_enriched <- top_genes %>%
  left_join(gene_mapping, by = c("feature" = "ensembl_gene_id"))
#
library(ggplot2)

# Keep top N genes and ensure symbols are not missing
top_symbols <- genes_enriched %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "NA" & hgnc_symbol != "") %>%
  arrange(desc(abs(value))) %>%
  slice(1:20)

ggplot(top_symbols, aes(x = reorder(hgnc_symbol, abs(value)), y = value)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top Genes for Factor 3 (by Symbol)",
       x = "Gene Symbol",
       y = "MOFA Weight (Factor 3)") +
  theme_minimal()
view(top_symbols)
#
ex <- expr_common


