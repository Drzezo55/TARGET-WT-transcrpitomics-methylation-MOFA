library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(here)
library(DESeq2)
library("biomartr")
library("clusterProfiler")
library("tidyverse")
library("enrichplot")
suppressPackageStartupMessages(library("org.At.tair.db"))
library("biomaRt")  # only use to remove cache bug
rna_data <- readRDS("TARGET_WT_RNAseq_data.rds")
row_data <- assay(rna_data, "unstranded")
tpm_data <- assay(rna_data, "tpm_unstrand")
metadata <- colData(rna_data)
metadata <- as.data.frame(metadata)
# assign groups based on our hypothesis 

metadata <- metadata %>%
  mutate(
    ten_year_survival = case_when(
      is.na(vital_status) | (vital_status == "Dead" & is.na(days_to_death)) ~ NA_character_,
      vital_status == "Dead" & days_to_death <= 10*365 ~ "dead",
      vital_status == "Dead" & days_to_death > 10*365 ~ "alive",
      vital_status == "Alive" ~ "alive",
      TRUE ~ NA_character_
    ),
    ten_year_survival = factor(ten_year_survival, levels = c("alive", "dead"))
  )
metadata$ten_year_survival
ggplot(metadata, aes(x = ten_year_survival)) +
  geom_bar(position = "dodge") +
  labs(
    title = "ten_year_survival",
    x = "ten_year_survival",
    y = "Number of Patients"
  ) +
  theme_minimal()
sum(is.na(metadata$ten_year_survival))
metadata <- metadata[!is.na(metadata$ten_year_survival),]
ids <- readRDS("data/barcode.rds")
matched <- match(ids, metadata$sample)
sum(is.na(matched))
# Check which ones matched

valid <- !is.na(matched)
sum(valid)
matched_metadata <- metadata[matched[valid], ]
dim(matched_metadata)
sum(is.na(matched_metadata$sample))
coldata <- matched_metadata
length(coldata)
# Subset the data (example for rna_data) for MOI
rna_matched <- rna_data[, colnames(rna_data) %in% matched_metadata$barcode]
# SAVE THE FINAL data
row <-  assay(rna_matched, "unstranded")
tpm <- assay(rna_matched, "tpm_unstrand")
fpkm <- assay(rna_matched, "fpkm_unstrand")
View(row)
#
genes_clean_row <- sub("\\..*", "", row.names(row))  # Remove version
genes_clean_tpm <- sub("\\..*", "", row.names(tpm))
genes_clean_fpkm <- sub("\\..*", "", row.names(fpkm))

row.names(row) <- genes_clean_row
row.names(tpm) <- genes_clean_tpm
row.names(fpkm) <- genes_clean_fpkm
#
expr_barcodes <- colnames(row)  # or colnames(assay(rna_data))
metadata_ordered <- metadata[match(expr_barcodes, metadata$barcode), ]
all(expr_barcodes == metadata_ordered$sample_id)
coldata <- metadata_ordered
dds <- DESeqDataSetFromMatrix(countData=row, 
                              colData=coldata, 
                              design=~ten_year_survival)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
res <- results(dds)
dds <- estimateSizeFactors(dds)
plotMA(res, main = "DESeq2 MA Plot", ylim = c(-5, 5))

#
library(biomaRt)
genes <- row.names(dds)
# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)
#
library(ggrepel)
gene_mapping$ensembl_gene_id
# Safely filter top DE genes
top_genes <- res[!is.na(res$padj) & !is.na(res$log2FoldChange) &
                 res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
top_genes <- top_genes[order(top_genes$padj), ][1:30, ]
top_genes <- as.data.frame(top_genes)
top_genes <- top_genes %>%
  mutate(names = row.names(.)) %>%
  left_join(gene_mapping, by= c("names" = "ensembl_gene_id"))
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.8) +
  geom_text_repel(data = top_genes, aes(label = top_genes$hgnc_symbol), size = 3) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot with Top 30 Genes", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme_minimal()
# extract up, down genes 
# Convert to data frame and remove NAs
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)
res_df <- res_df %>%
  mutate(names = row.names(res_df)) %>%
  left_join(gene_mapping, by= c("names" = "ensembl_gene_id"))


# Upregulated: padj < 0.05 and log2FC > 1
up_genes <- res_df[res_df$padj < 0.05 & res_df$log2FoldChange > 1, ]

# Downregulated: padj < 0.05 and log2FC < -1
down_genes <- res_df[res_df$padj < 0.05 & res_df$log2FoldChange < -1, ]

# Save to files (optional)
write.csv(up_genes, "results/upregulated_genes.csv", row.names = TRUE)
write.csv(down_genes, "results/downregulated_genes.csv", row.names = TRUE)
# Already filtered
deg <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
write.csv(deg, "results/significant_DEGs.csv", row.names = TRUE)
#
library(clusterProfiler)
library(org.Hs.eg.db)
ego <- enrichGO(
  gene = deg$entrezgene_id,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)
dotplot(ego, showCategory = 20)
ora_analysis_bp <- pairwise_termsim(ego, method = "JC")
emapplot(ora_analysis_bp, color = "qvalue")
universe <- as.character(gene_mapping$entrezgene_id)
#kegg
# Use enrichKEGG for pathway analysis
ora_analysis_kegg_pathways <- enrichKEGG(
  gene          = as.character(deg$entrezgene_id),
  universe      = as.character(gene_mapping$entrezgene_id),
  organism      = "hsa",
  keyType       = "ncbi-geneid",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.01
)

# This should now work and produce a plot
dotplot(ora_analysis_kegg_pathways, 
    color = "qvalue", 
    showCategory = 10, 
    size = "Count")
#
library(DOSE)

do_enrich <- enrichDO(
  gene          = deg$entrezgene_id,
  ont           = "HDO",            # other option: DOLite
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)
dotplot(do_enrich, showCategory = 20, title = "Disease Ontology Enrichment")
#
install.packages("enrichR")
library(enrichR)
library(tmod)
dbs <- listEnrichrDbs()
View(dbs)  
genes <- deg$hgnc_symbol
enriched <- enrichr(genes, databases = c("MSigDB_Hallmark_2020"))
#
library(ggplot2)

df <- enriched[["MSigDB_Hallmark_2020"]]

summary(drug)
view(df)
top_terms <- df[1:10, ]
# Split the "Overlap" column into numerator and denominator
overlap_parts <- strsplit(top_terms$Overlap, "/")
# Extract as numeric vectors
numerator <- as.numeric(sapply(overlap_parts, `[`, 1))
denominator <- as.numeric(sapply(overlap_parts, `[`, 2))
# Calculate percentage overlap
top_terms$percent_overlap <- (numerator / denominator) * 100
ggplot(top_terms, aes(x = reorder(Term, percent_overlap), 
                      y = -log10(Adjusted.P.value), 
                      fill = percent_overlap)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Oncogens",
       x = "Pathway",
       y = "-log10(Adj P-value)",
       fill = "% Overlap") +
  theme_minimal(base_size = 13)
#
genes_of_interest <- c("MYOG", "CKM", "CAV3", "ACTN2", "HSPB8", "TNNC2", "MYBPH",
                       "COX6A2", "MYL4", "ACTA1", "MYH2", "MYH3", "ACTC1", "DES",
                       "SLN", "SGCA", "MYL2", "MYF6", "CASQ2", "MYH8", "TNNI2",
                       "APOD", "CASQ1", "CACNG1")
# Ensure gene names match the rownames (case-sensitive!)
tpm_df <- as.data.frame(tpm)
tpm_df <- tpm_df %>%
  mutate(names = row.names(.)) %>%
  left_join(gene_mapping, by = c("names" = "ensembl_gene_id"))
common_genes <- intersect(genes_of_interest, tpm_df$hgnc_symbol)
length(common_genes)
# Subset the TPM data
tpm_subset <- tpm_df[match(common_genes, tpm_df$hgnc_symbol), ]
# Step 1: Remove row names (if any)
rownames(tpm_subset) <- NULL
# Step 2: Drop unwanted columns
tpm_clean <- dplyr::select(tpm_subset, -entrezgene_id, -names)
# Step 3: Convert 'hgnc_symbol' column to rownames
tpm_matrix <- tibble::column_to_rownames(tpm_clean, var = "hgnc_symbol")
dim(tpm_matrix)
#
# z-score by row (gene)
tpm_scaled <- t(apply(tpm_matrix, 1, scale))
colnames(tpm_scaled) <- colnames(tpm_matrix)
library(pheatmap)
matched_metadata$barcode
# Order columns by condition or survival
# Arrange metadata and extract properly ordered barcodes
matched_metadata_ordered <- matched_metadata %>%
  arrange(ten_year_survival)
ordered_samples <- matched_metadata_ordered$barcode
# Subset TPM data to match ordered barcodes
tpm_scaled <- tpm_scaled[, ordered_samples]
# Now annotation_col is ordered correctly
annotation_col <- matched_metadata_ordered[, "ten_year_survival", drop = FALSE]
rownames(annotation_col) <- matched_metadata_ordered$barcode
group_vector <- annotation_col$ten_year_survival
gap_index <- which(diff(as.numeric(factor(group_vector))) != 0)

pheatmap(tpm_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         gaps_col = gap_index,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Muscle-Related Genes Expression (TPM)",
         fontsize = 11,
         border_color = NA)
