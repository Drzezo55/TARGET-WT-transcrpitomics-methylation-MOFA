library(TCGAbiolinks)

# Step 1: Query RNA-seq gene expression data
query_rna <- GDCquery(
  project = "TARGET-WT",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Step 2: Download the RNA-seq data
GDCdownload(query_rna, method = "api", files.per.chunk = 5)

# Step 3: Prepare the RNA-seq data into a SummarizedExperiment object
rna_data <- GDCprepare(query_rna)

# Optional: Save data if needed
saveRDS(rna_data, "TARGET_WT_RNAseq_data.rds")
