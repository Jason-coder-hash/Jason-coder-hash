
library(TCGAbiolinks)
library(survival)
library(tidyverse)
library(DESeq2)

# Define cancer type (e.g., "TCGA-BRCA" for breast cancer)
cancer_type <- "TCGA-BRCA"

# Download RNA-seq gene expression data (raw counts)
query <- GDCquery(
  project = cancer_type,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)
expr_data=data
# Download clinical data
clinical_data <- GDCquery_clinic(project = cancer_type, type = "clinical")

# Extract counts matrix and normalize (e.g., DESeq2 normalization)
count_matrix <- assay(expr_data)
metadata <- colData(expr_data)

# Optional: Filter lowly expressed genes (e.g., keep genes with ≥10 counts in ≥10% of samples)
keep <- rowSums(count_matrix >= 10) >= 0.1*ncol(count_matrix)
count_matrix <- count_matrix[keep, ]

# Normalize using DESeq2 (or convert to log2(FPKM+1) if using FPKM)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~1)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

view(clinical_data)
# Subset relevant columns (adjust based on your dataset)
clinical_clean <- clinical_data %>%
  select(
    patient_id = submitter_id,
    survival_time = days_to_last_follow_up,
    vital_status,
    drug_name = treatments_pharmaceutical_therapeutic_agents,  # Verify column name
    age = age_at_index,
    stage = ajcc_pathologic_stage
  ) %>%
  mutate(
    survival_time = as.numeric(survival_time),
    status = ifelse(vital_status == "Dead", 1, 0)  # 1=event, 0=censored
  )
is.na(clinical_clean$drug_name)
clinical_clean_with_drug=clinical_clean[!is.na(clinical_clean$drug_name),]
dim(clinical_clean_with_drug)
dim(clinical_clean)
View(clinical_clean_with_drug)


drug_of_interest="Clodronate Disodium"
treated_patients <- clinical_clean %>%
  filter(drug_name == drug_of_interest) %>%
  drop_na(survival_time, status)

# Match gene expression data to treated patients
common_patients <- intersect(colnames(normalized_counts), treated_patients$patient_id)
expr_treated <- normalized_counts[, common_patients]
clinical_treated <- treated_patients %>% filter(patient_id %in% common_patients)

clinical_treated
