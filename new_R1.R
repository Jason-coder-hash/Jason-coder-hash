setwd("C:/GNBF6010")
getwd()
library(Seurat)
library(Matrix)
library(cowplot)
library(EnsDb.Mmusculus.v79)
library(Signac)
library(dplyr)
library(hdf5r)
library(presto)
library(SeuratWrappers)
library(ggplot2)
library(tidyr)
library(MAST)
set.seed(1234)

#Trial
inputdata.10x.LNCaP_DMSO=Read10X_h5("LNCaP-DMSO/filtered_feature_bc_matrix.h5")
rna.LNCaP_DMSO <- CreateSeuratObject(counts = inputdata.10x.LNCaP_DMSO,project = "LNCaP_DMSO")

inputdata.10x.LNCaP_ENZ48=Read10X_h5("LNCaP-ENZ48/filtered_feature_bc_matrix.h5")
rna.LNCaP_ENZ48 <- CreateSeuratObject(counts = inputdata.10x.LNCaP_ENZ48,project = "LNCaP_ENZ48")

inputdata.10x.LNCaP_RESA=Read10X_h5("LNCaP-RESA/filtered_feature_bc_matrix.h5")
rna.LNCaP_RESA <- CreateSeuratObject(counts = inputdata.10x.LNCaP_RESA,project = "LNCaP_RESA")

inputdata.10x.LNCaP_RESB=Read10X_h5("LNCaP-RESB/filtered_feature_bc_matrix.h5")
rna.LNCaP_RESB <- CreateSeuratObject(counts = inputdata.10x.LNCaP_RESB,project = "LNCaP_RESB")

rna.combined.trial=merge(rna.LNCaP_DMSO,y=c(rna.LNCaP_ENZ48, rna.LNCaP_RESA, rna.LNCaP_RESB),add.cell.ids=c("LNCaP_DMSO","LNCaP_ENZ48","LNCaP_RESA","LNCaP_RESB"), project='rna.combined.trial')
saveRDS(rna.combined.trial,"rna.combined.trial.rds")


# Load combined data
rna.combined.trial <- readRDS("rna.combined.trial.rds")

# Split into individual samples
sample_list <- SplitObject(rna.combined.trial, split.by = "orig.ident")

# Define filtering function
filter_samples <- function(seurat_obj) {
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Apply sample-specific filters
  sample_id <- unique(seurat_obj$orig.ident)
  
  if(sample_id == "LNCaP_DMSO") {
    cat("Processing LNCaP (DMSO control)\n")
    subset(seurat_obj,
           nFeature_RNA > 3000 & nFeature_RNA < 7000 &
             nCount_RNA > 16000 & nCount_RNA < 50000 &
             percent.mt < 15)
  } else if(sample_id == "LNCaP_ENZ48") {
    cat("Processing LNCaP-ENZ48\n")
    subset(seurat_obj,
           nFeature_RNA > 1500 & nFeature_RNA < 5000 &
             nCount_RNA > 5000 & nCount_RNA < 25000 &
             percent.mt < 15)
  } else if(sample_id == "LNCaP_RESA") {
    cat("Processing RES-A\n")
    subset(seurat_obj,
           nFeature_RNA > 1500 & nFeature_RNA < 5000 &
             nCount_RNA > 5000 & nCount_RNA < 25000 &
             percent.mt < 17)
  } else if(sample_id == "LNCaP_RESB") {
    cat("Processing RES-B\n")
    subset(seurat_obj,
           nFeature_RNA > 1500 & nFeature_RNA < 5000 &
             nCount_RNA > 5000 & nCount_RNA < 25000 &
             percent.mt < 20)
  }
}

# Apply filters
filtered_samples <- lapply(sample_list, filter_samples)

# Normalization and variable features
filtered_samples <- lapply(filtered_samples, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Integration steps
features <- SelectIntegrationFeatures(object.list = filtered_samples)
anchors <- FindIntegrationAnchors(object.list = filtered_samples, 
                                  anchor.features = features)
combined.integrated <- IntegrateData(anchorset = anchors)

combined.integrated

# Post-integration processing
DefaultAssay(combined.integrated) <- "integrated"
set.seed(123)

combined.integrated.finale <- ScaleData(combined.integrated, verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 1.3)

p1 <- DimPlot(combined.integrated.finale, 
              reduction = "umap",
              label = TRUE,
              label.size = 4,
              pt.size = 0.5) +
  ggtitle("Seurat Clusters 0-15") +
  theme(plot.title = element_text(hjust = 0.5))
p1
ggsave("Clusters in DimPlot.png")

markers <- c("KLK3", "AR", "FOLH1", "HOXC6", "NKX3-1")
FeaturePlot(combined.integrated.finale,
            features = markers,
            reduction = "umap",
            ncol = 3,
            order = TRUE)
ggsave("Verify the presence of genes.png",width=12,height=9, units='in',dpi = 600)
library(MAST)

#Check the cell distribution
# Extract metadata and calculate cell counts
meta_data <- FetchData(combined.integrated.finale, vars = c("seurat_clusters", "orig.ident"))
count_data <- meta_data %>%
  group_by(seurat_clusters, orig.ident) %>%
  tally() %>%
  ungroup()

# Create grouped bar plot
library(ggplot2)
ggplot(count_data, aes(x = seurat_clusters, y = n, fill = orig.ident)) +
  geom_col(position = position_dodge(preserve = "single")) +
  geom_text(aes(label = n), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Seurat Cluster", y = "Cell Count", 
       title = "Cell Type Distribution Across Clusters",
       fill = "Sample") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cell_distribution.png")
#Induced: 5,6; Persistent:1,2,3,7,8,10,11,12,13,14,15; Resistant: 0,4,9


# 1. Cell Cycle Scoring
cc_genes <- cc.genes
s_genes <- cc_genes$s.genes
g2m_genes <- cc_genes$g2m.genes

seurat_obj_scRNA_4_com <- CellCycleScoring(
  object = combined.integrated.finale,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)
saveRDS(seurat_obj_scRNA_4_com, "seurat_obj_scRNA_4_com.rds")


scRNA_obj <- seurat_obj_scRNA_4_com
dna_replication_genes<-s_genes
mitotic_genes<-g2m_genes

#Appendix: Find Cell cycle 
# Set active assay to RNA
DefaultAssay(scRNA_obj) <- "RNA"
# Calculate expression profiles using alternative methods
# Calculate expression profiles using alternative methods
replication_expression <- Seurat::AverageExpression(
  object = scRNA_obj,
  assays = "RNA",
  features = dna_replication_genes,
  slot = "counts"  # Using count data for raw expression comparison
)[[1]]

mitosis_expression <- Seurat::AverageExpression(
  object = scRNA_obj,
  assays = "RNA",
  features = mitotic_genes,
  slot = "counts"
)[[1]]

# Create analysis matrix with distinct column names
phase_comparison_matrix <- data.frame(
  ClusterID = colnames(replication_expression),
  DNA_Replication = colMeans(replication_expression),
  Mitosis = colMeans(mitosis_expression)
)

phase_comparison_matrix$CellCycleState <- ifelse(
  # Mitotic Phase condition
  (phase_comparison_matrix$Mitosis / phase_comparison_matrix$DNA_Replication > 1) & 
    (phase_comparison_matrix$Mitosis > 1),
  "Mitotic_Phase",
  ifelse(
    # Replication Phase condition
    (phase_comparison_matrix$DNA_Replication / phase_comparison_matrix$Mitosis > 1) & 
      (phase_comparison_matrix$DNA_Replication > 1),
    "Replication_Phase",
    "Undetermined"
  )
)
phase_comparison_matrix

library(ggplot2)
library(stringr)  # For string manipulation

# Extract numeric values from ClusterID for ordering
phase_comparison_matrix <- phase_comparison_matrix %>%
  mutate(
    ClusterNumber = as.numeric(str_remove(ClusterID, "g")),
    ClusterID = factor(ClusterID, levels = ClusterID[order(ClusterNumber)])
  )

# Create the plot with properly ordered gX labels
cell_cycle_plot <- ggplot(phase_comparison_matrix, 
                          aes(x = ClusterID)) +
  geom_point(aes(y = DNA_Replication, color = "S Phase"), size = 4) +
  geom_point(aes(y = Mitosis, color = "G2/M Phase"), size = 4) +
  scale_color_manual(values = c("S Phase" = "#FFA500", "G2/M Phase" = "#800080")) +
  labs(title = "Cell Cycle Signature Comparison",
       x = "Cluster Identifier",
       y = "Average Expression Level",
       color = "Cell Cycle Phase") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(cell_cycle_plot)
cell_cycle_plot

