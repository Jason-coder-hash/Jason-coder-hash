setwd("C:/GNBF6010")
getwd()
library(ggpubr)
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
library(Matrix)
library(CytoTRACE2)
# Load required packages
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggsignif)
library(forcats)
library(velocyto.R)
library(SeuratDisk)
set.seed(123)

# Load required libraries
library(Signac)
library(EnsDb.Hsapiens.v86)
data("blacklist_hg38_unified")

LNCaP_ENZ168_dge_data <- read.delim("LNCaP-ENZ168/LNCaP-ENZ168_gene_exon_tagged.dge.all.txt", 
                       header = TRUE, row.names = 1, check.names = FALSE)

LNCaP_ENZ168_counts <- as(as.matrix(LNCaP_ENZ168_dge_data), "dgCMatrix")
rna.LNCaP_ENZ168 <- CreateSeuratObject(
  counts = LNCaP_ENZ168_counts,
  project = "LNCaP_ENZ168",  # Custom project name
  min.cells = 3,             # Keep genes detected in ≥3 cells
  min.features = 200         # Keep cells with ≥200 detected genes
)
rna.LNCaP_ENZ168
saveRDS(rna.LNCaP_ENZ168,"rna.LNCaP_ENZ168.rds")


inputdata.10x.LNCaP_DMSO=Read10X_h5("LNCaP-DMSO/filtered_feature_bc_matrix.h5")
rna.LNCaP_DMSO <- CreateSeuratObject(counts = inputdata.10x.LNCaP_DMSO,project = "LNCaP_DMSO",min.cells = 3,min.features = 200)

inputdata.10x.LNCaP_ENZ48=Read10X_h5("LNCaP-ENZ48/filtered_feature_bc_matrix.h5")
rna.LNCaP_ENZ48 <- CreateSeuratObject(counts = inputdata.10x.LNCaP_ENZ48,project = "LNCaP_ENZ48",min.cells = 3,min.features = 200)

inputdata.10x.LNCaP_RESA=Read10X_h5("LNCaP-RESA/filtered_feature_bc_matrix.h5")
rna.LNCaP_RESA <- CreateSeuratObject(counts = inputdata.10x.LNCaP_RESA,project = "LNCaP_RESA",min.cells = 3,min.features = 200)

inputdata.10x.LNCaP_RESB=Read10X_h5("LNCaP-RESB/filtered_feature_bc_matrix.h5")
rna.LNCaP_RESB <- CreateSeuratObject(counts = inputdata.10x.LNCaP_RESB,project = "LNCaP_RESB",min.cells = 3,min.features = 200)


#Quality control

#LNCaP-ENZ168
rna.LNCaP_ENZ168[["percent.mt"]] <- PercentageFeatureSet(
  rna.LNCaP_ENZ168, 
  pattern = "^MT-"
)

# View distributions
summary(rna.LNCaP_ENZ168$nFeature_RNA)  # Genes per cell
summary(rna.LNCaP_ENZ168$nCount_RNA)    # UMIs per cell
summary(rna.LNCaP_ENZ168$percent.mt)    # Mitochondrial %

#LNCaP-DMSO
rna.LNCaP_DMSO[["percent.mt"]] <- PercentageFeatureSet(
  rna.LNCaP_DMSO, 
  pattern = "^MT-"
)

# View distributions
summary(rna.LNCaP_DMSO$nFeature_RNA)  # Genes per cell
summary(rna.LNCaP_DMSO$nCount_RNA)    # UMIs per cell
summary(rna.LNCaP_DMSO$percent.mt)    # Mitochondrial %

#LNCaP-ENZ48

rna.LNCaP_ENZ48[["percent.mt"]] <- PercentageFeatureSet(
  rna.LNCaP_ENZ48, 
  pattern = "^MT-"
)

# View distributions
summary(rna.LNCaP_ENZ48$nFeature_RNA)  # Genes per cell
summary(rna.LNCaP_ENZ48$nCount_RNA)    # UMIs per cell
summary(rna.LNCaP_ENZ48$percent.mt)    # Mitochondrial %

#LNCaP-RESA

rna.LNCaP_RESA[["percent.mt"]] <- PercentageFeatureSet(
  rna.LNCaP_RESA, 
  pattern = "^MT-"
)

# View distributions
summary(rna.LNCaP_RESA$nFeature_RNA)  # Genes per cell
summary(rna.LNCaP_RESA$nCount_RNA)    # UMIs per cell
summary(rna.LNCaP_RESA$percent.mt)    # Mitochondrial %

#LNCaP-RESB

rna.LNCaP_RESB[["percent.mt"]] <- PercentageFeatureSet(
  rna.LNCaP_RESB, 
  pattern = "^MT-"
)

# View distributions
summary(rna.LNCaP_RESB$nFeature_RNA)  # Genes per cell
summary(rna.LNCaP_RESB$nCount_RNA)    # UMIs per cell
summary(rna.LNCaP_RESB$percent.mt)    # Mitochondrial %

rna.combined.research=merge(rna.LNCaP_DMSO,y=c(rna.LNCaP_ENZ48, rna.LNCaP_ENZ168, rna.LNCaP_RESA, rna.LNCaP_RESB),add.cell.ids=c("LNCaP_DMSO","LNCaP_ENZ48","LNCaP_ENZ168", "LNCaP_RESA","LNCaP_RESB"), project='rna.combined.real')
saveRDS(rna.combined.research, "rna.combined.research.rds")
rna.combined.research=readRDS("rna.combined.research.rds")
sample_list <- SplitObject(rna.combined.research, split.by = "orig.ident")
# Modified filtering function
filter_samples <- function(seurat_obj) {
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Get sample ID
  sample_id <- unique(seurat_obj$orig.ident)
  
  # Apply sample-specific filters
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
  } else if(sample_id == "LNCaP_ENZ168") {
    cat("Processing LNCaP-ENZ168\n")
    subset(seurat_obj,
           nFeature_RNA > 600 & nFeature_RNA < 5000 &
             nCount_RNA > 700 & nCount_RNA < 10000 &
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

# Apply filters (same as before)
filtered_samples <- lapply(sample_list, filter_samples)
filtered_samples
raw.filtered_samples=filtered_samples
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
DefaultAssay(combined.integrated)='integrated'

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
ggsave("Clusters in DimPlot research.png")
#Induced:7,9,14; Persistent: 1-4,8,10-13,15; Resistent:5,6; Intermediate: 0
seurat_object=readRDS("seurat_object.rds")
#Check the cell distribution
# Extract metadata and calculate cell counts
meta_data <- FetchData(seurat_object, vars = c("seurat_clusters", "orig.ident"))
count_data <- meta_data %>%
  group_by(seurat_clusters, orig.ident) %>%
  tally() %>%
  ungroup()

# Ensure seurat_clusters in count_data uses the correct order (0-15)
count_data$seurat_clusters <- factor(
  count_data$seurat_clusters, 
  levels = c(0:15)  # Explicitly set levels here
)
# Convert orig.ident to a factor with custom order
count_data$orig.ident <- factor(count_data$orig.ident, 
                                levels = c("LNCaP_DMSO", "LNCaP_ENZ48", "LNCaP_ENZ168", 
                                           "LNCaP_RESA", "LNCaP_RESB"))

# Create grouped bar plot
library(ggplot2)
# Create grouped bar plot (modified)
ggplot(count_data, aes(x = seurat_clusters, y = n, fill = orig.ident)) +
  geom_col(position = position_dodge(preserve = "single")) +  # Removed geom_text
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Seurat Cluster", y = "Cell Count", 
       title = "Cell Type Distribution Across Clusters",
       fill = "Sample") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(drop = FALSE)  # Ensures all clusters (0-15) show even if empty

ggsave("cell_distribution_research.png")


#Cell Cycle Correction

# Step 1: Load cell cycle gene lists
# These are built-in lists of canonical cell cycle markers
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

combined.integrated.finale_wth_cc=CellCycleScoring(combined.integrated.finale, 
                                                   s.features = s.genes, 
                                                   g2m.features = g2m.genes, 
                                                   set.ident = FALSE)

# 2. Visualize S.Score and G2M.Score by cluster
FeaturePlot(combined.integrated.finale_wth_cc,
            reduction = "umap",
            features = c("S.Score", "G2M.Score"),
            pt.size = 0.4,
            order = TRUE,
            label = TRUE)

ggsave("CellCycleScoring research.png")


# 3. Calculate average expression of marker genes by cluster
seurat_object=combined.integrated.finale_wth_cc
Idents(seurat_object) <- "seurat_clusters"
s_markers_avg <- AverageExpression(seurat_object, features = s.genes, slot = "data")
g2m_markers_avg <- AverageExpression(seurat_object, features = g2m.genes, slot = "data")
View(s_markers_avg)

# 4. Determine dominant phase for each cluster
s_score_by_cluster <- colMeans(s_markers_avg$RNA)
g2m_score_by_cluster <- colMeans(g2m_markers_avg$RNA)
s_score_by_cluster
g2m_score_by_cluster

# Create dataframe of results
cell_cycle_dominance <- data.frame(
  Cluster = 0:15,  # Assuming 16 clusters numbered 0-15
  S_Score = s_score_by_cluster,
  G2M_Score = g2m_score_by_cluster
)

# Add dominant phase column
cell_cycle_dominance$Dominant_Phase <- ifelse(
  cell_cycle_dominance$S_Score > cell_cycle_dominance$G2M_Score, 
  "S", "G2M"
)
# View results
print(cell_cycle_dominance)
#From the results: S is dominant in 11,13, while G2M dominant in 8, 10 and 12


# 5. Visualize specific marker genes across clusters
DotPlot(seurat_object, 
        features = c("MCM6", "CCNE1", "CCNA2", "CCNB1", "AURKA", "NUSAP1"), 
        cols = c("lightgrey", "blue")) + 
  RotatedAxis()
saveRDS(seurat_object,"seurat_object.rds")

#Previously on only 4 ATAC and scRNA sequencing:
library(Signac)
library(EnsDb.Hsapiens.v86)  # Human genome annotation
library(biovizBase)
#ATAC:
input_ATAC_LNCaP_DMSO=Read10X_h5("LNCaP-DMSO/LNCaP-DMSO_filtered_peak_bc_matrix.h5")
input_ATAC_LNCaP_ENZ48=Read10X_h5("LNCaP-ENZ48/LNCaP-ENZ48_filtered_peak_bc_matrix.h5")
input_ATAC_LNCaP_RESA=Read10X_h5("LNCaP-RESA/LNCaP-RESA_filtered_peak_bc_matrix.h5")
input_ATAC_LNCaP_RESB=Read10X_h5("LNCaP-RESB/LNCaP-RESB_filtered_peak_bc_matrix.h5")

frag.file.LNCaP_DMSO="LNCaP-DMSO/LNCaP-DMSO_fragments.tsv.gz"
frag.file.LNCaP_ENZ48="LNCaP-ENZ48/LNCaP-ENZ48_fragments.tsv.gz"
frag.file.LNCaP_RESA="LNCaP-RESA/LNCaP-RESA_fragments.tsv.gz"
frag.file.LNCaP_RESB="LNCaP-RESB/LNCaP-RESB_fragments.tsv.gz"

metadata_LNCaP_DMSO=read.csv("LNCaP-DMSO/LNCaP-DMSO_singlecell.csv",header=TRUE,row.names=1)
metadata_LNCaP_ENZ48=read.csv("LNCaP-ENZ48/LNCaP-ENZ48_singlecell.csv",header=TRUE,row.names=1)
metadata_LNCaP_RESA=read.csv("LNCaP-RESA/LNCaP-RESA_singlecell.csv",header=TRUE,row.names=1)
metadata_LNCaP_RESB=read.csv("LNCaP-RESB/LNCaP-RESB_singlecell.csv",header=TRUE,row.names=1)

# Create ChromatinAssay objects for each sample
create_atac_assay <- function(count_matrix, frag_path, sample_name) {
  granges_counts <- StringToGRanges(rownames(count_matrix), sep = c(":", "-"))
  granges_use <- seqnames(granges_counts) %in% standardChromosomes(granges_counts)
  count_matrix <- count_matrix[as.vector(granges_use), ]
  
  # Get ENCODE-compatible annotations
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"
  
  CreateChromatinAssay(
    counts = count_matrix,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = frag_path,
    min.cells = 10,
    annotation = annotations
  )
}

# Create assays for each sample
assay_lncap_dmso <- create_atac_assay(input_ATAC_LNCaP_DMSO, frag.file.LNCaP_DMSO)
assay_lncap_enz <- create_atac_assay(input_ATAC_LNCaP_ENZ48, frag.file.LNCaP_ENZ48)
assay_resa <- create_atac_assay(input_ATAC_LNCaP_RESA, frag.file.LNCaP_RESA)
assay_resb <- create_atac_assay(input_ATAC_LNCaP_RESB, frag.file.LNCaP_RESB)

# Create Seurat objects and add QC metrics
atac_seurat <- list(
  "LNCaP_DMSO" = CreateSeuratObject(assay_lncap_dmso, assay = "ATAC",meta.data=metadata_LNCaP_DMSO,project="LNCaP-DMSO"),
  "LNCaP_ENZ48" = CreateSeuratObject(assay_lncap_enz, assay = "ATAC",meta.data=metadata_LNCaP_ENZ48,project="LNCaP-ENZ48"),
  "RESA" = CreateSeuratObject(assay_resa, assay = "ATAC",meta.data=metadata_LNCaP_RESA, project="RESA"),
  "RESB" = CreateSeuratObject(assay_resb, assay = "ATAC",meta.data=metadata_LNCaP_RESB, project="RESB")
)



# Define QC thresholds matrix
qc_thresholds <- data.frame(
  Sample = c("LNCaP", "LNCaP-ENZ48", "RES-A", "RES-B"),
  total_frags_min = c(2000, 1000, 1000, 1000),
  total_frags_max = c(20000, 20000, 20000, 25000),
  frac_peaks_min = c(30, 30, 40, 30),
  blacklist_max = c(0.01, 0.01, 0.01, 0.01),
  nucleosome_max = c(9, 9, 8, 8),
  tss_min = c(2, 2, 2, 2)
)

# Create mapping between Seurat objects and QC thresholds
sample_mapping <- c(
  "LNCaP_DMSO" = "LNCaP",
  "LNCaP_ENZ48" = "LNCaP-ENZ48",
  "RESA" = "RES-A",
  "RESB" = "RES-B"
)

# Add QC metrics to each sample
for (sample_name in names(atac_seurat)) {
  obj <- atac_seurat[[sample_name]]
  
  # Calculate QC metrics
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj, fast = FALSE)
  
  # Calculate additional metrics
  obj$pct_reads_in_peaks <- obj$nCount_ATAC / obj$passed_filters * 100
  obj$blacklist_ratio <- FractionCountsInRegion(
    object = obj,
    assay = "ATAC",
    regions = blacklist_hg38_unified
  )
  
  atac_seurat[[sample_name]] <- obj
}

atac_seurat_filtered <- lapply(names(atac_seurat), function(sample_name) {
  obj <- atac_seurat[[sample_name]]
  qc <- qc_thresholds[qc_thresholds$Sample == sample_mapping[sample_name], ]
  
  subset(
    x = obj,
    subset = nCount_ATAC > qc$total_frags_min &
      nCount_ATAC < qc$total_frags_max &
      pct_reads_in_peaks > qc$frac_peaks_min &
      blacklist_ratio < qc$blacklist_max &
      nucleosome_signal < qc$nucleosome_max &
      TSS.enrichment > qc$tss_min
  )
})

atac_seurat_filtered

saveRDS(atac_seurat_filtered,"atac_seurat_filtered.rds")

atac_seurat_filtered=readRDS("atac_seurat_filtered.rds")

atac_combined=merge(atac_seurat_filtered[[1]],y=list(atac_seurat_filtered[[2]],atac_seurat_filtered[[3]],atac_seurat_filtered[[4]]))


atac_combined$dataset <- atac_combined$orig.ident

# Processing pipeline
atac_combined <- FindTopFeatures(atac_combined, min.cutoff = 10)
atac_combined <- RunTFIDF(atac_combined)
atac_combined <- RunSVD(atac_combined)
atac_combined <- RunUMAP(atac_combined, reduction = "lsi", dims = 2:30)

# Visualization
DimPlot(atac_combined, group.by = "dataset") + 
  ggtitle("Integrated ATAC-Seq Datasets") +
  theme_minimal()

saveRDS(atac_combined, "atac_combined_processed.rds")

library(harmony)
atac_combined=readRDS("atac_combined_processed.rds")

atac_combined <- RunHarmony(
  object = atac_combined,
  group.by.vars = "orig.ident", # Replace with your batch variable if different
  reduction.use = 'lsi',        # Explicitly use 'reduction.use' instead of 'reduction'
  assay.use = 'ATAC',
  project.dim = FALSE
)

atac_combined <- RunUMAP(
  object = atac_combined,
  reduction = 'harmony',
  dims = 1:50
)

atac_combined <- FindNeighbors(
  object = atac_combined,
  reduction = 'harmony',
  dims = 1:50
)

atac_combined <- FindClusters(
  object = atac_combined,
  algorithm = 4, # SLM algorithm
  resolution = 0.4
)
saveRDS(atac_combined,"atac_combined.rds")
atac_combined=readRDS("atac_combined.rds")
atac_combined$seurat_clusters
# Visualize results
DimPlot(atac_combined, reduction = "umap", group.by = "seurat_clusters",  
        label = TRUE,
        label.size = 4,
        pt.size = 0.5) + 
  ggtitle("Integrated ATAC-Seq Datasets")

# Fetch cluster and sample data from the ATAC-seq object
meta_data <- FetchData(atac_combined, vars = c("seurat_clusters", "orig.ident"))

# Count cells by cluster and sample
count_data <- meta_data %>%
  group_by(seurat_clusters, orig.ident) %>%
  tally() %>%
  ungroup()

# Make sure the sample identities are in the right order
count_data$orig.ident <- factor(count_data$orig.ident, 
                                levels = c("LNCaP-DMSO", "LNCaP-ENZ48", "RESA", "RESB"))

# Create grouped bar plot for raw counts
ggplot(count_data, aes(x = seurat_clusters, y = n, fill = orig.ident)) +
  geom_col(position = position_dodge(preserve = "single")) +
  geom_text(aes(label = n), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "ATAC Cluster", y = "Number of Cells", 
       title = "Cell Type Distribution Across ATAC Clusters",
       fill = "Sample") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Cell Distribution ATAC.png")




# For seurat_obj_scRNA_4_com
cells_to_check <- c("LNCaP_RESA_CACCTGTTCCAGAATC-1_3", "LNCaP_RESA_GCCTACTAGAACCCGA-1_3")
existing_cells <- Cells(seurat_obj_scRNA_4_com)[Cells(seurat_obj_scRNA_4_com) %in% cells_to_check]
existing_cells

existing_cells_alt <- Cells(seurat_object)[Cells(seurat_object) %in% cells_to_check]
existing_cells_alt
dim(atac_combined@meta.data)
dim(seurat_obj_scRNA_4_com@meta.data)

##Check the ATAC cell DAR #Rerun ATAC data
atac_combined=readRDS("atac_combined.rds")


atac_combined@meta.data[atac_combined@meta.data$seurat_clusters==8,]

DefaultAssay(atac_combined) <- "ATAC"
Idents(atac_combined) <- "seurat_clusters" # Verify cluster column name

dar_list <- list()
for (cluster in levels(atac_combined)) {
  dar_list[[cluster]] <- FindMarkers(
    object = atac_combined,
    ident.1 = cluster,
    test.use = "LR", # Logistic regression with LRT
    latent.vars = "nCount_ATAC", # Total peaks as latent variable
    min.pct = 0.1, # Minimum % cells with accessibility
    logfc.threshold = 0.25, # Minimum log-fold change
    only.pos = TRUE # Identify cluster-specific accessible regions
  )
}
dar_list$`1`
saveRDS(dar_list,"dar_list.rds")

gene_list_atac=GeneActivity(atac_combined)
saveRDS(gene_list_atac,"gene_list_atac.rds")

atac_combined[["RNA"]]=CreateAssayObject(counts = gene_list_atac)
atac_combined <- NormalizeData(
  object = atac_combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac_combined$nCount_RNA)
)

atac_combined[["RNA"]]

#Integration with RNA data
# Convert to character for manipulation
seurat_obj_scRNA_4_com$orig.ident <- as.character(seurat_obj_scRNA_4_com$orig.ident)
seurat_obj_scRNA_4_com$orig.ident[seurat_obj_scRNA_4_com$orig.ident=='RESA']
# Rename entries
seurat_obj_scRNA_4_com$orig.ident <- ifelse(
  seurat_obj_scRNA_4_com$orig.ident %in% c("LNCaP_RESA", "LNCaP_RESB"),
  gsub("LNCaP_", "", seurat_obj_scRNA_4_com$orig.ident),  # Remove prefix for RESA/RESB
  ifelse(
    seurat_obj_scRNA_4_com$orig.ident %in% c("LNCaP_DMSO", "LNCaP_ENZ48"),
    gsub("_", "-", seurat_obj_scRNA_4_com$orig.ident),  # Replace _ with - for DMSO/ENZ48
    seurat_obj_scRNA_4_com$orig.ident  # Leave others unchanged
  )
)

# Optional: Convert back to factor if needed
seurat_obj_scRNA_4_com$orig.ident <- factor(seurat_obj_scRNA_4_com$orig.ident)
seurat_obj_scRNA_4_com$orig.ident
saveRDS(seurat_obj_scRNA_4_com,"seurat_obj_scRNA_4_com.rds")
atac_combined$orig.ident=factor(atac_combined$orig.ident)

DefaultAssay(atac_combined) <- 'RNA'
seurat_obj_scRNA_4_com <- UpdateSeuratObject(seurat_obj_scRNA_4_com)
transfer.anchors <- FindTransferAnchors(
  reference = seurat_obj_scRNA_4_com,
  query = atac_combined,
  reduction = 'cca'
)

scRNA_ATAC_predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = seurat_obj_scRNA_4_com$seurat_clusters,
  weight.reduction = atac_combined[['lsi']],
  dims = 2:30
)

scRNA_ATAC_predicted.labels

library(pheatmap)

# Add the predictions to your scATAC-seq object
atac_combined <- AddMetaData(object = atac_combined, metadata = scRNA_ATAC_predicted.labels)
atac_combined@meta.data
# Create a cross-tabulation of the predicted labels vs. actual scATAC clusters
prediction.table <- table(
  predicted = atac_combined$predicted.id,
  actual = atac_combined$seurat_clusters
)

# Convert to proportions - normalize by column (scATAC cluster) totals
prop.table <- prop.table(prediction.table, margin = 2)

prop.table
# Convert row and column labels to numeric values
numeric_rows <- as.numeric(rownames(prop.table))
numeric_cols <- as.numeric(colnames(prop.table))

# Reorder rows and columns using numeric values
prop.table <- prop.table[order(numeric_rows), order(numeric_cols)]
prop.table

# Convert to matrix for plotting if needed
prop.matrix <- as.matrix(prop.table)

# Generate the heatmap
pheatmap(
  prop.matrix,
  color = colorRampPalette(c("white", "orange", "red", "darkred"))(100),
  main = "Proportion of scATAC cells in each cluster",
  fontsize = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  breaks = seq(0, 1, by = 0.01)
)


#ATAC
#Induced: 3; Persistent: 1, 2,4,5; Resistant: 6,7

#scRNA:
#Induced: 5,6; Persistent:1,2,3,7,8,10,11,12,13,14,15; Resistant: 0,4,9


#Alternative: We have seurat_object which should have 5 clusters of cells
#Because cluster 0 contains most of LNCaP-ENZ168, we can check which clusters it does concentrate on!
seurat_object=readRDS("seurat_object.rds")
transfer.anchors.seurat_object <- FindTransferAnchors(
  reference = seurat_obj_scRNA_4_com,
  query = seurat_object,
  reduction = 'cca'
)
seurat_object_predicted.labels <- TransferData(
  anchorset = transfer.anchors.seurat_object,
  refdata = seurat_obj_scRNA_4_com$seurat_clusters,
  weight.reduction = seurat_object[['pca']],
  dims = 1:30
)

library(pheatmap)

# Add the predictions to your scATAC-seq object
seurat_object <- AddMetaData(object = seurat_object, metadata = seurat_object_predicted.labels)
seurat_object@meta.data

# Create a cross-tabulation of the predicted labels vs. actual scATAC clusters
prediction.table.seurat_object <- table(
  predicted = seurat_object$predicted.id,
  actual = seurat_object$seurat_clusters
)

# Convert to proportions - normalize by column (scATAC cluster) totals
prop.table.seurat_object <- prop.table(prediction.table.seurat_object, margin = 2)

prop.table.seurat_object

# Convert row and column labels to numeric values
numeric_rows.seurat_object <- as.numeric(rownames(prop.table.seurat_object))
numeric_cols.seurat_object <- as.numeric(colnames(prop.table.seurat_object))

# Reorder rows and columns using numeric values
prop.table.seurat_object <- prop.table.seurat_object[order(numeric_rows.seurat_object), order(numeric_cols.seurat_object)]
prop.table.seurat_object

# Convert to matrix for plotting if needed
prop.matrix.seurat_object <- as.matrix(prop.table.seurat_object)

library(pheatmap)
library(grid)

# Create the pheatmap with enough margin space
p <- pheatmap(
  prop.matrix.seurat_object,
  color = colorRampPalette(c("white", "orange", "red", "darkred"))(100),
  fontsize = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  breaks = seq(0, 1, by = 0.01),
  silent = TRUE,  # Don't display yet
  treeheight_row = 0,  # Remove dendrogram space
  treeheight_col = 0
)

# Start with clean plotting space with more room for labels
grid.newpage()

# Create a specific viewport for the heatmap 
# with ample margins for subtitles
pushViewport(viewport(
  x = 0.5, y = 0.55,  # Center horizontally, push up vertically
  width = 0.75, height = 0.75,  # Reduce size to create space
  name = "heatmap_area"
))

# Draw the heatmap in the adjusted viewport
grid.draw(p$gtable)
popViewport()

# Add x-axis subtitle with much more space
grid.text("scRNA clusters of DMSO, ENZ48, ENZ168, RESA and RESB",
          x = 0.5,  # Center horizontally
          y = 0.08,  # Position well below the heatmap
          gp = gpar(fontsize = 11, fontface = "bold"),
          just = "center")

# Add y-axis subtitle with ample space
grid.text("Predicted scRNA clusters of DMSO, ENZ48, RESA, RESB",
          x = 0.07,  # Position well to the left
          y = 0.55,  # Center vertically
          rot = 90,  # Rotate text
          gp = gpar(fontsize = 11, fontface = "bold"),
          just = "center")

#As expected, clusters 0 and 4 are predicted on the 0 clusters of scRNA where most of ENZ168 embedded.

#ATAC-RNA clusters showed cluster 6 has a strong correlation with the RNA clusters

#We will first analysize the gene expression in cluster 0.
saveRDS(seurat_object, "seurat_object.rds")
saveRDS(seurat_object_predicted.labels,"seurat_object_predicted.labels.rds")
seurat_object



# Set cluster identities
Idents(seurat_object) <- "seurat_clusters"

# Find markers for cluster 0
de_cluster0 <- FindMarkers(
  object = seurat_object,
  ident.1 = 0,
  logfc.threshold = 0.25,    # Adjust based on biological relevance
  min.pct = 0.25,            # Require detection in ≥25% of cells
  test.use = "wilcox"        # Default Wilcoxon rank-sum test
)


# For cluster 0 in old object
Idents(seurat_obj_scRNA_4_com) <- "seurat_clusters"
old_cluster0_genes <- rownames(FindMarkers(seurat_obj_scRNA_4_com, ident.1 = 0))

# For cluster 4 in old object
old_cluster4_genes <- rownames(FindMarkers(seurat_obj_scRNA_4_com, ident.1 = 4))

#For cluster 0 in new object
new_cluster0_genes<-rownames(de_cluster0)
overlap_genes <- intersect(
  new_cluster0_genes,
  intersect(old_cluster0_genes, old_cluster4_genes)
)
overlap_genes
de_cluster0_overlap=de_cluster0[overlap_genes,]

# Filter genes with adjusted p-value < 1e-10
filtered_genes <- de_cluster0_overlap %>%
  filter(p_val_adj < 1e-10)

# Get top 20 upregulated genes
top20_up <- filtered_genes %>%
  arrange(desc(avg_log2FC)) %>%
  head(20) %>%
  mutate(Direction = "Upregulated")

# Get bottom 20 downregulated genes
bottom20_down <- filtered_genes %>%
  arrange(avg_log2FC) %>%
  head(20) %>%
  mutate(Direction = "Downregulated")

top20_up
bottom20_down
saveRDS(de_cluster0_overlap,"de_cluster0_overlap.rds")
saveRDS(top20_up,"top20_up.rds")
saveRDS(bottom20_down,"bottom20_down.rds")

seurat_object=readRDS("seurat_object.rds")
# Split by cell line

cell_lines <- SplitObject(seurat_object, split.by = "orig.ident")

#Separate different cell lines
#LNCaP-DMSO expression matrix and phenotype
cell_lines
# Extract the counts layer matrix (sparse dgCMatrix)
DMSO_counts_matrix <- LayerData(cell_lines$LNCaP_DMSO[["RNA"]], 
                           layer = "counts.LNCaP_DMSO.1")
DMSO_dense_matrix <- as.matrix(DMSO_counts_matrix)
DMSO_annotation_table <- data.frame(
  seurat_clusters = cell_lines$LNCaP_DMSO$seurat_clusters,
  row.names = names(cell_lines$LNCaP_DMSO$seurat_clusters)
)
DMSO_cytotrace=cytotrace2(DMSO_dense_matrix,species="human")
DMSO_cytotrace$seurat_clusters=DMSO_annotation_table$seurat_clusters
saveRDS(DMSO_cytotrace,'DMSO_cytotrace.rds')
DMSO_cytotrace=readRDS("DMSO_cytotrace.rds")
DMSO_extracted_data=DMSO_cytotrace[, c("CytoTRACE2_Score","seurat_clusters")]

# Extract the necessary data from your Seurat object
DMSO_cyto_clusters_data <- data.frame(
  CytoTRACE2_Score = DMSO_extracted_data$CytoTRACE2_Score,
  seurat_clusters = DMSO_extracted_data$seurat_clusters
)
DMSO_cyto_clusters_data

library(forcats)
# Calculate median scores and reorder clusters
DMSO_cyto_clusters_data <- DMSO_cyto_clusters_data %>%
  mutate(seurat_clusters = fct_reorder(seurat_clusters, 
                                       CytoTRACE2_Score, 
                                       .fun = median, 
                                       .desc = TRUE))
saveRDS(DMSO_cyto_clusters_data,"DMSO_cyto_clusters_data.rds")
# Generate ordered boxplot
ggplot(DMSO_cyto_clusters_data, 
       aes(x = seurat_clusters, 
           y = CytoTRACE2_Score, 
           fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1) # Improve cluster label readability
  ) +
  guides(fill = "none")



# Perform Wilcoxon test between clusters 13 and 12
clusters_to_compare <- c("13", "8")
comparison_data <- DMSO_cyto_clusters_data %>%
  filter(seurat_clusters %in% clusters_to_compare)

wilcox_result <- wilcox.test(
  CytoTRACE2_Score ~ seurat_clusters, 
  data = comparison_data,
  exact = FALSE
)

# Format p-value for display
p_value <- wilcox_result$p.value
p_formatted <- formatC(p_value, format = "e", digits = 2)
signif_symbol <- paste("p =", p_formatted)

# Create the boxplot with significance annotation positioned at the top
ggplot(DMSO_cyto_clusters_data, aes(x = seurat_clusters, y = CytoTRACE2_Score, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  guides(fill = "none") +
  geom_signif(
    comparisons = list(c("13", "8")),
    annotations = signif_symbol,
    y_position = 0.75,  # Positioned at the top of the graph
    tip_length = 0.01,
    vjust = -0.2,
    textsize = 3,
    fontface = "bold"
  )
ggsave('Cytotrace LNCaP-DMSO.png',height=9,width=12,dpi=72)
# Extract the counts layer matrix (sparse dgCMatrix)
ENZ48_counts_matrix <- LayerData(cell_lines$LNCaP_ENZ48[["RNA"]], 
                                 layer = "counts.LNCaP_ENZ48.2")
ENZ48_dense_matrix <- as.matrix(ENZ48_counts_matrix)
ENZ48_annotation_table <- data.frame(
  seurat_clusters = cell_lines$LNCaP_ENZ48$seurat_clusters,
  row.names = names(cell_lines$LNCaP_ENZ48$seurat_clusters)
)
ENZ48_cytotrace=cytotrace2(ENZ48_dense_matrix,species="human")
ENZ48_cytotrace$seurat_clusters=ENZ48_annotation_table$seurat_clusters
saveRDS(ENZ48_cytotrace,'ENZ48_cytotrace.rds')
ENZ48_extracted_data=ENZ48_cytotrace[,c("CytoTRACE2_Score","seurat_clusters")]

# Extract the necessary data from your Seurat object
ENZ48_cyto_clusters_data <- data.frame(
  CytoTRACE2_Score = ENZ48_extracted_data$CytoTRACE2_Score,
  seurat_clusters = ENZ48_extracted_data$seurat_clusters
)
# Calculate median scores and reorder clusters
ENZ48_cyto_clusters_data <- ENZ48_cyto_clusters_data %>%
  mutate(seurat_clusters = fct_reorder(seurat_clusters, 
                                       CytoTRACE2_Score, 
                                       .fun = median, 
                                       .desc = TRUE))

# Generate ordered boxplot
ggplot(ENZ48_cyto_clusters_data, 
       aes(x = seurat_clusters, 
           y = CytoTRACE2_Score, 
           fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1) # Improve cluster label readability
  ) +
  guides(fill = "none")

# Perform Wilcoxon test between clusters 13 and 12
clusters_to_compare <- c("13", "12")
comparison_data <- ENZ48_cyto_clusters_data %>%
  filter(seurat_clusters %in% clusters_to_compare)

wilcox_result <- wilcox.test(
  CytoTRACE2_Score ~ seurat_clusters, 
  data = comparison_data,
  exact = FALSE
)

# Format p-value for display
p_value <- wilcox_result$p.value
p_formatted <- formatC(p_value, format = "e", digits = 2)
signif_symbol <- paste("p =", p_formatted)

# Create the boxplot with significance annotation positioned at the top
ggplot(ENZ48_cyto_clusters_data, aes(x = seurat_clusters, y = CytoTRACE2_Score, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  guides(fill = "none") +
  geom_signif(
    comparisons = list(c("13", "12")),
    annotations = signif_symbol,
    y_position = 0.75,  # Positioned at the top of the graph
    tip_length = 0.01,
    vjust = -0.2,
    textsize = 3,
    fontface = "bold"
  )
ggsave('Cytotrace LNCaP-ENZ48.png',height=9,width=12,dpi=72)


#LNCaP-ENZ168
# Extract the counts layer matrix (sparse dgCMatrix)
ENZ168_counts_matrix <- LayerData(cell_lines$LNCaP_ENZ168[["RNA"]], 
                                  layer = "counts.LNCaP_ENZ168.3")
ENZ168_dense_matrix <- as.matrix(ENZ168_counts_matrix)
ENZ168_annotation_table <- data.frame(
  seurat_clusters = cell_lines$LNCaP_ENZ168$seurat_clusters,
  row.names = names(cell_lines$LNCaP_ENZ168$seurat_clusters)
)
ENZ168_cytotrace=cytotrace2(ENZ168_dense_matrix,species="human")
ENZ168_cytotrace$seurat_clusters=ENZ168_annotation_table$seurat_clusters
saveRDS(ENZ168_cytotrace,'ENZ168_cytotrace.rds')
ENZ168_cytotrace=readRDS("ENZ168_cytotrace.rds")
ENZ168_extracted_data=ENZ168_cytotrace[, c("CytoTRACE2_Score","seurat_clusters")]

# Extract the necessary data from your Seurat object
ENZ168_cyto_clusters_data <- data.frame(
  CytoTRACE2_Score = ENZ168_extracted_data$CytoTRACE2_Score,
  seurat_clusters = ENZ168_extracted_data$seurat_clusters
)
library(forcats)
# Calculate median scores and reorder clusters
ENZ168_cyto_clusters_data <- ENZ168_cyto_clusters_data %>%
  mutate(seurat_clusters = fct_reorder(seurat_clusters, 
                                       CytoTRACE2_Score, 
                                       .fun = median, 
                                       .desc = TRUE))
saveRDS(ENZ168_cyto_clusters_data,"ENZ168_cyto_clusters_data.rds")
# Generate ordered boxplot
ggplot(ENZ168_cyto_clusters_data, 
       aes(x = seurat_clusters, 
           y = CytoTRACE2_Score, 
           fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1) # Improve cluster label readability
  ) +
  guides(fill = "none")

# Perform Wilcoxon test between clusters 13 and 8
clusters_to_compare <- c("13", "8")
comparison_data <- ENZ168_cyto_clusters_data %>%
  filter(seurat_clusters %in% clusters_to_compare)

wilcox_result <- wilcox.test(
  CytoTRACE2_Score ~ seurat_clusters, 
  data = comparison_data,
  exact = FALSE
)

# Format p-value for display
p_value <- wilcox_result$p.value
p_formatted <- formatC(p_value, format = "e", digits = 2)
signif_symbol <- paste("p =", p_formatted)

# Create the boxplot with significance annotation positioned at the top
ggplot(ENZ168_cyto_clusters_data, aes(x = seurat_clusters, y = CytoTRACE2_Score, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  guides(fill = "none") +
  geom_signif(
    comparisons = list(c("13", "8")),
    annotations = signif_symbol,
    y_position = 0.75,  # Positioned at the top of the graph
    tip_length = 0.01,
    vjust = -0.2,
    textsize = 3,
    fontface = "bold"
  )


#LNCaP-RESA
# Extract the counts layer matrix (sparse dgCMatrix)
RESA_counts_matrix <- LayerData(cell_lines$LNCaP_RESA[["RNA"]], 
                                layer = "counts.LNCaP_RESA.4")
RESA_dense_matrix <- as.matrix(RESA_counts_matrix)
RESA_annotation_table <- data.frame(
  seurat_clusters = cell_lines$LNCaP_RESA$seurat_clusters,
  row.names = names(cell_lines$LNCaP_RESA$seurat_clusters)
)
RESA_cytotrace=cytotrace2(RESA_dense_matrix,species="human")
RESA_cytotrace$seurat_clusters=RESA_annotation_table$seurat_clusters
saveRDS(RESA_cytotrace,'RESA_cytotrace.rds')
RESA_cytotrace=readRDS("RESA_cytotrace.rds")
RESA_extracted_data=RESA_cytotrace[, c("CytoTRACE2_Score","seurat_clusters")]

# Extract the necessary data from your Seurat object
RESA_cyto_clusters_data <- data.frame(
  CytoTRACE2_Score = RESA_extracted_data$CytoTRACE2_Score,
  seurat_clusters = RESA_extracted_data$seurat_clusters
)
library(forcats)
# Calculate median scores and reorder clusters
RESA_cyto_clusters_data <- RESA_cyto_clusters_data %>%
  mutate(seurat_clusters = fct_reorder(seurat_clusters, 
                                       CytoTRACE2_Score, 
                                       .fun = median, 
                                       .desc = TRUE))
saveRDS(RESA_cyto_clusters_data,"RESA_cyto_clusters_data.rds")
# Generate ordered boxplot
ggplot(RESA_cyto_clusters_data, 
       aes(x = seurat_clusters, 
           y = CytoTRACE2_Score, 
           fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1) # Improve cluster label readability
  ) +
  guides(fill = "none")


# Perform Wilcoxon test between clusters 13 and 8
clusters_to_compare <- c("13", "8")
comparison_data <- RESA_cyto_clusters_data %>%
  filter(seurat_clusters %in% clusters_to_compare)

wilcox_result <- wilcox.test(
  CytoTRACE2_Score ~ seurat_clusters, 
  data = comparison_data,
  exact = FALSE
)

# Format p-value for display
p_value <- wilcox_result$p.value
p_formatted <- formatC(p_value, format = "e", digits = 2)
signif_symbol <- paste("p =", p_formatted)

# Create the boxplot with significance annotation positioned at the top
ggplot(RESA_cyto_clusters_data, aes(x = seurat_clusters, y = CytoTRACE2_Score, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  guides(fill = "none") +
  geom_signif(
    comparisons = list(c("13", "8")),
    annotations = signif_symbol,
    y_position = 0.75,  # Positioned at the top of the graph
    tip_length = 0.01,
    vjust = -0.2,
    textsize = 3,
    fontface = "bold"
  )
ggsave('Cytotrace LNCaP-RESA.png',height=9,width=12,dpi=72)

# Extract the counts layer matrix (sparse dgCMatrix)
RESB_counts_matrix <- LayerData(cell_lines$LNCaP_RESB[["RNA"]], 
                                layer = "counts.LNCaP_RESB.5")
RESB_dense_matrix <- as.matrix(RESB_counts_matrix)
RESB_annotation_table <- data.frame(
  seurat_clusters = cell_lines$LNCaP_RESB$seurat_clusters,
  row.names = names(cell_lines$LNCaP_RESB$seurat_clusters)
)
RESB_cytotrace=cytotrace2(RESB_dense_matrix,species="human")
RESB_cytotrace$seurat_clusters=RESB_annotation_table$seurat_clusters
saveRDS(RESB_cytotrace,'RESB_cytotrace.rds')
RESB_cytotrace=readRDS("RESB_cytotrace.rds")
RESB_extracted_data=RESB_cytotrace[, c("CytoTRACE2_Score","seurat_clusters")]

# Extract the necessary data from your Seurat object
RESB_cyto_clusters_data <- data.frame(
  CytoTRACE2_Score = RESB_extracted_data$CytoTRACE2_Score,
  seurat_clusters = RESB_extracted_data$seurat_clusters
)
library(forcats)
# Calculate median scores and reorder clusters
RESB_cyto_clusters_data <- RESB_cyto_clusters_data %>%
  mutate(seurat_clusters = fct_reorder(seurat_clusters, 
                                       CytoTRACE2_Score, 
                                       .fun = median, 
                                       .desc = TRUE))
saveRDS(RESB_cyto_clusters_data,"RESB_cyto_clusters_data.rds")
# Generate ordered boxplot
ggplot(RESB_cyto_clusters_data, 
       aes(x = seurat_clusters, 
           y = CytoTRACE2_Score, 
           fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1) # Improve cluster label readability
  ) +
  guides(fill = "none")

# Perform Wilcoxon test between clusters 13 and 12
clusters_to_compare <- c("13", "12")
comparison_data <- RESB_cyto_clusters_data %>%
  filter(seurat_clusters %in% clusters_to_compare)

wilcox_result <- wilcox.test(
  CytoTRACE2_Score ~ seurat_clusters, 
  data = comparison_data,
  exact = FALSE
)

# Format p-value for display
p_value <- wilcox_result$p.value
p_formatted <- formatC(p_value, format = "e", digits = 2)
signif_symbol <- paste("p =", p_formatted)

# Create the boxplot with significance annotation positioned at the top
ggplot(RESB_cyto_clusters_data, aes(x = seurat_clusters, y = CytoTRACE2_Score, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(
    title = "CytoTRACE2_Score Distribution by Seurat Clusters",
    x = "Seurat Clusters (Ranked by Median Score)",
    y = "CytoTRACE2_Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  guides(fill = "none") +
  geom_signif(
    comparisons = list(c("13", "12")),
    annotations = signif_symbol,
    y_position = 0.75,  # Positioned at the top of the graph
    tip_length = 0.01,
    vjust = -0.2,
    textsize = 3,
    fontface = "bold"
  )



#Appendix
# Create UMAP visualization for S-phase
p1 <- FeaturePlot(seurat_object, 
                  features = "S.Score", 
                  reduction = "umap",
                  cols = c("lightgrey", "yellow", "red"),
                  order = TRUE,
                  label=T) +
  ggtitle("S-phase") +
  labs(color = "Average expression score")

# Create UMAP visualization for G2M-phase
p2 <- FeaturePlot(seurat_object, 
                  features = "G2M.Score", 
                  reduction = "umap",
                  cols = c("lightgrey", "yellow", "red"),
                  order = TRUE,
                  label=T) +
  ggtitle("G2M-phase") +
  labs(color = "Average expression score")
p1+p2
ggsave('S-phase, G2M-phase research.png',height=9,width=12,dpi=72)

# Get all the current identity levels
current_levels <- levels(seurat_object$orig.ident)

# Create new order with LNCaP_ENZ48 before LNCaP_ENZ168
# Make sure to include ALL conditions to avoid data loss
new_levels <- c("LNCaP_DMSO", "LNCaP_ENZ48", "LNCaP_ENZ168", "LNCaP_RESA", "LNCaP_RESB")

# Reorder factor levels (this won't lose any cells)
seurat_object$orig.ident <- factor(seurat_object$orig.ident, levels = new_levels)

# Recreate your plot with the new order
p3 <- DimPlot(seurat_object, 
              split.by = 'orig.ident',
              reduction = "umap",
              label = TRUE,
              label.size = 4,
              pt.size = 0.5) +
  ggtitle("Seurat Clusters 0-15 (each cell line)") +
  theme(plot.title = element_text(hjust = 0.5))
p3

#RNA velocity (Problematic)

#LNCaP_DMSO
cell_lines$LNCaP_DMSO
DefaultAssay(cell_lines$LNCaP_DMSO)='RNA'
trial <- RunVelocity(
  object = cell_lines$LNCaP_DMSO,
  spliced = "RNA",  # Name of assay containing spliced counts
  unspliced = "RNA",  # Name of assay containing unspliced counts
  deltaT = 1,
  kCells = 25,
  fit.quantile = 0.02
)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = cell_lines$LNCaP_DMSO)))
names(x = ident.colors) <- levels(x = cell_lines$LNCaP_DMSO)
cell.colors <- ident.colors[Idents(object = cell_lines$LNCaP_DMSO)]
names(x = cell.colors) <- colnames(x = cell_lines$LNCaP_DMSO)

cell_lines$LNCaP_DMSO <- SetAssayData(
  object = cell_lines$LNCaP_DMSO,
  assay = "RNA",
  layer = "counts",
  new.data = GetAssayData(cell_lines$LNCaP_DMSO, assay = "RNA", layer = "counts.LNCaP_DMSO.1")
)

# Check UMAP embeddings exist
head(Embeddings(cell_lines$LNCaP_DMSO, "umap"))

show.velocity.on.embedding.cor(emb = Embeddings(object = trial, reduction = "umap"), vel = Tool(object = trial,
                                                                                                                slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
                               cex = 0.8, arrow.scale = 4, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                               do.par = FALSE, cell.border.alpha = 0.1)
?Tool
cell_lines

#For pseudotime analysis (Switch to python)
cell_lines
# Convert list of Seurat objects (preserving original names)
lapply(names(cell_lines), function(obj_name) {
  obj <- cell_lines[[obj_name]]  # Access list element directly
  
  # Preserve original cell line identifier
  filename <- paste0(obj_name, ".h5Seurat")
  
  SaveH5Seurat(obj, filename = filename)
  Convert(filename, dest = "h5ad")
})

#Set memory on R

#Check the S/G2M genes in seurat_object_new.rds
seurat_object_new=readRDS("seurat_object_new.rds")
colnames(seurat_object_new@meta.data)


# Convert the cluster numbers from factor to numeric
clusters <- as.numeric(as.character(seurat_object_new$seurat_clusters))

# Create the 'cell_group' column based on cluster assignments
seurat_object_new$cell_group <- dplyr::case_when(
  clusters %in% c(7, 9, 14) ~ "Initial",
  clusters == 0 ~ "ENZ168",
  clusters %in% c(5, 6) ~ "Resistance",
  TRUE ~ "Persistence"
)
seurat_object_new$cell_group=factor(seurat_object_new$cell_group)
# Verify the new column and its distribution
table(seurat_object_new$cell_group, seurat_object_new$seurat_clusters)


# Extract metadata and filter if needed
metadata <- seurat_object_new@meta.data

# Define comparisons (ENZ168 vs others)
comparisons <- list(
  c("ENZ168", "Initial"),
  c("ENZ168", "Resistance"),
  c("ENZ168", "Persistence")
)

# Calculate p-values
pval_df <- data.frame(
  group1 = c("ENZ168", "ENZ168", "ENZ168"),
  group2 = c("Initial", "Resistance", "Persistence"),
  p.value = sapply(comparisons, function(pair) {
    wilcox.test(G2M.Score ~ cell_group, 
                data = metadata %>% filter(cell_group %in% pair))$p.value
  })
)

# Format p-values to 2 decimal places (scientific notation for small values)
pval_df <- pval_df %>%
  mutate(
    p.label = ifelse(p.value < 0.001,
                     formatC(p.value, format = "e", digits = 2),
                     sprintf("%.2f", round(p.value, 2)))
  )
pval_df
# Create boxplot with formatted p-values
ggboxplot(metadata, x = "cell_group", y = "G2M.Score", 
          fill = "cell_group", width = 0.6,
          palette = "npg", legend = "none") +
  stat_pvalue_manual(
    pval_df,
    y.position = max(metadata$S.Score) * seq(1.15, 1.45, length.out = 3),
    label = "p = {p.label}",
    tip.length = 0.005,
    size = 4
  ) +
  labs(title = "G2/M Score Distribution Across Cell Groups",
       x = "Cell Group", y = "G/2M Score") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(size = 12)
  )
library(dplyr)
median_scores <- seurat_object_new@meta.data %>%
  dplyr::group_by(cell_group) %>%
  dplyr::summarise(median_G2M.Score = median(G2M.Score, na.rm = TRUE)) %>%
  arrange(desc(median_G2M.Score))
print(median_scores)

seurat_object_new$cell_group <- case_when(
  as.numeric(as.character(seurat_object_new$seurat_clusters)) %in% c(7,9,14) ~ "Initial",
  as.numeric(as.character(seurat_object_new$seurat_clusters)) == 0 ~ "ENZ168",
  as.numeric(as.character(seurat_object_new$seurat_clusters)) %in% c(5,6) ~ "Resistance",
  TRUE ~ "Persistence"
)
r
unique(seurat_object_new$seurat_clusters)

seurat_obj_WGCNA_meta=readRDS("seurat_obj_WGCNA_meta.rds")

rna.combined.research
Layers(rna.combined.research[["RNA"]])
# Function to get number of cells in a layer
get_ncells_layer <- function(seurat_obj, assay = "RNA", layer) {
  mat <- LayerData(seurat_obj, assay = assay, layer = layer)
  ncol(mat)
}

# Get number of cells for each layer
ncells_per_layer <- sapply(layer_names, function(l) get_ncells_layer(rna.combined.research, layer = l))
ncells_per_layer

seurat_obj_WGCNA_meta
table(seurat_obj_WGCNA_meta$orig.ident)


library(dplyr)

# Extract metadata
meta <- seurat_obj_WGCNA_meta@meta.data

# Summarize per original object
summary_stats <- meta %>%
  group_by(orig.ident) %>%
  dplyr::summarize(
    n_cells = n(),
    median_UMI = median(nCount_RNA),
    median_genes = median(nFeature_RNA)
  )
print(summary_stats)
