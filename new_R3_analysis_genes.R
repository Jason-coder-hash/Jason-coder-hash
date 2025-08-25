setwd("C:/GNBF6010")
getwd()
library(msigdbdf)
library(tidyverse)
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
library(readxl)

# Load required libraries
library(Signac)
library(EnsDb.Hsapiens.v86)
data("blacklist_hg38_unified")



library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(AnnotationDbi)
library(KEGGREST)

library(msigdb)
library(limma)


seurat_object=readRDS('seurat_object.rds')
seurat_obj_scRNA_4_com=readRDS("seurat_obj_scRNA_4_com.rds")

#We revise the genes in clusters 0 in new object.
# Find markers for cluster 0
de_cluster0 <- FindMarkers(
  object = seurat_object,
  ident.1 = 0,
  logfc.threshold = 0.25,    # Adjust based on biological relevance
  min.pct = 0.25,            # Require detection in ≥25% of cells
  test.use = "wilcox"        # Default Wilcoxon rank-sum test
)

de_cluster0_filtered <-de_cluster0 %>%
  filter(p_val_adj < 0.05)
de_cluster0_filtered<-de_cluster0_filtered[order(-de_cluster0_filtered$avg_log2FC),]
saveRDS(de_cluster0_filtered,"de_cluster0_filtered.rds")
#Here, we report that hypoxia (HIF), CCN3 enhances AR transactivation in the presence of 0.05 and 0.1nM DHT in LNCaP prostate cancer cells.
#Herein, We include AR, MYC, MYO6, SLC4A4, etc for analysis
seurat_object_new
#Ok. Now Find all markers in every regions.
gene_expressed<-FindAllMarkers(
  object = seurat_object_new,
  only.pos = FALSE,          # Include both positive/negative markers
  logfc.threshold = 0.25,    # Minimum log fold change (adjust as needed)
  min.pct = 0.1,             # Minimum expression percentage difference
  p_val_adj = 0.05          # Adjusted p-value threshold
)
saveRDS(gene_expressed,"gene_expressed_seurat_object.rds")
gene_expressed=readRDS("gene_expressed_seurat_object.rds")
gene_list <-gene_expressed$gene
gene_list
excel_file<- read_excel("aay0267_Table_S8.xlsx", sheet = "Epi_Luminal_2Psca")
excel_file=excel_file[excel_file$Sig_filter=="TRUE",]
ref_gene<-c(excel_file$GeneSumbol)
ref_gene <- toupper(ref_gene)
ref_gene
gene_intersect=intersect(ref_gene,gene_list)
#gene_final=c(gene_intersect,"AR","MYC","MYO6","SLC4A4","CXXC5","CCN3")
gene_final=c(gene_intersect)
gene_final

saveRDS(gene_final,"gene_final.rds")
gene_list_new <- c("TSPAN8", "KRT19", "CLU", "GSTA4", "SLC6A3", "FAM83A", 
                   "KLF5", "LRRC26", "ATP1B1", "HEY1", "MALL", "MAP4K4", 
                   "MME", "TFF2", "S100A11", "CLDN2", "ANXA1", "KRT7", 
                   "CYB5A", "LINGO1", "PHLDA2", "S100A4", "ALCAM", "PARM1", 
                   "CHCH10", "BACP", "ID1", "CD2AP", "CLDN4", "MECOM", 
                   "NPOM1", "CRIPT1")
# Create the dotplot with the reordered clusters
the_dot_plot=DotPlot(seurat_object, 
        features = gene_list_new,
        cols = c("lightgreen", "red"),  # You can customize colors
        dot.scale = 8) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
colnames(the_dot_plot$data)

# Compute rankings for avg.exp.scaled
avg_rank <- the_dot_plot$data %>%
  group_by(features.plot) %>%
  arrange(desc(avg.exp.scaled), .by_group = TRUE) %>%
  mutate(rank_avg = row_number()) %>%
  ungroup() %>%
  select(features.plot, id, rank_avg)

library(tidyr)

avg_rank_wide <- avg_rank %>%
  
  pivot_wider(
    
    names_from = id,
    
    values_from = rank_avg
    
  )
avg_rank_wide
# Compute rankings for pct.exp
pct_rank <- the_dot_plot$data %>%
  group_by(features.plot) %>%
  arrange(desc(pct.exp), .by_group = TRUE) %>%
  mutate(rank_pct = row_number()) %>%
  ungroup() %>%
  select(features.plot, id, rank_pct)

pct_rank_wide <- pct_rank %>%
  
  pivot_wider(
    
    names_from = id,
    
    values_from = rank_pct
    
  )
pct_rank_wide

# Define the scoring function
calculate_score <- function(rank) {
  score <- 1 - (1/16) * (rank - 1)
  return(score)
}

# Apply scoring to avg.exp.scaled ranks
avg_score_wide <- avg_rank %>%
  mutate(score_avg = calculate_score(rank_avg)) %>%
  select(-rank_avg) %>%
  pivot_wider(
    names_from = id,
    values_from = score_avg
  )

# Apply scoring to pct.exp ranks
pct_score_wide <- pct_rank %>%
  mutate(score_pct = calculate_score(rank_pct)) %>%
  select(-rank_pct) %>%
  pivot_wider(
    names_from = id,
    values_from = score_pct
  )

View(avg_score_wide)

# 1. Calculate cluster averages for avg.exp.scaled scores
avg_cluster_means <- avg_score_wide %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  pivot_longer(
    everything(),
    names_to = "cluster",
    values_to = "avg_exp_score"
  )
# 2. Calculate cluster averages for pct.exp scores
pct_cluster_means <- pct_score_wide %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  pivot_longer(
    everything(),
    names_to = "cluster",
    values_to = "pct_exp_score"
  )

# 3. Combine both results into one table
final_cluster_scores <- avg_cluster_means %>%
  inner_join(pct_cluster_means, by = "cluster") %>%
  arrange(cluster)


# Convert scores to log scale, sum them, and rank clusters
final_cluster_ranking <- final_cluster_scores %>%
  mutate(
    # Natural log transformation (use log10() if preferred)
    log_avg = log(avg_exp_score),
    log_pct = log(pct_exp_score),
    # Sum the log-transformed scores
    total_log_score = log_avg + log_pct
  ) %>%
  # Rank clusters by descending total score
  arrange(desc(total_log_score)) %>%
  # Select relevant columns
  select(cluster, total_log_score, log_avg, log_pct, avg_exp_score, pct_exp_score)

# View the final ranked clusters
print(final_cluster_ranking)


# 1. Prepare data in long format for plotting
plot_data <- final_cluster_ranking %>%
  select(cluster, total_log_score, log_avg, log_pct) %>%
  pivot_longer(
    cols = -cluster,
    names_to = "metric",
    values_to = "score"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("total_log_score", "log_avg", "log_pct"),
                    labels = c("Total (Log)", "Avg.Exp (Log)", "Pct.Exp (Log)"))
  )

# 2. Create the plot
ggplot(plot_data, aes(x = reorder(cluster, -score), y = score, fill = metric)) +
  geom_col(position = "dodge") +
  labs(title = "Cluster Comparison: Log-Transformed Scores in terms of Expression Avg and Pct",
       x = "Cluster ID",
       y = "Log-Transformed Score",
       fill = "Metric") +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "bottom"
  )
#Representable pathways:
library(viridis)

cell_lines <- factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

integrated_features <- rownames(GetAssayData(seurat_object, assay = "integrated", slot = "data"))
integrated_features
# Create the pseudo_bulk matrix with proper dimensions and column names
pseudo_bulk <- matrix(0, nrow = length(integrated_features), ncol = length(cell_lines))
rownames(pseudo_bulk) <- integrated_features
colnames(pseudo_bulk) <- cell_lines
pseudo_bulk
seurat_object$seurat_clusters

# Now iterate through cell lines
for (i in seq_along(cell_lines)) {
  cell_line <- cell_lines[i]
  cells_in_line <- WhichCells(seurat_object, expression = seurat_clusters == cell_line)
  
  # Get integrated data for these cells
  counts_subset <- GetAssayData(seurat_object, assay = "integrated", slot = "data")[, cells_in_line]
  
  # Sum by row and assign to the correct column
  # Using both numeric and name indexing for clarity
  pseudo_bulk[, i] <- rowSums(counts_subset)
}
pseudo_bulk


#Msigdb
library(msigdbr)
# Obtain gene sets from MSigDb 
# You can choose different collections: H (hallmark), C2 (curated), C5 (GO), etc.
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)
hallmark_gene_sets

genes_in_expr <- rownames(pseudo_bulk)
genes_in_expr
hallmark_gene_sets <- lapply(hallmark_gene_sets, function(gene_set) {
  intersect(gene_set, genes_in_expr)
})

# Filter gene sets by size (typical thresholds)
min_genes <- 10
max_genes <- 500
hallmark_gene_sets <- hallmark_gene_sets[sapply(hallmark_gene_sets, length) >= min_genes & 
                                           sapply(hallmark_gene_sets, length) <= max_genes]
hallmark_gene_sets

#Annotation of pathway
# First create a parameter object with gsvaParam()
gsvaPar <- gsvaParam(expr = pseudo_bulk,
                     geneSets = hallmark_gene_sets,
                     kcdf = "Gaussian",
                     maxDiff = TRUE)
gsvaPar
# Then run gsva() with this parameter object
gsva_hallmark <- gsva(gsvaPar, verbose = TRUE)
gsva_hallmark_copy=gsva_hallmark
saveRDS(gsva_hallmark_copy,"gsva_hallmark.rds")

#Final
# Define your desired order
cell_line_order <- cell_lines

# Using ComplexHeatmap for better control
library(ComplexHeatmap)
library(circlize)

# Make sure your data is ordered correctly
gsva_hallmark <- gsva_hallmark[, cell_line_order]
gsva_hallmark
# Create color function
col_fun = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue", "#8888FF", "white", "#FF8888", "red"))


cell_line_colors <- viridis(16)
names(cell_line_colors) <- 0:15  # Match factor levels

ha <- HeatmapAnnotation(
  CellLine = factor(cell_line_order,levels = cell_line_order),
  col = list(CellLine = cell_line_colors),
  show_legend = TRUE
)
# # Set up annotation with correct colors
# ha = HeatmapAnnotation(
#   CellLine = factor(cell_line_order, levels = cell_line_order),
#   col = list(CellLine = setNames(
#     c("#4DAF4A", "#984EA3", "#FFFF33", "#377EB8", "#E41A1C"),
#     cell_line_order
#   )),
#   show_legend = TRUE
# )

# Create heatmap with forced column order
ht = Heatmap(gsva_hallmark,
             name = "GSVA score",
             col = col_fun,
             cluster_rows = TRUE,
             cluster_columns = FALSE,  # Turn off column clustering
             column_order = cell_line_order,  # Force exact column order
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 10),
             column_title = "Hallmark Pathway Enrichment Scores",
             top_annotation = ha,
             heatmap_legend_param = list(
               title = "GSVA score",
               at = c(-1, -0.5, 0, 0.5, 1),
               labels = c("-1", "-0.5", "0", "0.5", "1"),
               legend_direction = "horizontal",
               title_position = "topcenter"
             )
)
ht
