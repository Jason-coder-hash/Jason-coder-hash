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
library(reshape2)
library(SeuratData)
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

library(Rgraphviz)
library(CelliD)
library(supraHex)
library(dnet)
library(parallelly)


library(gsdensity)
library(future)
library(future.apply)
library(msigdbr)

# Install the package from source
library(genix)
library(igraph)


seurat_object=readRDS('seurat_object.rds')
saveRDS(seurat_object_new,'seurat_object_new.rds')
gsva_hallmark_copy=readRDS("gsva_hallmark.rds")
seurat_object_new=JoinLayers(seurat_object)
seurat_object_new
seurat_rna_ce=compute.mca(object = seurat_object_new)
saveRDS(seurat_rna_ce, 'seurat_rna_ce.rds')
seurat_rna_ce
seurat_rna_ce=readRDS('seurat_rna_ce.rds')

#Real pathway

#Load msigdb data
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
length(mdb_c5)
View(mdb_c5)
# If we just want to do biological process:
mdb_c5_bp <- mdb_c5[mdb_c5$gs_subcat == "GO:BP", ]
mdb_c5_bp$gs_name
# convert msigdbr gene sets to a list good for the input
gene.set.list <- list()
for (gene.set.name in unique(mdb_c5_bp$gs_name)) {
  gene.set.list[[gene.set.name]] <- mdb_c5_bp[mdb_c5_bp$gs_name %in%
                                                gene.set.name, ]$gene_symbol
}


genes <- sapply(gene.set.list, function(x) paste(x, collapse = ", "))
gene.set.list.df <- cbind(gene.set = names(gene.set.list), genes = genes)
rownames(gene.set.list.df) <- 1:nrow(gene.set.list.df)
head(gene.set.list.df)

gs.names<-names(gene.set.list)
# compute the deviation
res <- compute.kld(
  coembed = seurat_rna_ce,
  genes.use = intersect(rownames(seurat_rna_ce), rownames(seurat_object_new)),
  n.grids = 419,
  gene.set.list = gene.set.list[gs.names],
  gene.set.cutoff = 3,
  n.times = 100
)
saveRDS(res,'res.rds')
res=readRDS('res.rds')
View(res)
gene.set.deviated <- res[res$p.adj < 0.05, ]$gene.set
gene.set.deviated #567 enriched pathways
saveRDS(gene.set.deviated,"gene.set.deviated.rds")
gene.set.deviated=readRDS("gene.set.deviated.rds")
seurat_cells<-colnames(seurat_object_new)
seurat_el <- compute.nn.edges(coembed = seurat_rna_ce, nn.use = 300)
View(seurat_el)

library(future.apply)
plan(multisession)  # Use available cores
library(bigmemory)
options(future.globals.maxSize = 700 * 1024^2)  # 700 MiB
# Split gene sets into batches
batches <- split(gene.set.list[gene.set.deviated], ceiling(seq_along(gene.set.list[gene.set.deviated])/100))
batches

# Retry parallel processing
seurat_cv.df_list <- future_lapply(
  batches,
  function(batch) {
    run.rwr.list(
      el = seurat_el,
      gene_set_list = batch,
      cells = seurat_cells
    )
  },
  future.seed = TRUE
)
saveRDS(seurat_cv.df_list,"seurat_cv.df_list.rds")
# Split data into chunks for parallel processing
# Increase memory limit for large objects (if needed)
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
library(future.apply)


# Compute cell labels in parallel
seurat_labels_list <- future_lapply(
  seurat_cv.df_list,
  function(cell_df) {
    compute.cell.label.df(cell_df)  # S3 method dispatch for 'cell.label.df' class
  },
  future.seed = TRUE  # Maintain reproducibility if needed
)

saveRDS(seurat_labels_list,'seurat_labels_list.rds')
head(gene.set.deviated)


#We will now produce the true heatmap

DefaultAssay(seurat_object_new)<-"RNA"
seurat_object_new

# Load required packages
library(msigdbr)
library(GSVA)
library(ComplexHeatmap)
library(viridis)

cell_lines <- factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

# Switch to the RNA assay (contains all genes)
DefaultAssay(seurat_object_new) <- "RNA"
#Ok. Now Find all markers in every regions.
gene_expressed<-FindAllMarkers(
  object = seurat_object_new,
  only.pos = FALSE,          # Include both positive/negative markers
  logfc.threshold = 0.25,    # Minimum log fold change (adjust as needed)
  min.pct = 0.1,             # Minimum expression percentage difference
  p_val_adj = 0.05          # Adjusted p-value threshold
)
saveRDS(gene_expressed,"gene_expressed_new.rds")
gene_list_enriched=gene_expressed$gene
saveRDS(gene_list_enriched,"gene_list_enriched.rds")
gene_list_enriched
# Rebuild pseudo_bulk using RNA assay
pseudo_bulk <- matrix(0, nrow = nrow(seurat_object_new), ncol = 16)
rownames(pseudo_bulk) <- rownames(seurat_object_new)
colnames(pseudo_bulk) <- 0:15

for (i in 0:15) {
  cells <- WhichCells(seurat_object_new, expression = seurat_clusters == i)
  counts_subset <- GetAssayData(seurat_object_new, assay = "RNA", slot = "data")[, cells]
  pseudo_bulk[, i + 1] <- rowMeans(counts_subset)  # Use rowMeans for stability
}


#Using the enriched gene list
gene_list_enriched


prostate_adeno=readRDS('prostate_adeno.rds')
# 2. Filter gene sets by overlap with your data and size
# ------------------------------------------------------
genes_in_expr <- rownames(pseudo_bulk)  # Ensure this matches your integrated data
genes_in_expr
# Intersect gene sets with your dataset's genes
c5_gene_sets <- lapply(gene_list_enriched, function(gene_set) {
  intersect(gene_set, intersect(intersect(genes_in_expr,prostate_adeno)))
})

length(c5_gene_sets)
c5_gene_sets
# Filter by size (adjust min/max if needed)
min_genes <- 10
max_genes <- 500
c5_gene_sets <- c5_gene_sets[sapply(c5_gene_sets, length) >= min_genes & 
                               sapply(c5_gene_sets, length) <= max_genes]
length(c5_gene_sets)

# 3. Run GSVA with the C5 gene sets
# -----------------------------------
# Use the same pseudo_bulk matrix as before (ensure assay = "integrated")
gsvaPar <- gsvaParam(
  expr = pseudo_bulk,
  geneSets = c5_gene_sets,
  kcdf = "Gaussian",  # Suitable for normalized log-transformed data
  maxDiff = TRUE
)

gsva_c5 <- gsva(gsvaPar, verbose = TRUE)
View(gsva_c5)
# 4. Generate the Heatmap (same as before)
# ----------------------------------------
# Define cluster order (0-15)
cell_line_order <- 0:15  # Ensure this matches colnames(pseudo_bulk)

# Annotation colors for clusters
cell_line_colors <- viridis(16)
names(cell_line_colors) <- 0:15

ha <- HeatmapAnnotation(
  Cluster = factor(cell_line_order, levels = cell_line_order),
  col = list(Cluster = cell_line_colors),
  show_legend = TRUE
)

# Heatmap color scale
col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue", "#8888FF", "white", "#FF8888", "red"))
gsva_c5[1:30,]
# Create heatmap
ht <- Heatmap(
  gsva_c5[1:30,],  # Enforce column order
  name = "GSVA score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # Keep clusters in 0-15 order
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),  # Smaller column font
  width = unit(12, "cm"),  # Overall heatmap width (adjust as needed)
  column_title = "C5 GO Pathway Enrichment Scores",
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
# Draw the heatmap
draw(ht, heatmap_legend_side = "bottom",annotation_legend_side = "left")




# Create HeatmapAnnotation with a legend for clusters
ha <- HeatmapAnnotation(
  Cluster = factor(colnames(gsva_c5), levels = 0:15),
  col = list(Cluster = cell_line_colors),
  show_legend = TRUE,
  annotation_legend_param = list(
    Cluster = list(
      title = "Cluster",
      title_position = "topcenter",
      legend_direction = "vertical"
    )
  )
)

# Create the heatmap
ht <- Heatmap(
  gsva_c5[1:60,],  # Adjust rows as needed
  name = "GSVA score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),  # Smaller column font
  width = unit(12, "cm"),  # Overall heatmap width (adjust as needed)
  column_title = "C5 GO Pathway Enrichment Scores",
  top_annotation = ha,
  heatmap_legend_param = list(
    title = "GSVA score",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    legend_direction = "horizontal",
    title_position = "topcenter"
  )
)
# Draw the heatmap with cluster legend on the right and GSVA legend at the bottom
draw(ht, 
     heatmap_legend_side = "bottom",  # GSVA score legend
     annotation_legend_side = "left"  # Cluster color legend
)

# Create the heatmap
ht2 <- Heatmap(
  gsva_c5[61:120,],  # Adjust rows as needed
  name = "GSVA score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),  # Smaller column font
  width = unit(12, "cm"),  # Overall heatmap width (adjust as needed)
  column_title = "C5 GO Pathway Enrichment Scores",
  top_annotation = ha,
  heatmap_legend_param = list(
    title = "GSVA score",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    legend_direction = "horizontal",
    title_position = "topcenter"
  )
)
# Draw the heatmap with cluster legend on the right and GSVA legend at the bottom
draw(ht2, 
     heatmap_legend_side = "bottom",  # GSVA score legend
     annotation_legend_side = "left"  # Cluster color legend
)

# Create the heatmap
ht3 <- Heatmap(
  gsva_c5[121:180,],  # Adjust rows as needed
  name = "GSVA score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),  # Smaller column font
  width = unit(12, "cm"),  # Overall heatmap width (adjust as needed)
  column_title = "C5 GO Pathway Enrichment Scores",
  top_annotation = ha,
  heatmap_legend_param = list(
    title = "GSVA score",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    legend_direction = "horizontal",
    title_position = "topcenter"
  )
)
# Draw the heatmap with cluster legend on the right and GSVA legend at the bottom
draw(ht3, 
     heatmap_legend_side = "bottom",  # GSVA score legend
     annotation_legend_side = "left"  # Cluster color legend
)
saveRDS(gsva_c5, 'gsva_c5.rds')







