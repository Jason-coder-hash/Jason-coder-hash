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

seurat_object_new=readRDS('seurat_object_new.rds')

#Gene expression
DefaultAssay(seurat_object_new)<-"RNA"
gene_expressed_new<-FindAllMarkers(
  object = seurat_object_new,
  only.pos = FALSE,          # Include both positive/negative markers
  logfc.threshold = 0.25,    # Minimum log fold change (adjust as needed)
  min.pct = 0.1,             # Minimum expression percentage difference
  p_val_adj = 0.05          # Adjusted p-value threshold
)
saveRDS(gene_expressed_new,"gene_expressed_new.rds")
gene_list_new=gene_expressed_new$gene
saveRDS(gene_list_new,"gene_list_new.rds")

#Load other gene sets and H gene set
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
length(mdb_c5)
View(mdb_c5)
gene.set.deviated
C5_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP') %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)
class(C5_gene_sets)
class(gene.set.deviated)



####
# Get indices of matching gene sets
matched_sets <- names(C5_gene_sets)[names(C5_gene_sets) %in% gene.set.deviated]

# Subset C5_gene_sets
C5_new_gene_sets <- C5_gene_sets[matched_sets]

length(C5_new_gene_sets)

####

prostate_adeno=readRDS('prostate_adeno.rds')


#Prepare matrix
clusters <- factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

rna_features <- rownames(GetAssayData(seurat_object_new, assay = "RNA", slot = "counts"))
rna_features
# Create the pseudo_bulk matrix with proper dimensions and column names
pseudo_bulk <- matrix(0, nrow = length(rna_features), ncol = length(clusters))
rownames(pseudo_bulk) <- rna_features
colnames(pseudo_bulk) <- clusters
pseudo_bulk


# Now iterate through cell lines
for (i in seq_along(clusters)) {
  cell_line <- clusters[i]
  cells_in_line <- WhichCells(seurat_object_new, expression = seurat_clusters == cell_line)
  
  # Get integrated data for these cells
  counts_subset <- GetAssayData(seurat_object_new, assay = "RNA", slot = "counts")[, cells_in_line]
  
  # Sum by row and assign to the correct column
  # Using both numeric and name indexing for clarity
  pseudo_bulk[, i] <- rowSums(counts_subset)
}
pseudo_bulk

genes_in_expr <- rownames(pseudo_bulk)

C5_gene_sets_filter<-lapply(C5_new_gene_sets, function(gene_set) {
  intersect(gene_set, intersect(genes_in_expr,prostate_adeno))
})

#Filter gene sets by size (typical thresholds)
min_genes <- 10
max_genes <- 500
C5_gene_sets_filter_2 <- C5_gene_sets_filter[sapply(C5_gene_sets_filter, length) >= min_genes & 
                                           sapply(C5_gene_sets_filter, length) <= max_genes]
length(C5_gene_sets_filter_2)
gsvaPar <- gsvaParam(expr = pseudo_bulk,
                     geneSets = C5_gene_sets_filter_2,
                     kcdf = "Gaussian",
                     maxDiff = TRUE)
gsvaPar
# Then run gsva() with this parameter object
gsva_C5 <- gsva(gsvaPar, verbose = TRUE)
gsva_C5_copy=gsva_C5
gsva_C5
saveRDS(gsva_H_copy,"gsva_H_copy.rds")
saveRDS(H_gene_sets_filter_2,"H_gene_sets_filter_2.rds")
saveRDS(gsva_C5_copy,"gsva_C5_copy.rds")
saveRDS(C5_gene_sets_filter_2,"C5_gene_sets_filter_2.rds")
# Using ComplexHeatmap for better control
library(ComplexHeatmap)
library(circlize)

# Annotation colors for clusters
cell_line_colors <- viridis(16)
names(cell_line_colors) <- 0:15

ha <- HeatmapAnnotation(
  Cluster = factor(clusters, levels = cell_line_order),
  col = list(Cluster = cell_line_colors),
  show_legend = TRUE
)

# Heatmap color scale
col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue", "#8888FF", "white", "#FF8888", "red"))


#For supplementary
gsva_H=gsva_H_copy
# Create heatmap
ht <- Heatmap(
  gsva_H,  # Enforce column order
  name = "GSVA score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # Keep clusters in 0-15 order
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),  # Smaller column font
  width = unit(12, "cm"),  # Overall heatmap width (adjust as needed)
  column_title = "Hallmark Enrichment Scores",
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


##Hallmark (Real)
seurat_rna_ce=readRDS('seurat_rna_ce.rds')
mdb_H <- msigdbr(species = "Homo sapiens", category = "H")


#H section
# convert msigdbr gene sets to a list good for the input
H.gene.set.list <- list()
for (H.gene.set.name in unique(mdb_H$gs_name)) {
  H.gene.set.list[[H.gene.set.name]] <- mdb_H[mdb_H$gs_name %in%
                                                H.gene.set.name, ]$gene_symbol
}


H.genes <- sapply(H.gene.set.list, function(x) paste(x, collapse = ", "))
H.gene.set.list.df <- cbind(H.gene.set = names(H.gene.set.list), genes = H.genes)
rownames(H.gene.set.list.df) <- 1:nrow(H.gene.set.list.df)
H.gs.names<-names(H.gene.set.list)
# compute the deviation
H.res <- compute.kld(
  coembed = seurat_rna_ce,
  genes.use = intersect(rownames(seurat_rna_ce), rownames(seurat_object_new)),
  n.grids = 419,
  gene.set.list = H.gene.set.list[H.gs.names],
  gene.set.cutoff = 3,
  n.times = 100
)
H.res
saveRDS(H.res,'H.res.rds')
H.res=readRDS('H.res.rds')
View(H.res)
H.gene.set.deviated <- H.res[H.res$p.adj < 0.05, ]$gene.set
H.gene.set.deviated #25 enriched pathways
saveRDS(H.gene.set.deviated,"H.gene.set.deviated.rds")
H.gene.set.deviated=readRDS("H.gene.set.deviated.rds")


H_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)


####
# Get indices of matching gene sets
H.matched_sets <- names(H_gene_sets)[names(H_gene_sets) %in% H.gene.set.deviated]

# Subset H_gene_sets
H_new_gene_sets <- H_gene_sets[H.matched_sets]
H_new_gene_sets

#Filter and plot

H_gene_sets_filter<-lapply(H_new_gene_sets, function(gene_set) {
  intersect(gene_set, intersect(genes_in_expr,prostate_adeno))
})

#Filter gene sets by size (typical thresholds)
min_genes <- 10
max_genes <- 500
H_gene_sets_filter_2 <- H_gene_sets_filter[sapply(H_gene_sets_filter, length) >= min_genes & 
                                               sapply(H_gene_sets_filter, length) <= max_genes]
length(H_gene_sets_filter_2)
gsvaPar <- gsvaParam(expr = pseudo_bulk,
                     geneSets = H_gene_sets_filter_2,
                     kcdf = "Gaussian",
                     maxDiff = TRUE)
gsvaPar
# Then run gsva() with this parameter object
gsva_H <- gsva(gsvaPar, verbose = TRUE)
gsva_H_copy=gsva_H
gsva_H

# Using ComplexHeatmap for better control
library(ComplexHeatmap)
library(circlize)

# Annotation colors for clusters
cell_line_colors <- viridis(16)
names(cell_line_colors) <- 0:15

ha <- HeatmapAnnotation(
  Cluster = factor(clusters, levels = cell_line_order),
  col = list(Cluster = cell_line_colors),
  show_legend = TRUE
)

# Heatmap color scale
col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue", "#8888FF", "white", "#FF8888", "red"))


#For supplementary
gsva_H=gsva_H_copy
# Create heatmap
ht <- Heatmap(
  gsva_H,  # Enforce column order
  name = "GSVA score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # Keep clusters in 0-15 order
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),  # Smaller column font
  width = unit(12, "cm"),  # Overall heatmap width (adjust as needed)
  column_title = "Hallmark Enrichment Scores",
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


#ActivePathways
library(ActivePathways)
library(clusterProfiler)
library(GSA)
?ActivePathways

fname_scores <- system.file("extdata", "Adenocarcinoma_scores_subset.tsv",
                            package = "ActivePathways")
fname_scores
fname_GMT = system.file("extdata", "hsapiens_REAC_subset.gmt",
                        package = "ActivePathways")
dat <- as.matrix(read.table(fname_scores, header = TRUE, row.names = 'Gene'))
dat[is.na(dat)] <- 1
dat
class(fname_GMT)
gmt_data <- read.gmt(fname_GMT)
class(gmt_data)
ActivePathways(dat,fname_GMT)


#Demo finished.
gene_expressed_new=readRDS("gene_expressed_new.rds")


library(tidyverse)

# Convert NA to 1 (non-significant p-value)
pval_matrix <- gene_expressed_new %>%
  rownames_to_column("Gene") %>%
  select(Gene, cluster, p_val_adj) %>%
  mutate(cluster = factor(cluster, levels = 0:15)) %>%
  pivot_wider(
    names_from = cluster,
    values_from = p_val_adj,
    names_prefix = "cluster_",
    values_fill = 1  # Replace NA with 1
  ) %>%
  select(Gene, paste0("cluster_", 0:15))  # Ensure column order

write_tsv(pval_matrix, "activepathway_pval_matrix.tsv")

#Run the system
sin_scores=system.file("extdata", "activepathway_pval_matrix.tsv",
                       package = "ActivePathways")

sin_gmt=system.file("extdata", "filtered_gmt.gmt",
                       package = "ActivePathways")
sin_dat <- as.matrix(read.table(sin_scores, header = TRUE, row.names = 'Gene'))




ap_result=ActivePathways(sin_dat,sin_gmt)
ap_result
saveRDS(ap_result,'ap_result.rds')
ap_result[1:3,]
ap_result[,c('term_name','adjusted_p_val')]


#Important
library(future)
library(future.apply)
plan(multisession, workers = availableCores())  # Use available cores
library(bigmemory)


#Expression data for directional input
gene_expressed_new=readRDS("gene_expressed_new.rds")
# Convert NA to 1 (non-significant p-value)
dir_matrix <- gene_expressed_new %>%
  rownames_to_column("Gene") %>%
  select(Gene, cluster, avg_log2FC) %>%
  mutate(cluster = factor(cluster, levels = 0:15)) %>%
  pivot_wider(
    names_from = cluster,
    values_from = avg_log2FC,
    names_prefix = "cluster_"
  ) %>%
  select(Gene, paste0("cluster_", 0:15))%>%
  as.data.frame()# Ensure column order

row.names(dir_matrix)=dir_matrix$Gene
dir_matrix
dir_matrix=dir_matrix[-c(1)]
dir_matrix=as.matrix(dir_matrix)
dir_matrix[is.na(dir_matrix)]=0
dir_matrix=sign(dir_matrix)
saveRDS(dir_matrix,'dir_matrix.rds')
View(dir_matrix)
constraints_vector=c(1,1,1,1,1,1,1,-1,1,-1,1,1,1,1,-1,1)
constraints_vector
colnames(sin_dat)
colnames(dir_matrix)

directional_merged_pvals=merge_p_values(sin_dat,method='DPM',dir_matrix,constraints_vector)
directional_merged_pvals
merged_pvals<-merge_p_values(sin_dat,method="Brown")

lineplot_df <- data.frame(original = -log10(merged_pvals),
                          modified = -log10(directional_merged_pvals))
ggplot(lineplot_df) +
  geom_point(size = 2.4, shape = 19,
             aes(original, modified,
                 color = ifelse(original <= -log10(0.05),"gray",
                                ifelse(modified > -log10(0.05),"#1F449C","#F05039")))) +
  labs(title = "",
       x ="Merged -log10(P)",
       y = "Directional Merged -log10(P)") + 
  geom_hline(yintercept = 1.301, linetype = "dashed",
             col = 'black', size = 0.5) +
  geom_vline(xintercept = 1.301, linetype = "dashed",
             col = "black", size = 0.5) + 
  geom_abline(size = 0.5, slope = 1,intercept = 0) +
  scale_color_identity()

options(future.globals.maxSize = 1000 * 1024^2)  # 700 MiB


filtered_gmt=read.GMT(sin_gmt)
filterd_gmt=Filter(function(term) length(term$genes)>=10,filtered_gmt)
filterd_gmt=Filter(function(term) length(term$genes)<=500,filtered_gmt)


#Background
prostate_adeno=readRDS('prostate_adeno.rds')
class(prostate_adeno)
background=makeBackground(filtered_gmt)
filtered_background=subset(background,(prostate_adeno %in% background))
length(background)-length(filtered_background)
enriched_pathways_directional <- ActivePathways(
  sin_dat, gmt = filtered_gmt, cytoscape_file_tag = "Directional_",
  merge_method = "DPM", scores_direction = dir_matrix, constraints_vector = constraints_vector,background=filtered_background)
saveRDS(enriched_pathways,'enriched_pathways.rds')
cyto_result

cyto_result<-ActivePathways(sin_dat,filtered_gmt, cytoscape_file_tag="enrichment__",merge_method = "DPM", scores_direction = dir_matrix, constraints_vector = constraints_vector,background=filtered_background)
saveRDS(cyto_result,"cyto_result.rds")
enriched_pathways=readRDS('enriched_pathways.rds')
enriched_pathways_directional=readRDS('enriched_pathways_directional.rds')
pathways_lost_in_directional_integration = 
  setdiff(enriched_pathways$term_id, enriched_pathways_directional$term_id)
?ActivePathways


files <- c(system.file('extra', 'enrichmentMap__pathways.txt', package='ActivePathways'),
           system.file('extra', 'enrichmentMap__subgroups.txt', package='ActivePathways'),
           system.file('extra', 'enrichmentMap__pathways.gmt', package='ActivePathways'),
           system.file('extra', 'enrichmentMap__legend.pdf', package='ActivePathways'))

readLines(files[1])

gene_enriched=read.gmt(files[3])
as.vector(gene_enriched$gene)
write.table(as.vector(gene_enriched$gene), file='gene_alpha.txt', row.names=FALSE)

demo=cyto_result[1:10,c('term_id','overlap')]

strsplit(demo$overlap,',')

# Extract unique genes from all pathways
unique_genes <- unique(unlist(demo$overlap))
unique_genes

# Calculate recurrence across pathways
gene_counts <- table(unlist(demo$overlap))
recurrent_genes <- names(gene_counts[gene_counts >= 3])  # Genes in ≥3 pathways

library(rJava)
library(kmeRtone)
prostate_cancer_mutation=read_tsv("prostate_adeno_mutation.tsv")
mutated_gene=prostate_cancer_mutation$`Gene name`
intersect(mutated_gene,unique_genes)


gene_list_for_seurat=intersect(prostate_adeno, gene_expressed_new$gene)
seurat_object_new_with_module=AddModuleScore(seurat_object_new,features = gene_list_for_seurat, seed=123,name='prostate_cancer_sig')
saveRDS(seurat_object_new_with_module,"seurat_object_new_with_module")
View(seurat_object_new_with_module@meta.data)
