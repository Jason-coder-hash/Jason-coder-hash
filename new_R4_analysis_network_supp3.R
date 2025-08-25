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


#Load WGCNA library

library(WGCNA)
library(hdWGCNA)
library(patchwork)

library(ggraph)
# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)
seurat_object_new=readRDS('seurat_object_new.rds')


seurat_obj_WGCNA <- SetupForWGCNA(
  seurat_object_new,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "lncap" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj_WGCNA_meta <- MetacellsByGroups(
  seurat_obj = seurat_obj_WGCNA,
  group.by = c("seurat_clusters"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = "seurat_clusters" # set the Idents of the metacell seurat object
)
seurat_obj_WGCNA_meta
saveRDS(seurat_obj_WGCNA_meta,'seurat_obj_WGCNA_meta.rds')

#########################################################
seurat_obj_WGCNA_meta=readRDS("seurat_obj_WGCNA_meta.rds")

# normalize metacell expression matrix:
seurat_obj_WGCNA_meta <- NormalizeMetacells(seurat_obj_WGCNA_meta)
as.vector(unique(seurat_obj_WGCNA_meta$orig.ident))
metacell_obj <- GetMetacellObject(seurat_obj_WGCNA_meta)
unique(metacell_obj$seurat_clusters)
seurat_obj_WGCNA_meta <- SetDatExpr(
  seurat_obj_WGCNA_meta,
  group_name = as.vector(unique(metacell_obj@meta.data$seurat_clusters)), 
  group.by='seurat_clusters', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

seurat_obj_WGCNA_meta <- TestSoftPowers(
  seurat_obj_WGCNA_meta,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
saveRDS(seurat_obj_WGCNA_meta,"seurat_obj_WGCNA_meta.rds")

##############################################################################
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj_WGCNA_meta)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)


###############################################################################
# construct co-expression network:
seurat_obj_WGCNA_meta=readRDS("seurat_obj_WGCNA_meta.rds")
seurat_obj_WGCNA_meta <- ConstructNetwork(
  seurat_obj_WGCNA_meta,
  tom_name = 'LNCaP_network', # name of the topoligical overlap matrix written to disk
  overwrite_tom=TRUE)
saveRDS(seurat_obj_WGCNA_meta,"seurat_obj_WGCNA_meta.rds")
PlotDendrogram(seurat_obj_WGCNA_meta, main='LNCaP hdWGCNA Dendrogram')


# compute all MEs in the full single-cell dataset
seurat_obj_WGCNA_meta <- ModuleEigengenes(
  seurat_obj_WGCNA_meta,
  group.by.vars="seurat_clusters"
)
saveRDS(seurat_obj_WGCNA_meta,'seurat_obj_WGCNA_meta.rds')
seurat_obj_WGCNA_meta

# Get module assignments for all genes
modules <- GetModules(seurat_obj_WGCNA_meta)

seurat_obj_WGCNA_meta <- ModuleConnectivity(
  seurat_obj_WGCNA_meta,
  group.by = 'seurat_clusters', group_name = as.vector(unique(seurat_obj_WGCNA_meta@meta.data$seurat_clusters))
)
# rename the modules
seurat_obj_WGCNA_meta <- ResetModuleNames(
  seurat_obj_WGCNA_meta,
  new_name = "LNCaP_grp"
)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj_WGCNA_meta, ncol=5)

# get hub genes
hub_df <- GetHubGenes(seurat_obj_WGCNA_meta, n_hubs = 50)
hub_df
head(hub_df)
saveRDS(hub_df,"hub_df.rds")
#We will remove grey modules later on!
#Gene score
library(UCell)
seurat_obj_WGCNA_meta <- ModuleExprScore(
  seurat_obj_WGCNA_meta,
  n_genes = 25,
  method='UCell'
)
saveRDS(seurat_obj_WGCNA_meta,'seurat_obj_WGCNA_meta.rds')

#TFs:
seurat_obj_WGCNA_meta=readRDS('seurat_obj_WGCNA_meta.rds')


library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
library(RSQLite)

JASPAR2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

pfm_core <- TFBSTools::getMatrixSet(
  x = sq24,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
pfm_core

# run the motif scan
seurat_obj_WGCNA_meta <- MotifScan(
  seurat_obj_WGCNA_meta,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)
motif_df <- GetMotifs(seurat_obj_WGCNA_meta)

# keep all TFs, and then remove all genes from the grey module
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(seurat_obj_WGCNA_meta)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)
genes_use
saveRDS(genes_use,'genes_use_non_grey')
# update the gene list and re-run SetDatExpr
seurat_obj_WGCNA_meta <- SetWGCNAGenes(seurat_obj_WGCNA_meta, genes_use)
seurat_obj_WGCNA_meta <- SetDatExpr(
  seurat_obj_WGCNA_meta,
  group_name = as.vector(unique(seurat_obj_WGCNA_meta@meta.data$seurat_clusters)), 
  group.by='seurat_clusters', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

# define model params:
model_params <- list(
  objective = 'reg:squarederror',
  max_depth = 1,
  eta = 0.1,
  nthread=16,
  alpha=0.5
)

# construct the TF network
seurat_obj_WGCNA_meta <- ConstructTFNetwork(seurat_obj_WGCNA_meta, model_params)
saveRDS(seurat_obj_WGCNA_meta,"seurat_obj_WGCNA_meta.rds")
seurat_obj_WGCNA_meta=readRDS("seurat_obj_WGCNA_meta.rds")
results <- GetTFNetwork(seurat_obj_WGCNA_meta)
head(results)


seurat_obj_WGCNA_meta <- AssignTFRegulons(
  seurat_obj_WGCNA_meta,
  strategy = "C",
  reg_thresh = 0.1
)
results <- GetTFNetwork(seurat_obj_WGCNA_meta)
head(results)
tf_list=length(unique(results$tf))

saveRDS(results,"tf_results.rds")
saveRDS(tf_list,"tf_list.rds")

#Note that the visualization is done only after the filtering processes.
tf_regulons <- GetTFRegulons(seurat_obj_WGCNA_meta)
# get hub genes and subset by TFs. Beware: Before doing the analysis, we will run three more programmes.
hub_df <- GetHubGenes(seurat_obj_WGCNA_meta, n_hubs=Inf) %>%
  subset(gene_name %in% tf_regulons$tf)
hub_tf_gene=unique(hub_df$gene_name)
length(hub_tf_gene)
saveRDS(hub_tf_gene,'hub_tf_gene.rds')

hub_tf_gene=readRDS('hub_tf_gene.rds')
hub_tf_gene

#GSDensity needs to redo because of non filter process
#Read GMT
#mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = 'GO:BP')
mdb_h<-msigdbr(species='Homo sapiens', category='H')
mdb_reactome=msigdbr(species="Homo sapiens", category='C2', subcategory = "CP:REACTOME")

#mdb_c5$gs_name
# Combine all three datasets
combined_mdb <- rbind(mdb_h,mdb_reactome)
combined_mdb
combined_mdb$gs_name
gene.set.list=list()
# Create gene set list
for (gene.set.name in unique(combined_mdb$gs_name)) {
  gene.set.list[[gene.set.name]] <- combined_mdb %>% 
    filter(gs_name == gene.set.name) %>% 
    pull(gene_symbol)
}
View(gene.set.list)
# Subset pathways where 10 < number of genes < 500
pathway_lengths <- lengths(gene.set.list)
gene.set.list <- gene.set.list[pathway_lengths > 10 & pathway_lengths < 500]
head(gene.set.list)
genes <- sapply(gene.set.list, function(x) paste(x, collapse = ", "))
gene.set.list.df <- cbind(gene.set = names(gene.set.list), genes = genes)
rownames(gene.set.list.df) <- 1:nrow(gene.set.list.df)
View(gene.set.list.df)
saveRDS(gene.set.list.df,'gene.set.list.df.rds')
saveRDS(gene.set.list,'gene.set.list.rds')
seurat_object_new=readRDS('seurat_object_new.rds')
gs.names=names(gene.set.list)
seurat_rna_ce=readRDS("seurat_rna_ce.rds")
ce=seurat_rna_ce
ce
rownames(seurat_object_new)
res <- compute.kld(
  coembed = ce,
  genes.use = intersect(rownames(ce), rownames(seurat_object_new)),
  n.grids = 419,
  gene.set.list = gene.set.list[gs.names],
  gene.set.cutoff = 3,
  n.times = 100
)
View(res)
res=res[res$p.adj<=0.05,]
View(res)
saveRDS(res,"res.rds")
res=readRDS("res.rds")
res$gene.set
true.gene.sets=res$gene.set
true.gene.sets
shared_genesets <- intersect(true.gene.sets, names(gene.set.list))
extracted_genes_path_gs <- gene.set.list[shared_genesets]
head(extracted_genes_path_gs)
saveRDS(extracted_genes_path_gs,"extracted_genes_path_gs.rds")
extracted_genes_path_gs=readRDS("extracted_genes_path_gs.rds")
names(extracted_genes_path_gs) #608
all_genes_gs <- unlist(extracted_genes_path_gs, use.names = FALSE)
gene_list_new=readRDS("gene_list_new.rds")

#Where
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
all_genes_gs<-intersect(all_genes_gs,gene_list_new)

saveRDS(all_genes_gs,"all_genes_gs.rds")

all_genes_gs=readRDS("all_genes_gs.rds")

# Filter genes in each pathway to those in all_genes_gs
filtered_pathways <- lapply(extracted_genes_path_gs, function(genes) {
  genes[genes %in% all_genes_gs]
})

# Remove pathways with zero genes after filtering (optional)
filtered_pathways_gs <- filtered_pathways[sapply(filtered_pathways, length) > 0]
length(filtered_pathways_gs) 
saveRDS(filtered_pathways_gs,"filtered_pathways_gs.rds")

filtered_pathways_gs=readRDS("filtered_pathways_gs.rds")
seurat_object_new=readRDS("seurat_object_new.rds")
clusters <- factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
rna_features <- rownames(GetAssayData(seurat_object_new, assay = "RNA", layer = "counts"))
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
gene.set.list=readRDS('gene.set.list.rds')
names(gene.set.list)
#Finally, GSVA
prostate_adeno=readRDS('prostate_adeno.rds')
gene_list_new=readRDS("gene_list_new.rds")
gene_list_new
gsva_gene_sets_filter<-lapply(gene.set.list, function(gene_set) {
  intersect(gene_set, intersect(gene_list_new,prostate_adeno))
})
View(gsva_gene_sets_filter) #1373 pathways



pathway_lengths=lengths(gsva_gene_sets_filter)
gsva_gene_sets_filter <- gsva_gene_sets_filter[pathway_lengths > 10 & pathway_lengths < 500]
View(gsva_gene_sets_filter)
names(gsva_gene_sets_filter)
gsvaPar <- gsvaParam(expr = pseudo_bulk,
                     geneSets = gsva_gene_sets_filter,
                     kcdf = "Gaussian",
                     maxDiff = TRUE)
gsvaPar
gsva_msig_c5_bp <- gsva(gsvaPar, verbose = TRUE)
View(gsva_msig_c5_bp)
saveRDS(gsva_msig_c5_bp,"gsva_msig_c5_bp.rds")


# Using ComplexHeatmap for better control
library(ComplexHeatmap)
library(circlize)



#Write CSV

gsva_msig_c5_bp %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  readr::write_tsv(file.path(
    "./",
    "research_gsva_results.tsv"
  ))

# Using ComplexHeatmap for better control
library(ComplexHeatmap)
library(circlize)
library(viridis)
# Annotation colors for clusters
gsva_msig_c5_bp_df <- gsva_msig_c5_bp %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway")
saveRDS(gsva_msig_c5_bp_df,"gsva_msig_c5_bp_df.rds")
cell_line_colors <- viridis(16)
names(cell_line_colors) <- 0:15

clusters <- factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
ha <- HeatmapAnnotation(
  Cluster = factor(clusters),
  col = list(Cluster = cell_line_colors),
  show_legend = TRUE
)
# Heatmap color scale
col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue", "#8888FF", "white", "#FF8888", "red"))


gsva_pathway=gsva_msig_c5_bp_df$pathway
saveRDS(gsva_pathway,"gsva_pathway.rds")
gsva_pathway=readRDS("gsva_pathway.rds")
View(filtered_pathways_gs)
true_pathway=intersect(names(filtered_pathways_gs),gsva_pathway) #84 pathways
#9 hallmark and 75 reactome pathways
saveRDS(true_pathway,"true_pathway.rds")
length(true_pathway)

prostate_adeno=readRDS('prostate_adeno.rds')
pathway_genes=unique(unlist(gene.set.list[true_pathway]))
prostate_pathway_gene=intersect(prostate_adeno,pathway_genes)
saveRDS(prostate_pathway_gene,"prostate_pathway_gene.rds")
prostate_pathway_gene=readRDS('prostate_pathway_gene.rds')
length(prostate_pathway_gene) #526 genes

# Create the pseudo_bulk matrix with proper dimensions and column names
clusters <- factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
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
gsvaPar_path <- gsva_msig_c5_bp[true_pathway,]
View(gsvaPar_path)

clusters <- factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
ha <- HeatmapAnnotation(
  Cluster = factor(clusters),
  col = list(Cluster = cell_line_colors),
  show_legend = TRUE
)
# Heatmap color scale
col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("blue", "#8888FF", "white", "#FF8888", "red"))



#All the hub and TF genes will be revealed!!!
all_hub_tf_genes=intersect(hub_tf_gene,prostate_pathway_gene)
all_hub_tf_genes
length(all_hub_tf_genes)
saveRDS(all_hub_tf_genes,"all_hub_tf_genes.rds")
all_hub_tf_genes=readRDS("all_hub_tf_genes.rds")

#Learn the turqoise genes
# Add cell group annotations based on clusters
seurat_obj_WGCNA_meta@meta.data$cell_group <- case_when(
  seurat_obj_WGCNA_meta@meta.data$seurat_clusters %in% c(7, 9, 14) ~ "Initial",
  seurat_obj_WGCNA_meta@meta.data$seurat_clusters == 0 ~ "ENZ168",
  seurat_obj_WGCNA_meta@meta.data$seurat_clusters %in% c(5, 6) ~ "Resistance",
  TRUE ~ "Persistence"  # All other clusters
)

# Convert to factor for proper ordering
seurat_obj_WGCNA_meta@meta.data$cell_group <- factor(
  seurat_obj_WGCNA_meta@meta.data$cell_group,
  levels = c("Initial", "ENZ168", "Resistance", "Persistence")
)

# Get module assignments from hdWGCNA
modules <- GetModules(seurat_obj_WGCNA_meta)

# Extract genes in the turquoise module
turquoise_genes <- subset(modules, color == "turquoise")$gene_name


# Add module score to the Seurat object
seurat_obj_WGCNA_meta <- AddModuleScore(
  object = seurat_obj_WGCNA_meta,
  features = list(turquoise_module = turquoise_genes),
  name = "turquoise_score"
)

seurat_obj_WGCNA_meta$turquoise_score=seurat_obj_WGCNA_meta$turquoise_score1
# 2. Extract scores and groups into a data frame
plot_data <- data.frame(
  turquoise_score = seurat_obj_WGCNA_meta$turquoise_score,
  cell_group = seurat_obj_WGCNA_meta@meta.data$cell_group
)


# Subset data for ENZ168 vs. others
groups <- unique(seurat_obj_WGCNA_meta$cell_group)
p_values <- list()

for (group in groups[groups != "ENZ168"]) {
  test_result <- wilcox.test(
    seurat_obj_WGCNA_meta$turquoise_score[seurat_obj_WGCNA_meta$cell_group == "ENZ168"],
    seurat_obj_WGCNA_meta$turquoise_score[seurat_obj_WGCNA_meta$cell_group == group]
  )
  p_values[[paste0("ENZ168_vs_", group)]] <- test_result$p.value
}

# Adjust p-values for multiple comparisons (Bonferroni)
adjusted_p <- p.adjust(unlist(p_values), method = "bonferroni")

library(ggpubr) # For stat_compare_means

ggplot(seurat_obj_WGCNA_meta@meta.data, aes(x = cell_group, y = turquoise_score, fill = cell_group)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = list(c("ENZ168", "Initial"), c("ENZ168", "Resistance"), c("ENZ168", "Persistence")),
    method = "wilcox.test",
    label = "p.label" # Add significance stars (*, **, ***)
  ) +
  labs(
    title = "Turquoise Module Scores Across Cell Groups",
    x = "Cell Group",
    y = "Module Score"
  ) +
  theme_minimal()



# Violin plot of module scores by cluster
VlnPlot(
  seurat_obj_WGCNA_meta,
  features = "turquoise_score",
  pt.size = 0.1,
  group.by = "seurat_clusters"  # Replace with your cluster column name
) + 
  geom_boxplot(width = 0.1, fill = "white") +
  ggtitle("Turquoise Module Scores by Cluster")

# Extract module scores and cluster IDs
scores <- seurat_obj_WGCNA_meta$turquoise_score
clusters <- seurat_obj_WGCNA_meta$seurat_clusters  # Replace with your cluster metadata column name

# Split scores into cluster 0 vs. others
cluster0_scores <- scores[clusters == 0]
other_scores <- scores[clusters != 0]

# Subset TFs present in the Turquoise module
turquoise_tfs <- intersect(turquoise_genes, all_hub_tf_genes)
print(paste("Number of TFs in Turquoise module:", length(turquoise_tfs)))


# Find markers for Cluster 0 vs. all other clusters
cluster0_markers <- FindMarkers(
  seurat_obj_WGCNA_meta,
  ident.1 = 0,                       # Cluster 0
  ident.2 = NULL,                     # Compare to all other clusters
  features = turquoise_tfs,           # Test only Turquoise TFs
  test.use = "wilcox",                # Default test
  logfc.threshold = 0,                # Include all genes (no threshold)
  min.pct = 0.1                       # Minimum expression percentage
)

cluster0_markers
# Filter significant TFs (adjust thresholds as needed)
sig_tfs <- subset(cluster0_markers, p_val_adj < 0.05 & avg_log2FC > 0)
sig_tfs
saveRDS(sig_tfs,"sig_tfs.rds")
#Related pathway
head(gene.set.list[true_pathway])
# Filter pathway genes to include only hub TFs
filtered_true_pathways <- lapply(gene.set.list[true_pathway], function(pathway_genes) {
  pathway_genes[pathway_genes %in% all_hub_tf_genes]
})
length(filtered_true_pathways)
# Remove pathways with no overlapping hub genes
filtered_true_pathways <- filtered_true_pathways[lengths(filtered_true_pathways) > 0]
names(filtered_true_pathways)
# View results
the_42_pathways=gene.set.list[names(filtered_true_pathways)]
saveRDS(the_42_pathways,"the_42_pathways.rds")
the_42_pathways=readRDS("the_42_pathways.rds")
the_42_pathways=lapply(the_42_pathways,function(pathway_genes){
  pathway_genes[pathway_genes %in% gene_list_new]
})
head(the_42_pathways)

View(the_42_pathways)
gsva_msig_c5_bp=readRDS("gsva_msig_c5_bp.rds")
gsva_choice=gsva_msig_c5_bp[names(the_42_pathways),]
View(gsva_choice)
# Create heatmap
ht <- Heatmap(
  gsva_choice,
  name = "GSVA score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  width = unit(12, "cm"),
  column_title = "Enrichment Scores for Hallmark and Reactome",
  top_annotation = ha,
  row_names_max_width = max_text_width(rownames(gsva_choice), gp = gpar(fontsize = 6)),
  heatmap_legend_param = list(
    title = "GSVA score",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    legend_direction = "horizontal",
    title_position = "topcenter"
  )
)
ht

cluster_0=FindMarkers(seurat_object_new,'0')


#Last ht analysis

# Load required packages
library(msigdbr)
library(dplyr)
library(ggplot2)

# 1. Get pathway gene sets from MSigDB (Reactome in this case)
pathway_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
pathway_list <- split(pathway_df$gene_symbol, pathway_df$gs_name)
pathway_list



your_pathways <-names(the_42_pathways)
pathway_genes <- unique(unlist(pathway_list[your_pathways]))

cluster_0_select=cluster_0[intersect(rownames(cluster_0),pathway_genes),]
cluster_0_select_down=cluster_0_select[as.numeric(cluster_0_select$avg_log2FC)<as.numeric(-1),]
cluster_0_select_down

#This will be a big one.
seurat_obj_WGCNA_meta=readRDS("seurat_obj_WGCNA_meta.rds")
results <- GetTFNetwork(seurat_obj_WGCNA_meta)
head(results)

# positive regulons
seurat_obj_WGCNA_meta <- RegulonScores(
  seurat_obj_WGCNA_meta,
  target_type = 'positive',
  ncores=8
)
saveRDS(pos_regulon_scores, "pos_regulon_scores.rds")
# negative regulons
seurat_obj_WGCNA_meta <- RegulonScores(
  seurat_obj_WGCNA_meta,
  target_type = 'negative',
  ncores=8
)
# access the results:

pos_regulon_scores <- GetRegulonScores(seurat_obj_WGCNA_meta, target_type='positive')
neg_regulon_scores <- GetRegulonScores(seurat_obj_WGCNA_meta, target_type='negative')

pos_regulon_scores=readRDS("pos_regulon_scores.rds")
neg_regulon_scores=readRDS("neg_regulon_scores.rds")
all_hub_tf_genes=unique(c(intersect(all_hub_tf_genes,colnames(pos_regulon_scores)),intersect(all_hub_tf_genes, colnames(neg_regulon_scores))))
all_hub_tf_genes
saveRDS(all_hub_tf_genes,"all_hub_tf_genes.rds")

all_hub_tf_genes=readRDS('all_hub_tf_genes.rds')
# Generate plots for all TFs
plot_list <- lapply(all_hub_tf_genes, function(tf) {
  # Create the FeaturePlot with dynamic title and feature
  FeaturePlot(seurat_obj_WGCNA_meta, 
              features = tf,  # Dynamic feature name
              label=T) + 
    umap_theme() + 
    labs(title = tf)  # Title = TF name
})
combined_plots<-wrap_plots(plot_list, ncol = 4, nrow = 5)
combined_plots
# Generate plots for all TFs
plot_list <- lapply(all_hub_tf_genes, function(tf) {
  # Add the TF-specific regulon score to metadata
  seurat_obj_WGCNA_meta[[paste0("pos_regulon_score_", tf)]] <- pos_regulon_scores[, tf]
  
  # Create the FeaturePlot with dynamic title and feature
  FeaturePlot(seurat_obj_WGCNA_meta, 
              features = paste0("pos_regulon_score_", tf),  # Dynamic feature name
              cols = c('lightgrey', 'red'),label=T) + 
    umap_theme() + 
    labs(title = tf)  # Title = TF name
})
# Combine plots into a 4x4 grid
positive_regulon_combined_plots <- wrap_plots(plot_list, ncol = 5, nrow = 4)
positive_regulon_combined_plots

# Generate plots for all TFs
plot_list <- lapply(intersect(all_hub_tf_genes, colnames(neg_regulon_scores)), function(tf) {
  # Add the TF-specific regulon score to metadata
  seurat_obj_WGCNA_meta[[paste0("neg_regulon_score_", tf)]] <- neg_regulon_scores[, tf]
  
  # Create the FeaturePlot with dynamic title and feature
  FeaturePlot(seurat_obj_WGCNA_meta, 
              features = paste0("neg_regulon_score_", tf),  # Dynamic feature name
              cols = c('lightgrey', "purple"),label=T) + 
    umap_theme() + 
    labs(title = tf)  # Title = TF name
})
# Combine plots into a 4x4 grid
negative_regulon_combined_plots <- wrap_plots(plot_list, ncol = 4, nrow = 3)
negative_regulon_combined_plots



seurat_obj_WGCNA_meta <- RunModuleUMAP(
  seurat_obj_WGCNA_meta,
  n_hubs = 10,
  n_neighbors=16, 
  min_dist=0.1
)

# select TF
cur_tf <- c("JUN","CTCF","RXRA", "SOX9","MYC")
# get the modules table
modules <- GetModules(seurat_obj_WGCNA_meta)
umap_df <- GetModuleUMAP(seurat_obj_WGCNA_meta)
View(modules)
# get module color scheme
mods <- levels(modules$module)
mod_colors <- dplyr::select(modules, c(module, color)) %>%
  distinct %>% arrange(module) %>% .$color
cp <- mod_colors; names(cp) <- mods
# get top 10 hub genes per module:
hub_df <- GetHubGenes(seurat_obj_WGCNA_meta, n_hubs=10)

# get TF regulons
tf_net <- GetTFNetwork(seurat_obj_WGCNA_meta)
tf_regulons <- GetTFRegulons(seurat_obj_WGCNA_meta) %>% 
  subset(gene %in% umap_df$gene & tf %in% umap_df$gene)
all(tf_regulons$gene %in% umap_df$gene)
tf_subnetworks$JUN[tf_subnetworks$JUN$color=='black',]
black_target=tf_subnetworks$JUN[tf_subnetworks$JUN$color=='black',]$target
write.table(black_target,"black_target_JUN.txt",row.names = F,col.names = F, quote=F)
# get the target genes
cur_network <- GetTFTargetGenes(
  seurat_obj_WGCNA_meta,
  selected_tfs=cur_tf, 
  depth=2, 
  target_type='both'
) %>% subset(gene %in% umap_df$gene & tf %in% umap_df$gene)

# get the max depth of each gene
gene_depths <- cur_network %>% 
  group_by(gene) %>% 
  slice_min(n=1, order_by=depth) %>% 
  select(c(gene, depth)) %>% distinct()


# rename columns
cur_network <- cur_network %>%
  dplyr::rename(c(source=tf, target=gene)) 

# only include connections between TFs:
cur_network <- subset(cur_network, target %in% unique(tf_net$tf) | target %in% hub_df$gene_name)
cur_network

# make a tidygraph object
graph <- tidygraph::as_tbl_graph(cur_network) %>% 
  tidygraph::activate(nodes) %>% 
  mutate(degree = centrality_degree())  

# compute the degree for each TF:
tf_degrees <- table(tf_regulons$tf)
tmp <- tf_degrees[names(V(graph))]; tmp <- tmp[!is.na(tmp)]
V(graph)[names(tmp)]$degree <- as.numeric(tmp)



# specify the selected TFs vs TFs vs genes
V(graph)$gene_type <- ifelse(names(V(graph)) %in% unique(tf_regulons$tf), 'TF', 'Gene')
V(graph)$gene_type <- ifelse(names(V(graph)) == cur_tf, 'selected', V(graph)$gene_type)

# make the layout table using the umap coords:
umap_layout <- umap_df[names(V(graph)),] %>% dplyr::rename(c(x=UMAP1, y = UMAP2, name=gene))
rownames(umap_layout) <- 1:nrow(umap_layout)
lay <- create_layout(graph, umap_layout)

# add the depth info:
gene_depths <- subset(gene_depths, gene %in% lay$name)
tmp <- dplyr::left_join(lay, gene_depths, by = c('name' = 'gene'))
lay$depth <- tmp$depth
lay$depth <- ifelse(lay$name %in% cur_tf, 0, lay$depth)
lay$depth <- factor(lay$depth, levels=0:max(as.numeric(lay$depth)))

# shape layout:
cur_shapes <- c(23, 24, 25); names(cur_shapes) <- levels(lay$depth)

# set up plotting labels
label_tfs <- subset(cur_network, target %in% tf_regulons$tf) %>% .$target %>% unique
lay$lab <- ifelse(lay$name %in% c(cur_tf, label_tfs), lay$name, NA)


p <- ggraph(lay) 

# 1: the full module umap showing all genes
p <- p + geom_point(inherit.aes=FALSE, data=umap_df, aes(x=UMAP1, y=UMAP2), color=umap_df$color, alpha=0.3, size=2)


# 2: Network edges
p <- p + geom_edge_fan(
  aes(color=Cor, alpha=abs(Cor)),
  arrow = arrow(length = unit(2, 'mm'), type='closed'), 
  end_cap = circle(3, 'mm')
) 

# 3: Network nodes (hub genes)
p <- p + geom_node_point(
  data=subset(lay, gene_type == 'Gene'), aes(fill=module), shape=21, color='black', size=2
)

# 4: Network nodes (TFs)
p <- p + geom_node_point(
  data=subset(lay, gene_type == 'TF'),
  aes(fill=module, size=degree, shape=depth), color='black'
) 

# 5: add labels
p <- p + geom_node_label(
  aes(label=lab), repel=TRUE, max.overlaps=Inf, 
  fontface='italic', color='black'
) 

# 6: set colors, shapes, clean up legends
p <- p +  scale_edge_colour_gradient2(high='orange2', mid='white', low='dodgerblue')  + 
  scale_colour_manual(values=cp) + 
  scale_fill_manual(values=cp) + 
  scale_shape_manual(values=cur_shapes) + 
  guides(
    edge_alpha="none", 
    size = "none",
    shape = "none",
    fill = "none"
  ) 

p

# make a featureplot of hMEs for each module
plot_lister <- ModuleFeaturePlot(
  seurat_obj_WGCNA_meta,
  features='hMEs', # plot the hMEs
  order=TRUE, # order so the points with highest hMEs are on top
  label_legend =TRUE
)
# plot module correlagram
ModuleCorrelogram(seurat_obj_WGCNA_meta)
# stitch together with patchwork
wrap_plots(plot_lister, ncol=2,nrow=4)
umap_df['ELF3',]



# Update selected TFs to include SOX9,JUN,CTCF and ETV1
cur_tf <- c("JUN","CTCF","RXRA", "SOX9","MYC","MXI1")
cur_network_2=cur_network
# Get TF-target network including all three TFs
cur_network_2 <- GetTFTargetGenes(
  seurat_obj_WGCNA_meta,
  selected_tfs = cur_tf, 
  depth = 2, 
  target_type = 'both'
) %>% 
  subset(gene %in% umap_df$gene & tf %in% umap_df$gene) %>%
  # Fix column renaming syntax:
  dplyr::rename(source = tf, target = gene) %>%
  # Ensure 'source' is character (not factor):
  dplyr::mutate(source = as.character(source))

# Split network into individual TF subnetworks
tf_subnetworks <- lapply(cur_tf, function(tf) {
  sub <- subset(cur_network_2, source == tf)
  # Include indirect connections through SOX9 for JUN network
  if(tf == "JUN") {
    sub <- subset(cur_network_2, source == tf | (source == "SOX9" & depth > 0))
  }
  return(sub)
}) 
names(tf_subnetworks) <- cur_tf

# After creating tf_subnetworks, add color based on target's module
tf_subnetworks <- lapply(tf_subnetworks, function(subnet) {
  # Join with umap_df to get the color for each target gene
  subnet <- dplyr::left_join(subnet, 
                             umap_df %>% dplyr::select(gene, color), 
                             by = c("target" = "gene"))
  return(subnet)
})
tf_subnetworks
colnames(tf_subnetworks$JUN)
# Filter each subnetwork to depth == 2
tf_subnetworks <- lapply(tf_subnetworks, function(subnet) {
  dplyr::filter(subnet, depth == 2)  # Keep only depth = 2 rows
})

# Verify (e.g., for JUN):
colnames(tf_subnetworks$JUN)  # Columns should remain the same
nrow(tf_subnetworks$JUN)      # Number of rows after filtering

#JUN targets 1260 genes.
View(tf_subnetworks$JUN)
saveRDS(tf_subnetworks,"tf_subnetworks.rds")

tf_subnetworks=readRDS("tf_subnetworks.rds")

#Back to basics.
#Initial clusters:7,9,14

cluster_7=FindMarkers(seurat_obj_WGCNA_meta,ident.1=7)
cluster_7_markers=cluster_7[cluster_7$p_val_adj<0.05,]
cluster_7_markers
cluster_9=FindMarkers(seurat_obj_WGCNA_meta,ident.1=9)
cluster_9_markers=cluster_9[cluster_7$p_val_adj<0.05,]
cluster_9_markers
cluster_14=FindMarkers(seurat_obj_WGCNA_meta,ident.1=14)
cluster_14_markers=cluster_14[cluster_14$p_val_adj<0.05,]
cluster_14_markers

#Trial: Find intersection with ETV1's regulated genes
initial_markers=c(rownames(cluster_7_markers),rownames(cluster_14_markers),rownames(cluster_9_markers))
ETV1_initial=tf_subnetworks$ETV1[tf_subnetworks$ETV1$target %in% initial_markers,]
head(ETV1_initial)
dim(ETV1_initial)
"ETV1" %in% tf_subnetworks$ETV1$target
library(ggrepel)
library(dplyr)
# Sort data by Correlation (Cor) descending
ETV1_initial_sorted <- ETV1_initial %>%
  arrange(desc(Cor)) %>%
  mutate(target = factor(target, levels = rev(unique(target))))  # Reverse factor levels for proper y-axis ordering

# Identify top and bottom 5 targets + SOX9
label_data <- ETV1_initial_sorted %>%
  # Get original top/bottom 5
  dplyr::slice(c(1:5, (n()-4):n())) %>%  
  # Add SOX9 specifically
  bind_rows(filter(ETV1_initial_sorted, target == "SOX9")) %>%  
  # Remove potential duplicates if SOX9 was already in top/bottom 5
  distinct(target, .keep_all = TRUE)
label_data
ggplot(ETV1_initial_sorted, aes(x = Cor, y = target)) +
  geom_point(aes(color = color), size = 3) +
  geom_vline(xintercept = mean(ETV1_initial$Cor),
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_label_repel(
    data = label_data,
    aes(label = target, color = color),
    fill = "white",          # White background
    label.size = 0,          # Remove border
    size = 3.5,
    direction = "y",
    box.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "grey50" # Add connecting line color
  ) +
  scale_color_identity() +
  labs(x = "Correlation Coefficient (Cor)",
       title = "Correlation Coefficients of Target Genes of ETV1 on Initial clusters") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.ticks.y = element_blank()
  )
#Trial
ETV1_initial_for_cytoscape=c(tf_subnetworks$ETV1$target, "ETV1",TP53_initial$target,"TP53")
tf_subnetworks$ETV1[tf_subnetworks$ETV1$target=="TP53",]
write.table(ETV1_initial_for_cytoscape,'ETV1_initial_for_cytoscape.txt',row.names =F,col.names=F,quote=F)
#Now investigate the relation with cluster 0: ENZ168 cluster
cluster_0=FindMarkers(seurat_obj_WGCNA_meta,ident.1=0)
enz168_markers=cluster_0[cluster_0$p_val_adj<0.05,]
ETV1_enz168=tf_subnetworks$ETV1[tf_subnetworks$ETV1$target %in% rownames(enz168_markers),]
ETV1_enz168
# Sort data by Correlation (Cor) descending
ETV1_enz168_sorted <- ETV1_enz168 %>%
  arrange(desc(Cor)) %>%
  mutate(target = factor(target, levels = rev(unique(target))))  # Reverse factor levels for proper y-axis ordering

# Identify top and bottom 5 targets + SOX9
label_data <- ETV1_enz168_sorted %>%
  # Get original top/bottom 5
  slice(c(1:6, (n()-5):n())) %>%  
  # Add SOX9 specifically
  bind_rows(filter(ETV1_enz168_sorted, target == "SOX9")) %>%  
  # Remove potential duplicates if SOX9 was already in top/bottom 5
  distinct(target, .keep_all = TRUE)
label_data
ggplot(ETV1_enz168_sorted, aes(x = Cor, y = target)) +
  geom_point(aes(color = color), size = 3) +
  geom_vline(xintercept = mean(ETV1_enz168$Cor),
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_label_repel(
    data = label_data,
    aes(label = target, color = color),
    fill = "white",          # White background
    label.size = 0,          # Remove border
    size = 3.5,
    direction = "y",
    box.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "grey50" # Add connecting line color
  ) +
  scale_color_identity() +
  labs(x = "Correlation Coefficient (Cor)",
       title = "Correlation Coefficients of Target Genes of ETV1 on ENZ168 clusters") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.ticks.y = element_blank()
  )
ETV1_enz168_table_data=unique(c(as.vector(ETV1_enz168_sorted$target),"ETV1"))
write.table(ETV1_enz168_table_data,'ETV1_enz168_table_data',row.names =F,col.names=F,quote=F)


#Investigate JUN and Cluster 0
JUN_enz168=tf_subnetworks$JUN[tf_subnetworks$JUN$target %in% rownames(enz168_markers) & tf_subnetworks$JUN$depth==2,]
JUN_enz168

# Sort data by Correlation (Cor) descending
JUN_enz168_sorted <- JUN_enz168 %>%
  arrange(desc(Cor)) %>%
  mutate(target = factor(target, levels = rev(unique(target))))  # Reverse factor levels for proper y-axis ordering

# Identify top and bottom 5 targets + SOX9
label_data <- JUN_enz168_sorted %>%
  # Get original top/bottom 5
  dplyr::slice(c(1:5), (n()-4):n()) %>%  
  # Add SOX9 specifically
  bind_rows(filter(JUN_enz168_sorted, target == "SOX9")) %>%  
  # Remove potential duplicates if SOX9 was already in top/bottom 5
  distinct(target, .keep_all = TRUE)
label_data
ggplot(JUN_enz168_sorted, aes(x = Cor, y = target)) +
  geom_point(aes(color = color), size = 3) +
  geom_vline(xintercept = mean(JUN_enz168$Cor),
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_label_repel(
    data = label_data,
    aes(label = target, color = color),
    fill = "white",          # White background
    label.size = 0,          # Remove border
    size = 3.5,
    direction = "y",
    box.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "grey50" # Add connecting line color
  ) +
  scale_color_identity() +
  labs(x = "Correlation Coefficient (Cor)",
       title = "Correlation Coefficients of Target Genes of JUN on ENZ168 clusters") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.ticks.y = element_blank()
  )
JUN_enz168_table_data=unique(c(as.vector(JUN_enz168_sorted$target),"JUN"))
write.table(JUN_enz168_table_data,'JUN_enz168_table_data',row.names =F,col.names=F,quote=F)
saveRDS(JUN_enz168_table_data,"JUN_enz168_table_data.rds")


#Reanalysis
saveRDS(cur_network,"cur_network.rds")
cur_network=readRDS("cur_network.rds")
cur_network

#Back to basics.
#Initial clusters:7,9,14

cluster_7=FindMarkers(seurat_obj_WGCNA_meta,ident.1=7)
cluster_7_markers=cluster_7[cluster_7$p_val_adj<0.05,]
cluster_7_markers
cluster_9=FindMarkers(seurat_obj_WGCNA_meta,ident.1=9)
cluster_9_markers=cluster_9[cluster_7$p_val_adj<0.05,]
cluster_9_markers
cluster_14=FindMarkers(seurat_obj_WGCNA_meta,ident.1=14)
cluster_14_markers=cluster_14[cluster_14$p_val_adj<0.05,]
cluster_14_markers

#Trial: Find intersection with ETV1's regulated genes
initial_markers=c(rownames(cluster_7_markers),rownames(cluster_14_markers),rownames(cluster_9_markers))
ETV1_initial=cur_network[cur_network$target %in% initial_markers & cur_network$source=='ETV1' & cur_network$depth==2,]
ETV1_initial
ETV1_initial_sorted <- ETV1_initial %>%
  arrange(desc(Cor))
ETV1_initial_sorted
#Trial
ETV1_initial_data=c(ETV1_initial_sorted$target,"ETV1")
ETV1_initial_data
write.table(ETV1_initial_data,'ETV1_initial_data',row.names =F,col.names=F,quote=F)



#Trial: Find intersection with TP53 regulated genes

TP53_initial=tf_subnetworks$TP53[tf_subnetworks$TP53$target %in% initial_markers,]
TP53_initial
#Supp
cluster_0_markers=cluster_0[cluster_0$p_val_adj<0.05,]
cluster_0_markers
enz168_markers=rownames(cluster_0_markers)
tf_target=cur_network[cur_network$target %in% enz168_markers,"target"]
JUN_enz168_tfs=tf_subnetworks$JUN[tf_subnetworks$JUN$target %in% tf_target & tf_subnetworks$JUN$depth==2,]
JUN_enz168_tfs_sorted <- JUN_enz168_tfs %>%
  arrange(desc(Cor))
JUN_enz168_tfs_sorted
# Identify top and bottom 5 targets + SOX9
label_data <- JUN_enz168_tfs_sorted %>%
  # Get original top/bottom 5
  slice(c(1:5), (n()-4):n()) %>%  
  # Add SOX9 specifically
  bind_rows(filter(JUN_enz168_sorted, target == "SOX9")) %>%  
  # Remove potential duplicates if SOX9 was already in top/bottom 5
  distinct(target, .keep_all = TRUE)
ggplot(JUN_enz168_tfs_sorted, aes(x = Cor, y = target)) +
  geom_point(aes(color = color), size = 3) +
  geom_vline(xintercept = mean(JUN_enz168_tfs_sorted$Cor),
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_label_repel(
    data = label_data,
    aes(label = target, color = color),
    fill = "white",          # White background
    label.size = 0,          # Remove border
    size = 3.5,
    direction = "y",
    box.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "grey50" # Add connecting line color
  ) +
  scale_color_identity() +
  labs(x = "Correlation Coefficient (Cor)",
       title = "Correlation Coefficients of TF Genes associated with JUN on ENZ168 clusters") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.ticks.y = element_blank()
  )

# Update selected TFs to include SOX9,JUN,CTCF and ETV1
cur_tf <- c("JUN","CTCF","RXRA", "SOX9","MXI1")
#Investigate ELF3 and Cluster 0
cluster_0=FindMarkers(seurat_obj_WGCNA_meta,ident.1=0)
cluster_0_markers=cluster_0[cluster_0$p_val_adj<0.05,]
enz168_markers=rownames(cluster_0_markers)
CTCF_enz168_tfs=tf_subnetworks$CTCF[tf_subnetworks$CTCF$target %in% rownames(enz168_markers) & tf_subnetworks$ELF3$depth==2,]

CTCF_enz168_tfs_sorted <- CTCF_enz168_tfs %>%
  arrange(desc(Cor))
CTCF_enz168_tfs_sorted
# Identify top and bottom 5 targets + SOX9
label_data <- CTCF_enz168_tfs_sorted %>%
  # Get original top/bottom 5
  dplyr::slice(c(1:5), (n()-4):n()) %>%  
  # Add SOX9 specifically
  bind_rows(filter(CTCF_enz168_sorted, target == "SOX9")) %>%  
  # Remove potential duplicates if SOX9 was already in top/bottom 5
  distinct(target, .keep_all = TRUE)
ggplot(CTCF_enz168_tfs_sorted, aes(x = Cor, y = target)) +
  geom_point(aes(color = color), size = 3) +
  geom_vline(xintercept = mean(CTCF_enz168_tfs_sorted$Cor),
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_label_repel(
    data = label_data,
    aes(label = target, color = color),
    fill = "white",          # White background
    label.size = 0,          # Remove border
    size = 3.5,
    direction = "y",
    box.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "grey50" # Add connecting line color
  ) +
  scale_color_identity() +
  labs(x = "Correlation Coefficient (Cor)",
       title = "Correlation Coefficients of TF Genes associated with CTCF on ENZ168 clusters") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.ticks.y = element_blank()
  )


#Investigate JUN and Cluster 2 and 4
cluster_2=FindMarkers(seurat_obj_WGCNA_meta,ident.1=2)
cluster_4=FindMarkers(seurat_obj_WGCNA_meta,ident.1=4)
cluster_2_markers=cluster_2[cluster_2$p_val_adj<0.05,]
cluster_4_markers=cluster_4[cluster_4$p_val_adj<0.05,]
persist_marker=c(rownames(cluster_2_markers),rownames(cluster_4_markers))
persist_marker
"MYC" %in%rownames(cluster_2_markers)
JUN_persist_tfs=tf_subnetworks$MXI1[tf_subnetworks$MXI1$target %in% rownames(persist_marker) & tf_subnetworks$JUN$depth==2,]

new_network=tf_subnetworks$JUN[tf_subnetworks$JUN$target=="MYC",]

library(enrichR)
saveRDS(seurat_obj_WGCNA_meta,"seurat_obj_WGCNA_meta.rds")
seurat_obj_WGCNA_meta <- RunEnrichrRegulons(seurat_obj_WGCNA_meta, wait_time=1)
# get the enrichr results table 
enrich_df <- GetEnrichrRegulonTable(seurat_obj_WGCNA_meta)
# select a TF to plot 
cur_tf <- 'JUN'
plot_df <- subset(enrich_df, tf == cur_tf & P.value < 0.05)
table(plot_df$target_type)
new_plot=plot_df[plot_df$Adjusted.P.value<0.05 & plot_df$db=="GO_Biological_Process_2021",]
table(new_plot$target_type)
new_plot
# barplot for positively-correlated target gene enrichment
p2 <- new_plot %>% 
  subset(target_type == 'positive') %>%
  slice_max(n=10, order_by=Combined.Score) %>%
  mutate(Term = stringr::str_replace(Term, " \\s*\\([^\\)]+\\)", "")) %>% head(20) %>%
  ggplot(aes(x=log(Combined.Score), y=reorder(Term, Combined.Score)))+
  geom_bar(stat='identity', position='identity', fill='lightgrey') +
  geom_text(aes(label=Term), x=.1, color='black', size=3.5, hjust='left') +
  xlab('log(Enrichment)') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ggtitle('Positively correlated target genes') +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    axis.line.y=element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank()
  )
p2
ploter <-new_plot %>% 
  subset(target_type == 'positive') %>%
  mutate(Term = stringr::str_replace(Term, " \\s*\\([^\\)]+\\)", ""))
View(ploter)

neg_ploter <-new_plot %>%
  subset(target_type == 'negative') %>%
  mutate(Term = stringr::str_replace(Term, " \\s*\\([^\\)]+\\)", ""))
neg_ploter

new_plotter <-new_plot %>%  mutate(Term = stringr::str_replace(Term, " \\s*\\([^\\)]+\\)", ""))
View(new_plotter)
colnames(new_plotter)
head(new_plotter$Genes)
saveRDS(new_plotter,"new_plotter.rds")


# Identify rows where PAM16 appears in the Genes column
new_plotter$PAM16_present <- sapply(
  strsplit(new_plotter$Genes, ";"), 
  function(gene_list) "PAM16" %in% gene_list
)

# Filter pathways containing PAM16
PAM16_data <- new_plotter[new_plotter$PAM16_present, ]

# Stop if no PAM16 found
if (nrow(PAM16_data) == 0) stop("PAM16 is not present in any pathway.")


# Separate pathways into positive/negative based on target_type
positive_pathways <- PAM16_data[PAM16_data$target_type == "positive", ]
negative_pathways <- PAM16_data[PAM16_data$target_type == "negative", ]

# Combine all pathways (no top/bottom filtering)
plot_data <- rbind(positive_pathways, negative_pathways)

# Order pathways by Odds.Ratio within their target_type groups
plot_data <- plot_data %>%
  arrange(target_type, desc(Odds.Ratio)) %>%
  mutate(Term = factor(Term, levels = unique(Term)))


library(ggplot2)

ggplot(plot_data, aes(x = Term, y = Odds.Ratio, fill = target_type)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +  # Reference line
  geom_text(
    aes(label = round(Odds.Ratio, 2)), 
    hjust = ifelse(plot_data$Odds.Ratio > 1, -0.1, 1.1),  # Adjust label position
    size = 3.5
  ) +
  coord_flip() +  # Horizontal bars
  scale_fill_manual(
    values = c("positive" = "#2c7bb6", "negative" = "#d7191c"),
    name = "Pathway Type"
  ) +
  labs(
    x = "Pathway Term",
    y = "Odds Ratio",
    title = "Pathways Containing PAM16"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))  # Add space for labels





#For overall
# Get top 5 positive pathways (highest Odds.Ratio)
positive_top5 <- new_plotter %>%
  filter(target_type == "positive") %>%
  arrange(desc(Odds.Ratio)) %>%
  slice_head(n = 5)
View(positive_top5)
# Get top 5 negative pathways (lowest Odds.Ratio)
negative_top5 <- new_plotter %>%
  filter(target_type == "negative") %>%
  arrange(Odds.Ratio) %>%
  slice_tail(n = 5)

# Combine data and create directional Odds Ratio
plot_data <- bind_rows(positive_top5, negative_top5) %>%
  mutate(
    Directional_OR = ifelse(target_type == "positive", Odds.Ratio, -Odds.Ratio)
  )

# Order terms: positives first (sorted high to low), negatives next (sorted low to high)
term_order <- c(
  positive_top5$Term[order(positive_top5$Odds.Ratio, decreasing = TRUE)],
  negative_top5$Term[order(negative_top5$Odds.Ratio)]
)
plot_data$Term <- factor(plot_data$Term, levels = term_order)

library(ggplot2)

ggplot(plot_data, aes(x = Term, y = Directional_OR, fill = target_type)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +  # Divider line
  coord_flip() +  # Horizontal orientation
  scale_fill_manual(
    values = c("positive" = "#377eb8", "negative" = "#e41a1c"),
    name = "Pathway Type"
  ) +
  scale_y_continuous(
    labels = abs,  # Show absolute values for Odds Ratio
    breaks = scales::pretty_breaks()
  ) +
  labs(
    x = "Pathway Term",
    y = "Odds Ratio",
    title = "Top 5 Positive and Negative Pathways",
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, hjust = 0),
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )




#Last check pathways
# Get Reactome pathways (C2:CP:REACTOME)
reactome_df <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
)

# Get Hallmark gene sets (H)
hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

# Combine both data frames
pathway_df <- dplyr::bind_rows(reactome_df, hallmark_df)



#Last effort!
head(new_plotter$Genes)
# Split each semicolon-separated string into individual genes
genes_list <- strsplit(new_plotter$Genes, ";")
head(new_plotter)
# Combine all genes into a single vector and remove duplicates
unique_genes <- unique(unlist(genes_list))
cluster_0=FindMarkers(seurat_obj_WGCNA_meta,ident.1=0,p_val_adj=0.05)
cluster_0_gene_express=cluster_0[rownames(cluster_0) %in% unique_genes,]
dim(cluster_0_gene_express)
new_cluster_0=data.frame('Gene'=rownames(cluster_0_gene_express),'avg_log2FC'=cluster_0_gene_express$avg_log2FC)

head(new_cluster_0)
new_unique_genes=new_cluster_0$Gene
# Write the unique genes to a text file (one gene per line)
write.table(
  data.frame(Genes = new_unique_genes),  # Convert to data frame
  "cluster_0_genes.txt",                       # Output file name
  row.names = FALSE,                 # Exclude row numbers
  col.names = FALSE,                 # Exclude column header
  quote = FALSE                      # Avoid quoting gene names
)
write.csv(
  new_cluster_0,              # Your data frame
  "new_cluster_0_expression_table.csv",     # Output file name
  row.names = FALSE,          # Exclude row numbers
  quote = FALSE               # Avoid quoting gene names
)


# Split the 'Genes' column into lists of gene vectors and name them with the 'Term' column
gene_list_cluster_0 <- strsplit(as.character(new_plotter$Genes), ";")
names(gene_list_cluster_0) <- paste(new_plotter$Term, new_plotter$target_type, sep = "_")
gene_list_cluster_0


# Create formatted term labels with "(positive/negative)"
term_labels <- paste0(new_plotter$Term, ": (", new_plotter$target_type, ")")

# Split genes and collapse with spaces
gene_vectors <- strsplit(as.character(new_plotter$Genes), ";")
gene_strings <- sapply(gene_vectors, paste, collapse = " ")

# Combine into final output
formatted_output <- data.frame(
  Term_Target = term_labels,
  Genes = gene_strings,
  stringsAsFactors = FALSE
)

# Write to file without quotes
write.table(formatted_output, "gene_terms.txt",
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE,
            col.names = FALSE)
new_plotter$Term



# Load required libraries
library(tm)
library(dplyr)
new_plotter_new <- data.frame(
  Term = new_plotter$Term,
  stringsAsFactors = FALSE
)
# Create text corpus
corpus <- Corpus(VectorSource(new_plotter_new$Term))

# Text preprocessing pipeline
corpus <- corpus %>%
  tm_map(content_transformer(tolower)) %>%
  tm_map(removePunctuation) %>%
  tm_map(removeNumbers) %>%
  tm_map(removeWords, stopwords("en")) %>%
  tm_map(stripWhitespace) %>%
  tm_map(stemDocument)

# Create document-term matrix
dtm <- DocumentTermMatrix(corpus)
dtm <- removeSparseTerms(dtm, 0.95)  # Keep most frequent terms

# Convert to matrix and calculate distances
dtm_matrix <- as.matrix(dtm)
dist_matrix <- dist(dtm_matrix, method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")

# Cut tree into 7 clusters
new_plotter$similar_pathways <- as.factor(cutree(hc, k = 7))

# View clustered results
new_plotter %>% arrange(similar_pathways)

View(new_plotter)

# Convert numeric factors to named categories using your specified grouping
new_plotter_new <- new_plotter %>%
  mutate(similar_pathways = factor(similar_pathways,
                                   levels = 1:7,
                                   labels = c("Protein Modification", 
                                              "Assembly",
                                              "Regulation",
                                              "biosynthetic process",
                                              "Others",
                                              "Mitochondrial",
                                              "Protein targeting")))

# Verify the mapping
enrich_table=table(new_plotter_new$similar_pathways)
saveRDS(enrich_table,"enrich_table.rds")
# View the full mapping
new_plotter_new_arranged<-new_plotter_new %>% 
  arrange(similar_pathways)

colnames(new_plotter_new_arranged)
cluster_0_gene_express_new=data.frame("Genes"=rownames(cluster_0_gene_express), "avg_log2FC"=cluster_0_gene_express$avg_log2FC)
colnames(cluster_0_gene_express_new)



library(pheatmap)
library(dplyr)
library(stringr)

# 1. Filter genes by expression
filtered_genes <- cluster_0_gene_express_new

# 2. Set factor levels for similar_pathways for desired order
new_plotter_new_arranged$similar_pathways <- factor(
  new_plotter_new_arranged$similar_pathways,
  levels = c(
    "Protein Modification", "Assembly", "Regulation",
    "biosynthetic process", "Others", "Mitochondrial", "Protein targeting"
  )
)

# 3. Order pathways (columns) by similar_pathways
new_plotter_new_arranged <- new_plotter_new_arranged %>%
  arrange(similar_pathways)

saveRDS(new_plotter_new_arranged, "new_plotter_new_arranged.rds")
# 4. Build the matrix
heatmap_matrix <- matrix(
  0,
  nrow = nrow(filtered_genes),
  ncol = nrow(new_plotter_new_arranged)
)
rownames(heatmap_matrix) <- filtered_genes$Genes
colnames(heatmap_matrix) <- new_plotter_new_arranged$Term

for(i in 1:nrow(new_plotter_new_arranged)) {
  pathway_genes <- unlist(strsplit(new_plotter_new_arranged$Genes[i], ";"))
  for(j in 1:nrow(filtered_genes)) {
    gene <- filtered_genes$Genes[j]
    if(gene %in% pathway_genes) {
      heatmap_matrix[j, i] <- ifelse(filtered_genes$avg_log2FC[j] > 1, 1, -1)
    }
  }
}

# 5. Annotation for pathways
anno_col <- data.frame(similar_pathways = new_plotter_new_arranged$similar_pathways)
rownames(anno_col) <- new_plotter_new_arranged$Term

# 6. Define categorical color palette
ann_colors <- list(
  similar_pathways = setNames(
    c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3", "#FFFF33", "#A65628"),
    levels(new_plotter_new_arranged$similar_pathways)
  )
)

# 7. Plot with pheatmap, with small x-axis and legend fonts
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE, # Do not cluster columns, keep your grouping!
  color = c("blue", "white", "red"),
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("Negatively expressed", "Not found", "Positively expressed"),
  annotation_col = anno_col,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 3.5,
  fontsize_col = 2.5, # x-axis label size
  fontsize_row = 1.5,
  main = "Differential Gene Expression of ENZ168 group against pathways"
)

View(new_plotter_new_arranged)
#For cytoscape
# 1. Filter genes by expression
filtered_genes <- cluster_0_gene_express_new
colnames(filtered_genes)
write.csv(filtered_genes, "filtered_gene_expression.csv", row.names = FALSE)
write.table(filtered_genes$Genes,"gene_list.txt", col.names=F, row.names=F, quote = F)
#571 genes

library(RCy3)
cytoscapePing()
# Example: Map from HGNC symbol (gene symbol) to UniProt or Ensembl
mapTableColumn(
  column = "shared name",      # or "name" or whichever column holds node IDs
  species = "Human",
  map.from = "display name",           # or "Gene Symbol"
  map.to = "display name"           # or "UniProt-TrEMBL", depending on your network
)

node_table=read.csv("STRING network default node.csv")
View(node_table)
merged_genes <- merge(
  filtered_genes,
  node_table[, c("display.name", "name")],
  by.x = "Genes",
  by.y = "display.name",
  all.x = TRUE
)
merged_genes=merged_genes[,c("avg_log2FC", "name")]
write.csv(merged_genes, "filtered_gene_expression.csv", row.names = FALSE)



target_group <- "Cell cycle"  # change as needed
subset_pathways <- new_plotter_new_arranged[new_plotter_new_arranged$similar_pathways == target_group, ]

all_genes_in_group <- unique(unlist(strsplit(paste(subset_pathways$Genes, collapse = ";"), ";")))
overlap_genes <- filtered_genes$Genes[filtered_genes$Genes %in% all_genes_in_group]
write.table(overlap_genes,"Cell_cycle_genes.txt",col.names = F,row.names=F, quote=F)
View(new_plotter_new_arranged)

insert_gene=cluster_0_gene_express_new[abs(cluster_0_gene_express_new$avg_log2FC)>1,]$Genes
sum(insert_gene %in% tf_subnetworks$JUN$target)
cluster_0_gene_express_new[cluster_0_gene_express_new$Genes=="PPP6C",]

cluster_0[rownames(cluster_0)=="NAMPT",]

sig_tfs=readRDS("sig_tfs.rds")
rownames(sig_tfs) %in% tf_subnetworks$JUN$target

