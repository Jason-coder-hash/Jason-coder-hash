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


library(SeuratDisk)
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

#SCPA
library(SCPA)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(magrittr)
library(gsdensity)


# Load required libraries
library(Seurat)
library(monocle3)  # For pseudotime analysis
library(AUCell)    # For gene set scoring
library(msigdbr)   # For pathway gene sets
library(cluster)   # For pam clustering
library(tidyverse)

library(monocle3)
library(AUCell)
library(future)
library(future.apply)
plan(multisession, workers = availableCores())  # Use available cores  # Use available cores
plan(multisession, workers = 13) 
library(bigmemory)
library(furrr)
options(future.globals.maxSize = 10000 * 1024^2)  # 700 MiB


?AUCell_buildRankings
seurat_rna_ce=readRDS("seurat_rna_ce.rds")
ce=seurat_rna_ce
?compute.kld

#GSDensity needs to redo because of non filter process
#Read GMT
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = 'GO:BP')
mdb_h<-msigdbr(species='Homo sapiens', category='H')
mdb_reactome=msigdbr(species="Homo sapiens", category='C2', subcategory = "CP:REACTOME")
mdb_c5
# Combine all three datasets
combined_mdb <- rbind(mdb_c5, mdb_h, mdb_reactome)

# Create gene set list
for (gene.set.name in unique(combined_mdb$gs_name)) {
  gene.set.list[[gene.set.name]] <- combined_mdb %>% 
    filter(gs_name == gene.set.name) %>% 
    pull(gene_symbol)
}

head(gene.set.list)
genes <- sapply(gene.set.list, function(x) paste(x, collapse = ", "))
gene.set.list.df <- cbind(gene.set = names(gene.set.list), genes = genes)
rownames(gene.set.list.df) <- 1:nrow(gene.set.list.df)
View(gene.set.list.df)
saveRDS(gene.set.list.df,'gene.set.list.df.rds')
saveRDS(gene.set.list,'gene.set.list.rds')
seurat_object_new=readRDS('seurat_object_new.rds')
gs.names=names(gene.set.list)
res <- compute.kld(
    coembed = ce,
    genes.use = intersect(rownames(ce), rownames(seurat_object_new)),
    n.grids = 419,
    gene.set.list = gene.set.list[gs.names],
    gene.set.cutoff = 3,
    n.times = 100
  )
saveRDS(res,'res.rds')
res=res[res$p.adj<=0.05,]
View(res)
res=readRDS('res.rds')

#AUC
true.gene.sets=res$gene.set
true.gene.sets

# 2. Pseudotime analysis with Monocle3
cds <- as.cell_data_set(seurat_object_new)
cds <- cluster_cells(cds)
colData(cds)
# find all possible partitions
all_partitions <- unique(cds@clusters$UMAP$partitions)
all_partitions <- all_partitions[all_partitions != "1"]

# set all partitions to 1
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions %in% all_partitions] <- "1"
cds <- learn_graph(cds, use_partition = FALSE)

get_earliest_principal_node <- function(cds, ident=9){
  cell_ids <- which(colData(cds)[, "ident"] == ident)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
cds
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


seurat_object_new$pseudotime <- monocle3::pseudotime(cds)
seurat_object_new@meta.data

#Run the AUC
cell
expr_matrix=GetAssayData(object = seurat_object_new, assay = "RNA", layer = "data")
cells_rankings <- AUCell_buildRankings(expr_matrix)
# Check for overlapping gene set names
shared_genesets <- intersect(true.gene.sets, names(gene.set.list))

# Extract matching entries #trial
extracted_genes <- gene.set.list[shared_genesets]
extracted_genes
pathway_AUC <- AUCell_calcAUC(extracted_genes, cells_rankings)
saveRDS(pathway_AUC,'pathway_AUC.rds')
pathway_AUC@NAMES





#SCPA (Failed)
gene.set.list=readRDS('gene.set.list.rds')
head(gene.set.list)
filtered_gene_sets <- gene.set.list[lengths(gene.set.list) >= 10 & lengths(gene.set.list) <= 500]
gene.set.list=filtered_gene_sets
saveRDS(gene.set.list,'gene.set.list.rds')


expr_matrix=seurat_extract(seurat_object_new)
View(head(expr_matrix))
gmt_files <- list.files(path='.',pattern = 'filtered_gene_sets.gmt',full.names = TRUE)
scpa_path=compare_pathways(expr_matrix,gmt_files, min_genes = 10, max_genes=500,parallel = TRUE, cores=13)
scpa_path

SaveH5Seurat(seurat_object_new, filename="seurat_filtered.h5seurat", overwrite=TRUE)
Convert("seurat_filtered.h5seurat",dest='h5ad', overwrite=TRUE)
