# Immune Cell Subset Reanalysis Pipeline
# Focus: Myeloid and Lymphoid cell populations from integrated dataset

# ------------------------------------------------------
# 1. Extract Immune Cell Subsets
# ------------------------------------------------------
# Subset myeloid and lymphoid cells from the main dataset
immune_sc <- subset(mice_PD1_sc, subset = names_4 %in% c("Myeloid cell", "Lymphoid cell"))

# ------------------------------------------------------
# 2. Preprocessing for Immune Subset
# ------------------------------------------------------
# Normalization and feature selection
immune_sc <- NormalizeData(immune_sc, normalization.method = "LogNormalize", scale.factor = 10000)
immune_sc <- FindVariableFeatures(immune_sc)
immune_sc <- ScaleData(immune_sc)

# Dimensionality reduction with PCA
immune_sc <- RunPCA(immune_sc)
ElbowPlot(immune_sc, ndims = 50)  # Determine optimal number of PCs

# ------------------------------------------------------
# 3. Pre-integration Analysis
# ------------------------------------------------------
# Cluster unintegrated data
immune_sc <- FindNeighbors(immune_sc, dims = 1:40, reduction = "pca")
immune_sc <- FindClusters(immune_sc, resolution = 2, cluster.name = "immune_unintegrated_clusters")

# Generate UMAP for unintegrated data
immune_sc <- RunUMAP(immune_sc, dims = 1:40, reduction = "pca", reduction.name = "immune_umap.unintegrated")

# Visualize unintegrated clusters and batch effects
DimPlot(immune_sc, reduction = "immune_umap.unintegrated", group.by = "immune_unintegrated_clusters")
DimPlot(immune_sc, reduction = "immune_umap.unintegrated", split.by = "orig.ident")  # Sample distribution
DimPlot(immune_sc, reduction = "immune_umap.unintegrated", label = TRUE)

# Check marker gene expression in unintegrated data
FeaturePlot(immune_sc, features = c("Lyz2", "Cd3e", "Cd3d", "Ctla4", "Trac", "Trdc", "Foxp3", "Cd4", "Cd8a"), 
            reduction = "immune_umap.unintegrated")

# Cluster-sample distribution table
table(immune_sc$immune_unintegrated_clusters, immune_sc$orig.ident)
DimPlot(immune_sc, reduction = "immune_umap.unintegrated", group.by = "ident_column", label = TRUE)

# ------------------------------------------------------
# 4. Integrate Immune Subset
# ------------------------------------------------------
# Integrate using RPCA (Seurat v5)
immune_sc <- IntegrateLayers(
  object = immune_sc, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "immune_integrated.rpca",
  verbose = FALSE
)

# Merge layers (Seurat v5 specific step)
immune_sc <- JoinLayers(immune_sc)

# ------------------------------------------------------
# 5. Post-integration Clustering and Visualization
# ------------------------------------------------------
# Cluster integrated data
immune_sc <- FindNeighbors(immune_sc, reduction = "immune_integrated.rpca", dims = 1:35)
immune_sc <- FindClusters(immune_sc, resolution = 0.6, cluster.name = "rpca_clusters")

# Generate UMAP for integrated data
immune_sc <- RunUMAP(immune_sc, dims = 1:35, reduction = "immune_integrated.rpca")

# Visualize integrated results
DimPlot(immune_sc, reduction = "umap", split.by = "orig.ident")  # Check batch correction
DimPlot(immune_sc, reduction = "umap", group.by = "rpca_clusters", label = TRUE)
DimPlot(immune_sc, reduction = "umap", group.by = "ident_column", label = TRUE)

# ------------------------------------------------------
# 6. Marker Gene Validation
# ------------------------------------------------------
# T cell markers
FeaturePlot(immune_sc, features = c("Cd8a", "Cd4", "Trbc2", "Mki67"), reduction = "umap")
FeaturePlot(immune_sc, features = c("Cd3e", "Ctla4", "Trac", "Foxp3", "Cd4", "Cd8a"), reduction = "umap")

# Innate lymphoid and myeloid markers
FeaturePlot(immune_sc, features = c("Siglecf", "Il5ra", "F5", "Klf2"), reduction = "umap")
FeaturePlot(immune_sc, features = c("Ccr6", "Il17a", "Ifng"), reduction = "umap")
FeaturePlot(immune_sc, features = c("Mrc1", "Lyz2", "Cybb", "Gata3"), reduction = "umap")
FeaturePlot(immune_sc, features = c("Irf4", "Lyz2", "Cybb", "Gata3"), reduction = "umap")

# Cluster-sample distribution after integration
table(immune_sc$seurat_clusters, immune_sc$orig.ident)

# ------------------------------------------------------
# 7. Immune Cell Type Annotation
# ------------------------------------------------------
# Assign cell type labels to clusters
immune_cell_names_3 <- c(
  "Th cell", "CD8+ T cell", "mono-macro", "mono-macro", "γδT cell", 
  "mono-macro", "mono-macro", "Mast cell", "ILC2", "Cycling αβT cell", 
  "NK cell", "Basophil", "cDC2", "Treg cell", "Langerhans cell",
  "γδT cell", "mregDC", "cDC1", "mono-macro", "mono-macro",
  "Neutrophil", "NK cell", "Eosinophil", "mono-macro"
)

# Apply annotations to metadata
immune_sc@meta.data$Immune_cell_names_3 <- immune_sc@meta.data$seurat_clusters
levels(immune_sc@meta.data$Immune_cell_names_3) <- immune_cell_names_3

# Organize cell type order for visualization
immune_sc$Immune_cell_names_3 <- factor(immune_sc$Immune_cell_names_3, 
                                        levels = c("CD8+ T cell", "Th cell", "Treg cell", "Cycling αβT cell", "γδT cell", 
                                                   "ILC2", "NK cell", "mono-macro", "cDC1", "cDC2", "mregDC", 
                                                   "Langerhans cell", "Mast cell", "Neutrophil", "Eosinophil", "Basophil")
)

# ------------------------------------------------------
# 8. Visualization with Custom Colors
# ------------------------------------------------------
# Define color palette
immune_colors <- c(
  "#D51F26",    # CD8+ T cell
  "#C06CAB",    # Th cell
  "#89288F",    # Treg cell
  "pink",       # Cycling αβT cell
  "#717171",    # γδT cell
  "#785A5B",    # ILC2
  "#90D5E4",    # NK cell
  "#8A9FD1",    # mono-macro
  "#357DB9",    # cDC1
  "#29306B",    # cDC2
  "#41888a",    # mregDC
  "#208A42",    # Langerhans cell
  "#97CE72",    # Mast cell
  "#FEE501",    # Neutrophil
  "#D8A767",    # Eosinophil
  "#F79B5C"     # Basophil
)

# Plot annotated clusters
DimPlot(immune_sc, reduction = "umap", group.by = "Immune_cell_names_3", 
        label = TRUE, label.size = 5, pt.size = 0.2) + 
  scale_color_manual(values = immune_colors)

# ------------------------------------------------------
# 9. Marker Gene Expression Analysis
# ------------------------------------------------------
# Load CD45 marker genes from Excel file
cd45_markers <- read_xlsx("D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/CD45marker.xlsx")
cd45_markers <- cd45_markers[[2]]  # Extract second column as gene list

# Generate dot plot of marker expression
DotPlot(immune_sc, features = cd45_markers, assay = 'RNA', group.by = "Immune_cell_names_3") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
  ) +
  labs(x = NULL, y = NULL, title = "CD45 Marker Gene Expression Across Cell Types") +
  guides(size = guide_legend("Percent Expression")) +
  scale_color_gradientn(colours = c("white", 'red')) +
  scale_y_discrete(limits = rev(levels(immune_sc$Immune_cell_names_3)))

# Generate and save violin plots for each marker
pdf("CD45markers.vln2.pdf", 4, 3)
lapply(cd45_markers, function(x) {
  VlnPlot(immune_sc, features = x, group.by = "Immune_cell_names_3", pt.size = 0) + 
    labs(x = "") +
    NoLegend() +
    scale_fill_manual(values = immune_colors)
})
dev.off()

# Generate and save feature plots for each marker
pdf("Immune_sc_featureplot.pdf", 4, 3.5)
lapply(cd45_markers, function(x) {
  FeaturePlot(immune_sc, features = x, reduction = "umap", pt.size = 0.2) +
    scale_color_gradient(low = "lightgrey", high = "darkred", na.value = "white") +
    labs(title = x)
})
dev.off()



# ------------------------------------------------------
# 10. Differential Expression Analysis (PBS vs DT)
# ------------------------------------------------------
# Set Seurat object and define cell type identities
seurat_obj <- immune_sc
clusters <- unique(immune_sc$Immune_cell_names_3)
Idents(seurat_obj) <- seurat_obj$Immune_cell_names_3

# Load required library for Excel handling
library(openxlsx)

# Create workbook to store DEG results
wb_deg <- createWorkbook()

# Calculate DEGs for each immune cell type (DT vs PBS) without p-value filtering
for (cluster in clusters) {
  print(paste("Processing cell type:", cluster))
  
  # Subset to current cell type
  cluster_cells <- subset(seurat_obj, idents = cluster)
  
  # Identify differentially expressed genes (DT vs PBS)
  # Parameters: minimum 25% expression in either group, logFC threshold of 0.25
  markers <- FindMarkers(
    object = cluster_cells,
    ident.1 = "DT",          # Reference group
    ident.2 = "PBS",         # Comparison group
    group.by = "group",      # Metadata column defining groups
    min.pct = 0.25,          # Minimum expression percentage in either group
    logfc.threshold = 0.25   # Minimum log2 fold change
  )                          # No p-value filtering applied
  
  # Add results to Excel workbook
  addWorksheet(wb_deg, sheetName = cluster)
  writeData(wb_deg, sheet = cluster, x = markers, rowNames = TRUE)
}

# Save DEG results
saveWorkbook(wb_deg, file = "DEG_Results_NoPvalFilter.xlsx", overwrite = TRUE)

print("Differential expression analysis completed. Results saved in 'DEG_Results_NoPvalFilter.xlsx'")



