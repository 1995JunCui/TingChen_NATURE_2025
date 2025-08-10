
# ------------------------------------------------------
# Keratinocyte Subset Analysis Pipeline
# ------------------------------------------------------

# ------------------------------------------------------
# 1. Memory Allocation and Data Preparation
# ------------------------------------------------------
# Increase memory limit for large dataset processing
options(future.globals.maxSize = 20 * 1024^3)  # 20GB memory allocation

# Split RNA assay by sample origin for integration
KC_sc[["RNA"]] <- split(KC_sc[["RNA"]], f = KC_sc$orig.ident)

# ------------------------------------------------------
# 2. Preprocessing
# ------------------------------------------------------
# Normalization and variable feature selection
KC_sc <- NormalizeData(KC_sc, normalization.method = "LogNormalize", scale.factor = 10000)
KC_sc <- FindVariableFeatures(KC_sc)

# Scale data with regression of potential confounders
KC_sc <- ScaleData(
  KC_sc, 
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")  # Correct for mitochondrial content and cell cycle
)

# ------------------------------------------------------
# 3. Dimensionality Reduction and Integration
# ------------------------------------------------------
# Perform PCA
KC_sc <- RunPCA(KC_sc)
ElbowPlot(KC_sc, ndims = 50)  # Determine optimal number of PCs

# Integrate data using Harmony to correct batch effects
KC_sc <- IntegrateLayers(
  object = KC_sc, 
  method = HarmonyIntegration, 
  orig.reduction = "pca", 
  new.reduction = "harmony", 
  verbose = FALSE
)

# Merge layers (Seurat v5 specific step)
KC_sc[["RNA"]] <- JoinLayers(KC_sc[["RNA"]])

# ------------------------------------------------------
# 4. Clustering and Visualization
# ------------------------------------------------------
# Cluster integrated data
KC_sc <- FindNeighbors(KC_sc, reduction = "harmony", dims = 1:13)
KC_sc <- FindClusters(KC_sc, resolution = 0.8, cluster.name = "harmony_clusters")

# Generate UMAP and t-SNE visualizations
KC_sc <- RunUMAP(KC_sc, dims = 1:13, reduction = "harmony", reduction.name = "harmony_umap")
KC_sc <- RunTSNE(KC_sc, dims = 1:13, reduction = "harmony", reduction.name = "harmony_tsne")

# Visualize clusters
DimPlot(KC_sc, reduction = "harmony_umap", group.by = "harmony_clusters", label = TRUE, label.size = 5, pt.size = 0.2)
DimPlot(KC_sc, reduction = "harmony_tsne", group.by = "harmony_clusters", label = TRUE, label.size = 5, pt.size = 0.2)
DimPlot(KC_sc, reduction = "harmony_umap", group.by = "rpca_clusters", label = TRUE, label.size = 5, pt.size = 0.2)

# Check sample distribution across clusters
table(KC_sc$orig.ident, KC_sc$harmony_clusters)
DimPlot(KC_sc, reduction = "harmony_umap", split.by = "orig.ident")

# ------------------------------------------------------
# 5. Cell Type Annotation
# ------------------------------------------------------
# Define cluster identities
cluster_ids <- KC_sc$harmony_clusters
new_cluster_names <- rep(length(cluster_ids))

new_cluster_names[cluster_ids %in% c(14)] <- "Epi_supra"
new_cluster_names[cluster_ids %in% c(7)] <- "Epi_basal"
new_cluster_names[cluster_ids %in% c(1, 15)] <- "IFD_basal"
new_cluster_names[cluster_ids %in% c(3)] <- "IFD_supra"
new_cluster_names[cluster_ids %in% c(17)] <- "SG"
new_cluster_names[cluster_ids %in% c(8)] <- "Istmus"
new_cluster_names[cluster_ids %in% c(0)] <- "HFSC_1"
new_cluster_names[cluster_ids %in% c(12)] <- "HFSC_2"
new_cluster_names[cluster_ids %in% c(10)] <- "Niche"
new_cluster_names[cluster_ids %in% c(4, 5)] <- "ORS"
new_cluster_names[cluster_ids %in% c(6)] <- "CL"
new_cluster_names[cluster_ids %in% c(13)] <- "IRS"
new_cluster_names[cluster_ids %in% c(9, 2)] <- "Matrix"
new_cluster_names[cluster_ids %in% c(11)] <- "CX"
new_cluster_names[cluster_ids %in% c(16)] <- "Medulla"

# Assign and order cell types
KC_sc$harmony_KC_types <- new_cluster_names
desired_order <- c(
  "Epi_supra", "Epi_basal", "IFD_basal", "IFD_supra", "SG", "Istmus",
  "HFSC_1", "HFSC_2", "Niche", "ORS", "CL", "IRS", "Matrix", "CX", "Medulla"
)
KC_sc$harmony_KC_types <- factor(KC_sc$harmony_KC_types, levels = desired_order)

# Visualize annotated cell types
DimPlot(KC_sc, reduction = "harmony_umap", group.by = "harmony_KC_types", label = TRUE)
DimPlot(KC_sc, reduction = "harmony_tsne", group.by = "harmony_KC_types", label = TRUE)

# ------------------------------------------------------
# 6. Custom Visualization with Color Coding
# ------------------------------------------------------
# Define color palette
immunecol <- c(
  "#F79B5C", "#D8A767", "#fbda41", "#c3d94e", "#5a9f35", "#53b76f",
  "pink", "#D51F26", "#6E318E", "#A67EB7", "#90D5E4", "#29306B",
  "#357DB9", "#785A5B", "#979797"
)

# Plot annotated clusters with custom colors
DimPlot(KC_sc, reduction = "harmony_umap", group.by = "harmony_KC_types", 
        label = TRUE, label.size = 5, pt.size = 0.2) + 
  scale_color_manual(values = immunecol)

# ------------------------------------------------------
# 7. Marker Gene Validation
# ------------------------------------------------------
# Load marker genes from Excel file
markers <- read_xlsx("D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/Markers.xlsx")
markers <- markers[[6]]  # Extract 6th column

# Remove NA values from marker list
if (any(is.na(markers))) {
  markers <- markers[!is.na(markers)]
}

# Generate dot plot of marker expression
DotPlot(KC_sc, features = markers, assay = 'RNA', group.by = "harmony_KC_types") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
  ) +
  labs(x = NULL, y = NULL, title = "Marker Gene Expression Across Cell Types") +
  guides(size = guide_legend("Percent Expression")) +
  scale_color_gradientn(colours = c("white", 'red')) +
  scale_y_discrete(limits = rev(levels(KC_sc$harmony_KC_types)))

# Save t-SNE feature plots
pdf("KC_sc_Tsne_featureplot.pdf", 4, 3.5)
lapply(markers, function(x) {
  FeaturePlot(KC_sc, features = x, reduction = "harmony_tsne", pt.size = 0.2) +
    scale_color_gradient(low = "lightgrey", high = "darkred", na.value = "white") +
    labs(title = x)
})
dev.off()

# Save violin plots
pdf("KC_sc.vln.pdf", 4, 3)
lapply(markers, function(x) {
  VlnPlot(KC_sc, features = x, group.by = "harmony_KC_types", pt.size = 0) + 
    labs(x = "") +
    NoLegend() +
    scale_fill_manual(values = immunecol)
})
dev.off()

# ------------------------------------------------------
# 8. Differential Expression Analysis (PBS vs DT)
# ------------------------------------------------------
# Set up Seurat object and define cell type identities
seurat_obj <- KC_sc
clusters <- unique(KC_sc$KC_types_merged)
Idents(seurat_obj) <- seurat_obj$KC_types_merged

# Load required library
library(openxlsx)

# Create workbook for DEG results
wb_deg <- createWorkbook()

# Calculate DEGs for each keratinocyte subtype (DT vs PBS) without logFC filtering
for (cluster in clusters) {
  print(paste("Processing cell type:", cluster))
  
  # Subset to current cell type
  cluster_cells <- subset(seurat_obj, idents = cluster)
  
  # Identify differentially expressed genes (DT vs PBS)
  markers <- FindMarkers(
    object = cluster_cells,
    ident.1 = "DT",
    ident.2 = "PBS",
    group.by = "group",
    min.pct = 0.25,
    logfc.threshold = 0  # No log fold change filtering
  )
  
  # Add results to Excel workbook
  addWorksheet(wb_deg, sheetName = cluster)
  writeData(wb_deg, sheet = cluster, x = markers, rowNames = TRUE)
}

# Save DEG results
saveWorkbook(wb_deg, file = "DEG_Results_KC_sc_KC_types_merged.xlsx", overwrite = TRUE)
print("DEG analysis completed. Results saved in 'DEG_Results_KC_sc_KC_types_merged.xlsx'")








