# Single-Cell RNA Sequencing Data Analysis Pipeline for Murine Samples
# This script processes scRNA-seq data from 4 murine samples (PBS_1, PBS_2, DT_1, DT_2), including quality control, doublet removal, integration, clustering, and cell type annotation.
# All steps follow standard best practices for scRNA-seq analysis using Seurat (v5) and DoubletFinder.


# ------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------
library(dplyr)          # For data manipulation
library(Seurat)         # Core single-cell analysis package
library(patchwork)      # For combining plots
library(tidyverse)      # For data wrangling and visualization
library(DoubletFinder)  # For doublet detection
library(readxl)         # For reading Excel files (if needed)
library(ggplot2)        # For custom visualization

# ------------------------------------------------------
# 2. Single Sample Processing
# ------------------------------------------------------
# Each sample is processed independently to ensure quality control and doublet removal before integration.

## 2.1 Process Sample PBS_1
# Load 10X Genomics data (filtered feature-barcode matrix)
PBS_1.data <- Read10X(data.dir = "F:/20241119/single cell seq/80-1652712019/Report/Result/EXP/01_Cellranger/PBS-1/filtered_feature_bc_matrix/")

# Create Seurat object with basic filtering:
# - Retain cells with ≥200 detected features
# - Retain features present in ≥3 cells
PBS_1 <- CreateSeuratObject(counts = PBS_1.data, project = "PBS_1", min.cells = 3, min.features = 200)

# Quality Control (QC) metrics:
# Calculate percentage of mitochondrial genes (using "^mt-" pattern for murine mitochondrial genes)
PBS_1[["percent.mt"]] <- PercentageFeatureSet(PBS_1, pattern = "^mt-")

# Filter low-quality cells based on QC metrics:
# - nFeature_RNA: 200–5000 (removes cells with too few features or potential multiplets)
# - percent.mt < 20% (removes cells with high mitochondrial content, indicative of cell stress/death)
PBS_1 <- subset(PBS_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

# Normalization and feature selection:
# Log-normalization (scales UMI counts to 10,000 per cell)
PBS_1 <- NormalizeData(PBS_1, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify 2000 most variable features using variance stabilizing transformation (vst)
PBS_1 <- FindVariableFeatures(PBS_1, selection.method = "vst", nfeatures = 2000)

# Scale data (centers and scales gene expression across cells)
all.genes <- rownames(PBS_1)
PBS_1 <- ScaleData(PBS_1, features = all.genes)

# Dimensionality reduction and clustering:
# Perform PCA using variable features
PBS_1 <- RunPCA(PBS_1, features = VariableFeatures(object = PBS_1))

# Construct nearest-neighbor graph using top 10 PCs
PBS_1 <- FindNeighbors(PBS_1, dims = 1:10)

# Cluster cells with resolution = 0.8 (balances cluster granularity)
PBS_1 <- FindClusters(PBS_1, resolution = 0.8)

# Generate UMAP visualization for clustering validation
PBS_1 <- RunUMAP(PBS_1, dims = 1:10)
DimPlot(PBS_1, reduction = "umap", label = TRUE) + ggtitle("PBS_1 Clusters (Pre-Doublet Removal)")

# Doublet removal using DoubletFinder:
# Step 1: Parameter optimization (sweep for optimal pK)
sweep.res.list_PBS_1 <- paramSweep(PBS_1, PCs = 1:10, sct = FALSE)
sweep.stats_PBS_1 <- summarizeSweep(sweep.res.list_PBS_1, GT = FALSE)
bcmvn_PBS_1 <- find.pK(sweep.stats_PBS_1)
pK_bcmvn_PBS_1 <- bcmvn_PBS_1$pK[which.max(bcmvn_PBS_1$BCmetric)] %>% as.character() %>% as.numeric()

# Step 2: Estimate expected doublet number
DoubletRate <- ncol(PBS_1) * 8 * 1e-6  # Estimated doublet rate (8% per 1000 cells, standard for 10X)
homotypic.prop <- modelHomotypic(PBS_1$seurat_clusters)  # Adjust for homotypic doublets
nExp_poi <- round(DoubletRate * ncol(PBS_1))  # Expected doublets
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))  # Adjusted for homotypic doublets

# Step 3: Identify doublets
PBS_1 <- doubletFinder(PBS_1, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
PBS_1 <- doubletFinder(PBS_1, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_1564", sct = FALSE)

# Validate doublet calls with UMAP and QC metrics
DimPlot(PBS_1, reduction = "umap", group.by = "DF.classifications_0.25_0.005_1471") + ggtitle("PBS_1 Doublets vs. Singlets")
VlnPlot(PBS_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "DF.classifications_0.25_0.005_1471")

# Save processed object
saveRDS(PBS_1, file = "D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/PBS_1.rds")

## 2.2 Process Sample DT_1 (Following the same workflow as PBS_1)
DT_1.data <- Read10X(data.dir = "F:/20241119/single cell seq/80-1652712019/Report/Result/EXP/01_Cellranger/DT-1/filtered_feature_bc_matrix/")
DT_1 <- CreateSeuratObject(counts = DT_1.data, project = "DT_1", min.cells = 3, min.features = 200)
DT_1[["percent.mt"]] <- PercentageFeatureSet(DT_1, pattern = "^mt-")
DT_1 <- subset(DT_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
DT_1 <- NormalizeData(DT_1)
DT_1 <- FindVariableFeatures(DT_1, selection.method = "vst", nfeatures = 2000)
DT_1 <- ScaleData(DT_1, features = rownames(DT_1))
DT_1 <- RunPCA(DT_1, features = VariableFeatures(object = DT_1))
DT_1 <- FindNeighbors(DT_1, dims = 1:36)  # Using 36 PCs based on variance explained
DT_1 <- FindClusters(DT_1, resolution = 0.8)
DT_1 <- RunUMAP(DT_1, dims = 1:36)
DimPlot(DT_1, reduction = "umap", label = TRUE) + ggtitle("DT_1 Clusters (Pre-Doublet Removal)")

# Doublet removal for DT_1
sweep.res.list_DT_1 <- paramSweep(DT_1, PCs = 1:36, sct = FALSE)
sweep.stats_DT_1 <- summarizeSweep(sweep.res.list_DT_1, GT = FALSE)
bcmvn_DT_1 <- find.pK(sweep.stats_DT_1)
pK_bcmvn_DT_1 <- bcmvn_DT_1$pK[which.max(bcmvn_DT_1$BCmetric)] %>% as.character() %>% as.numeric()

DoubletRate <- ncol(DT_1) * 8 * 1e-6 
homotypic.prop <- modelHomotypic(DT_1$seurat_clusters)
nExp_poi <- round(DoubletRate * ncol(DT_1)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

DT_1 <- doubletFinder(DT_1, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DT_1 <- doubletFinder(DT_1, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_579", sct = FALSE)
saveRDS(DT_1, file = "D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/DT_1.rds")

## 2.3 Process Sample DT_2 (Following the same workflow)
DT_2.data <- Read10X(data.dir = "F:/20241119/single cell seq/80-1652712019/Report/Result/EXP/01_Cellranger/DT-2/filtered_feature_bc_matrix/")
DT_2 <- CreateSeuratObject(counts = DT_2.data, project = "DT_2", min.cells = 3, min.features = 200)
DT_2[["percent.mt"]] <- PercentageFeatureSet(DT_2, pattern = "^mt-")
DT_2 <- subset(DT_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
DT_2 <- NormalizeData(DT_2)
DT_2 <- FindVariableFeatures(DT_2)
DT_2 <- ScaleData(DT_2)
DT_2 <- RunPCA(DT_2)
DT_2 <- FindNeighbors(DT_2, dims = 1:36)
DT_2 <- FindClusters(DT_2, resolution = 0.8)
DT_2 <- RunUMAP(DT_2, dims = 1:36)
DimPlot(DT_2, reduction = "umap", label = TRUE) + ggtitle("DT_2 Clusters (Pre-Doublet Removal)")

# Doublet removal for DT_2
sweep.res.list_DT_2 <- paramSweep(DT_2, PCs = 1:36, sct = FALSE)
sweep.stats_DT_2 <- summarizeSweep(sweep.res.list_DT_2, GT = FALSE)
bcmvn_DT_2 <- find.pK(sweep.stats_DT_2)
pK_bcmvn_DT_2 <- bcmvn_DT_2$pK[which.max(bcmvn_DT_2$BCmetric)] %>% as.character() %>% as.numeric()

DoubletRate <- ncol(DT_2) * 8 * 1e-6 
homotypic.prop <- modelHomotypic(DT_2$seurat_clusters)
nExp_poi <- round(DoubletRate * ncol(DT_2)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

DT_2 <- doubletFinder(DT_2, PCs = 1:10, pN = 0.25, pK = 0.2, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DT_2 <- doubletFinder(DT_2, PCs = 1:10, pN = 0.25, pK = 0.2, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.2_764", sct = FALSE)
saveRDS(DT_2, file = "D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/DT_2.rds")

## 2.4 Process Sample PBS_2 (Following the same workflow)
PBS_2.data <- Read10X(data.dir = "F:/20241119/single cell seq/80-1652712019/Report/Result/EXP/01_Cellranger/PBS-2/filtered_feature_bc_matrix/")
PBS_2 <- CreateSeuratObject(counts = PBS_2.data, project = "PBS_2", min.cells = 3, min.features = 200)
PBS_2[["percent.mt"]] <- PercentageFeatureSet(PBS_2, pattern = "^mt-")
PBS_2 <- subset(PBS_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
PBS_2 <- NormalizeData(PBS_2)
PBS_2 <- FindVariableFeatures(PBS_2)
PBS_2 <- ScaleData(PBS_2)
PBS_2 <- RunPCA(PBS_2)
PBS_2 <- FindNeighbors(PBS_2, dims = 1:30)  # Using 30 PCs based on variance explained
PBS_2 <- FindClusters(PBS_2, resolution = 0.8)
PBS_2 <- RunUMAP(PBS_2, dims = 1:30)
DimPlot(PBS_2, reduction = "umap", label = TRUE) + ggtitle("PBS_2 Clusters (Pre-Doublet Removal)")

# Doublet removal for PBS_2
sweep.res.list_PBS_2 <- paramSweep(PBS_2, PCs = 1:30, sct = FALSE)
sweep.stats_PBS_2 <- summarizeSweep(sweep.res.list_PBS_2, GT = FALSE)
bcmvn_PBS_2 <- find.pK(sweep.stats_PBS_2)
pK_bcmvn_PBS_2 <- bcmvn_PBS_2$pK[which.max(bcmvn_PBS_2$BCmetric)] %>% as.character() %>% as.numeric()

DoubletRate <- ncol(PBS_2) * 8 * 1e-6 
homotypic.prop <- modelHomotypic(PBS_2$seurat_clusters)
nExp_poi <- round(DoubletRate * ncol(PBS_2)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

PBS_2 <- doubletFinder(PBS_2, PCs = 1:30, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
PBS_2 <- doubletFinder(PBS_2, PCs = 1:30, pN = 0.25, pK = 0.04, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.04_1321", sct = FALSE)
saveRDS(PBS_2, file = "D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/PBS_2.rds")

# ------------------------------------------------------
# 3. Sample Integration
# ------------------------------------------------------
# Integrate all 4 samples to correct for batch effects and enable cross-sample comparisons.

# 3.1 Load processed single-sample objects
PBS_1 <- readRDS("D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/PBS_1.rds")
PBS_2 <- readRDS("D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/PBS_2.rds")
DT_1 <- readRDS("D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/DT_1.rds")
DT_2 <- readRDS("D:/OneDrive/数据/immune privilege/生信分析/20241114_mice sc-seq/single cell analysis_by myself/DT_2.rds")

# 3.2 Extract singlets (remove doublets identified by DoubletFinder)
PBS_1_singlet <- subset(PBS_1, subset = DF.classifications_0.25_0.005_1471 == "Singlet")
PBS_2_singlet <- subset(PBS_2, subset = DF.classifications_0.25_0.04_1258 == "Singlet")
DT_1_singlet <- subset(DT_1, subset = DF.classifications_0.25_0.005_532 == "Singlet")
DT_2_singlet <- subset(DT_2, subset = DF.classifications_0.25_0.2_714 == "Singlet")

# 3.3 Merge samples into a single Seurat object
# Add unique cell IDs per sample and set project name
mice_PD1 <- merge(PBS_1_singlet, y = c(PBS_2_singlet, DT_1_singlet, DT_2_singlet), 
                  add.cell.ids = c("PBS_1", "PBS_2", "DT_1", "DT_2"), project = "mice_PD1")
mice_PD1 <- subset(mice_PD1, subset = percent.mt < 10)
# Verify cell counts per sample
table(mice_PD1$orig.ident)
dim(mice_PD1)  # Total number of cells after merging

# 3.4 Pre-integration analysis to assess batch effects
# Normalize, select variable features, and scale merged data
mice_PD1 <- NormalizeData(mice_PD1, normalization.method = "LogNormalize", scale.factor = 10000)
mice_PD1 <- FindVariableFeatures(mice_PD1)
mice_PD1 <- ScaleData(mice_PD1)
mice_PD1 <- RunPCA(mice_PD1)

# Generate UMAP of unintegrated data to visualize batch effects
mice_PD1 <- FindNeighbors(mice_PD1, dims = 1:30, reduction = "pca")
mice_PD1 <- FindClusters(mice_PD1, resolution = 2, cluster.name = "unintegrated_clusters")  # High resolution for batch check
mice_PD1 <- RunUMAP(mice_PD1, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Plot unintegrated clusters and sample distribution
DimPlot(mice_PD1, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = TRUE) + ggtitle("Unintegrated Clusters")
DimPlot(mice_PD1, reduction = "umap.unintegrated", group.by = "orig.ident") + ggtitle("Sample Distribution (Pre-Integration)")

# 3.5 Integrate samples using Canonical Correlation Analysis (CCA)
# Corrects batch effects while preserving biological variation (Seurat v5)
mice_PD1 <- IntegrateLayers(object = mice_PD1, method = CCAIntegration, 
                            orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
mice_PD1[["RNA"]] <- JoinLayers(mice_PD1[["RNA"]])  # Merge layers post-integration

# ------------------------------------------------------
# 4. Cell Type Annotation and Refinement
# ------------------------------------------------------
# Annotate cell types based on canonical markers and refine clusters.

# 4.1 Coarse clustering on integrated data
mice_PD1 <- FindNeighbors(mice_PD1, reduction = "integrated.cca", dims = 1:10)  # Use integrated CCA reduction
mice_PD1 <- FindClusters(mice_PD1, resolution = 0.3)  # Coarse resolution for broad cell types
mice_PD1 <- RunUMAP(mice_PD1, dims = 1:10, reduction = "integrated.cca")  # UMAP on integrated data
DimPlot(mice_PD1, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.2) + ggtitle("Integrated Clusters (Coarse)")

# 4.2 Manual cell type annotation (coarse)
# Assign broad cell type labels based on canonical markers (e.g., myeloid: Itgam, lymphocytes: Cd3e, keratinocytes: Krt14)
names_3 <- c(
  "Myeloid cell", "Keratinocyte", "Fibroblast", "Lymphoid cell", "Lymphoid cell", 
  "Myeloid cell", "Lymphoid cell","Myeloid cell", "Myeloid cell", "Myeloid cell", 
  "Endothelial cell", "Melanocyte", "Myeloid cell")
mice_PD1@meta.data$names_3 <- mice_PD1@meta.data$seurat_clusters
levels(mice_PD1@meta.data$names_3) <- names_3

# Visualize coarse annotations
cluster_colors3 <- c(
  "Myeloid cell" = "#E6DF84", "Keratinocyte" = "#AEC7E8", "Lymphoid cell" = "#C5B0D5", 
  "Fibroblast" = "#a4d38e", "Endothelial cell" = "#E377C2", "Melanocyte" = "#8C564B")
DimPlot(mice_PD1, reduction = "umap", group.by = "names_3", pt.size = 0.2) + 
  scale_color_manual(values = cluster_colors3) + ggtitle("Coarse Cell Type Annotations")

# 4.3 Refine clusters by removing artifactual cells
# Manually select and remove ambiguous/artifact cells using interactive CellSelector
DimPlot(mice_PD1, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
plot <- DimPlot(mice_PD1, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
Doublet <- CellSelector(plot = plot)  # Interactively select cells to remove
Doublet2 <- CellSelector(plot = plot)
Doublet3 <- CellSelector(plot = plot)
Doublets <- c(Doublet, Doublet2, Doublet3)
mice_PD1_sc <- subset(mice_PD1, cells = setdiff(Cells(mice_PD1), Doublets))  # Retain clean cells

# 4.4 Re-integrate and fine-cluster refined dataset
mice_PD1_sc[["RNA"]] <- split(mice_PD1_sc[["RNA"]], f = mice_PD1_sc$orig.ident)
mice_PD1_sc <- IntegrateLayers(
  object = mice_PD1_sc, method = RPCAIntegration,  # RPCA for improved integration of refined data
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

# Fine clustering with higher resolution
mice_PD1_sc <- FindNeighbors(mice_PD1_sc, reduction = "integrated.rpca", dims = 1:45)
mice_PD1_sc <- FindClusters(mice_PD1_sc, resolution = 1)  # Higher resolution for subpopulations
mice_PD1_sc <- RunUMAP(mice_PD1_sc, dims = 1:45, reduction = "integrated.rpca")
DimPlot(mice_PD1_sc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Fine Clusters (Refined)")

# 4.5 Fine-grained cell type annotation
# Assign specific cell types based on marker genes (e.g., Schwann cells: Mpz, melanocytes: Mitf)
names_4 <- c( "Lymphoid cell", "Myeloid cell", "Fibroblast", "Lymphoid cell", "Lymphoid cell", 
              "Keratinocyte", "Myeloid cell","Myeloid cell", "Myeloid cell", "Keratinocyte", 
              "Fibroblast",  "Keratinocyte", "Myeloid cell","Keratinocyte","Lymphoid cell",
              "Lymphoid cell","Keratinocyte","Keratinocyte","Keratinocyte","Keratinocyte",
              "Lymphoid cell","Myeloid cell","Fibroblast","Keratinocyte","Myeloid cell",
              "Lymphoid cell","Fibroblast","Keratinocyte","Myeloid cell","Keratinocyte",
              "Myeloid cell","Keratinocyte","Fibroblast","Myeloid cell","Myeloid cell",
              "Schwann cell","Melanocyte","Endothelial cell")
mice_PD1_sc@meta.data$names_4 <- mice_PD1_sc@meta.data$seurat_clusters
levels(mice_PD1_sc@meta.data$names_4) <- names_4

# Validate annotations with marker genes
FeaturePlot(mice_PD1_sc, features = c("H2-Aa", "Itgam"), reduction = "umap") + ggtitle("Marker Genes: H2-Aa (Myeloid) and Itgam (Immune)")

# Final visualization of refined annotations
cluster_colors3 <- c(
  "Myeloid cell" = "#E6DF84", "Keratinocyte" = "#AEC7E8", "Lymphoid cell" = "#C5B0D5", 
  "Fibroblast" = "#a4d38e", "Endothelial cell" = "#D62728", "Melanocyte" ="#8C564B", "Schwann cell"="#1F77B4")
DimPlot(mice_PD1_sc, reduction = "umap", group.by = "names_4", label = FALSE, pt.size = 0.2) + 
  scale_color_manual(values = cluster_colors3) + ggtitle("Final Cell Type Annotations")




