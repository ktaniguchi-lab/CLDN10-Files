library(Seurat)
d1 <- CreateSeuratObject(Read10X_h5(filename = "Gel3D_D1_filtered_feature_bc_matrix.h5",
                                              use.names = TRUE,
                       unique.features = TRUE), project = "D1")
d2 <- CreateSeuratObject(Read10X_h5(filename = "Gel3D_D2_filtered_feature_bc_matrix.h5",
                                              use.names = TRUE,
                       unique.features = TRUE), project = "D2")
d3 <- CreateSeuratObject(Read10X_h5(filename = "Gel3D_D3_filtered_feature_bc_matrix.h5",
                                           use.names = TRUE,
                                           unique.features = TRUE), project = "D3")
d4 <- CreateSeuratObject(Read10X_h5(filename = "Gel3D_D4_filtered_feature_bc_matrix.h5",
                                           use.names = TRUE,
                                           unique.features = TRUE), project = "D4")
a = merge(d1, list(d2, d3, d4))
a[["percent.mt"]] <- PercentageFeatureSet(object = a, pattern = "^MT-")
VlnPlot(object = a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
a <- subset(x = a, subset = nFeature_RNA > 1500 & percent.mt > 3 & percent.mt < 15 & nCount_RNA < 30000)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
a <- CellCycleScoring(a, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
a <- NormalizeData(a)
a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(a)
a <- ScaleData(a, features = all.genes)
a <- ScaleData(a, vars.to.regress = c("S.Score", "G2M.Score"))
a <- RunPCA(a, features = VariableFeatures(object = a))
a <- FindNeighbors(a, dims = 1:10)
a <- FindClusters(a, resolution = 0.2)
a <- RunUMAP(a, dims = 1:15)