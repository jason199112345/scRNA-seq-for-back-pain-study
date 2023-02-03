devtools::install_github("kevinblighe/EnhancedVolcano")
install.packages("remotes")
remotes::install_github("dynverse/scvelo")

library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)
library("BiocFileCache")
library("writexl")
library(viridis)
library(RColorBrewer)
library(EnhancedVolcano)
library(plotly)
library(scvelo)
library(monocle3)

#[Part A]===============Reading Data and Create Seurat Objects ==================
#Each "data.dir" points to a folder containing the matrix file, the barcode file, and the feature file.

NPC68m=Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)
bpNPC1s <- CreateSeuratObject(counts = NPC68m, project = "bpNPC1", min.cells = 3, min.features = 200)

NPC69m=Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)
bpNPC2s <- CreateSeuratObject(counts = NPC69m, project = "bpNPC2", min.cells = 3, min.features = 200)

NPC85m=Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)
daNPC1s <- CreateSeuratObject(counts = NPC85m, project = "daNPC1", min.cells = 3, min.features = 200)

NPC86m=Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)
bpNPC3s <- CreateSeuratObject(counts = NPC86m, project = "bpNPC3", min.cells = 3, min.features = 200)

NPCXm=Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)
daNPC2s <- CreateSeuratObject(counts = NPCXm, project = "daNPC2", min.cells = 3, min.features = 200)

daNPC3m <- Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

daNPC3s = CreateSeuratObject(counts = daNPC3m, project = "daNPC3")

daNPC4m <- Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

daNPC4s = CreateSeuratObject(counts = daNPC4m, project = "daNPC4")


# Quality control
bpNPC1s[["percent.mt"]] <- PercentageFeatureSet(object = bpNPC1s, pattern = "^MT-")
bpNPC2s[["percent.mt"]] <- PercentageFeatureSet(object = bpNPC2s, pattern = "^MT-")
bpNPC3s[["percent.mt"]] <- PercentageFeatureSet(object = bpNPC3s, pattern = "^MT-")
daNPC1s[["percent.mt"]] <- PercentageFeatureSet(object = daNPC1s, pattern = "^MT-")
daNPC2s[["percent.mt"]] <- PercentageFeatureSet(object = daNPC2s, pattern = "^MT-")
daNPC3s[["percent.mt"]] <- PercentageFeatureSet(object = daNPC3s, pattern = "^MT-")
daNPC4s[["percent.mt"]] <- PercentageFeatureSet(object = daNPC4s, pattern = "^MT-")


VlnPlot(bpNPC1s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
bpNPC1s.qc <- subset(bpNPC1s, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
bpNPC1s.qc <- NormalizeData(bpNPC1s)

VlnPlot(bpNPC2s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
bpNPC2s.qc <- subset(bpNPC2s, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
bpNPC2s.qc <- NormalizeData(bpNPC2s)


VlnPlot(bpNPC3s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
bpNPC3s.qc <- subset(bpNPC3s, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
bpNPC3s.qc <- NormalizeData(bpNPC3s)

VlnPlot(daNPC1s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
daNPC1s.qc <- subset(daNPC1s, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
daNPC1s.qc <- NormalizeData(daNPC1s)


VlnPlot(daNPC2s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
daNPC2s.qc <- subset(daNPC2s, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
daNPC2s.qc <- NormalizeData(daNPC2s)


VlnPlot(daNPC3s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
daNPC3s.qc <- subset(daNPC3s, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
daNPC3s.qc <- NormalizeData(daNPC3s)


VlnPlot(daNPC4s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
daNPC4s.qc <- subset(daNPC4s, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
daNPC4s.qc <- NormalizeData(daNPC4s)



DataAll <- c(bpNPC1s.qc, bpNPC2s.qc, bpNPC3s.qc, daNPC1s.qc, daNPC2s.qc, daNPC3s.qc, daNPC4s.qc)


for (i in 1:length(DataAll)) {
  DataAll[[i]] <- NormalizeData(DataAll[[i]], verbose = FALSE)
  DataAll[[i]] <- FindVariableFeatures(DataAll[[i]], selection.method = "vst", 
                                       nfeatures = 6000, verbose = FALSE)
}

IVD.anchors <- FindIntegrationAnchors(object.list = DataAll, dims = 1:30)
IVD.integrated <- IntegrateData(anchorset = IVD.anchors)
DefaultAssay(IVD.integrated) <- "integrated"

IVD.integrated <- ScaleData(IVD.integrated, verbose = FALSE)
IVD.integrated <- RunPCA(IVD.integrated, npcs = 30, verbose = FALSE)
IVD.integrated <- RunUMAP(IVD.integrated, reduction = "pca", dims = 1:30)
p01 <- DimPlot(IVD.integrated, reduction = "umap", pt.size = 1, label.size = 20, cols = c("red", "red", "red", "gray80", "gray80", "gray80", "gray80")) + theme(aspect.ratio = 1) 

IVD.integrated3 <- FindNeighbors(IVD.integrated, reduction = "pca", dims = 1:30)
IVD.integrated3 <- FindClusters(IVD.integrated3, resolution = 0.5)
DefaultAssay(IVD.integrated3) <- "RNA"
p02a <- DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 10) + coord_fixed(ratio=1)
p02b <- DimPlot(IVD.integrated3, reduction = "umap", pt.size = 0.3, label = FALSE, label.size = 10) + coord_fixed(ratio=1)

IVD.integrated3$ClusterNumber <- Idents(IVD.integrated3)
Idents(IVD.integrated3) <- "orig.ident"
IVD.integrated3 <- RenameIdents(object = IVD.integrated3, "bpNPC1" = "bpNPC", "bpNPC2" = "bpNPC", 'bpNPC3' = 'bpNPC', 'daNPC1' = 'daNPC', 'daNPC2' = 'daNPC', 'daNPC3' = 'daNPC', 'daNPC4' = 'daNPC')
IVD.integrated3$SampleType <- Idents(IVD.integrated3)
table(IVD.integrated3$SampleType)
Idents(IVD.integrated3) <- "ClusterNumber"


#[Part B]=============Analyzing NPC clusters===============

#Identifying cluster identity
DotPlot(IVD.integrated3, features = c("MMP3", "SPP1", "CRTAC1", "MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", 
                                      "ACAN", "COL2A1", "SOX9",
                                      "COL1A1", "CALR", "HSPA6",
                                      "MEG3", "CD44", "CD14", "HBB", "HBA1"), scale.by = "radius", dot.scale = 6, col.min = -2, col.max = 2) + 
  theme(axis.text.x = element_text(angle = 90)) 

#Subset NPC, NPC-like cells
IVD.integrated3.NPC <- subset(x = IVD.integrated3, idents = c("0", "1", "4", "6", "7", "8", "9") )
levels(IVD.integrated3.NPC)

IVD.integrated3.NPC.C0 <- subset(x = IVD.integrated3.NPC, idents = "0")
IVD.integrated3.NPC.C1 <- subset(x = IVD.integrated3.NPC, idents = "1")
IVD.integrated3.NPC.C4 <- subset(x = IVD.integrated3.NPC, idents = "4")
IVD.integrated3.NPC.C6 <- subset(x = IVD.integrated3.NPC, idents = "6")
IVD.integrated3.NPC.C7 <- subset(x = IVD.integrated3.NPC, idents = "7")
IVD.integrated3.NPC.C8 <- subset(x = IVD.integrated3.NPC, idents = "8")
IVD.integrated3.NPC.C9 <- subset(x = IVD.integrated3.NPC, idents = "9")

Idents(IVD.integrated3.NPC) <- "SampleType"
IVD.integrated3.NPC.bp <- subset(x = IVD.integrated3.NPC, idents = "bpNPC")
IVD.integrated3.NPC.da <- subset(x = IVD.integrated3.NPC, idents = "daNPC")
Idents(IVD.integrated3.NPC) <- "ClusterNumber"
Idents(IVD.integrated3.NPC.bp) <- "ClusterNumber"
Idents(IVD.integrated3.NPC.da) <- "ClusterNumber"


levels(IVD.integrated3.NPC)

#Find in silico top markers for each cluster
Markers.integrated3.NPC.C0 <- FindMarkers(IVD.integrated3.NPC, ident.1 = "0",  min.pct = 0.1)
Markers.integrated3.NPC.C1 <- FindMarkers(IVD.integrated3.NPC, ident.1 = "1",  min.pct = 0.1)
Markers.integrated3.NPC.C4 <- FindMarkers(IVD.integrated3.NPC, ident.1 = "4",  min.pct = 0.1)
Markers.integrated3.NPC.C6 <- FindMarkers(IVD.integrated3.NPC, ident.1 = "6",  min.pct = 0.1)
Markers.integrated3.NPC.C7 <- FindMarkers(IVD.integrated3.NPC, ident.1 = "7",  min.pct = 0.1)
Markers.integrated3.NPC.C8 <- FindMarkers(IVD.integrated3.NPC, ident.1 = "8",  min.pct = 0.1)
Markers.integrated3.NPC.C9 <- FindMarkers(IVD.integrated3.NPC, ident.1 = "9",  min.pct = 0.1)
Markers.integrated3.C5 <- FindMarkers(IVD.integrated3, ident.1 = "5",  min.pct = 0.1)
Idents(IVD.integrated3) <- "SampleType"
Markers.integrated3.bp <- FindMarkers(IVD.integrated3, ident.1 = "bpNPC",  min.pct = 0.1)

#Cell Count
table(IVD.integrated3.NPC$SampleType, IVD.integrated3.NPC$ClusterNumber)
table(IVD.integrated3$SampleType, IVD.integrated3$ClusterNumber)


p03a <- DimPlot(IVD.integrated3.NPC,  reduction = "umap", 
                cols = c("gray75", "gray75",  "cyan4",
                         "steelblue3", 
                         "firebrick3", "indianred2",  "cyan3"),
                pt.size = 0.1,  label = FALSE, label.size = 5) + coord_fixed(ratio=1)

p03b <- DimPlot(IVD.integrated3.NPC,  reduction = "umap", 
                cols = c("gray75", "gray75",  "cyan4",
                         "steelblue3", 
                         "firebrick3", "indianred2",  "cyan3"),
                pt.size = 0.1,  label = TRUE, label.size = 5) + coord_fixed(ratio=1)



p03c <- FeaturePlot(IVD.integrated3.NPC, features = c("BDNF", "NGF", "HIF1A"), 
            label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.3, coord.fixed = TRUE)   

#Markers for cluster of interests in their prevalent IVDs
Markers.integrated3.NPC.bp.C7 <- FindMarkers(IVD.integrated3.NPC.bp, ident.1 = "7",  min.pct = 0.1)
Markers.integrated3.NPC.da.C9 <- FindMarkers(IVD.integrated3.NPC.bp, ident.1 = "9",  min.pct = 0.1)
Markers.integrated3.NPC.da.C6 <- FindMarkers(IVD.integrated3.NPC.bp, ident.1 = "6",  min.pct = 0.1)


#Highlight C9
p04 <- DimPlot(IVD.integrated3.NPC,  reduction = "umap", split.by = "SampleType",
               cols = c("gray75", "gray75",  "gray75",
                        "gray75", 
                        "gray75", "gray75",  "cyan3"),
               pt.size = 0.7,  label = FALSE, label.size = 5)  + coord_fixed(ratio=1)

#Analyzing C7
Idents(IVD.integrated3.NPC) <- "ClusterNumber"
IVD.integrated3.NPC.C7 <- subset(x = IVD.integrated3.NPC, idents = "7" )

IVD.integrated3.NPC.C7 <- FindVariableFeatures(IVD.integrated3.NPC.C7, selection.method = "vst")
IVD.integrated3.NPC.C7 <- ScaleData(IVD.integrated3.NPC.C7, verbose = FALSE)
IVD.integrated3.NPC.C7 <- RunPCA(IVD.integrated3.NPC.C7, npcs = 30, verbose = FALSE)
IVD.integrated3.NPC.C7 <- RunUMAP(IVD.integrated3.NPC.C7, reduction = "pca", dims = 1:30)
IVD.integrated3.NPC.C7 <- FindNeighbors(IVD.integrated3.NPC.C7, reduction = "pca", dims = 1:30)
IVD.integrated3.NPC.C7 <- FindClusters(IVD.integrated3.NPC.C7, resolution = 0.5)

Idents(IVD.integrated3.NPC.C7) <- "SampleType"
p04a1 <- DimPlot(IVD.integrated3.NPC.C7, reduction = "umap", pt.size = 0.3, label = FALSE, label.size = 10) + coord_fixed(ratio=1)

IVD.integrated3.NPC.C7.allmarkers <- FindAllMarkers(IVD.integrated3.NPC.C7, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
IVD.integrated3.NPC.C7.allmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

p04a2 <- FeaturePlot(IVD.integrated3.NPC.C7, features = c("MMP3", "SERPINA1", "MGST1", "RGCC"),
                     label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.3, coord.fixed = TRUE)   
p04a3 <-VlnPlot(IVD.integrated3.NPC.C7, features = c("MMP3", "SERPINA1", "MGST1", "RGCC"), ncol = 2) 
FeaturePlot(IVD.integrated3.NPC.C7, features = c("VEGFA", "NGF", "BDNF"),
            label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.3, coord.fixed = TRUE)

#Analyzing C6
Idents(IVD.integrated3.NPC) <- "ClusterNumber"
IVD.integrated3.NPC.C6 <- subset(x = IVD.integrated3.NPC, idents = "6" )

IVD.integrated3.NPC.C6 <- FindVariableFeatures(IVD.integrated3.NPC.C6, selection.method = "vst")
IVD.integrated3.NPC.C6 <- ScaleData(IVD.integrated3.NPC.C6, verbose = FALSE)
IVD.integrated3.NPC.C6 <- RunPCA(IVD.integrated3.NPC.C6, npcs = 30, verbose = FALSE)
IVD.integrated3.NPC.C6 <- RunUMAP(IVD.integrated3.NPC.C6, reduction = "pca", dims = 1:30)
IVD.integrated3.NPC.C6 <- FindNeighbors(IVD.integrated3.NPC.C6, reduction = "pca", dims = 1:30)
IVD.integrated3.NPC.C6 <- FindClusters(IVD.integrated3.NPC.C6, resolution = 0.5)

Idents(IVD.integrated3.NPC.C6) <- "SampleType"
p04b1 <- DimPlot(IVD.integrated3.NPC.C6, reduction = "umap", pt.size = 0.3, label = FALSE, label.size = 10) + coord_fixed(ratio=1)


p04b2 <- FeaturePlot(IVD.integrated3.NPC.C6, features = c("XIST", "COMP", "MT-ND4L", "MT-CO3"),
                     label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.3, coord.fixed = TRUE)   

p04b3 <-VlnPlot(IVD.integrated3.NPC.C6, features = c("XIST", "COMP", "MT-ND4L", "MT-CO3"), ncol = 2) 


#Analyzing C9
Idents(IVD.integrated3.NPC) <- "ClusterNumber"
IVD.integrated3.NPC.C9 <- subset(x = IVD.integrated3.NPC, idents = "9" )

IVD.integrated3.NPC.C9 <- FindVariableFeatures(IVD.integrated3.NPC.C9, selection.method = "vst")
IVD.integrated3.NPC.C9 <- ScaleData(IVD.integrated3.NPC.C9, verbose = FALSE)
IVD.integrated3.NPC.C9 <- RunPCA(IVD.integrated3.NPC.C9, npcs = 30, verbose = FALSE)
IVD.integrated3.NPC.C9 <- RunUMAP(IVD.integrated3.NPC.C9, reduction = "pca", dims = 1:30)
IVD.integrated3.NPC.C9 <- FindNeighbors(IVD.integrated3.NPC.C9, reduction = "pca", dims = 1:30)
IVD.integrated3.NPC.C9 <- FindClusters(IVD.integrated3.NPC.C9, resolution = 0.5)

Idents(IVD.integrated3.NPC.C9) <- "SampleType"
p04d1 <- DimPlot(IVD.integrated3.NPC.C9, reduction = "umap", pt.size = 0.3, label = FALSE, label.size = 10) + coord_fixed(ratio=1)


p04d2 <- FeaturePlot(IVD.integrated3.NPC.C9, features = c("CILP", "MTRNR2L8", "AKR1C1", "RGCC"),
                     label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.3, coord.fixed = TRUE)   

p04d3 <-VlnPlot(IVD.integrated3.NPC.C9, features = c("CILP", "MTRNR2L8", "AKR1C1", "RGCC"), ncol = 2) 




#Volcano

Idents(IVD.integrated3.NPC.C7) <- "SampleType"
Markers.integrated3.NPC.C7.bpvsda <- FindMarkers(IVD.integrated3.NPC.C7, ident.1 = "bpNPC",  min.pct = 0.1)
p05a1 <- EnhancedVolcano(Markers.integrated3.NPC.C7.bpvsda,
                        lab = rownames(Markers.integrated3.NPC.C7.bpvsda),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = "In C7, bpNPC vs daNPC",
                        labSize = 4,
                        pointSize = 2.5,
                        pCutoff = 10e-36,
                        FCcutoff = 1,
                        col = c("grey30", "grey30", "grey30", "red2"),
                        colAlpha = 3/5,
                        legendPosition = 'none',
                        drawConnectors = TRUE,
                        max.overlaps = "inf",
                        maxoverlapsConnectors = 60,
                        widthConnectors = 0.5,    
                        colConnectors = 'black',
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE
) + theme(aspect.ratio = 1)
p05a2 <- EnhancedVolcano(Markers.integrated3.NPC.C7.bpvsda,
                         lab = rownames(Markers.integrated3.NPC.C7.bpvsda),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = "In C7, bpNPC vs daNPC",
                         labSize = 4,
                         pointSize = 2.5,
                         pCutoff = 10e-36,
                         FCcutoff = 1,
                         col = c("grey30", "grey30", "grey30", "red2"),
                         colAlpha = 3/5,
                         legendPosition = 'none',
                         drawConnectors = TRUE,
                         max.overlaps = "inf",
                         selectLab = c(""),
                         maxoverlapsConnectors = 60,
                         widthConnectors = 0.5,    
                         colConnectors = 'black',
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE
) + theme(aspect.ratio = 1)


Idents(IVD.integrated3.NPC.C6) <- "SampleType"
Markers.integrated3.NPC.C6.davsbp <- FindMarkers(IVD.integrated3.NPC.C6, ident.1 = "daNPC",  min.pct = 0.1)
p05b1 <- EnhancedVolcano(Markers.integrated3.NPC.C6.davsbp,
                        lab = rownames(Markers.integrated3.NPC.C6.davsbp),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = "In C6, daNPC vs bpNPC",
                        labSize = 4,
                        pointSize = 2.5,
                        pCutoff = 10e-4,
                        FCcutoff = 1,
                        col = c("grey30", "grey30", "grey30", "red2"),
                        colAlpha = 3/5,
                        legendPosition = 'none',
                        drawConnectors = TRUE,
                        max.overlaps = "inf",
                        maxoverlapsConnectors = 60,
                        widthConnectors = 0.5,    
                        colConnectors = 'black',
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE
) + theme(aspect.ratio = 1)
p05b2 <- EnhancedVolcano(Markers.integrated3.NPC.C6.davsbp,
                         lab = rownames(Markers.integrated3.NPC.C6.davsbp),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = "In C6, daNPC vs bpNPC",
                         labSize = 4,
                         pointSize = 2.5,
                         pCutoff = 10e-4,
                         FCcutoff = 1,
                         col = c("grey30", "grey30", "grey30", "red2"),
                         colAlpha = 3/5,
                         legendPosition = 'none',
                         drawConnectors = FALSE,
                         selectLab = c(""),
                         maxoverlapsConnectors = 60,
                         widthConnectors = 0.5,    
                         colConnectors = 'black',
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE
) + theme(aspect.ratio = 1)



Idents(IVD.integrated3.NPC.C9) <- "SampleType"
Markers.integrated3.NPC.C9.davsbp <- FindMarkers(IVD.integrated3.NPC.C9, ident.1 = "daNPC",  min.pct = 0.1)
p05c1 <- EnhancedVolcano(Markers.integrated3.NPC.C9.davsbp,
                        lab = rownames(Markers.integrated3.NPC.C9.davsbp),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = "In C9, daNPC vs bpNPC",
                        labSize = 4,
                        pointSize = 2.5,
                        pCutoff = 10e-4,
                        FCcutoff = 1,
                        xlim = c(-3,3),
                        ylim = c(0,6),
                        col = c("grey30", "grey30", "grey30", "red2"),
                        colAlpha = 3/5,
                        legendPosition = 'none',
                        drawConnectors = TRUE,
                        max.overlaps = "inf",
                        maxoverlapsConnectors = 60,
                        widthConnectors = 0.5,    
                        colConnectors = 'black',
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE
) + theme(aspect.ratio = 1)
p05c2 <- EnhancedVolcano(Markers.integrated3.NPC.C9.davsbp,
                        lab = rownames(Markers.integrated3.NPC.C9.davsbp),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        title = "In C9, daNPC vs bpNPC",
                        labSize = 4,
                        pointSize = 2.5,
                        pCutoff = 10e-4,
                        FCcutoff = 1,
                        xlim = c(-3,3),
                        ylim = c(0,6),
                        col = c("grey30", "grey30", "grey30", "red2"),
                        colAlpha = 3/5,
                        legendPosition = 'none',
                        drawConnectors = FALSE,
                        selectLab = c(""),
                        maxoverlapsConnectors = 60,
                        widthConnectors = 0.5,    
                        colConnectors = 'black',
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE
) + theme(aspect.ratio = 1)



#Dotplot
levels(x = IVD.integrated3.NPC) <- c("7", "8", "0", "1", "4", "9", "6") 

p06a <- DotPlot(IVD.integrated3.NPC, features = c(
                                    "MT-ATP6",
                                    "PRG4", 
                                    "HSPA6",
                                    "CILP",
                                    "S100A2",
                                    "SPP1",
                                    "MMP3"
                                    
                                    ), scale.by = "radius", dot.scale = 9) +coord_flip() 

p06b <- DotPlot(IVD.integrated3.NPC, features = c(  "COL1A1",   "SOX9",
                                                    "ACAN", 
                                                    "COL2A1"
), scale.by = "radius", dot.scale = 9) +coord_flip() 





#Dual markers

p07 <- FeaturePlot(IVD.integrated3.NPC, features = c("SLC7A2", "TM4SF1"),  blend = TRUE, blend.threshold = 0.1,
                    label = FALSE,  order = TRUE, cols = c("black", "red", "blue"), split.by = "SampleType",  pt.size = 1, coord.fixed = TRUE)    + coord_fixed(ratio=1)

#[Part C]============Analyzing gene of interest

GoI <- c( "MMP3", 
          "NGF", "BDNF", "CRCP", "CGRP", "TRPV4", "VEGFA",
          "CXCL8", "IL6", "IL1B", 
          "TNF",  "TNFA", "TNFAIP1", "TNFSF4", 
          "TRPC6",
          "NOS2","PRDX1", "HIF1A", "BNIP3", "PRG4",
          "ASIC3",  "SEMA3A",  "SHH", "ADAMTS5", "FGF2" ,
          "CD9", "CD109",
          "KRT9", "KRT19", "MIR155HG")

#Calculate average expression of gene expression


AverageExpression(
  IVD.integrated3.NPC.da,
  assays = "RNA",
  features = GoI,
  return.seurat = FALSE,
  group.by = "ClusterNumber",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)

AverageExpression(
  IVD.integrated3.NPC.bp,
  assays = "RNA",
  features = GoI,
  return.seurat = FALSE,
  group.by = "ClusterNumber",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)

AverageExpression(
  IVD.integrated3.NPC,
  assays = "RNA",
  features = GoI,
  return.seurat = FALSE,
  group.by = "SampleType",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)

Idents(IVD.integrated3.NPC) <- "SampleType"
p08a <- VlnPlot(IVD.integrated3.NPC, features = c("MMP3", "NGF", "CXCL8", "PRDX1", "HIF1A", "NOS2", "BNIP3", "FGF2", "VEGFA", "MIR155HG"), ncol = 5, pt.size = 0.01, log = FALSE) 
p08b <- FeaturePlot(IVD.integrated3.NPC, features = c("MMP3", "NGF", "HIF1A"), label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.3, coord.fixed = TRUE)   
p08c <- FeaturePlot(IVD.integrated3.NPC, features = c("NOS2", "BNIP3", "PRDX1"), label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.3, coord.fixed = TRUE)   
p08d <- FeaturePlot(IVD.integrated3.NPC, features = c("FGF2", "VEGFA", "MIR155HG"), label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.3, coord.fixed = TRUE)   
Idents(IVD.integrated3.NPC) <- "ClusterNumber"



#Export for Qiangen IPA Pathway Analysis
MarkersToExcel.Markers.integrated3.NPC.C7.bpvsda <- cbind("Gene ID"=rownames(Markers.integrated3.NPC.C7.bpvsda), Markers.integrated3.NPC.C7.bpvsda)
write_xlsx(MarkersToExcel.Markers.integrated3.NPC.C7.bpvsda,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 3/IPA pathways/C7.bp.xlsx")

MarkersToExcel.Markers.integrated3.NPC.C9.davsbp <- cbind("Gene ID"=rownames(Markers.integrated3.NPC.C9.davsbp), Markers.integrated3.NPC.C9.davsbp)
write_xlsx(MarkersToExcel.Markers.integrated3.NPC.C9.davsbp,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 3/IPA pathways/C9.da.xlsx")

MarkersToExcel.Markers.integrated3.NPC.C6.davsbp <- cbind("Gene ID"=rownames(Markers.integrated3.NPC.C6.davsbp), Markers.integrated3.NPC.C6.davsbp)
write_xlsx(MarkersToExcel.Markers.integrated3.NPC.C6.davsbp,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 3/IPA pathways/C6.da.xlsx")


#[Part C]====Monocle 3 Trajectory=======
#Install and prepare Monocle 3 package
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
install.packages("remotes")
remotes::install_github("satijalab/seurat-wrappers")

#Subset C6, C7, C9
IVD.integrated3.NPC.traj <- subset(x = IVD.integrated3.NPC, idents = c("4", "6", "7", "8", "9") )
levels(IVD.integrated3.NPC.traj)
#Load Seurat into Monocle3
data <- IVD.integrated3.NPC.traj@assays$RNA@data
cell_metadata <- new('AnnotatedDataFrame', data = IVD.integrated3.NPC.traj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
gene_metadata <- new('AnnotatedDataFrame', data = fData)

gene_metadata <- as(gene_metadata, "data.frame")
cell_metadata <- as(cell_metadata, "data.frame")

cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_metadata
)

#-----------Data processing and pre-analysis------------
#Pre-processing
cds <- preprocess_cds(cds, num_dim = 100)
T01 <- plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

#Remove batch effect
cds = align_cds(cds, num_dim = 100, alignment_group = "orig.ident")
#Run the align_cds code to show this citation information:  Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). 'Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.' Nat. Biotechnol., 36(5), 421-427. doi: 10.1038/nbt.4091

cds = reduce_dimension(cds)
T02a <- plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE)
T02b <- plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=TRUE, group_label_size = 7) 

#Clustering
cds = cluster_cells(cds, resolution=1e-5)
T03a <- plot_cells(cds, group_label_size = 7, label_cell_groups = FALSE)
T03b <- plot_cells(cds, group_label_size = 7, label_cell_groups=TRUE) 

#Identifying top markers by SampleType
marker_test_res_SampleType <- top_markers(cds, group_cells_by="SampleType", 
                               reference_cells=1000, cores=8)


top_specific_markers_SampleType <- marker_test_res_SampleType %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

top_specific_markers_SampleType_ids <- unique(top_specific_markers_SampleType %>% pull(gene_id))


#Identifying top markers by ClusterNumber
marker_test_res_ClusterNumber <- top_markers(cds, group_cells_by="ClusterNumber", 
                                             reference_cells=1000, cores=8)

top_specific_markers_ClusterNumber <- marker_test_res_ClusterNumber %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_markers_ClusterNumber_ids <- unique(top_specific_markers_ClusterNumber %>% pull(gene_id))


#Plot gene expression to identify cell types

plot_cells(cds, genes="MMP3",  
           show_trajectory_graph=TRUE, 
           label_cell_groups=TRUE)

plot_cells(cds, genes="SLC7A2", 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)
plot_cells(cds, genes="TM4SF1", 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)
plot_cells(cds, genes="NEAT1", 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)


#Subset (if needed)
cds_subset <- choose_cells(cds)

#---------Constructing Single-Cell Trajectory---------
cds <- learn_graph(cds)
T04a <- plot_cells(cds, cell_size = 0.6,
           color_cells_by = "SampleType",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_cell_groups=FALSE)  + coord_fixed(ratio=1)


#Pseudo-time
T04b <- plot_cells(cds, cell_size = 0.6,
           color_cells_by = "ClusterNumber",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_cell_groups=FALSE) + coord_fixed(ratio=1)


cds <- order_cells(cds)

T05a <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_cell_groups=FALSE) + coord_fixed(ratio=1)

T05a2 <- plot_cells(cds,
                   color_cells_by = "pseudotime",
                   label_groups_by_cluster=TRUE,
                   label_leaves=TRUE,
                   label_branch_points=TRUE,
                   graph_label_size=4,
                   label_cell_groups=FALSE) + coord_fixed(ratio=1)

T05a3 <- plot_cells(cds,
                    color_cells_by = "SampleType",
                    label_groups_by_cluster=TRUE,
                    label_leaves=TRUE,
                    label_branch_points=TRUE,
                    graph_label_size=4,
                    label_cell_groups=FALSE) + coord_fixed(ratio=1)

plot_cells(cds,
           genes="MMP3",
           cell_size = 1,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE) + coord_fixed(ratio=1)


#Identifying genes that change with time
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=64)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 10e-16))
T07a <- plot_cells(cds,
           genes="NGF",
           cell_size = 1,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE) + coord_fixed(ratio=1)

T07b <-  plot_cells(cds,
                    genes="PRDX1",
                    cell_size = 1,
                    label_branch_points = FALSE,
                    label_roots = FALSE,
                    label_leaves = FALSE,
                    label_cell_groups=FALSE,
                    show_trajectory_graph=TRUE) + coord_fixed(ratio=1)

T07c <- plot_cells(cds,
           genes="VEGFA",
           cell_size = 1,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE) + coord_fixed(ratio=1)

T07d <- plot_cells(cds,
           genes="HIF1A",
           cell_size = 1,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE) + coord_fixed(ratio=1)

T07e <- plot_cells(cds,
           genes="CXCL8",
           cell_size = 1,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE) + coord_fixed(ratio=1)

T07f <- plot_cells(cds,
           genes="NOS2",
           cell_size = 1,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE) + coord_fixed(ratio=1)


T07g <- plot_cells(cds,
           genes="MMP3",
           cell_size = 1,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE) + coord_fixed(ratio=1)

#---------- Plot genes in pseudotime----------
ciliated_genes <- c("PRDX1", "FGF2", 
                    "NGF", "VEGFA", "HIF1A")


cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

T08 <- plot_genes_in_pseudotime(cds_subset, cell_size = 0, ncol = 3, vertical_jitter = TRUE, horizontal_jitter = FALSE,  trend_formula = "~ splines::ns(pseudotime, df=5)",
                         color_cells_by="SampleType",
                         min_expr=NULL)
