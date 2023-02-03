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
#NPCs1:sNPC
#NPCs2:nsNPC

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

NPCs1m <- Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

NPCs1s = CreateSeuratObject(counts = NPCs1m, project = "NPCs1")

NPCs2m <- Read10X(
  data.dir = "/data directory/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

NPCs2s = CreateSeuratObject(counts = NPCs2m, project = "NPCs2")


# Quality control
bpNPC1s[["percent.mt"]] <- PercentageFeatureSet(object = bpNPC1s, pattern = "^MT-")
bpNPC2s[["percent.mt"]] <- PercentageFeatureSet(object = bpNPC2s, pattern = "^MT-")
bpNPC3s[["percent.mt"]] <- PercentageFeatureSet(object = bpNPC3s, pattern = "^MT-")
daNPC1s[["percent.mt"]] <- PercentageFeatureSet(object = daNPC1s, pattern = "^MT-")
daNPC2s[["percent.mt"]] <- PercentageFeatureSet(object = daNPC2s, pattern = "^MT-")
daNPC3s[["percent.mt"]] <- PercentageFeatureSet(object = daNPC3s, pattern = "^MT-")
daNPC4s[["percent.mt"]] <- PercentageFeatureSet(object = daNPC4s, pattern = "^MT-")
NPCs1s[["percent.mt"]] <- PercentageFeatureSet(object = NPCs1s, pattern = "^MT-")
NPCs2s[["percent.mt"]] <- PercentageFeatureSet(object = NPCs2s, pattern = "^MT-")


VlnPlot(bpNPC1s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
bpNPC1s.qc <- subset(bpNPC1s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
bpNPC1s.qc <- NormalizeData(bpNPC1s)

VlnPlot(bpNPC2s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
bpNPC2s.qc <- subset(bpNPC2s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
bpNPC2s.qc <- NormalizeData(bpNPC2s)


VlnPlot(bpNPC3s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
bpNPC3s.qc <- subset(bpNPC3s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
bpNPC3s.qc <- NormalizeData(bpNPC3s)

VlnPlot(daNPC1s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
daNPC1s.qc <- subset(daNPC1s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
daNPC1s.qc <- NormalizeData(daNPC1s)


VlnPlot(daNPC2s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
daNPC2s.qc <- subset(daNPC2s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
daNPC2s.qc <- NormalizeData(daNPC2s)


VlnPlot(daNPC3s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
daNPC3s.qc <- subset(daNPC3s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
daNPC3s.qc <- NormalizeData(daNPC3s)


VlnPlot(daNPC4s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
daNPC4s.qc <- subset(daNPC4s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
daNPC4s.qc <- NormalizeData(daNPC4s)

VlnPlot(NPCs1s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
NPCs1s.qc <- subset(NPCs1s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
NPCs1s.qc <- NormalizeData(NPCs1s)

VlnPlot(NPCs2s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
NPCs2s.qc <- subset(NPCs2s, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
NPCs2s.qc <- NormalizeData(NPCs2s)

#=============A1. All sample analysis===============
data1 <- c(bpNPC1s.qc, bpNPC2s.qc, bpNPC3s.qc, daNPC1s.qc, daNPC2s.qc, daNPC3s.qc, daNPC4s.qc, NPCs1s.qc, NPCs2s.qc)

for (i in 1:length(data1)) {
  data1[[i]] <- NormalizeData(data1[[i]], verbose = FALSE)
  data1[[i]] <- FindVariableFeatures(data1[[i]], selection.method = "vst", 
                                     nfeatures = 6000, verbose = FALSE)
}

IVD.anchors <- FindIntegrationAnchors(object.list = data1, dims = 1:30)
IVD.integrated <- IntegrateData(anchorset = IVD.anchors)
DefaultAssay(IVD.integrated) <- "integrated"

IVD.integrated <- ScaleData(IVD.integrated, verbose = FALSE)
IVD.integrated <- RunPCA(IVD.integrated, npcs = 30, verbose = FALSE)
IVD.integrated <- RunUMAP(IVD.integrated, reduction = "pca", dims = 1:30)

IVD.integrated <- FindNeighbors(IVD.integrated, reduction = "pca", dims = 1:30)
IVD.integrated <- FindClusters(IVD.integrated, resolution = 0.5)

DefaultAssay(IVD.integrated) <- "RNA"
IVD.integrated$ClusterNumber <- Idents(IVD.integrated)
Idents(IVD.integrated) <- "orig.ident"
IVD.integrated <- RenameIdents(object = IVD.integrated, "bpNPC1" = "bpNPC", "bpNPC2" = "bpNPC", 'bpNPC3' = 'bpNPC','daNPC1' = 'aNPC', 'daNPC2' = 'aNPC', 'daNPC3' = 'aNPC', 'daNPC4' = 'aNPC', 'NPCs1' = 'NPCs1', 'NPCs2' = 'NPCs2')
IVD.integrated$SampleType <- Idents(IVD.integrated)
table(IVD.integrated$SampleType)
Idents(IVD.integrated) <- "ClusterNumber"
Idents(IVD.integrated) <- "SampleType"


#Dot plot and cluster markers to assigne cell type
DotPlot(IVD.integrated, features = c("MMP3", "SPP1", "CRTAC1", "MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", 
                                            "ACAN", "COL2A1", "SOX9",
                                            "COL1A1", "CALR", "HSPA6",
                                            "MEG3", "CD44", "CD14", "HBB", "HBA1"), scale.by = "radius", dot.scale = 6, col.min = -2, col.max = 2) + 
  theme(axis.text.x = element_text(angle = 90)) 


Markers.IVD.integrated.C0 <- FindMarkers(IVD.integrated, ident.1 = "0",  min.pct = 0.1)
Markers.IVD.integrated.C1 <- FindMarkers(IVD.integrated, ident.1 = "1",  min.pct = 0.1)
Markers.IVD.integrated.C2 <- FindMarkers(IVD.integrated, ident.1 = "2",  min.pct = 0.1)
Markers.IVD.integrated.C3 <- FindMarkers(IVD.integrated, ident.1 = "3",  min.pct = 0.1)
Markers.IVD.integrated.C4 <- FindMarkers(IVD.integrated, ident.1 = "4",  min.pct = 0.1)
Markers.IVD.integrated.C5 <- FindMarkers(IVD.integrated, ident.1 = "5",  min.pct = 0.1)
Markers.IVD.integrated.C6 <- FindMarkers(IVD.integrated, ident.1 = "6",  min.pct = 0.1)
Markers.IVD.integrated.C7 <- FindMarkers(IVD.integrated, ident.1 = "7",  min.pct = 0.1)
Markers.IVD.integrated.C8 <- FindMarkers(IVD.integrated, ident.1 = "8",  min.pct = 0.1)
Markers.IVD.integrated.C9 <- FindMarkers(IVD.integrated, ident.1 = "9",  min.pct = 0.1)
Markers.IVD.integrated.C10 <- FindMarkers(IVD.integrated, ident.1 = "10",  min.pct = 0.1)
Markers.IVD.integrated.C11 <- FindMarkers(IVD.integrated, ident.1 = "11",  min.pct = 0.1)
Markers.IVD.integrated.C12 <- FindMarkers(IVD.integrated, ident.1 = "12",  min.pct = 0.1)
Markers.IVD.integrated.C13 <- FindMarkers(IVD.integrated, ident.1 = "13",  min.pct = 0.1)


#We should subset C0,1,4, 5, 10 for NPC analysis
IVD.integrated.NPC <- subset(x = IVD.integrated, idents =c("0", "1", "4", "5", "10") )


IVD.integrated.NPC.C5 <- subset(x = IVD.integrated.NPC, idents ="5" )
Idents(IVD.integrated.NPC) <- "SampleType"
IVD.integrated.NPC.bp <- subset(x = IVD.integrated.NPC, idents ="bpNPC" )
IVD.integrated.NPC.s1 <- subset(x = IVD.integrated.NPC, idents ="NPCs1" )
IVD.integrated.NPC.a <- subset(x = IVD.integrated.NPC, idents ="aNPC" )
IVD.integrated.NPC.s2 <- subset(x = IVD.integrated.NPC, idents ="NPCs2" )
IVD.integrated.NPC.bpands1 <- subset(x = IVD.integrated.NPC, idents =c("bpNPC", "NPCs1") )

IVD.integrated.NPC.noC5 <- subset(x = IVD.integrated, idents =c("0", "1", "4", "10") )


#-----cell count----
table(IVD.integrated.NPC$SampleType)
table(IVD.integrated.NPC$ClusterNumber)
table(IVD.integrated.NPC$SampleType, IVD.integrated.NPC$ClusterNumber)


#----Count marker+ NPC-----
Idents(IVD.integrated.NPC.bp) <- "ClusterNumber"
Idents(IVD.integrated.NPC.a) <- "ClusterNumber"
Idents(IVD.integrated.NPC.s1) <- "ClusterNumber"
Idents(IVD.integrated.NPC.s2) <- "ClusterNumber"


p <- VlnPlot(object = IVD.integrated.NPC.bp, features =c("MMP3"))
p$data %>% group_by(ident) %>% summarize(counts = sum(MMP3, na.rm = TRUE))


p01 <- DimPlot(IVD.integrated.NPC, reduction = "umap", pt.size = 0.3, split.by = "SampleType", label = TRUE, label.size = 10) + coord_fixed(ratio=1)


p02 <- DimPlot(IVD.integrated.NPC, reduction = "umap", pt.size = 0.3, split.by = "SampleType", label = FALSE, label.size = 10, cols = c("gray80", "gray80", "gray80", "red",  "gray80")) + coord_fixed(ratio=1)


#marker feature plot, single marker and dual marker
p03 <- FeaturePlot(IVD.integrated.NPC, features = c("MMP3", "NGF", "FGF2",  "CXCL8", "PRDX1", "NOS2"), 
            label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.2, coord.fixed = TRUE)   
p03b <- FeaturePlot(IVD.integrated.NPC, features = c("MMP3",  "CXCL8"),  label = FALSE,  order = TRUE, cols = brewer.pal(11, "PRGn"), split.by = "SampleType",  pt.size = 0.2, coord.fixed = TRUE)
                   

p04 <- FeaturePlot(IVD.integrated.NPC, features = c("SLC7A2", "TM4SF1"),  blend = TRUE, blend.threshold = 0.1,
            label = FALSE,  order = TRUE, cols = c("black", "red", "blue"), split.by = "SampleType",  pt.size = 1, coord.fixed = TRUE)    + coord_fixed(ratio=1)

AverageExpression(
  IVD.integrated.NPC,
  assays = "RNA",
  features = c("NGF", "NOS2", "IL6", "CCN2", "CXCL8", "MMP3", "FGF2", "PRDX1", "TM4SF1", "SLC7A2" ),
  return.seurat = FALSE,
  group.by = "SampleType",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)

Idents(IVD.integrated.NPC.C5) <- "SampleType"

AverageExpression(
  IVD.integrated.NPC.C5,
  assays = "RNA",
  features = c("NGF", "NOS2", "IL6", "CCN2", "CXCL8", "MMP3", "FGF2", "PRDX1", "TM4SF1", "SLC7A2", "CD36" ),
  return.seurat = FALSE,
  group.by = "SampleType",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)

#Count cell by gene expression
sum(GetAssayData(object = IVD.integrated.NPC.bp, slot = "data")["MMP3",]>0)

Idents(IVD.integrated.NPC.C5) <- "SampleType"
IVD.integrated.NPC.C5.bp <- subset(x = IVD.integrated.NPC.C5, idents ="bpNPC" )
sum(GetAssayData(object = IVD.integrated.NPC.C5.bp, slot = "data")["MMP3",]>0)
sum(GetAssayData(object = IVD.integrated.NPC.C5.bp, slot = "data")["CXCL8",]>0)
sum(GetAssayData(object = IVD.integrated.NPC.C5.bp, slot = "data")["NGF",]>0)
sum(GetAssayData(object = IVD.integrated.NPC.C5.bp, slot = "data")["FGF2",]>0)


#-------Stress changes marker, which are they------

Markers.IVD.integrated.NPC.bpvsa <- FindMarkers(IVD.integrated.NPC, ident.1 = "bpNPC", ident.2 = "aNPC",  min.pct = 0.1)
Markers.IVD.integrated.NPC.s1vss2 <- FindMarkers(IVD.integrated.NPC, ident.1 = "NPCs1", ident.2 = "NPCs2",  min.pct = 0.1)

IVD.integrated3.bpNPC <- subset(x = IVD.integrated3, idents ="bpNPC" )
DimPlot(IVD.integrated3.bpNPC, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 10) + coord_fixed(ratio=1)
Idents(IVD.integrated3.bpNPC) <- "ClusterNumber"


#-----Export to Qiagen IPA-----
Markers.IVD.integrated.NPC.bp <- FindMarkers(IVD.integrated.NPC, ident.1 = "bpNPC",  min.pct = 0.1)
Markers.IVD.integrated.NPC.a <- FindMarkers(IVD.integrated.NPC, ident.1 = "aNPC",  min.pct = 0.1)
Markers.IVD.integrated.NPC.s1 <- FindMarkers(IVD.integrated.NPC, ident.1 = "NPCs1",  min.pct = 0.1)
Markers.IVD.integrated.NPC.s2 <- FindMarkers(IVD.integrated.NPC, ident.1 = "NPCs2",  min.pct = 0.1)
Markers.IVD.integrated.NPC.C5.bp <- FindMarkers(IVD.integrated.NPC.C5, ident.1 = "bpNPC",  min.pct = 0.1)
Markers.IVD.integrated.NPC.C5.a <- FindMarkers(IVD.integrated.NPC.C5, ident.1 = "aNPC",  min.pct = 0.1)
Markers.IVD.integrated.NPC.C5.s1 <- FindMarkers(IVD.integrated.NPC.C5, ident.1 = "NPCs1",  min.pct = 0.1)
Markers.IVD.integrated.NPC.C5.s2 <- FindMarkers(IVD.integrated.NPC.C5, ident.1 = "NPCs2",  min.pct = 0.1)

Markers.IVD.integrated.NPC.C5 <- FindMarkers(IVD.integrated.NPC, ident.1 = "5",  min.pct = 0.1)

MarkersToExcel.Markers.IVD.integrated.NPC.bp <- cbind("Gene ID"=rownames(Markers.IVD.integrated.NPC.bp), Markers.IVD.integrated.NPC.bp)
write_xlsx(MarkersToExcel.Markers.IVD.integrated.NPC.bp,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/all.bp.xlsx")

MarkersToExcel.Markers.IVD.integrated.NPC.a <- cbind("Gene ID"=rownames(Markers.IVD.integrated.NPC.a), Markers.IVD.integrated.NPC.a)
write_xlsx(MarkersToExcel.Markers.IVD.integrated.NPC.a,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/all.a.xlsx")

MarkersToExcel.Markers.IVD.integrated.NPC.s1 <- cbind("Gene ID"=rownames(Markers.IVD.integrated.NPC.s1), Markers.IVD.integrated.NPC.s1)
write_xlsx(MarkersToExcel.Markers.IVD.integrated.NPC.s1,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/all.s1.xlsx")

MarkersToExcel.Markers.IVD.integrated.NPC.s2 <- cbind("Gene ID"=rownames(Markers.IVD.integrated.NPC.s2), Markers.IVD.integrated.NPC.s2)
write_xlsx(MarkersToExcel.Markers.IVD.integrated.NPC.s2,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/all.s2.xlsx")

MarkersToExcel.Markers.IVD.integrated.NPC.C5.bp <- cbind("Gene ID"=rownames(Markers.IVD.integrated.NPC.C5.bp), Markers.IVD.integrated.NPC.C5.bp)
write_xlsx(MarkersToExcel.Markers.IVD.integrated.NPC.C5.bp,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/C5.bp.xlsx")

MarkersToExcel.Markers.IVD.integrated.NPC.C5.a <- cbind("Gene ID"=rownames(Markers.IVD.integrated.NPC.C5.a), Markers.IVD.integrated.NPC.C5.a)
write_xlsx(MarkersToExcel.Markers.IVD.integrated.NPC.C5.a,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/C5.a.xlsx")

MarkersToExcel.Markers.IVD.integrated.NPC.C5.s1 <- cbind("Gene ID"=rownames(Markers.IVD.integrated.NPC.C5.s1), Markers.IVD.integrated.NPC.C5.s1)
write_xlsx(MarkersToExcel.Markers.IVD.integrated.NPC.C5.s1,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/C5.s1.xlsx")

MarkersToExcel.Markers.IVD.integrated.NPC.C5.s2 <- cbind("Gene ID"=rownames(Markers.IVD.integrated.NPC.C5.s2), Markers.IVD.integrated.NPC.C5.s2)
write_xlsx(MarkersToExcel.Markers.IVD.integrated.NPC.C5.s2,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/C5.s2.xlsx")



#========= compare aNPC, nsNPC, sNPC and pain cluser bpNPC1==========


data2 <- c(IVD.integrated.NPC.s1, IVD.integrated.NPC.C5.bp, IVD.integrated.NPC.s2, IVD.integrated.NPC.a)

for (i in 1:length(data2)) {
  data2[[i]] <- NormalizeData(data2[[i]], verbose = FALSE)
  data2[[i]] <- FindVariableFeatures(data2[[i]], selection.method = "vst", 
                                     nfeatures = 6000, verbose = FALSE)
}

IVD.anchors2 <- FindIntegrationAnchors(object.list = data2, dims = 1:30)
IVD.integrated2 <- IntegrateData(anchorset = IVD.anchors2)
DefaultAssay(IVD.integrated2) <- "integrated"

IVD.integrated2 <- ScaleData(IVD.integrated2, verbose = FALSE)
IVD.integrated2 <- RunPCA(IVD.integrated2, npcs = 30, verbose = FALSE)
IVD.integrated2 <- RunUMAP(IVD.integrated2, reduction = "pca", dims = 1:30)

IVD.integrated2 <- FindNeighbors(IVD.integrated2, reduction = "pca", dims = 1:30)
IVD.integrated2 <- FindClusters(IVD.integrated2, resolution = 0.5)


DefaultAssay(IVD.integrated2) <- "RNA"
IVD.integrated2$ClusterNumber <- Idents(IVD.integrated2)
Idents(IVD.integrated2) <- "orig.ident"
IVD.integrated2 <- RenameIdents(object = IVD.integrated2, "bpNPC1" = "bpNPC_C5", "bpNPC2" = "bpNPC_C5", "bpNPC3" = "bpNPC_C5", "NPCs1" = "NPCs1",  "NPCs2" = "NPCs2", "daNPC1" = "aNPC",  "daNPC2" = "aNPC",  "daNPC3" = "aNPC",  "daNPC4" = "aNPC")
IVD.integrated2$SampleType <- Idents(IVD.integrated2)
table(IVD.integrated2$SampleType)
Idents(IVD.integrated2) <- "ClusterNumber"
Idents(IVD.integrated2) <- "SampleType"

DotPlot(IVD.integrated2, features = c("MMP3", "SPP1", "CRTAC1", "MAP1B", "SOX4","SOX17", "ITGA6", "BASP1", 
                                      "ACAN", "COL2A1", "SOX9",
                                      "COL1A1", "CALR", "HSPA6",
                                      "MEG3", "CD44", "CD14", "HBB", "HBA1"), scale.by = "radius", dot.scale = 6, col.min = -2, col.max = 2) + 
  theme(axis.text.x = element_text(angle = 90)) 

#C5 is the MMP3+ pain cluster

p05a <- DimPlot(IVD.integrated2, reduction = "umap", pt.size = 0.3, label = FALSE, split.by = "SampleType", cols = c("gray80", "gray80", "gray80", "gray80", "gray80", "red", "gray80",  "gray80",  "gray80",  "gray80"),  label.size = 10) + coord_fixed(ratio=1)


p05b <- DimPlot(IVD.integrated2, reduction = "umap", pt.size = 0.3, label = FALSE, split.by = "SampleType",  label.size = 10) + coord_fixed(ratio=1)

p06 <- FeaturePlot(IVD.integrated2, features = c("MMP3",  "CXCL8"),  label = FALSE,  order = TRUE, 
                   cols = brewer.pal(11, "PuRd"), split.by = "SampleType",  pt.size = 0.2, coord.fixed = TRUE)

p07 <- FeaturePlot(IVD.integrated2, features = c("CXCL8", "FGF2"),  blend = TRUE, blend.threshold = 0.1,
            label = FALSE,  order = TRUE, cols = c("black", "blue", "red"), split.by = "SampleType",  pt.size = 1, coord.fixed = TRUE)    + coord_fixed(ratio=1)



IVD.integrated2$seurat_clusters

GoI <- c( "MMP3", 
          "NGF", "BDNF", "CRCP", "CGRP", "TRPV4", "VEGFA",
          "CXCL8", "IL6", "IL1B", 
          "TNF",  "TNFA", "TNFAIP1", "TNFSF4", 
          "CCN2", "TM4SF1", "SLC7A2",
          "TRPC6",
          "NOS2","PRDX1", "HIF1A", "BNIP3", "PRG4",
          "ASIC3",  "SEMA3A",  "SHH", "ADAMTS5", "FGF2" ,
          "CD9", "CD109",
          "KRT9", "KRT19", "MIR155HG")

AverageExpression(
  IVD.integrated2,
  assays = "RNA",
  features = GoI,
  return.seurat = FALSE,
  group.by = "SampleType",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)




Markers.IVD.integrated2.bpC5 <- FindMarkers(IVD.integrated2, ident.1 = "bpNPC_C5",  min.pct = 0.1)
Markers.IVD.integrated2.s1 <- FindMarkers(IVD.integrated2, ident.1 = "NPCs1",  min.pct = 0.1)
Markers.IVD.integrated2.s2 <- FindMarkers(IVD.integrated2, ident.1 = "NPCs2",  min.pct = 0.1)
Markers.IVD.integrated2.a <- FindMarkers(IVD.integrated2, ident.1 = "aNPC",  min.pct = 0.1)


MarkersToExcel.Markers.IVD.integrated2.bpC5 <- cbind("Gene ID"=rownames(Markers.IVD.integrated2.bpC5), Markers.IVD.integrated2.bpC5)
write_xlsx(MarkersToExcel.Markers.IVD.integrated2.bpC5,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/A4.bpC5.xlsx")

MarkersToExcel.Markers.IVD.integrated2.s1 <- cbind("Gene ID"=rownames(Markers.IVD.integrated2.s1), Markers.IVD.integrated2.s1)
write_xlsx(MarkersToExcel.Markers.IVD.integrated2.s1,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/A4.s1.xlsx")

MarkersToExcel.Markers.IVD.integrated2.s2 <- cbind("Gene ID"=rownames(Markers.IVD.integrated2.S2), Markers.IVD.integrated2.s2)
write_xlsx(MarkersToExcel.Markers.IVD.integrated2.s2,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/A4.s2.xlsx")

MarkersToExcel.Markers.IVD.integrated2.a <- cbind("Gene ID"=rownames(Markers.IVD.integrated2.a), Markers.IVD.integrated2.a)
write_xlsx(MarkersToExcel.Markers.IVD.integrated2.a,  path = "/Users/JiangW2/Desktop/Wensen/Single cell_back pain_trajectory/Outputs 4_stress compare/IPA pathway/A4.a.xlsx")

