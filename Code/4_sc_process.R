#1. scRNA-seq data(Periodontitis) -------------------------------------------

library(harmony)
library(ggplot2)
library(dplyr)
library(hdf5r)
library(SeuratDisk)
library(reticulate)
library(GEOquery)

rm(list = ls())
gc()
path <- getwd()
setwd(path)

input <- "sc_data/sc_P1"
output <- "10_sc_Preprocessed/sc_P1"

if (!dir.exists(output)) {
  dir.create(output,recursive = TRUE)
}


# GSE152042/sc_P1 
# Load the PBMC dataset
gsm_896_h1 <- Read10X_h5(paste0(input,"/GSM4600896_filtered_feature_bc_matrix_h1.h5"))
gsm_897_h2 <- Read10X_h5(paste0(input,"/GSM4600897_filtered_feature_bc_matrix_h2.h5"))
gsm_898_pm <- Read10X_h5(paste0(input,"/GSM4600898_filtered_feature_bc_matrix_p1.h5"))
gsm_899_ps <- Read10X_h5(paste0(input,"/GSM4600899_filtered_feature_bc_p2.h5"))

ss_h1 <- CreateSeuratObject(counts = gsm_896_h1)
ss_h2 <- CreateSeuratObject(counts = gsm_897_h2)
ss_pm <- CreateSeuratObject(counts = gsm_898_pm)
ss_ps <- CreateSeuratObject(counts = gsm_899_ps)

ss_h1$group <- "Healthy"
ss_h2$group <- "Healthy"
ss_pm$group <- "Periodontitis"
ss_ps$group <- "Periodontitis"

ss_h1$label <- "Healthy"
ss_h2$label <- "Healthy"
ss_pm$label <- "Mild Periodontitis"
ss_ps$label <- "Severe Periodontitis"

ss_h1$batch <- "batch_h1"
ss_h2$batch <- "batch_h2"
ss_pm$batch <- "batch_pm"
ss_ps$batch <- "batch_ps"



ss <- merge(ss_h1, 
            y = c(ss_h2, ss_pm, ss_ps), 
            add.cell.ids = c("h1", "h2", "pm", "ps"), 
            project = "Combined_Project")



table(ss$group)
table(ss$batch)
table(ss$label)


##QC
ss[["percent.mt"]] <- PercentageFeatureSet(ss, pattern = "^MT-")
ss[["percent.HB"]] <- PercentageFeatureSet(ss, pattern = "^HB")
ss[["percent.Ribo"]] <- PercentageFeatureSet(ss, pattern = "^RPS|^RPL")

VlnPlot(ss, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.HB", "percent.Ribo"), 
        group.by = "batch",
        ncol = 3);
ggsave(filename = paste0(output,"/Vln_before_QC.png"), 
       width = 16, 
       height = 16)



FeatureScatter(ss, "nCount_RNA", "nFeature_RNA", group.by = "batch", pt.size = 0.5)
FeatureScatter(ss, "nCount_RNA", "percent.mt", group.by = "batch", pt.size = 0.5)

table(ss$nFeature_RNA > 200 & 
        ss$nFeature_RNA < 5000 &
        ss$nCount_RNA < 50000 & 
        ss$percent.mt < 15 & 
        ss$percent.HB < 5 & 
        ss$percent.Ribo < 40)


ss <- subset(ss, 
             subset = nFeature_RNA > 200 & 
               nFeature_RNA < 5000 & 
               nCount_RNA < 50000 & 
               percent.mt < 15 & 
               percent.HB < 5 & 
               percent.Ribo < 40)

length(Cells(ss))



p_qc <- VlnPlot(ss, 
                features = c("nFeature_RNA", "nCount_RNA",
                             "percent.mt", "percent.HB", "percent.Ribo"),  
                split.plot = FALSE, 
                group.by = "batch",
                ncol = 3)


ggsave(filename = paste0(output,"/Vln_after_QC.png"), 
       plot = p_qc, 
       width = 16, 
       height = 16)


##Normalize
ss <- NormalizeData(ss, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
ss <- FindVariableFeatures(ss, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ss <- ScaleData(ss, vars.to.regress = c("percent.Ribo", "percent.mt"), verbose = FALSE)

# PCA
ss <- RunPCA(ss, dims = 1:20)
ss <- RunUMAP(ss, dims = 1:20)

p1 <- DimPlot(ss, reduction = "umap", group.by = "batch")
ss_harmony <- RunHarmony(ss, group.by.vars = "batch", dims.use = 1:30)


ss_harmony <- FindNeighbors(ss_harmony,reduction = "harmony",dims = 1:20)
ss_harmony <- FindClusters(ss_harmony,resolution = 0.5)

ss_harmony <- RunUMAP(ss_harmony,reduction = "harmony",dims = 1:20,reduction.name = "umap")
DimPlot(ss_harmony,reduction = "umap",label = T)

p2 <- DimPlot(ss_harmony,reduction = "umap",group.by = "batch",label = F)

p = p1|p2
ggsave(paste0(output,'/umap_harmony.png'), 
       plot = p, device = 'png', height = 5, width = 10)


##vital umap
vital_genes <- c("CXCL5","FOSB","PTGR1")
p_feature <- FeaturePlot(ss_harmony, 
                         features = vital_genes,
                         cols = c("lightgrey", "#BC3C29FF"),
                         ncol = 3)
p_feature
FeaturePlot(ss_harmony, 
            features = c("CXCL5"),
            split.by = "group",
            cols = c("lightgrey", "#BC3C29FF"))
ggsave(filename = paste0(output,"/","UMAP_vital_genes_all.png"), 
       plot = p_feature, 
       width = 9, 
       height = 3)


##marker
markers_paper <- list(
  "T/NK" = c("CD3D","CD3G","NKG7"),#"CXCR6"
  "B" = c("CD79A","CD19"),
  "Follicular B" = c( "MZB1", "DERL3", "IGHG4"),
  "Memory B" = c("SELL","CD27","MS4A1", "CD37"),
  "Cycling B" = c("MKI67", "TOP2A", "HMGB2", "TUBA1B", "UBE2C"), 
  "IgA Plasma B" = c("IGHA1", "IGHA2"),
  "IgG Plasma B" = c("IGHG1", "IGHG2", "IGHG3", "IGHG4"),
  "Endothelial" = c("VWF","PECAM1"),
  "Epithelial" = c("HOPX", "IGFBP5", "LAMB3"),
  "DC" = c("CLEC9A", "IRF8"),
  "Macrophages" = c("LYZ", "AIF1"),
  "Mast" = c("TPSB2", "CPA3"),
  "Fibroblast" = c("COL1A1","COL1A2","LUM","CXCL14","ADAMDEC1"),
  "Myofibroblast" = c("RGS5","ACTA2")
)


marker_genes <- unlist(markers_paper)

p_marker <- DotPlot(ss_harmony,
                    features =unique(marker_genes) ,
                    cluster.idents = T)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_marker

ggsave(filename = paste0(output,"/","DOT_marker.png"), 
       plot = p_marker, 
       width = 12, 
       height = 10)


##annotation
ss_anno <- ss_harmony

ss_anno@meta.data <- ss_anno@meta.data %>%
  mutate(cell_type = case_when(
    seurat_clusters %in% c(0, 2, 12, 14) ~ "IgG Plasma B",   
    seurat_clusters %in% c(3, 11) ~ "Macrophages", 
    seurat_clusters %in% c(9, 17) ~ "Endothelial", 
    seurat_clusters %in% c(6, 10) ~ "Epithelial", 
    seurat_clusters %in% c(15,16)  ~ "DC",
    seurat_clusters == 5  ~ "Fibroblast",
    seurat_clusters == 8  ~ "Myofibroblast",
    seurat_clusters == 7  ~ "Memory B",
    seurat_clusters == 13 ~ "Cycling B",
    seurat_clusters == 1 ~ "Mast",
    seurat_clusters == 4  ~ "T/NK"))

table(ss_anno@meta.data$cell_type)
DimPlot(ss_anno, group.by = "cell_type", label = TRUE)

##umap
p_manual <- DimPlot(ss_anno,reduction = "umap",group.by = "cell_type",label = T,label.size = 6)
ggsave(filename = paste0(output,"/","UMAP_anno.png"), 
       plot = p_manual, 
       width = 10, 
       height = 10)


saveRDS(ss_anno, file = paste0(output, "/ss_anno.rds"))


# B Cells ---------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)

input12 <- "12_sc_analysis/P"
output <- "13_sc_cluster2/B/P"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss <- readRDS(file = paste0(input12, "/ss_anno.rds"))


table(ss$label)
table(ss$cell_type)
table(ss$anno)



Idents(ss) <- "anno"
ss_B <- subset(ss, idents = "B")
table(ss_B$seurat_clusters)
table(ss_B$label)
table(ss_B$group)


ss_B <- RunPCA(ss_B)
ss_B <- FindNeighbors(ss_B,reduction = "harmony",dims = 1:15)
ss_B <- FindClusters(ss_B, resolution = 0.3)
ss_B <- RunUMAP(ss_B,reduction = "harmony",dims = 1:15)


DimPlot(ss_B, reduction = "umap", group.by = "batch")
DimPlot(ss_B, reduction = "umap",label = T)



B_cell_markers <- list(
  "Follicular B" = c( "IGHD","TCL1A", "FCER2"),
  "Cycling B Cells" = c("MKI67"), 
  "Memory B Cells" = c("CD27"),
  "IgA Memory B" = c("IGHA1", "IGHA2"),
  "IgG Memory B" = c("IGHG1", "IGHG2", "IGHG3", "IGHG4")
)
all_marker <- unique(unlist(B_cell_markers))
DotPlot(ss_B,
        features = all_marker,
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# annotation
ss_anno <- ss_B
ss_anno@meta.data <- ss_anno@meta.data %>%
  mutate(cell_type = case_when(
    seurat_clusters %in% c(1) ~ "Memory B",   
    seurat_clusters %in% c(0) ~ "Follicular B",
    seurat_clusters %in% c(2) ~ "Cycling B"
  ))


table(ss_anno@meta.data$cell_type)
DimPlot(ss_anno, group.by = "cell_type", label = TRUE)

##umap
p_manual <- DimPlot(ss_anno,reduction = "umap",group.by = "cell_type",label = T,label.size = 6)
ggsave(filename = paste0(output,"/","UMAP_anno_B.png"), 
       plot = p_manual, 
       width = 10, 
       height = 10)


saveRDS(ss_anno, file = paste0(output, "/ss_anno_B.rds"))


# Fibroblast --------------------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)


input12 <- "12_sc_analysis/P"
output <- "13_sc_cluster/fibro"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss <- readRDS(file = paste0(input12, "/ss_anno.rds"))
table(ss$label)
table(ss$cell_type)


Idents(ss) <- "cell_type"
ss_fibro <- subset(ss, idents = "Fibroblast")
table(ss_fibro$seurat_clusters)
table(ss_fibro$label)
table(ss_fibro$group)
saveRDS(ss_fibro, file = paste0(output, "/ss_P_fibro.rds"))


ss_fibro <- NormalizeData(ss_fibro)
ss_fibro <- FindVariableFeatures(ss_fibro)
ss_fibro <- ScaleData(ss_fibro)
ss_fibro <- RunPCA(ss_fibro)
ss_fibro <- RunUMAP(ss_fibro, dims = 1:10)


DimPlot(ss_fibro, reduction = "umap", group.by = "label")

actived_markers <- c("IL11","IL24","OSM","OSMR","WNT2B","WNT5B",
                     "CXCL1","CXCL5","CXCL8","CCL2","CCL7","TNFSF14",
                     "CLU","IL33","IL6","LOX","GREM1","MMP3","IL1B","WNT5A","IL6")
DotPlot(ss_fibro, 
        features =  unique(actived_markers),
        group.by = "group") + 
  coord_flip()


saveRDS(ss_fibro, file = paste0(output, "/ss_P_fibro2.rds"))



# Epithelial --------------------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)


input12 <- "12_sc_analysis/P"
output <- "13_sc_cluster/Epi/P"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss <- readRDS(file = paste0(input12, "/ss_anno.rds"))
table(Idents(ss))
Idents(ss) <- "cell_type"
Epi <- subset(ss, idents = "Epithelial")

ss_Epi <- Epi


DimPlot(ss_Epi, reduction = "umap", label = TRUE)
DimPlot(ss_Epi, reduction = "umap", split.by  = "group",label = TRUE) 

saveRDS(ss_Epi, file = paste0(output, "/ss_P_Epi.rds"))








#2.scRNA-seq data (UC)------------------------------------------------------
library(Seurat)
library(ggplot2)
library(devtools)
library(harmony)
library(dplyr)
library(matrixStats)


rm(list = ls())
gc()
path <- getwd()
setwd(path)


input <- "10_sc_Preprocessed"
output <- "11_sc_merged"

if (!dir.exists(output)) {
  dir.create(output,recursive = TRUE)
}

rds_files <- list.files(input, pattern = "*.rds", full.names = TRUE)
seurat_list <- lapply(rds_files, readRDS)



# Merging Seurat Objects
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1], 
                         add.cell.ids = gsub(".rds", "", basename(rds_files)))

## QC
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")
seurat_combined[["percent.HB"]] <- PercentageFeatureSet(seurat_combined, pattern = "^HB")
seurat_combined[["percent.Ribo"]] <- PercentageFeatureSet(seurat_combined, pattern = "^RPS|^RPL")

VlnPlot(seurat_combined, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.HB", "percent.Ribo"), 
        group.by = "batch",
        ncol = 3);
ggsave(filename = paste0(output,"/Vln_before_QC.pdf"), 
       width = 16, 
       height = 16)



FeatureScatter(seurat_combined, "nCount_RNA", "nFeature_RNA", group.by = "batch", pt.size = 0.5)
FeatureScatter(seurat_combined, "nCount_RNA", "percent.mt", group.by = "batch", pt.size = 0.5)
seurat_qc <- subset(seurat_combined, 
                    subset = nFeature_RNA > 200 & 
                      nFeature_RNA < 4000 & 
                      nCount_RNA < 50000 & 
                      percent.mt < 30 & 
                      percent.HB < 5 & 
                      percent.Ribo < 60)

length(Cells(seurat_qc))

table(seurat_combined$nFeature_RNA > 200 & 
        seurat_combined$nFeature_RNA < 4000 &
        seurat_combined$nCount_RNA < 50000 & 
        seurat_combined$percent.mt < 30 & 
        seurat_combined$percent.HB < 5 & 
        seurat_combined$percent.Ribo < 60)


gene_counts_per_cell <- colSums(seurat_combined@assays$SCT@counts > 0)


average_gene_count <- mean(gene_counts_per_cell)
average_gene_count


p_qc <- VlnPlot(seurat_qc, 
                features = c("nFeature_RNA", "nCount_RNA",
                             "percent.mt", "percent.HB", "percent.Ribo"),  
                split.plot = FALSE, 
                group.by = "batch",
                ncol = 3)


ggsave(filename = paste0(output,"/Vln_after_QC.pdf"), 
       plot = p_qc, 
       width = 16, 
       height = 16)

#saveRDS(seurat_qc,file = paste0(output,'/ss_qc.Rds'))



# Normalization 

ss <- SCTransform(seurat_qc,vars.to.regress = "percent.mt",verbose = FALSE)

#PCA
ss <- RunPCA(object = ss,features = VariableFeatures(object = ss))
DimPlot(ss, reduction = "pca",group.by = "batch")

ElbowPlot(ss)
ss <- FindNeighbors(ss,reduction = "pca",dims = 1:20)
ss <- FindClusters(ss,resolution = 0.5)
ss <- RunTSNE(ss,dims = 1:20)
ss <- RunUMAP(ss,dims = 1:20)

colnames(ss@meta.data)
umap_p_1 <- DimPlot(ss,reduction = "umap",label = T,group.by = 'batch')
umap_p_2 <- DimPlot(ss,reduction = "umap",label = T,split.by = 'batch',ncol = 2)

ggsave(paste0(output,'/umap_before_harmony.png'), 
       plot = umap_p_1, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/umap_before_harmony_split.png'), 
       plot = umap_p_2, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')



# harmony /remove batch effect
ss_harmony <- RunHarmony(ss, group.by.vars = "batch",assay.use = "SCT",max_iter = 10)
ss_harmony <- RunUMAP(ss_harmony, reduction = "harmony", dims = 1:20,reduction.name = "umap")


# UMAP
ss_harmony <- FindNeighbors(ss_harmony,reduction = "harmony",dims = 1:20)
ss_harmony <- FindClusters(ss_harmony,resolution = 0.5)

ss_harmony <- RunUMAP(ss_harmony,reduction = "harmony",dims = 1:20,reduction.name = "umap")

umap_p_3 <- DimPlot(ss_harmony,reduction = "umap",group.by = "batch",label = T)
umap_p_4 <- DimPlot(ss_harmony,reduction = "umap",split.by = "batch",label = T,ncol = 2)


ggsave(paste0(output,'/umap_after_harmony.png'), 
       plot = umap_p_3, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/umap_after_harmony_split.png'), 
       plot = umap_p_4, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')


# marker genes
ss_harmony <- PrepSCTFindMarkers(ss_harmony)

all_marker <- FindAllMarkers(
  ss_harmony,
  assay = 'SCT',    
  slot = 'data',   
  only.pos = TRUE   
)
write.csv(all_marker,file = paste0(output,'/all_marker.csv'),fileEncoding = 'UTF-8',row.names = F)


filtered_markers <- all_marker

top10 <- filtered_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
write.csv(top10,file = paste0(output,'/cluster_top10_marker_0.5.csv'),fileEncoding = 'UTF-8',row.names = F)

saveRDS(ss_harmony,file = paste0(output,"/ss_after_harmony.Rds"))



table(Idents(ss_harmony))
markers_paper <- list(
  "T/NK" = c("CD3D","CD3G","NKG7"),
  "Plasma" = c("IGHA1","IGHA2"),
  "B" = c("CD79A","CD19"),
  "Endothelial" = c("VWF","PECAM1"),
  "Fibroblasts" = c("COL1A1","COL1A2","LUM"),
  "Gial" = c("S100B"),
  "Melyoid" = c("LYZ"), 
  "Epithelial" = c("EPCAM")
)
marker_genes <- unlist(markers_paper)

p_marker <- DotPlot(ss_harmony,
                    features =unique(marker_genes),
                    cluster.idents = T)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_marker
ggsave(paste0(output,'/dot_all_cluster_marker.png'), 
       plot = p_marker, device = 'png', height = 25, width = 40, limitsize = FALSE, units = 'cm')


##umap
p <- DimPlot(ss_harmony,reduction = "umap",label = T,label.size = 3)
ggsave(filename = paste0(output,"/","umap_all_cluster.pdf"), 
       plot = p, 
       width = 24, 
       height = 12)

##vital umap
vital_genes <- c("CXCL5","FOSB","PTGR1")
p_feature <- FeaturePlot(ss_harmony, 
                         features = vital_genes,
                         split.by  = "group")
p_feature
ggsave(filename = paste0(output,"/","UMAP_vital_genes_all.pdf"), 
       plot = p_feature, 
       width = 16, 
       height = 16)




# Removal of Low-Quality Clusters
ss_harmony <- readRDS(paste0(output,'/ss_after_harmony.Rds'))
filtered_markers <- read.csv(paste0(output,"/all_marker.csv"))
top10 <- read.csv(paste0(output,"/cluster_top10_marker_0.5.csv"))
table(ss_harmony$seurat_clusters)


markers_paper <- list(
  "T/NK" = c("CD3D","CD3G","NKG7"),
  "Plasma" = c("IGHA1","IGHA2"),
  "B" = c("CD79A","CD19"),
  "Endothelial" = c("VWF","PECAM1"),
  "Fibroblasts" = c("COL1A1","COL1A2","LUM"),
  "Gial" = c("S100B"),
  "Melyoid" = c("LYZ"), 
  "Epithelial" = c("EPCAM")
)
marker_genes <- unlist(markers_paper)

DotPlot(ss_harmony,
        features =unique(marker_genes) ,
        cluster.idents = T)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


cluster_n <- filtered_markers %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(cluster == 17) %>%
  arrange(desc(pct.1)) 

write.csv(cluster_n, 
          file = paste0(output, "/cluster_17.csv"), 
          row.names = FALSE)


cells_to_keep <- WhichCells(ss_harmony, 
                            idents = setdiff(unique(Idents(ss_harmony)), c("17","32","37")))


ss_filtered <- subset(ss_harmony, cells = cells_to_keep)


p1 <- DimPlot(ss_filtered,reduction = "umap",group.by = "seurat_clusters", label = TRUE)

ggsave(paste0(output,'/umap_after_filter.png'), 
       plot = p1, device = 'png', height = 15, width = 20, limitsize = FALSE,units = 'cm')



# Re-clustering
ss_filtered <- FindNeighbors(ss_filtered, reduction = "harmony",dims = 1:20)
ss_filtered <- FindClusters(ss_filtered, resolution = 0.3)
ss_filtered <- RunUMAP(ss_filtered, reduction = "harmony",dims = 1:20)


DimPlot(ss_filtered, reduction = "umap", label = TRUE)

umap_p_5 <- DimPlot(ss_filtered,reduction = "umap",group.by = "batch",label = T)
umap_p_6 <- DimPlot(ss_filtered,reduction = "umap",split.by = "batch",label = T,ncol = 2)
umap_p_7 <- DimPlot(ss_filtered,reduction = "umap",label = T)


ggsave(paste0(output,'/Umap_filtered.png'), 
       plot = umap_p_5, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/Umap_filtered_split.png'), 
       plot = umap_p_6, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/Umap_filtered_cluster.png'), 
       plot = umap_p_7, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')

markers_paper <- list(
  "T/NK" = c("CD3D","CD3G","NKG7"),
  "Plasma" = c("IGHA1","IGHA2"),
  "B" = c("CD79A","CD19"),
  "Endothelial" = c("VWF","PECAM1"),
  "Fibroblasts" = c("COL1A1","COL1A2","LUM"),
  "Gial" = c("S100B"),
  "Melyoid" = c("LYZ"), 
  "Epithelial" = c("EPCAM")
)
marker_genes <- unlist(markers_paper)



p_filtered_marker <- DotPlot(ss_filtered,
                             features =unique(marker_genes) ,
                             cluster.idents = T)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_filtered_marker

ggsave(paste0(output,'/dot_filtered_marker.png'), 
       plot = p_filtered_marker, height = 10, width = 15)



ss_filtered <- PrepSCTFindMarkers(ss_filtered)
all_marker <- FindAllMarkers(
  ss_filtered,
  assay = 'SCT',   
  slot = 'data',   
  only.pos = TRUE   
)
write.csv(all_marker,file = paste0(output,'/new_all_marker.csv'),fileEncoding = 'UTF-8',row.names = F)
filtered_markers <- all_marker

top10 <- filtered_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
write.csv(top10,file = paste0(output,'/new_cluster_top10_marker.csv'),fileEncoding = 'UTF-8',row.names = F)


cluster_n <- filtered_markers %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(cluster == 11) %>%
  arrange(desc(pct.1))  

write.csv(cluster_n, 
          file = paste0(output, "/new2_cluster_11.csv"), 
          row.names = FALSE)


saveRDS(ss_filtered,file = paste0(output,'/ss_filtered.Rds'))




# Second Round of Low-Quality Cluster Removal
rm(list = ls())
gc()
path <- getwd()
setwd(path)


output <- "11_sc_merged"

ss_filtered <- readRDS(paste0(output,'/ss_filtered.Rds'))



cells_to_keep2 <- WhichCells(ss_filtered, 
                             idents = setdiff(unique(Idents(ss_filtered)), c("11")))


ss_filtered2 <- subset(ss_filtered, cells = cells_to_keep2)

gc()

p <- DimPlot(ss_filtered2,reduction = "umap",group.by = "seurat_clusters", label = TRUE)
p
ggsave(paste0(output,'/umap_after_filter2.png'), 
       plot = p, device = 'png', height = 15, width = 20, limitsize = FALSE,units = 'cm')


ss_filtered2 <- FindNeighbors(ss_filtered2, reduction = "harmony",dims = 1:20)
ss_filtered2 <- FindClusters(ss_filtered2, resolution = 0.3)
ss_filtered2 <- RunUMAP(ss_filtered2, reduction = "harmony",dims = 1:20)


DimPlot(ss_filtered2, reduction = "umap", label = TRUE)

umap_p_5 <- DimPlot(ss_filtered2,reduction = "umap",group.by = "batch",label = T)
umap_p_6 <- DimPlot(ss_filtered2,reduction = "umap",split.by = "batch",label = T,ncol = 2)
umap_p_7 <- DimPlot(ss_filtered2,reduction = "umap",label = T)


ggsave(paste0(output,'/Umap_filtered2.png'), 
       plot = umap_p_5, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/Umap_filtered2_split.png'), 
       plot = umap_p_6, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/Umap_filtered2_cluster.png'), 
       plot = umap_p_7, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')

markers_paper <- list(
  "T/NK" = c("CD3D","CD3G","NKG7"),
  "Plasma" = c("IGHA1","IGHA2"),
  "B" = c("CD79A","CD19"),
  "Endothelial" = c("VWF","PECAM1"),
  "Fibroblasts" = c("COL1A1","COL1A2","LUM"),
  "Gial" = c("S100B"),
  "Melyoid" = c("LYZ"), 
  "Epithelial" = c("EPCAM")
)
marker_genes <- unlist(markers_paper)



p_filtered_marker <- DotPlot(ss_filtered2,
                             features =unique(marker_genes) ,
                             cluster.idents = T)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_filtered_marker

ggsave(paste0(output,'/dot_filtered2_marker.png'), 
       plot = p_filtered_marker, height = 10, width = 15)

table(Idents(ss_filtered2))  



# cell annotation
new.cluster.ids <- c("Epithelial",
                     #1-5
                     "Fibroblasts","B","T/NK","T/NK","Plasma",
                     #6-10
                     "Endothelial","Plasma","Plasma","T/NK","Myeloid",
                     #11-15
                     "Plasma","T/NK","Myofibroblast","Epithelial","Plasma",
                     #16-20
                     "T/NK","Fibroblasts","Epithelial","Gial","T/NK",
                     #21-25
                     "Plasma","Plasma","Plasma","Plasma","Plasma",
                     #26-29
                     "T/NK","T/NK","Plasma","Plasma")


ss_anno <- ss_filtered2
names(new.cluster.ids) <- levels(ss_anno)
ss_anno <- RenameIdents(ss_anno, new.cluster.ids)


ss_anno$anno <- Idents(ss_anno)


##umap
p_manual <- DimPlot(ss_anno,reduction = "umap",label = T,label.size = 6)
ggsave(filename = paste0(output,"/","UMAP_anno_manual2.pdf"), 
       plot = p_manual, 
       width = 10, 
       height = 10)


table(ss_anno$anno)



##DOT
dot_manual <- DotPlot(ss_anno,
                      features =unique(marker_genes),
                      
                      cluster.idents = T)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste0(output,"/","DOT_anno_manual2.pdf"), 
       plot = dot_manual, 
       width = 12, 
       height = 10)


ss <- ss_anno
table(ss$anno)
table(ss$label)
table(ss$histo3)

ss$label <- factor(ss$label, levels = c("healthy", "non inflamed UC", "inflamed UC"))
ss$histo3 <- factor(ss$histo3, levels = c("Healthy (N=14)", "Non inflamed UC (N=13)", "Inflamed UC (N=23)"))
saveRDS(ss, file = paste0(output, "/ss_final.rds"))


# B Cells -----------------------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)


input11 <- "11_sc_merged"
output <- "13_sc_cluster2/B/UC"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss <- readRDS(file = paste0(input11, "/ss_final.rds"))
table(Idents(ss))
table(ss$histo3)

Idents(ss) <- "anno"
B <- subset(ss, idents = "B")


DimPlot(B, reduction = "umap",group.by = "batch")


# PCA
ss_B <- RunPCA(B, features = VariableFeatures(B))
# UMAP
ss_B <- RunUMAP(ss_B, dims = 1:20)
ss_B <- FindNeighbors(ss_B, dims = 1:15)  
ss_B <- FindClusters(ss_B, resolution = 0.3)  

DimPlot(ss_B, reduction = "umap", group.by = "batch")
DimPlot(ss_B, reduction = "umap",label = T)


define_mark <- c("IGHA1","IGHA2","CD79A","CD19")
DotPlot(ss_B,
        features = define_mark,
        cluster.idents = TRUE,
        scale = T) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


b_cell_markers <- list(
  #"all" = c("CD79A","CD19"),
  #"Pre-B Cells" = c("CD38", "CD19"), 
  "Follicular B" = c("TCL1A","FCER2","IGHD"),
  "Naive B Cells" = c("IGHM", "IGHD","FCMR"),
  "Activated B Cells" = c("CD38","CD69","CD86","CD40"),
  "Memory B Cells" = c("SELL","CD27","MS4A1", "CD37"),
  "Brm" = c("CD69","CD44"),#tissue-resident memory B
  "Cycling B Cells" = c("MKI67", "HMGB2", "TUBA1B", "UBE2C"), 
  
  "Plasma" = c("ITGAE"),
  "IgA Producing B Cells" = c("IGHA1", "IGHA2"),
  "IgG Producing B Cells" = c("IGHG1", "IGHG2", "IGHG3", "IGHG4")
  # "IgE Producing B Cells" = c("IGHE")
)
all_marker <- unique(unlist(b_cell_markers))



DotPlot(ss_B,
        features = all_marker,
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Removal of Low-Quality Clusters
cells_to_keep <- WhichCells(ss_B, 
                            idents = setdiff(unique(Idents(ss_B)), c("5","9","10")))


ss_filtered <- subset(ss_B, cells = cells_to_keep)
ncol(ss_B)
ncol(ss_filtered) 

DimPlot(ss_filtered,reduction = "umap",group.by = "seurat_clusters", label = TRUE)


# Re-clustering
ss_filtered <- FindNeighbors(ss_filtered,reduction = "harmony",dims = 1:15)
ss_filtered <- FindClusters(ss_filtered, resolution = 0.3)
ss_filtered <- RunUMAP(ss_filtered,reduction = "harmony",dims = 1:15)


table(Idents(ss_filtered))


DimPlot(ss_filtered, reduction = "umap", label = TRUE)
DimPlot(ss_filtered, reduction = "umap",group.by = "batch")

DotPlot(ss_filtered,
        features = define_mark,
        scale = T,
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1 <- DimPlot(ss_filtered, reduction = "umap", label = TRUE) + 
  ggtitle("B")
p1
ggsave(filename = paste0(output,"/","UMAP_filtered_B.png"), 
       plot = p1, 
       width = 12, 
       height =12,
       limitsize = FALSE, 
       units = 'cm')

p2 <- DimPlot(ss_filtered, reduction = "umap", split.by  = "group",label = TRUE) + 
  ggtitle("B")
p2
ggsave(filename = paste0(output,"/UMAP_filtered_B_split.png"), 
       plot = p2, 
       width = 16, 
       height = 12,
       limitsize = FALSE, 
       units = 'cm')





b_cell_markers <- list(
  "all" = c("CD79A","CD79B"),
  "Follicular B" = c( "IGHD","TCL1A", "FCER2"),
  "Cycling B Cells" = c("MKI67"),
  "Memory B Cells" = c("CD27"),
  "IgA Memory B" = c("IGHA1", "IGHA2"),
  "IgG Memory B" = c("IGHG1", "IGHG2", "IGHG3", "IGHG4")
)
all_marker <- unique(unlist(b_cell_markers))

DotPlot(ss_filtered,
        features = all_marker,
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


p4 <- DotPlot(ss_filtered,
              features = all_marker,
              cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4
ggsave(paste0(output, '/Dot_filtered_B.png'), 
       plot = p4,
       width = 40,
       height = 20,  
       limitsize = FALSE, 
       units = 'cm')


# annotation
ss_anno <- ss_filtered
ss_anno@meta.data <- ss_anno@meta.data %>%
  mutate(cell_type = case_when(
    seurat_clusters %in% c(7,8) ~ "Cycling B",
    seurat_clusters %in% c(0,2,3,4,5,6) ~ "Memory B", 
    seurat_clusters %in% c(1) ~ "Follicular B"
  ))


table(ss_anno@meta.data$cell_type)
DimPlot(ss_anno, group.by = "cell_type", label = TRUE)

##umap
p_manual <- DimPlot(ss_anno,reduction = "umap",group.by = "cell_type",label = T,label.size = 6)
ggsave(filename = paste0(output,"/","UMAP_anno_B.png"), 
       plot = p_manual, 
       width = 10, 
       height = 10)

saveRDS(ss_anno, file = paste0(output, "/ss_anno_B.rds"))


# Fibroblast --------------------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)


input11 <- "11_sc_merged"
output <- "13_sc_cluster/fibro"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss <- readRDS(file = paste0(input11, "/ss_final.rds"))


table(ss$label)
table(ss$anno)


Idents(ss) <- "anno"
ss_fibro <- subset(ss, idents = "Fibroblast")
table(ss_fibro$seurat_clusters)
table(ss_fibro$label)


saveRDS(ss_fibro, file = paste0(output, "/ss_UC_fibro.rds"))



# Epithelial --------------------------------------------------------------
rm(list = ls())
gc()
path <- getwd()
setwd(path)


input11 <- "11_sc_merged"
output <- "13_sc_cluster/Epi/UC"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

ss <- readRDS(file = paste0(input11, "/ss_final.rds"))
table(Idents(ss))
table(ss$histo3)
 
Idents(ss) <- "anno"
Epi <- subset(ss, idents = "Epithelial")



Epi_counts <- GetAssayData(Epi, assay = "RNA", layer = "counts")



Epi_meta <- Epi@meta.data  
colnames(Epi_meta)


Epi_meta <- Epi_meta[, c("orig.ident","nCount_RNA","nFeature_RNA",
                         "project","group","Disease","Site","sample_id",
                         "Inflammation","label","batch","histo3")]  

rownames(Epi_meta) <- colnames(Epi_counts)


ss_Epi <- CreateSeuratObject(counts = Epi_counts, meta.data = Epi_meta)


# QC
ss_Epi[["percent.mt"]] <- PercentageFeatureSet(ss_Epi, pattern = "^MT-")
ss_Epi[["percent.HB"]] <- PercentageFeatureSet(ss_Epi, pattern = "^HB")
ss_Epi[["percent.Ribo"]] <- PercentageFeatureSet(ss_Epi, pattern = "^RPS|^RPL")
VlnPlot(ss_Epi, 
        features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.HEpi", "percent.Ribo"), 
        group.by = "batch",
        ncol = 3)


FeatureScatter(ss_Epi, "nCount_RNA", "percent.mt", group.by = "batch", pt.size = 0.5)

ss_Epi <- subset(ss_Epi,
                 subset = nFeature_RNA > 200 &
                   nFeature_RNA < 4000 &
                   nCount_RNA < 50000 &
                   percent.mt < 20 &
                   percent.HB < 5 &
                   percent.Ribo < 40)
length(Cells(ss_Epi))
p_ss_Epi <- VlnPlot(ss_Epi,
                    features = c("nFeature_RNA", "nCount_RNA",
                                 "percent.mt", "percent.HB", "percent.Ribo"),
                    split.plot = FALSE,
                    group.by = "batch",
                    ncol = 3)
ggsave(filename = paste0(output,"/Vln_after_QC.png"),
       plot = p_ss_Epi,
       width = 16,
       height = 16)



# Ribosomal and Mitochondrial Gene Removal
ribosomal_genes <- grep("^RPS|^RPL", rownames(ss_Epi), value = TRUE)
mitochondrial_genes <- grep("^MT-", rownames(ss_Epi), value = TRUE)

excluded_genes <- union(ribosomal_genes, mitochondrial_genes)
selected_genes <- setdiff(rownames(ss_Epi), excluded_genes)


ss_Epi <- subset(ss_Epi, features = selected_genes)


# Normalize
ss_Epi <- NormalizeData(ss_Epi)
ss_Epi <- FindVariableFeatures(ss_Epi, selection.method = "vst", nfeatures = 5000)
ss_Epi <- ScaleData(ss_Epi)


ss_Epi <- RunPCA(ss_Epi)
ElbowPlot(ss_Epi)
ss_Epi <- RunUMAP(ss_Epi, dims = 1:20)


ss_Epi <- FindNeighbors(ss_Epi, dims = 1:20)
ss_Epi <- FindClusters(ss_Epi, resolution = 0.3)

DimPlot(ss_Epi, reduction = "umap", label = TRUE)
DimPlot(ss_Epi, reduction = "umap", split.by  = "group",label = TRUE) 

umap_p_1 <- DimPlot(ss_Epi,reduction = "umap",label = T,group.by = 'batch')
umap_p_2 <- DimPlot(ss_Epi,reduction = "umap",label = T,split.by = 'batch',ncol = 2)

ggsave(paste0(output,'/umap_before_harmony.png'), 
       plot = umap_p_1, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/umap_before_harmony_split.png'), 
       plot = umap_p_2, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')


# Batch Effect Removal
ss_harmony <- RunHarmony(ss_Epi, group.by.vars = "batch",max_iter = 20)
ss_harmony <- RunUMAP(ss_harmony, reduction = "harmony", dims = 1:20,reduction.name = "umap")



ss_harmony <- FindNeighbors(ss_harmony,reduction = "harmony",dims = 1:20)
ss_harmony <- FindClusters(ss_harmony,resolution = 0.3)

ss_harmony <- RunUMAP(ss_harmony,reduction = "harmony",dims = 1:20,reduction.name = "umap")

DimPlot(ss_harmony, reduction = "umap", label = TRUE)
DimPlot(ss_harmony, reduction = "umap", split.by  = "group",label = TRUE) 


umap_p_3 <- DimPlot(ss_harmony,reduction = "umap",group.by = "batch",label = T)
umap_p_4 <- DimPlot(ss_harmony,reduction = "umap",split.by = "batch",label = T,ncol = 2)


ggsave(paste0(output,'/umap_after_harmony.png'), 
       plot = umap_p_3, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/umap_after_harmony_split.png'), 
       plot = umap_p_4, device = 'png', height = 150, width = 200, limitsize = FALSE, units = 'mm')



p1 <- DimPlot(ss_Epi, reduction = "umap", label = TRUE)
p2 <- DimPlot(ss_Epi, reduction = "umap", split.by  = "batch",label = TRUE,ncol = 2) 


p3 <- DimPlot(ss_harmony, reduction = "umap", label = TRUE)
p4 <- DimPlot(ss_harmony, reduction = "umap", split.by  = "batch",label = TRUE,ncol = 2) 


p13 <- p1|p3
p24 <- p2|p4

ggsave(paste0(output,'/batch_effect.png'), 
       plot = p13, device = 'png', height = 150, width = 300, limitsize = FALSE, units = 'mm')
ggsave(paste0(output,'/batch_effect_split.png'), 
       plot = p24, device = 'png', height = 150, width = 300, limitsize = FALSE, units = 'mm')




# Removal of Low-Quality Clusters
define_mark <- c("EPCAM")
DotPlot(ss_harmony,
        features = define_mark,
        cluster.idents = TRUE,
        scale = T) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Epi_cell_markers <- list(
  "Common Epithelial Cells" = list(
    "1" = c("EPCAM"),
    "Stem Cells" = c("LGR5", "ASCL2", "SMOC2", "RGMB", "OLFM4"),
    "Paneth Cells" = c("DEFA5", "DEFA6", "REG3A"),
    "Transit-Amplifying (TA)" = c("MKI67", "TOP2A", "PCNA"),
    "Goblet Cells" = c("CLCA1", "SPDEF", "FCGBP", "ZG16", "MUC2"),
    "BEST4 Enterocytes" = c("BEST4", "OTOP2", "CA7"),
    "Enterocytes" = c("RBP2", "ANPEP", "FABP2"),
    "time1" = c("FREM1","PCCA","DMBT1",
                "RBP2","APOA1","APOC3","APOA4","GUCA2B"),
    "Colonocytes" = c("CA2", "SLC26A2", "FABP1"),
    "time2" = c("B3GNT7","ABR","ADH1C","STEAP3","ATP5G1","PCNP", ##Colonocyte
                "AQP8","GUCA2A","CA4","CEACAM1"),
    "Enteroendocrine Cells" = c("CHGA", "CHGB", "NEUROD1")
    #"Microfold Cells" = c("SPIB", "CCL20", "GP2")
    #"Tuft Cells" = c("POU2F3", "LRMP", "TRPM5")
  )
)

all_marker <- unique(unlist(Epi_cell_markers))
DotPlot(ss_harmony,
        features = all_marker,
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


cells_to_keep <- WhichCells(ss_harmony, 
                            idents = setdiff(unique(Idents(ss_harmony)), c("11","13")))


ss_filtered <- subset(ss_harmony, cells = cells_to_keep)
ncol(ss_Epi)
ncol(ss_harmony)
ncol(ss_filtered) 

DimPlot(ss_filtered,reduction = "umap",group.by = "seurat_clusters", label = TRUE)


EPCAM_expression <- FetchData(ss_filtered, vars = "EPCAM")
ss_EPCAM <- ss_filtered[, EPCAM_expression > 0]
table(Idents(ss_filtered))
table(Idents(ss_EPCAM))
ncol(ss_filtered)
ncol(ss_EPCAM)
ss_filtered <- ss_EPCAM



# Re-clustering
ss_filtered <- FindNeighbors(ss_filtered,reduction = "harmony",dims = 1:20)
ss_filtered <- FindClusters(ss_filtered, resolution = 0.5)
ss_filtered <- RunUMAP(ss_filtered,reduction = "harmony",dims = 1:20)


table(Idents(ss_filtered))


DimPlot(ss_filtered, reduction = "umap", label = TRUE)
DimPlot(ss_filtered, reduction = "umap",group.by = "batch")

DotPlot(ss_filtered,
        features = define_mark,
        scale = T,
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1 <- DimPlot(ss_filtered, reduction = "umap", label = TRUE) + 
  ggtitle("Epi")
p1
ggsave(filename = paste0(output,"/","UMAP_filtered_Epi.png"), 
       plot = p1, 
       width = 8, 
       height =6)

p2 <- DimPlot(ss_filtered, reduction = "umap", split.by  = "group",label = TRUE) + 
  ggtitle("Epi")
p2
ggsave(filename = paste0(output,"/UMAP_filtered_Epi_split.png"), 
       plot = p2, 
       width = 8, 
       height = 6)


Epi_cell_markers <- list(
  "all" = c("EPCAM"),
  "Stem Cells" = c("LGR5", "ASCL2", "SMOC2", "RGMB", "OLFM4"),
  "Paneth Cells" = c("DEFA5", "DEFA6", "REG3A"),
  "Transit-Amplifying (TA)" = c("MKI67", "TOP2A", "PCNA"),
  "Goblet Cells" = c("CLCA1", "SPDEF", "FCGBP", "ZG16", "MUC2","TFF3"),
  "BEST4 Enterocytes" = c("BEST4", "OTOP2", "CA7"),
  # "Enterocytes" = c("RBP2", "ANPEP", "FABP2"),
  # "time1" = c("FREM1","PCCA","DMBT1",
  #             "RBP2","APOA1","APOC3","APOA4","GUCA2B"),
  "Colonocytes" = c("CA2", "SLC26A2", "FABP1"),
  "time2" = c("B3GNT7","ABR","ADH1C","STEAP3","ATP5G1","PCNP", ##Colonocyte
              "AQP8","GUCA2A","CA4","CEACAM1"),
  "Enteroendocrine Cells" = c("CHGA", "CHGB", "NEUROD1"),
  "Microfold Cells" = c("SPIB", "CCL20", "GP2"),
  "Tuft Cells" = c("POU2F3", "LRMP", "TRPM5"),
  "LND" = c( "LCN2","NOS2" ,"DUOX2" )

)

all_marker <- unique(unlist(Epi_cell_markers))



DotPlot(ss_filtered,
        features = all_marker,
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



p4 <- DotPlot(ss_filtered,
              features = all_marker,
              cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
p4
ggsave(paste0(output, '/Dot_filtered_Epi.png'), 
       plot = p4,
       width = 20,
       height = 40,  
       limitsize = FALSE, 
       units = 'cm')

##
all_marker <- FindAllMarkers(
  ss_filtered,
  slot = 'data',    
  only.pos = TRUE   
)
write.csv(all_marker,file = paste0(output,'/all_marker_Epi.csv'),fileEncoding = 'UTF-8',row.names = T)



filtered_markers <- all_marker[!grepl("^MT-|^RPS|^RPL", all_marker$gene), ]
write.csv(filtered_markers, file = paste0(output, '/filtered_markers_Epi.csv'), fileEncoding = 'UTF-8', row.names = TRUE)

top10 <- filtered_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
write.csv(top10,file = paste0(output,'/cluster_top10_marker_Epi.csv'),fileEncoding = 'UTF-8',row.names = T)


p5 <- DoHeatmap(ss_filtered, 
                slot = "data", #"counts","scale.data"
                features = top10$gene)+
  NoLegend()

ggsave(filename = paste0(output,"/Heatmap_top10_markers_Epi.png"), 
       plot = p5,
       width = 40, 
       height = 40, 
       limitsize = FALSE, 
       units = 'cm')



# annotation
ss_anno <- ss_filtered
ss_anno@meta.data <- ss_anno@meta.data %>%
  mutate(cell_type = case_when(
    seurat_clusters %in% c(3, 7) ~ "Early Colonocytes",   
    seurat_clusters %in% c(1, 5) ~ "Intermediate Colonocytes", 
    seurat_clusters %in% c(0) ~ "Mature Colonocytes", 
    seurat_clusters == 4  ~ "LND",
    seurat_clusters == 8  ~ "BEST4/OTOP2", 
    seurat_clusters == 6 ~ "Goblet",
    seurat_clusters == 9 ~ "undifferented cell",
    seurat_clusters == 2 ~ "TA",
    seurat_clusters == 10 ~ "Enteroendocrine",
    seurat_clusters == 11 ~ "Tuft"
  ))

table(ss_anno@meta.data$cell_type)
DimPlot(ss_anno, group.by = "cell_type", label = TRUE)

##umap
p_manual <- DimPlot(ss_anno,reduction = "umap",group.by = "cell_type",label = T,label.size = 6)
ggsave(filename = paste0(output,"/","UMAP_anno_Epi.png"), 
       plot = p_manual, 
       width = 10, 
       height = 10)




saveRDS(ss_anno, file = paste0(output, "/ss_anno2_Epi.rds"))

ss <- ss_anno
stem_genes <- c("LGR5","OLFM4")

expression_data <- FetchData(ss, vars = stem_genes)

ss$express_stem_genes <- ifelse(expression_data$LGR5 > 0 & expression_data$OLFM4 > 0, "Express", "Not Express")


table(ss$express_stem_genes)
p_stem <- FeaturePlot(ss, 
                      features = stem_genes,
                      reduction = "umap", 
                      cols = c("lightgrey", "#BC3C29FF")) +
  ggtitle("Cells expressing both LGR5 and OLFM4")


DimPlot(ss, group.by = "express_stem_genes", cols = c("#BC3C29FF", "grey")) + 
  ggtitle("Both LGR5 and OLFM4")

ggsave(filename = paste0(output, "/", "UMAP_stem.png"), 
       width = 6, height = 4, dpi = 300)

ss$cell_type <- ifelse(ss$express_stem_genes == "Express", "Stem", ss$cell_type)
table(ss$cell_type)

DimPlot(ss, group.by = "cell_type", cols = rainbow(length(unique(ss$cell_type)))) +
  ggtitle("Cell Types with Stem Label")

saveRDS(ss, file = paste0(output, "/ss_anno3_Epi.rds"))






