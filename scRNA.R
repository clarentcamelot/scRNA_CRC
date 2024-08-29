library(scImpute)
library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(ComplexHeatmap) 
library(dendextend)
library(circlize)
library(RColorBrewer)
library(infercnv)
library(rjags)
library(limma)
library(GSVA)
library(GSEABase)
library(ggh4x)
library(ggplot2)
library(gplots)
library(monocle)
library(ggsci)
library(CellChat)


scimpute(count_path = 'D:/Rwork/matrix/GSM4904235_Patient1_adenoma.rds',infile = 'rds',outfile = 'csv', out_dir = 'D:/Rwork/p1_ad',Kcluster = 15,ncores = 1)
scimpute(count_path = 'D:/Rwork/matrix/GSM4904236_Patient1_carcinoma.rds',infile = 'rds',outfile = 'csv', out_dir = 'D:/Rwork/p1_ca',Kcluster = 15,ncores = 1)
scimpute(count_path = 'D:/Rwork/matrix/GSM4904247_Patient1_normal.rds',infile = 'rds',outfile = 'csv', out_dir = 'D:/Rwork/p1_no',Kcluster = 15,ncores = 1)

filelist <- list.files(path = 'D:/Rwork/scimpute/',full.names = TRUE, recursive = FALSE)
filelist <- filelist[1:3]
for (file in filelist) {
  data <- read.csv(file = file, row.names = 1)
  file_name <- basename(file)
  project_name <- sub("\\.csv$", "", file_name)
  seurat_obj <- CreateSeuratObject(counts = data, project = project_name, min.cells = 3, min.features = 200)
  seurat_list[[project_name]] <- seurat_obj
}

seurat_merge <- merge(seurat_list[[1]],y = seurat_list[c(2:3)],add.cell.ids = names(seurat_list))

seurat_merge[["percent.mt"]] <- PercentageFeatureSet(seurat_merge, pattern = "^MT-")

seurat_merge@meta.data <- seurat_merge@meta.data  %>%
  mutate(orig.ident = case_when(
    orig.ident == "Patient1_adenoma" ~ "Adenoma",
    orig.ident == "Patient1_carcinoma" ~ "Carcinoma",
    orig.ident == "Patient1_normal" ~ "Normal",
    TRUE ~ orig.ident 
  ))
seurat_merge$orig.ident <- factor(seurat_merge$orig.ident, levels = c('Normal', 'Adenoma', 'Carcinoma'))

qc_vln2 <- VlnPlot(seurat_merge, features =c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3, group.by = 'orig.ident') + scale_fill_observable()
qc_vln2[[1]] <- qc_vln2[[1]] + ggtitle('Gene count')+ scale_fill_observable()
qc_vln2[[2]] <- qc_vln2[[2]] + ggtitle('UMI count')+ scale_fill_observable()
qc_vln2[[3]] <- qc_vln2[[3]] + ggtitle('Mitochondrial gene level')+ scale_fill_observable()

plot1 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident' )+
  labs(x = 'UMI count', y = 'Mitochondrial gene level')+scale_color_observable()
plot2 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')+
  labs(x = 'UMI count', y = 'Gene count') + scale_color_observable()
qc_sca2 <- plot1+plot2
ggsave(filename = 'qc_vln2.jpg',qc_vln2,width = 12,height = 7)
ggsave(filename = 'qc_sca2.jpg',qc_sca2,width = 12,height = 6)

rm(seurat_list)
rm(seurat_merge)
seurat_filt <- NormalizeData(seurat_filt, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_filt <- FindVariableFeatures(seurat_filt, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(seurat_filt), 20)
plot3 <- VariableFeaturePlot(seurat_filt)
plot4 <- LabelPoints(plot = plot3, points = top20, repel = TRUE)
plot3
plot4

all.genes <- rownames(seurat_filt)
seurat_filt <- ScaleData(seurat_filt, features = all.genes)

Idents(seurat_filt) <- 'orig.ident'
seurat_filt <- RunPCA(seurat_filt, features = VariableFeatures(object = seurat_filt))
VizDimLoadings(seurat_filt, dims = 1:2, reduction = "pca")
DimPlot(seurat_filt, reduction = "pca")
seurat_harmony <- RunHarmony(seurat_filt,group.by.vars = 'orig.ident')
DimPlot(seurat_harmony, reduction = "harmony")


ElbowPlot(seurat_harmony,ndims = 50,reduction = 'harmony')
seurat_harmony <- FindNeighbors(seurat_harmony,dims = 1:20,reduction = 'harmony')
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.6, cluster.name = 'harmony_clusters')
seurat_harmony <- RunUMAP(seurat_harmony,dims = 1:20,n.neighbors = 30, min.dist = 0.3,reduction = 'harmony',reduction.name = 'umap.harmony')
seurat_harmony <- RunTSNE(seurat_harmony,dims = 1:20,reduction = 'harmony',reduction.name = 'tsne.harmony')

seurat_harmony <- JoinLayers(seurat_harmony)
Idents(seurat_harmony) <- seurat_harmony$harmony_clusters
all.markers <- FindAllMarkers(seurat_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top50 <- all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

markers <- c(
  'GUCA2A', #enterocyte 
  "GUCA2B", "CA1", 'SLC26A3','MUC2', 
  "TFF3", #Goblet
  'EPCAM',  #epithelial
  'MS4A1', 'CD79A','CD79B',  #follicular B
  'MZB1',  # Plasma B
  'CD3D',   #T
  'CD8A','CD3E' ,  #CD8+ T
  'IL7R',   #central memory T 
  'KLRD1',  #natural killer T 
  'CD68', 'CD86', #M1 macrophage 
  'CD163','IL1B',  #M2 macrophage
  'CD14', 'FCGR3A' #monocyte
)
VlnPlot(seurat_harmony, features = markers, stack = TRUE, sort = TRUE)


epi_clusters <- c(14,19,15,13,16,5,4,6,2,7,18,21,11,20)

celltype <- data.frame(ClusterID=0:22, celltype='others')
celltype$celltype[celltype$ClusterID %in% epi_clusters] <- 'Epithelial Mix'
celltype$celltype[celltype$ClusterID %in% c(17)] <- 'Monocytes & Macrophages'
celltype$celltype[celltype$ClusterID %in% c(0,1,3,9,10)] <- 'T cells'
#celltype$celltype[celltype$ClusterID %in% c(10)] <- 'Tcm'
#celltype$celltype[celltype$ClusterID %in% c(9)] <- 'T?'
celltype$celltype[celltype$ClusterID %in% c(8)] <- 'Plasma B'
celltype$celltype[celltype$ClusterID %in% c(12)] <- 'Follicular B'
celltype$celltype[celltype$ClusterID %in% c(14,19)] <- 'Enterocytes'
celltype$celltype[celltype$ClusterID %in% c(20)] <- 'Goblet'
celltype$celltype[celltype$ClusterID %in% c(4,13,21)] <- 'Epithelial Carcinoma'
celltype$celltype[celltype$ClusterID %in% c(2,5)] <- 'Epithelial Adenoma'

sce <- seurat_harmony
sce@meta.data$celltype = 'NA'
for (i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$harmony_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}

sce$celltype <- factor(sce$celltype, levels = c('Enterocytes', 'Goblet', 'Epithelial Mix','Epithelial Adenoma','Epithelial Carcinoma',
                                                'Follicular B','Plasma B','T cells','Monocytes & Macrophages'))

p_tsne <- DimPlot(sce, reduction = "tsne.harmony", group.by = 'celltype', label = TRUE, label.size = 7)  +
  labs(title = 'Cell type')+
  xlim(-55,75)+ 
  ylim(-50,50)+
  theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  scale_color_observable()
p_umap <- DimPlot(sce, reduction = "umap.harmony", group.by = 'celltype', label = TRUE, label.size = 7) +
  labs(x = "UMAP_1", y = "UMAP_2",title = 'Cell type') +
  theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  scale_color_observable()
ggsave("p1_tsne.png", p, width = 26,height = 11)

p_tsne_t <- DimPlot(sce, reduction = "tsne.harmony", group.by = 'orig.ident')  +
  labs(title = 'Tissue distribution')+
  xlim(-55,75)+ 
  ylim(-50,50)+
  theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  scale_color_observable()
p_umap_t <- DimPlot(sce, reduction = "umap.harmony", group.by = 'orig.ident') +
  labs(x = "UMAP_1", y = "UMAP_2",title = 'Tissue distribution') +
  theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  scale_color_observable()



p <- DotPlot(sce, features = markers, group.by = 'celltype')
dat <- p$data
dot <- ggplot(dat, aes(id,features.plot,size=pct.exp, fill=avg.exp.scaled)) + 
  geom_point(shape = 21) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank(),
    axis.text.y = element_text(color='black',size=10),
    axis.text.x = element_text(color='black',size=10, angle = 45, hjust = 1, vjust = 1))+
  scale_fill_gradientn(colours = c("#FFFFFF","#FFF4F4", "#FF0101"), values = scales::rescale(c(-2,-1, 1)))

data_summary <- sce@meta.data %>%
  group_by(celltype, orig.ident) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(celltype) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup()


empty_rows <- data.frame(
  celltype = c(" ", "  ", "   "), 
  fraction = c(NA, NA, NA),
  orig.ident = c("Normal", "Adenoma", "Carcinoma"), 
  count=c(NA, NA, NA),
  stringsAsFactors = FALSE      
)

data_summary <- rbind(data_summary, empty_rows)

data_summary$celltype <- factor(data_summary$celltype, 
                                levels = rev(c('Enterocytes', 'Goblet', 'Epithelial Mix','Epithelial Adenoma','Epithelial Carcinoma',' ',
                                               'Follicular B','Plasma B','  ',
                                               'T cells','   ','Monocytes & Macrophages')))


data_summary$orig.ident <- factor(data_summary$orig.ident, 
                                  levels = c("Normal",'Adenoma','Carcinoma'))

bar <- ggplot(data_summary, aes(y = celltype, 
                                x = fraction, 
                                fill = orig.ident)) +
  geom_bar(stat = "identity", na.rm = TRUE) +  
  scale_x_reverse() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    strip.text.y = element_blank(),  # 去掉分面标签
    panel.spacing = unit(0.5, "lines"),
    axis.ticks.y = element_blank()  # 不显示空白的刻度
  ) +
  labs(
    y = "Cell Type",
    x = "Fraction of cells",
    fill = "Tissue") +
  scale_fill_observable()
ggsave(filename = 'barplot0.jpg',bar, width = 7,height = 7)


Idents(sce)= 'celltype'
data.input = sce@assays$RNA$data 
meta = sce@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
data.input <- sce[["RNA"]]$data
labels <- Idents(sce)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
meta$samples <- factor(c(rep("sample1", 10112) ))

cellchat <- createCellChat(object = sce, group.by = "celltype", assay = "RNA")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
netAnalysis_river(cellchat, pattern = "outgoing")
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)

saveRDS(file = 'cellchat.rds',cellchat)

cellorder <- c('Follicular B','Plasma B','Enterocytes', 'Epithelial Mix','Epithelial Adenoma','Epithelial Carcinoma', 'Goblet',
               'Monocytes & Macrophages','T cells')

netAnalysis_contribution(cellchat, signaling = 'MHC-I')
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = 'MHC-I', geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
netVisual_individual(cellchat, signaling = 'MHC-I', pairLR.use = LR.show, layout = "chord" )
pdf("cellchat.pdf", width = 12, height = 10)
netVisual_individual(cellchat, signaling = 'MHC-I', pairLR.use = LR.show, layout = "chord",cell.order =cellorder )


#netAnalysis_contribution(cellchat, signaling = 'CD200')
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = 'CD200', geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
netVisual_individual(cellchat, signaling = 'CD200', pairLR.use = LR.show, layout = "chord",cell.order =cellorder)


#netAnalysis_contribution(cellchat, signaling = 'FASLG')
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = 'FASLG', geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
netVisual_individual(cellchat, signaling = 'FASLG', pairLR.use = LR.show, layout = "chord",cell.order =cellorder)

dev.off()


epi_clusters <- c(14,19,15,13,16,5,4,6,2,7,18,21,11,20)

celltype <- data.frame(ClusterID=0:22, celltype='others')
celltype$celltype[celltype$ClusterID %in% epi_clusters] <- 'Epithelial Cells'

for (i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$harmony_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}

epi <- subset(sce, idents = c('Epithelial Cells'))
epi_obj = CreateSeuratObject(counts = GetAssayData(epi,assay = 'RNA',layer = 'counts'),meta.data = epi@meta.data)
rm(epi)
epi_obj <- NormalizeData(epi_obj, normalization.method = "LogNormalize", scale.factor = 10000)
epi_obj <- FindVariableFeatures(epi_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(epi_obj)
epi_obj <- ScaleData(epi_obj, features = all.genes)

epi_obj <- RunPCA(epi_obj, features = VariableFeatures(object = epi_obj))
epi_harmony <- RunHarmony(epi_obj,group.by.var = 'orig.ident' )

ElbowPlot(epi_obj,ndims = 50)
epi_obj <- FindNeighbors(epi_obj,dims = 1:20)
epi_obj <- FindClusters(epi_obj, resolution = 0.2)

ElbowPlot(epi_harmony,ndims = 50, reduction = 'harmony')
epi_harmony <- FindNeighbors(epi_harmony,dims = 1:25,reduction = 'pca')
epi_harmony <- FindClusters(epi_harmony, resolution = 0.15, cluster.name = 'noh_clusters')
DimPlot(epi_obj, reduction = "tsne",group.by = 'orig.ident')
epi_obj <- RunUMAP(epi_obj,dims = 1:20)
DimPlot(epi_obj, reduction = "umap",group.by = 'orig.ident')
DimPlot(epi_obj, reduction = "umap")
epi_harmony <- RunTSNE(epi_harmony,dims = 1:25)
DimPlot(epi_harmony, reduction = "tsne",group.by = 'orig.ident')
DimPlot(epi_harmony, reduction = "tsne")


epi_harmony <- RunUMAP(epi_harmony,dims = 1:25,n.neighbors = 30, min.dist = 0.3,reduction = 'harmony',reduction.name = 'umap.harmony')
epi_harmony <- RunTSNE(epi_harmony,dims = 1:25,reduction = 'harmony',reduction.name = 'tsne.harmony')

DimPlot(epi_harmony, reduction = "tsne.harmony")
DimPlot(epi_harmony, reduction = "tsne.harmony",group.by = 'orig.ident')
DimPlot(epi_harmony, reduction = "umap.harmony")
DimPlot(epi_harmony, reduction = "umap.harmony",group.by = 'orig.ident')

cl <- data.frame(ClusterID=0:9, celltype='others')
cl$celltype[celltype_fin$ClusterID %in% c(0)] <- 1
cl$celltype[celltype_fin$ClusterID %in% c(1)] <- 2
cl$celltype[celltype_fin$ClusterID %in% c(2)] <- 3
cl$celltype[celltype_fin$ClusterID %in% c(3)] <- 4
cl$celltype[celltype_fin$ClusterID %in% c(4)] <- 5
cl$celltype[celltype_fin$ClusterID %in% c(5)] <- 6
cl$celltype[celltype_fin$ClusterID %in% c(6)] <- 7
cl$celltype[celltype_fin$ClusterID %in% c(7)] <- 8
cl$celltype[celltype_fin$ClusterID %in% c(8)] <- 9
cl$celltype[celltype_fin$ClusterID %in% c(9)] <- 10

for (i in 1:nrow(cl)){
  epi_harmony@meta.data[which(epi_harmony@meta.data$harmony_clusters == cl$ClusterID[i]),'cl'] <- cl$celltype[i]
}
matrix_counts <- as.matrix(GetAssayData(epi_harmony,assay = 'RNA',layer = 'counts'))
epi_harmony$newcelltype <- paste0(epi_harmony$orig.ident,'_',epi_harmony$cl)
unique(epi_harmony$newcelltype)
write.table(epi_harmony$newcelltype,"celltype_label.txt",sep = '\t',quote = F, col.names = F)

library(readr)
gencode <- read_tsv('hg38_gencode_v27.txt',col_names = c('gene','chr','start','end'))
gencode <- gencode[!duplicated(gencode$gene),]
common_genes <- intersect(gencode$gene,rownames(matrix_counts))
groupnames <- sort(unique(epi_harmony$newcelltype))
infer_obj <- CreateInfercnvObject(raw_counts_matrix = matrix_counts[common_genes,],
                                  annotations_file = "D:/Rwork/celltype_label.txt",
                                  delim = '\t' ,
                                  gene_order_file = "D:/Rwork/hg38_gencode_v27.txt",
                                  ref_group_names = groupnames[19:26],
                                  chr_exclude = c('chrY','chrM'))


colnames <- names(infer_obj@reference_grouped_cell_indices)
numeric_order <- as.numeric(gsub("Normal_", "", colnames))
sorted_colnames <- colnames[order(numeric_order)]
infer_obj@reference_grouped_cell_indices <- infer_obj@reference_grouped_cell_indices[sorted_colnames]


colnames <- names(infer_obj@observation_grouped_cell_indices)
numeric_order1 <- as.numeric(gsub("Adenoma_", "", colnames[1:9]))
numeric_order2 <- as.numeric(gsub("Carcinoma_", "", colnames[10:18]))
cm1 <- colnames[1:9]
cm2 <- colnames[10:18]
sorted_colnames1 <- cm1[order(numeric_order1)]
sorted_colnames2 <- cm2[order(numeric_order2)]
sorted_colnames <- c(sorted_colnames1 ,sorted_colnames2)

infer_obj@observation_grouped_cell_indices <- infer_obj@observation_grouped_cell_indices[sorted_colnames]

infer_obj = infercnv::run(infercnv_obj = infer_obj,
                          cutoff = 0.1,
                          out_dir = './infercnv3',
                          no_prelim_plot = T,
                          cluster_by_groups = T,
                          denoise = TRUE,
                          HMM = F,
                          min_cells_per_gene = 10,
                          num_threads = 10,
                          write_expr_matrix = T
)


saveRDS(epi_harmony,file = 'epi_harmony.rds')
saveRDS(infer_obj,file = 'infer_obj.rds')


epi_harmony <- readRDS("D:/Rwork/epi_harmony.rds")
celltype_fin <- data.frame(ClusterID=0:9, celltype='others')
celltype_fin$celltype[celltype_fin$ClusterID %in% c(0)] <- '1_Benign'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(1)] <- '2_Malignant_1'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(2)] <- '3_Mix_1'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(3)] <- '4_Enterocytes'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(4)] <- '5_Mix_2'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(5)] <- '6_Mix_3'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(6)] <- '7_Mix_4'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(7)] <- '8_Malignant_2'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(8)] <- '9_Goblet'
celltype_fin$celltype[celltype_fin$ClusterID %in% c(9)] <- '10_Mix_5'

cell_type <- data.frame(ClusterID=0:9, celltype='Others')
cell_type$celltype[cell_type$ClusterID %in% c(3,8)] <- 'Physiological'
cell_type$celltype[cell_type$ClusterID %in% c(0)] <- 'Benign'
cell_type$celltype[cell_type$ClusterID %in% c(1,7)] <- 'Malignant'
cell_type$celltype[cell_type$ClusterID %in% c(5)] <- '6_Mix_3'
#cell_type$celltype[cell_type$ClusterID %in% c(9)] <- '10_Mix_5'
cell_type$celltype[cell_type$ClusterID %in% c(2)] <- '3_Mix_1'



epi_harmony@meta.data <- epi_harmony@meta.data  %>%
  mutate(orig.ident = case_when(
    orig.ident == "Patient1_adenoma" ~ "Adenoma",
    orig.ident == "Patient1_carcinoma" ~ "Carcinoma",
    orig.ident == "Patient1_normal" ~ "Normal",
    TRUE ~ orig.ident 
  ))


sce_epi <- epi_harmony
#sce_epi@meta.data$celltype_fin = 'NA'
for (i in 1:nrow(celltype_fin)){
  sce_epi@meta.data[which(sce_epi@meta.data$harmony_clusters == celltype_fin$ClusterID[i]),'celltype_fin'] <- celltype_fin$celltype[i]
}
for (i in 1:nrow(cell_type)){
  sce_epi@meta.data[which(sce_epi@meta.data$harmony_clusters == cell_type$ClusterID[i]),'cell_type'] <- cell_type$celltype[i]
}


Idents(sce_epi) <- sce_epi$celltype_fin
all.markers2 <- FindAllMarkers(sce_epi, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)

top50 <- all.markers2 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

sce_epi$celltype_fin <- factor(sce_epi$celltype_fin, levels = c('1_Benign','2_Malignant_1','3_Mix_1', '4_Enterocytes','5_Mix_2','6_Mix_3','7_Mix_4','8_Malignant_2','9_Goblet','10_Mix_5'))

p_tsne <- DimPlot(sce_epi, reduction = "tsne.harmony", group.by = 'celltype_fin', label = T, label.size = 7)  +
  labs(x = "tSNE_1", y = "tSNE_2",title = 'Cell type')+theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  scale_color_npg()
p_umap <- DimPlot(sce_epi, reduction = "umap.harmony", group.by = 'celltype_fin', label = T, label.size = 7) +
  labs(x = "UMAP_1", y = "UMAP_2",title = 'Cell type') +theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  scale_color_npg()
p <- p_tsne+p_umap
ggsave("epi_tsne.png", p, width = 22,height = 11)

sce_epi$orig.ident <- factor(sce_epi$orig.ident, 
                             levels = c("Normal",'Adenoma','Carcinoma'))


DimPlot(sce_epi, reduction = "tsne.harmony", group.by = 'orig.ident')  +
  labs(x = "tSNE_1", y = "tSNE_2",title = 'Tissue distribution')+
  scale_color_observable()
normal_cells <- subset(sce_epi, subset = orig.ident == "Normal")
Adenoma_cells <- subset(sce_epi, subset = orig.ident == "Adenoma")
Carcinoma_cells <- subset(sce_epi, subset = orig.ident == "Carcinoma")


p1<- DimPlot(sce_epi, reduction = "tsne.harmony", group.by = "orig.ident") +
  xlim(-45,45)+ 
  ylim(-60,54)+
  theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  ggtitle("All tissues")+scale_color_observable()+labs(x = "tSNE_1", y = "tSNE_2")
p2 <- DimPlot(normal_cells, reduction = "tsne.harmony", group.by = "orig.ident") +
  ggtitle("Normal") +
  xlim(-45,45)+ 
  ylim(-60,54)+
  theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  labs(x = "tSNE_1", y = "tSNE_2") +
  scale_color_manual(values = c("Normal" = '#4269D0FF'))

p3 <- DimPlot(Adenoma_cells, reduction = "tsne.harmony", group.by = "orig.ident") +
  ggtitle("Adenoma") +
  xlim(-45,45)+ 
  ylim(-60,54)+
  theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  labs(x = "tSNE_1", y = "tSNE_2") +
  scale_color_manual(values = c("Adenoma" = '#EFB118FF'))

p4 <- DimPlot(Carcinoma_cells, reduction = "tsne.harmony", group.by = "orig.ident") +
  ggtitle("Carcinoma") +
  xlim(-45,45)+ 
  ylim(-60,54)+
  theme(
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 20))+
  labs(x = "tSNE_1", y = "tSNE_2") +
  scale_color_manual(values = c("Carcinoma" = '#FF725CFF'))

p <- (p1 | p2) / (p3 | p4)
ggsave("epi_tsne_tissue.png", p, width = 18,height = 15)


sce_epi$celltype_fin <- factor(sce_epi$celltype_fin, levels = c( '4_Enterocytes',
                                                                 '9_Goblet',
                                                                 '1_Benign',
                                                                 '2_Malignant_1',
                                                                 '8_Malignant_2',
                                                                 '3_Mix_1',
                                                                 '5_Mix_2',
                                                                 '6_Mix_3',
                                                                 '7_Mix_4',
                                                                 '10_Mix_5'))


data_summary <- sce_epi@meta.data %>%
  group_by(celltype_fin, orig.ident) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(celltype_fin) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup()


empty_rows <- data.frame(
  celltype_fin = c(" ", "  ", "   "), 
  fraction = c(NA, NA, NA),
  orig.ident = c("Normal", "Adenoma", "Carcinoma"), 
  count=c(NA, NA, NA),
  stringsAsFactors = FALSE  
)

data_summary <- rbind(data_summary, empty_rows)

data_summary$celltype_fin <- factor(data_summary$celltype_fin, 
                                    levels = c('10_Mix_5','7_Mix_4', '6_Mix_3', '5_Mix_2', '3_Mix_1' ,' ',
                                               '8_Malignant_2','2_Malignant_1','  ',
                                               '1_Benign','   ',
                                               '9_Goblet', '4_Enterocytes'))


data_summary$orig.ident <- factor(data_summary$orig.ident, 
                                  levels = c("Normal",'Adenoma','Carcinoma'))

bar <- ggplot(data_summary, aes(y = celltype_fin, 
                                x = fraction, 
                                fill = orig.ident)) +
  geom_bar(stat = "identity", na.rm = TRUE) +  
  scale_x_reverse() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    strip.text.y = element_blank(),  
    panel.spacing = unit(0.5, "lines"),
    axis.ticks.y = element_blank()  
  ) +
  labs(
    y = "Cell Type",
    x = "Fraction of cells",
    fill = "Tissue") +
  scale_fill_observable()


markers2 <- c("BST2", "TNFRSF11B", "NDUFA4L2", "CYP2W1", "TGFB1", "OLFM4", "MYC", "SIVA1", "ASCL2", "DPEP1", "LEFTY1", "PCCA", "REG1A", "SMOC2", "APCDD1", "TFF3", "MUC2", "GUCA2B", "GUCA2A", "CA1", "SLC26A3")

Idents(sce_epi) = 'celltype_fin'
selected_groups <- c("4_Enterocytes", "9_Goblet", "1_Benign", "2_Malignant_1",'8_Malignant_2')
filtered_seurat_object <- subset(sce_epi, idents = selected_groups)

ca.markers <- FindAllMarkers(filtered_seurat_object, only.pos = TRUE, min.pct = 0.9, logfc.threshold = 1)

top50ca <- ca.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


markers2 <-rev(c( "GUCA2B", "GUCA2A", "CA1", "SLC26A3","TFF3", "MUC2",
                  "REG1A", "SMOC2", "APCDD1",'RASSF10',
                  "PCCA", "DPEP1", "ASCL2", "SIVA1", "OLFM4", "MYC",
                  "BST2", "TNFRSF11B", "CYP2W1", "TGFB1",
                  'TRIB3','CA9','WNT4','RGS16'
)) 

p <- DotPlot(filtered_seurat_object, features = markers2, group.by = 'celltype_fin')
dat <- p$data
markers_group_df <- do.call(rbind, lapply(names(markers_group), function(group) {
  data.frame(Group = group, Gene = markers_group[[group]], stringsAsFactors = FALSE)
}))
cluster.order <- c("4_Enterocytes", "9_Goblet", "1_Benign", "2_Malignant_1",'8_Malignant_2')
df <- left_join(dat, markers_group_df, by = c("features.plot" = "Gene"))
df$id <- factor(df$id, levels= cluster.order)
df$features.plot <- factor(df$features.plot, levels= rev(markers2))

dot <- ggplot(df, aes(id,interaction(features.plot, Group),size=pct.exp, fill=avg.exp.scaled)) + 
  geom_point(shape = 21) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title = element_blank(),
    axis.text.y = element_text(color='black',size=10),
    axis.text.x = element_text(color='black',size=10, angle = 45, hjust = 1, vjust = 1))+
  scale_fill_gradientn(colours = c("#FFFFFF","#FFF4F4", "#FF0101"), values = scales::rescale(c(-2,-1, 1)))+
  guides(y = "axis_nested") +
  theme(ggh4x.axis.nestline.y = element_line(size=4,
                                             colour = c('#E58606', '#5D69B1', '#52BCA3', '#99C945','#9467BDFF')))+ 
  geom_rect(aes(xmin = 0.8, xmax = 2.2, ymin = 0.6, ymax = 6.4), color = "red", fill = NA, size = 0.7) + 
  geom_rect(aes(xmin = 2.8, xmax = 3.2, ymin = 6.6, ymax = 10.4), color = "red", fill = NA, size = 0.7) + 
  geom_rect(aes(xmin = 2.8, xmax = 5.2, ymin = 10.6, ymax = 16.4), color = "red", fill = NA, size = 0.7) +
  geom_rect(aes(xmin = 3.8, xmax = 5.2, ymin = 16.6, ymax = 20.4), color = "red", fill = NA, size = 0.7) + 
  geom_rect(aes(xmin = 3.8, xmax = 5.2, ymin = 20.6, ymax = 24.4), color = "red", fill = NA, size = 0.7) 

ggsave(filename = 'dotplot2.jpg',dot, width = 7,height = 7)


plots <- FeaturePlot(
  object = sce_epi, 
  features = c("CA1", "REG1A", "BST2", 
               "GUCA2B", "SMOC2", "CA9", 
               "MUC2", "APCDD1", "RGS16"), 
  cols = c('grey', 'red'), 
  reduction = 'tsne.harmony', 
  min.cutoff = 'q10'
)

for (i in 1:9) {
  plots[[i]][["labels"]][["x"]] <-'tSNE_1'
  plots[[i]][["labels"]][["y"]] <-'tSNE_2'
}
ggsave(filename = 'feature.png',plots,width = 9, height = 9)


Idents(sce_epi) = 'cell_type'
selected_groups <- c('Physiological','Benign','Malignant')
selected_groups <- c('6_Mix_3')
sce_epi_flit <- subset(sce_epi, idents = selected_groups)

sce_epi_flit$cell_type <- factor(sce_epi_flit$cell_type,levels = c('Physiological','Benign','Malignant') )

gsva_data <- as.data.frame(sce_epi_flit@assays$RNA$data)
GeneSet1 <- getGmt("D:/Rwork/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
GeneSet2 <- getGmt("D:/Rwork/c5.go.bp.v2023.2.Hs.symbols.gmt")
GeneSet3 <- getGmt("D:/Rwork/c6.all.v2023.2.Hs.symbols.gmt")
GeneSet4 <- getGmt("D:/Rwork/c7.all.v2023.2.Hs.symbols.gmt")

GeneSet <- GeneSetCollection(c(GeneSet1, GeneSet2, GeneSet3, GeneSet4))

gsva_param <- gsvaParam(as.matrix(gsva_data),GeneSet)
gsva_result <- gsva(gsva_param, verbose = TRUE, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
saveRDS(gsva_result, file = 'gsva_result.rds')



table(sce_epi_flit@meta.data$cell_type)
cluster_bp <- sce_epi_flit@meta.data %>% dplyr::filter(cell_type %in% c('Physiological', "Benign")) %>% dplyr::arrange(cell_type)
use_gsva_bp <- as.data.frame(gsva_result) %>% dplyr::select(rownames(cluster_bp))
group_bp <- c(rep("Physiological", 555), rep("Benign", 1492)) %>% as.factor()
desigN_bp <- model.matrix(~ 0 + group_bp) 
colnames(desigN_bp) <- levels(group_bp)
fit_bp = lmFit(use_gsva_bp, desigN_bp)
cont.matrix <- makeContrasts(contrasts = c('Benign-Physiological'), levels = desigN_bp)
fit2_bp <- contrasts.fit(fit_bp, cont.matrix)
fit3_bp <- eBayes(fit2_bp)


cluster_mp <- sce_epi_flit@meta.data %>% dplyr::filter(cell_type %in% c('Physiological','Malignant')) %>% dplyr::arrange(cell_type)
use_gsva_mp <- as.data.frame(gsva_result) %>% dplyr::select(rownames(cluster_mp))
group_mp <- c(rep("Physiological", 555), rep("Malignant", 1396)) %>% as.factor()
desigN_mp <- model.matrix(~ 0 + group_mp) 
colnames(desigN_mp) <- levels(group_mp)
fit_mp = lmFit(use_gsva_mp, desigN_mp)
cont.matrix <- makeContrasts(contrasts = c('Malignant-Physiological'), levels = desigN_mp)
fit2_mp <- contrasts.fit(fit_mp, cont.matrix)
fit3_mp <- eBayes(fit2_mp)

cluster_bm <- sce_epi_flit@meta.data %>% dplyr::filter(cell_type %in% c('Benign','Malignant')) %>% dplyr::arrange(cell_type)
use_gsva_bm <- as.data.frame(gsva_result) %>% dplyr::select(rownames(cluster_bm))
group_bm <- c(rep("Benign", 1492), rep("Malignant", 1396)) %>% as.factor()
desigN_bm <- model.matrix(~ 0 + group_bm) 
colnames(desigN_bm) <- levels(group_bm)
fit_bm = lmFit(use_gsva_bm, desigN_bm)
cont.matrix <- makeContrasts(contrasts = c('Benign-Malignant'), levels = desigN_bm)
fit2_bm <- contrasts.fit(fit_bm, cont.matrix)
fit3_bm <- eBayes(fit2_bm)
summary(decideTests(fit3_bp,lfc=0.5, p.value=0.05)) #统计查看差异结果

diff_bp <- topTable(fit3_bp,adjust='fdr', coef=1, number=Inf)
cluster_diff_bp <- na.omit(diff_bp)
cluster_diff_bp$pathway <- rownames(cluster_diff_bp)

diff_mp <- topTable(fit3_mp,adjust='fdr', coef=1, number=Inf)
cluster_diff_mp <- na.omit(diff_mp)
cluster_diff_mp$pathway <- rownames(cluster_diff_mp)

diff_bm <- topTable(fit3_bm,adjust='fdr', coef=1, number=Inf)
cluster_diff_bm <- na.omit(diff_bm)
cluster_diff_bm$pathway <- rownames(cluster_diff_bm)

deg <- rbind(cluster_diff_bp, cluster_diff_mp,cluster_diff_bm)
sig_deg <- deg %>% dplyr::filter(abs(logFC) > 0.2 & adj.P.Val < 0.05)

old_pathways <- c(
  "GOBP_PHOSPHAGEN_METABOLIC_PROCESS",
  "GOBP_HYDROGEN_SULFIDE_METABOLIC_PROCESS",
  "GOBP_LEUKOTRIENE_D4_BIOSYNTHETIC_PROCESS",
  "GOBP_PHOSPHAGEN_METABOLIC_PROCESS",
  "GOBP_ARGININE_BIOSYNTHETIC_PROCESS",
  "GOBP_CELLULAR_RESPONSE_TO_ACIDIC_PH",
  "GOBP_CELLULAR_RESPONSE_TO_CADMIUM_ION",
  "GOBP_RIBOSOME_BIOGENESIS",
  "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION",
  "GOBP_RIBOSOMAL_SMALL_SUBUNIT_ASSEMBLY",
  "GOBP_POSITIVE_REGULATION_OF_NUCLEAR_TRANSCRIBED_MRNA_CATABOLIC_PROCESS_DEADENYLATION_DEPENDENT_DECAY",
  "GOBP_FORMATION_OF_TRANSLATION_PREINITIATION_COMPLEX",
  "GOBP_RIBOSOME_ASSEMBLY",
  "GOBP_RIBOSOMAL_LARGE_SUBUNIT_ASSEMBLY",
  "GOBP_GROWTH_HORMONE_RECEPTOR_SIGNALING_PATHWAY_VIA_JAK_STAT",
  "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY_INVOLVED_IN_CARDIAC_MUSCLE_CELL_FATE_COMMITMENT",
  "KEGG_MEDICUS_REFERENCE_EREG_EGFR_PI3K_SIGNALING_PATHWAY",
  "KEGG_MEDICUS_REFERENCE_PKR_EIF2ALPHA_SIGNALING_PATHWAY",
  "GOBP_NEGATIVE_REGULATION_OF_INTERLEUKIN_6_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_NEGATIVE_REGULATION_OF_MACROPHAGE_DIFFERENTIATION",
  'GOBP_ANTIBODY_DEPENDENT_CELLULAR_CYTOTOXICITY'
)

new_pathways <- c(
  "Amino or nucleotide sugar metabolism",
  "Sulfur metabolism",
  "Leukotriene D4 biosynthesis",
  "Phosphagen metabolism",
  "Arginine biosynthesis",
  "Cellular response to pH",
  "Cellular response to cadmium ion",
  "Ribosome biogenesis",
  "Translational initiation",
  "Ribosomal small subunit assembly",
  "Decay mRNA catabolic process",
  "Formation of translation preinitiation complex",
  "Ribosome assembly",
  "Ribosomal large subunit assembly",
  "Jak-Stat signaling pathway",
  "Wnt signaling pathway",
  "EGFR PI3K signaling pathway",
  "PKR EIF2α signaling pathway",
  "Negative regulation of macrophage differentation",
  "Negative regulation of IL6 mediated signaling pathway",
  'Antibody dependent cellular cytotoxicity'
)

matched_indices <- match(old_pathways, rownames(gsva_result))
extracted_results <- gsva_result[matched_indices, ]
rownames(extracted_results) <- new_pathways

sce_epi_flit$orig.ident <- factor(sce_epi_flit$orig.ident, levels = c('Normal', 'Adenoma', 'Carcinoma'))
cell_types <- sce_epi_flit@meta.data$cell_type
cell_tissue <- sce_epi_flit@meta.data$orig.ident
cell_ids <- colnames(sce_epi_flit@assays$RNA$data)
cell_info <- data.frame(Cell_ID = cell_ids, Cell_Type = cell_types, Tissue = cell_tissue)
cell_info$cell_types <- factor(cell_types, levels = c( "Benign",'6_Mix_3', 'Malignant'))
cell_info$Tissue <- factor(cell_tissue, levels = c('Normal', 'Adenoma', 'Carcinoma'))
cell_info_ordered <- cell_info[order(cell_info$cell_types), ]

common_columns <- intersect(colnames(gsva_result), cell_info$Cell_ID)
gsva_result_filtered <- gsva_result[, common_columns, drop = FALSE]

annotation_col <- data.frame(
  Tissue = cell_info_ordered$Tissue,CellType = cell_info_ordered$cell_types
)
rownames(annotation_col) <- ordered_cell_ids

col_annotation <- HeatmapAnnotation(
  df = annotation_col,
  col = list(
    Tissue = c("Normal" = "#4269D0FF", 
               "Adenoma" = "#EFB118FF", 
               "Carcinoma" = "#FF725cff"),
    CellType = c("Physiological" = "#0099B4FF", 
                 "Benign" = "#00468BFF", 
                 "Malignant" = "#FDAF91FF")
    
  ),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  annotation_legend_param = list(direction='horizontal')
)
col_annotation@anno_list[["CellType"]]@color_mapping@name <- 'Cell Type'
col_annotation@anno_list[["CellType"]]@label <- 'Cell Type'
col_annotation@anno_list[["CellType"]]@name <- 'Cell Type'

ht <- Heatmap(
  as.matrix(extracted_results_ordered),
  name = "GSVA Score",
  row_title = "Pathways",
  column_title = " ",
  top_annotation = col_annotation,
  show_row_names = TRUE,
  row_names_max_width = unit(10,'cm'),
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  use_raster = TRUE
)
pdf("heat_gsava.pdf", width = 10, height = 10)
draw(ht)
dev.off()


Idents(sce_epi) = 'cell_type'
selected_groups <- c("Benign", "Malignant",'6_Mix_3')
sce_epi_flit4 <- subset(sce_epi, idents = selected_groups)

sce_epi_flit4$cell_type <- factor(sce_epi_flit4$cell_type, levels = c("Benign",'6_Mix_3', "Malignant"))

exp <- as(as.matrix(sce_epi_flit4[["RNA"]]$counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sce_epi_flit4@meta.data)
fData <- data.frame(gene_short_name = row.names(exp), row.names = row.names(exp))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(exp,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily=VGAM::negbinomial.size(),
                      lowerDetectionLimit=1)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) 

pData(cds)$Cluster  <- as.factor(pData(cds)$cell_type)
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~Cluster")
head(diff_test_res)

ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:900]
cds <- setOrderingFilter(cds,ordering_genes = ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds,reverse = T)

p1 <- plot_cell_trajectory(cds, color_by = "Cluster",cell_size = 0.8)
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.8)
p1 <- p1 + theme(
  legend.text = element_text(size = 14), 
  legend.title = element_text(size = 16) 
)+ guides(color = guide_legend(override.aes = list(size = 6))) + scale_color_manual(breaks = c("Benign", "6_Mix_3", "Malignant"),  # 设置类群的顺序
                                                                                    values = c("#00468BFF", "#925e9fff", "#FDAF91FF"))
p <-p1 + p2
ggsave(filename = 'pst7.png',p, width = 16,height = 8)



Idents(sce_epi) = 'cell_type'
selected_groups <- c('6_Mix_3')
sce_epi_flit <- subset(sce_epi, idents = selected_groups)

gsva_data <- as.data.frame(sce_epi_flit@assays$RNA$data)


gsva_param <- gsvaParam(as.matrix(gsva_data),GeneSet)
gsva_result6 <- gsva(gsva_param, verbose = TRUE, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
Idents(sce_epi) = 'cell_type'
selected_groups <- c( "Benign",'6_Mix_3', 'Physiological')
sce_epi_flit <- subset(sce_epi, idents = selected_groups)

sce_epi_flit <- JoinLayers(sce_epi_flit)
Idents(sce_epi_flit) <- sce_epi_flit$celltype_fin
all.markers <- FindAllMarkers(sce_epi_flit, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top50 <- all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

sce_epi_flit$cell_type <- factor(sce_epi_flit$cell_type, levels =  c( "Benign",'6_Mix_3', 'Malignant') )

DoHeatmap(sce_epi_flit, features = top50$gene, group.by = 'cell_type') + NoLegend()


Idents(sce_epi) = 'cell_type'
selected_groups <- c( "Benign",'6_Mix_3', 'Malignant')
sce_epi_flit <- subset(sce_epi, idents = selected_groups)

sce_epi_flit$orig.ident <- factor(sce_epi_flit$orig.ident, levels = c('Normal', 'Adenoma', 'Carcinoma'))
cell_types <- sce_epi_flit@meta.data$cell_type
cell_tissue <- sce_epi_flit@meta.data$orig.ident
cell_ids <- colnames(sce_epi_flit@assays$RNA$data)
cell_info <- data.frame(Cell_ID = cell_ids, Cell_Type = cell_types, Tissue = cell_tissue)
cell_info$cell_types <- factor(cell_types, levels = c( "Benign",'6_Mix_3', 'Malignant'))
cell_info$Tissue <- factor(cell_tissue, levels = c('Normal', 'Adenoma', 'Carcinoma'))
cell_info_ordered <- cell_info[order(cell_info$cell_types), ]

common_columns <- intersect(colnames(gsva_result), cell_info$Cell_ID)
gsva_result_filtered <- gsva_result[, common_columns, drop = FALSE]

common_rows <- intersect(rownames(gsva_result_filtered), rownames(gsva_result6))

gsva_result_filtered_common <- gsva_result_filtered[common_rows, , drop = FALSE]
gsva_result6_common <- gsva_result6[common_rows, , drop = FALSE]

gsva_result2 <- cbind(gsva_result_filtered_common , gsva_result6_common) 
ordered_cell_ids <- cell_info_ordered$Cell_ID
gsva_result2_ordered <- gsva_result2[, ordered_cell_ids]



cluster_bp <- sce_epi_flit@meta.data %>% dplyr::filter(cell_type %in% c("Benign",'6_Mix_3')) %>% dplyr::arrange(cell_type)
use_gsva_bp <- as.data.frame(gsva_result2_ordered) %>% dplyr::select(rownames(cluster_bp))
group_bp <- c(rep('Mix', 420),rep("Benign", 1492)) %>% as.factor()
desigN_bp <- model.matrix(~ 0 + group_bp) 
colnames(desigN_bp) <- levels(group_bp)
fit_bp = lmFit(use_gsva_bp, desigN_bp)
cont.matrix <- makeContrasts(contrasts = c('Benign-Mix'), levels = desigN_bp)
fit2_bp <- contrasts.fit(fit_bp, cont.matrix)
fit3_bp <- eBayes(fit2_bp)


cluster_mp <- sce_epi_flit@meta.data %>% dplyr::filter(cell_type %in% c('6_Mix_3','Malignant')) %>% dplyr::arrange(cell_type)
use_gsva_mp <- as.data.frame(gsva_result2_ordered) %>% dplyr::select(rownames(cluster_mp))
group_mp <- c(rep('Mix', 420),rep("Malignant", 1396)) %>% as.factor()
desigN_mp <- model.matrix(~ 0 + group_mp) 
colnames(desigN_mp) <- levels(group_mp)
fit_mp = lmFit(use_gsva_mp, desigN_mp)
cont.matrix <- makeContrasts(contrasts = c('Malignant-Mix'), levels = desigN_mp)
fit2_mp <- contrasts.fit(fit_mp, cont.matrix)
fit3_mp <- eBayes(fit2_mp)

summary(decideTests(fit3_bp, p.value=0.05))
diff_bp <- topTable(fit3_bp,adjust='fdr', coef=1, number=Inf)
cluster_diff_bp <- na.omit(diff_bp)
cluster_diff_bp$pathway <- rownames(cluster_diff_bp)

diff_mp <- topTable(fit3_mp,adjust='fdr', coef=1, number=Inf)
cluster_diff_mp <- na.omit(diff_mp)
cluster_diff_mp$pathway <- rownames(cluster_diff_mp)

sig_mp <- cluster_diff_mp %>% dplyr::filter(abs(logFC) > 0.3 & adj.P.Val < 0.05)
sig_bp <- cluster_diff_bp %>% dplyr::filter(abs(logFC) > 0.3 & adj.P.Val < 0.05)

common_pathways <- intersect(sig_mp$pathway, sig_bp$pathway)

filtered_mp <- sig_mp %>% filter(pathway %in% common_pathways)
filtered_bp <- sig_bp %>% filter(pathway %in% common_pathways)

merged_data <- rbind(filtered_mp, filtered_bp)


old_pathways <- c('GOBP_NEGATIVE_REGULATION_OF_SERINE_TYPE_ENDOPEPTIDASE_ACTIVITY',
                  'GOBP_NEGATIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_DISASSEMBLY',
                  'GOBP_NEGATIVE_REGULATION_OF_MACROPHAGE_APOPTOTIC_PROCESS',
                  "GOBP_AMYLIN_RECEPTOR_SIGNALING_PATHWAY",
                  "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION",
                  "KEGG_MEDICUS_VARIANT_SCRAPIE_CONFORMATION_PRPSC_TO_PERK_ATF4_SIGNALING_PATHWAY","KEGG_MEDICUS_REFERENCE_REGULATION_OF_GF_RTK_RAS_ERK_SIGNALING_PTP",
                  "GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_BY_P53_CLASS_MEDIATOR",
                  "KEGG_MEDICUS_VARIANT_LOSS_OF_RASSF1_TO_RAS_RASSF1_SIGNALING_PATHWAY",
                  'KEGG_MEDICUS_REFERENCE_HIF_1_SIGNALING_PATHWAY'
                  
)
new_pathways <- c(
  "Regulation of serine type edopeptidase activity",
  "Regulation of extracellular matrix disassembly",
  "Regulation of macrophage apoptotic processing",
  "Amylin receptor signaling pathway",
  "PERK/ATF4 pathway","RTK-Ras-ERK signaling pathway",
  "Regulation of intrinsic apoptotic signaling pathway by p53 class medeiator",
  "Loss of RASSF1 to RAS-RASSF1 signaling signaling pathway",
  "HIF 1 signaling pathway"
  
)

matched_indices <- match(old_pathways, rownames(gsva_result))
extracted_results <- gsva_result2_ordered[matched_indices, ]
rownames(extracted_results) <- new_pathways
extracted_results_ordered <- extracted_results[, ordered_cell_ids]

annotation_col <- data.frame(
  Tissue = cell_info_ordered$Tissue,CellType = cell_info_ordered$cell_types
)
rownames(annotation_col) <- ordered_cell_ids

col_annotation <- HeatmapAnnotation(
  df = annotation_col,
  col = list(
    Tissue = c("Normal" = "#4269D0FF", 
               "Adenoma" = "#EFB118FF", 
               "Carcinoma" = "#FF725cff"),
    CellType = c( "Benign" = "#00468BFF",
                 "6_Mix_3" = "#925e9fff", 
                 "Malignant" = "#FDAF91FF")
    
  ),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  annotation_legend_param = list(direction='horizontal')
)
col_annotation@anno_list[["CellType"]]@color_mapping@name <- 'Cell Type'
col_annotation@anno_list[["CellType"]]@label <- 'Cell Type'
col_annotation@anno_list[["CellType"]]@name <- 'Cell Type'

ht <- Heatmap(
  as.matrix(extracted_results_ordered),
  name = "GSVA Score",
  row_title = "Pathways",
  column_title = " ",
  top_annotation = col_annotation,
  show_row_names = TRUE,
  row_names_max_width = unit(10,'cm'),
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  use_raster = TRUE # Set to TRUE to use raster graphics
)
pdf("heat_gsava6.pdf", width = 10, height = 10)
draw(ht)
dev.off()
