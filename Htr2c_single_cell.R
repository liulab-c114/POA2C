library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(plyr)
library(tidydr)
#Read 10X data 
preoptic.data <- Read10X(data.dir = "/data/xiongtianze/data/preoptic/")
dim(preoptic.data)
summary(colSums(preoptic.data))
at_least_one <- apply(preoptic.data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
hist(colSums(preoptic.data),
     breaks = 100,
     main = "Expression sum per cell",
     xlab = "Sum expression")

#Create Seuratobject
preoptic <- CreateSeuratObject(counts = preoptic.data, project = "preoptic3w", min.cells = 3, min.features = 200)
rm(preoptic.data)

#data feature
preoptic[["percent.mt"]] <- PercentageFeatureSet(preoptic, pattern = "^MT-")
VlnPlot(preoptic, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(preoptic, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(preoptic, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Data normalization
preoptic <- NormalizeData(preoptic, normalization.method = "LogNormalize", scale.factor = 10000)

#Find variable features
preoptic <- FindVariableFeatures(preoptic, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(preoptic), 10)
plot1 <- VariableFeaturePlot(preoptic)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Data scaling
all.genes <- rownames(preoptic)
preoptic <- ScaleData(preoptic, features = all.genes)

#PCA
preoptic <- RunPCA(preoptic, features = VariableFeatures(object = preoptic))
DimPlot(preoptic, reduction = "pca") + NoLegend()
DimHeatmap(preoptic, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(preoptic)

# mapping csv annotation to seurat object
cell_class <- read.csv("/data/xiongtianze/data/preoptic/cell_class.csv", skip = 1)
preoptic@meta.data$cell_type <- plyr::mapvalues(x = Cells(preoptic),
                                                from = cell_class$Cell.name,
                                                to = cell_class$Cell.class..determined.from.clustering.of.all.cells.)

#subset into neurons
sub_neurons <- subset(x = preoptic, cell_type %in% c("Excitatory", "Inhibitory")) #18553cells
sub_htr2c <- subset(x = sub_neurons, Htr2c > 0)#6196 cells
rm(sub_neurons)
rm(preoptic)
FeaturePlot(object = sub_htr2c, features = c('Htr2c'))

#PCA
sub_htr2c <- RunPCA(object = sub_htr2c, features = VariableFeatures(object = sub_htr2c))
ElbowPlot(sub_htr2c, ndims = 50)

#Clustering
sub_htr2c <- FindNeighbors(object = sub_htr2c, dims = 1:10)
sub_htr2c <- FindClusters(object = sub_htr2c, resolution = 0.2)
sub_htr2c <- RunUMAP(object = sub_htr2c, dims = 1:10)

#Visualization
DimPlot(object = sub_htr2c, reduction = "umap", label = T)  

#Find clustering markers
all_markers <- FindAllMarkers(sub_htr2c, only.pos = TRUE)
write.csv(all_markers, '/data/xiongtianze/data/preoptic/cluster_markers.csv')
marker_0 <- all_markers[all_markers$cluster == 0, ]$gene[50:75]
marker_1 <- all_markers[all_markers$cluster == 1, ]$gene[1:20]
marker_2 <- all_markers[all_markers$cluster == 2, ]$gene[1:20]
marker_3 <- all_markers[all_markers$cluster == 3, ]$gene[1:20]
marker_4 <- all_markers[all_markers$cluster == 4, ]$gene[1:20]
marker_5 <- all_markers[all_markers$cluster == 5, ]$gene[1:20]
VlnPlot(sub_htr2c, features = marker_5, stack = TRUE, pt.size = 0, flip = TRUE)

#Markers for clustering and markers for Htr family
total_markers <- c('Htr2c', 'Slc17a6', 'Slc32a1', 'Gad1', 'Gad2',
                   'Ptprt', 'Ppp1r1b', 'Nrn1', 'Hmx2', 'Sox6', 'St18',
                   'Nos1', 'Avp', 'Rxfp1', 'Pdyn', 'Cck', 'Grp', 'Ptprd', 'Ang', 'Nts'
                   )
Htr_features <- c('Htr1a', 'Htr1b', 'Htr1d', 'Htr1f', 'Htr2a', 'Htr2b', 'Htr3a', 'Htr3b', 'Htr4', 'Htr5a', 'Htr6', 'Htr7')

#Add new ids
new.cluster.ids <- c('i1:Htr2c/Ptprt', 'i2:Htr2c/Ppp1r1b', 'e3:Htrc2/Nrn1', 'i4:Htr2c/Hmx2', 'i5:Htr2c/Sox6', 'i6:Htr2c/St18')
names(new.cluster.ids) <- levels(sub_htr2c)
UMAP_color <- c('#b4403e', '#256ea2', '#ea841e', '#399335', '#603c87','#fff200')
names(UMAP_color) <- new.cluster.ids
sub_htr2c <- RenameIdents(sub_htr2c, new.cluster.ids)

#Violin plot
pdf('/data/xiongtianze/data/preoptic/Vlnplot_htr2c_revised.pdf', width = 10, height = 12)
VlnPlot(object = sub_htr2c, features = total_markers, pt = 0, stack = TRUE, flip = TRUE) +
  theme(legend.text = element_text(face = 'bold', size = 14),
        axis.title.x = element_text(face = 'bold', size = 14),
        axis.title.y = element_text(face = 'bold', size = 14))
dev.off()
#UMAP plot
p <- DimPlot(sub_htr2c, cols = UMAP_color,reduction = "umap", pt.size = 0.3, alpha = 0.8) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme_dr(xlength = 0.2,
           ylength = 0.3,
           arrow = arrow(angle = 15, length = unit(1.5, 'mm'))) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(hjust = 0.015, size = 8),
        axis.title.y = element_text(vjust = 0, size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))
p
ggsave('/data/xiongtianze/data/preoptic/UMAP_htr2c_revised.pdf', p, width = 5.5, height = 3, dpi = 500)
#Dotplot
pdf('/data/xiongtianze/data/preoptic/Dotplot_htr2c_revised.pdf', width = 5, height = 3.5)
p <- DotPlot(sub_htr2c, features = Htr_features, col.min = -1.5, col.max = 1.5) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#256ea2', '#95afc0','#b4403e')) #颜色
p
dev.off()
