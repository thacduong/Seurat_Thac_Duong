library(dplyr)
library(Seurat)
library(patchwork)

# Load the aml dataset
aml.data <- Read10X(data.dir = "/single_cell_analysis_GSE254282")
# Initialize the Seurat object with the raw (non-normalized data).
aml <- CreateSeuratObject(counts = aml.data, project = "PRJNA509910", min.cells = 3, min.features = 200)
aml

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
aml[["percent.mt"]] <- PercentageFeatureSet(aml, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(aml@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(aml, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(aml, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(aml, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Normalizing the data
aml <- NormalizeData(aml, normalization.method = "LogNormalize", scale.factor = 10000)
aml <- NormalizeData(aml)

# Identification of highly variable features
aml <- FindVariableFeatures(aml, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(aml), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(aml)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scailing the data
all.genes <- rownames(aml)
aml <- ScaleData(aml, features = all.genes)

# Linear dimensional reduction
aml <- RunPCA(aml, features = VariableFeatures(object = aml))
# Examine and visualize PCA results a few different ways
print(aml[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(aml, dims = 1:2, reduction = "pca")
DimPlot(aml, reduction = "pca") + NoLegend()
DimHeatmap(aml, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(aml, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(aml)

# Cluster the cells 
aml <- FindNeighbors(aml, dims = 1:10)
aml <- FindClusters(aml, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(aml), 5)

# Non-linear dimensional reduction 
aml <- RunUMAP(aml, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(aml, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(aml, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(aml, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
aml.markers <- FindAllMarkers(aml, only.pos = TRUE)
aml.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

##### Cluster annotation #####
# Automated annotation using SingleR

# reference by SingleR
reference_HPCA <- HumanPrimaryCellAtlasData()
reference_BED <- celldex::BlueprintEncodeData()
reference_DICE <- celldex::DatabaseImmuneCellExpressionData()

singleR_result_HPCA <- SingleR(test = aml@assays$RNA$data, ref = reference_HPCA, labels = reference$label.main)
singleR_result_BED <- SingleR(test = aml@assays$RNA$data, ref = reference_BED, labels = reference$label.main)
singleR_result_DICE <- SingleR(test = aml@assays$RNA$data, ref = reference_DICE, labels = reference$label.main)
singleR_result_Hematopoietic <- SingleR(test = aml@assays$RNA$data, ref = reference_Hematopoietic, labels = reference$label.main)
aml$SingleR.labels.HPCA <- singleR_result_HPCA$labels
aml$SingleR.labels.BED <- singleR_result_BED$labels
aml$SingleR.labels.DICE <- singleR_result_BED$labels
aml$SingleR.labels.MonacoImmune <- singleR_result_MonacoImmune$labels

Idents(aml) = "SingleR.labels.DICE"
DimPlot(aml, reduction = "umap")

##### find the percentage of each cell types in each clusters #####
# generate a df with cluster and cell type information (HPCA, BED, DICE)
Idents(aml) = "seurat_clusters"
cluster_celltype_df_DICE = data.frame(
  cluster = Idents(aml), 
  cell_type = aml$SingleR.labels.DICE
)
View(cluster_celltype_df_DICE)

# generate a table of cluster vs cell type counts 
cluster_celltype_count_DICE = table(cluster_celltype_df_DICE$cluster, cluster_celltype_df_DICE$cell_type)
View(cluster_celltype_count_DICE)

# convert counts to percentages 
cluster_celltype_percentage_DICE = prop.table(cluster_celltype_count_DICE, margin = 1) * 100
View(cluster_celltype_percentage_DICE)

# convert to df for easier viewing
cluster_celltype_df_DICE = as.data.frame(cluster_celltype_percentage_DICE)
colnames(cluster_celltype_df_DICE) = c('cluster', 'cell_type', 'percentage')
View(cluster_celltype_df_DICE)

# cell type percentage group by cluster
cluster_percentage_strings_DICE <- cluster_celltype_df_DICE %>% 
  group_by(cluster) %>% 
  arrange(desc(percentage)) %>% 
  summarise(cell_type_percentage_DICE = paste0(cell_type, "(", round(percentage, 2), ")", collapse = " "))
View(cluster_percentage_strings_DICE)

# merge the cluster annotation with the markers
df = merge(aml.markers, cluster_percentage_strings, by = "cluster")
df = merge(df, cluster_percentage_strings_DICE, by = "cluster")
View(df)




