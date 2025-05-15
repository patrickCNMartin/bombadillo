library(vesalius)

coord <- read.csv("data/spatial_data_test_coordinates.csv")
coord$x <- coord$x * 1000
coord$y <- coord$y * 1000

counts <- read.csv("data/spatial_data_test_counts.csv", row.names = 1)

cells <- coord$cell_labels
names(cells) <- coord$barcodes


ves <- build_vesalius_assay(coord[c("barcodes","x","y")], counts)
ves <- add_cells(ves, cells, add_name = "cell_labels")

ves <- generate_embeddings(ves, dim_reduction = "PCA")
image_plot(ves)

pcs <- vector("list", 25)
for (i in seq_along(pcs)){
    pcs[[i]] <- image_plot(ves, dimensions = i, embedding = "PCA")
}

p3 <- ggarrange(plotlist = pcs, ncol = 5, nrow = 5)

library(Seurat)

seu <- CreateSeuratObject(counts)
coord.df <- coord[,c("x","y")]
rownames(coord.df) <- coord$barcodes

seu@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )

seu <- seu %>% 
    NormalizeData() %>% 
    ScaleData() %>% 
    FindVariableFeatures(nfeatures = 2000) %>% 
    RunPCA() %>% 
    FindNeighbors() %>% 
    FindClusters(resolution = 0.1) %>% 
    RunUMAP(dims = 1:10)

g <- DimPlot(seu)
g1 <- territory_plot(ves)
new_coord <- coord
new_coord$seurat_clusters <- seu@meta.data$seurat_clusters
g2 <- ggplot(new_coord,aes(x,y,col = as.factor(seurat_clusters))) + geom_point()
g + g2 +g1