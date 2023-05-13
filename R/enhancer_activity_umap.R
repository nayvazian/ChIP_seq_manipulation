library(Seurat)

#https://rdrr.io/github/qenvio/dryhic/src/R/logpseudo.R
pseudo <- function(x) {
  x + min(x[x>0], na.rm = TRUE)/2
}

enhancer_df <- pseudo(read.csv(url("http://users.wenglab.org/moorej3/Yu-Project/Mouse-Enhancer-Matrix.txt"), sep = "\t", row.names=1))

enhancer_seurat <- CreateSeuratObject(counts = enhancer_df, 
                                 project = "UMAP of enhancer activity",
                                 meta.data = data.frame(CellID = colnames(enhancer_df)),
                                 assay = "enhancer",
                                 features = rownames(enhancer_df))

# Preprocess and normalize the data
#enhancer_seurat <- subset(enhancer_seurat, subset = !is.na(enhancer_seurat@assays$enhancer@counts))


enhancer_seurat <- NormalizeData(enhancer_seurat)
enhancer_seurat <- FindVariableFeatures(enhancer_seurat)
enhancer_seurat <- ScaleData(enhancer_seurat)

enhancer_seurat <- RunPCA(enhancer_seurat, features = VariableFeatures(object = enhancer_seurat), npcs = 30)

enhancer_seurat <- RunUMAP(enhancer_seurat, reduction = "pca", dims = 1:10)

UMAPPlot(enhancer_seurat, label = TRUE)




#Xseurat <- NormalizeData(Xseurat, normalization.method = "LogNormalize")

#features <- rownames(Xseurat)
#Xseurat <- ScaleData(Xseurat, features = features)


#Xseurat <- RunPCA(Xseurat, features = VariableFeatures(object = Xseurat))
#Xseurat <- RunUMAP(Xseurat, dims = 1:10)
#https://willwerscheid.github.io/scFLASH/pseudocount.html
#logpseudo <- function(x) {
#  log(x/0.5+1)
#}

#hist(xlog$C57BL.6_embryonic_facial_prominence_embryo_11.5_days)


#Xumap <- umap(Xnorm)
#Xumap
