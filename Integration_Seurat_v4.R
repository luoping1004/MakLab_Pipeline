# Seurat v4 pipeline for scRNAseq data integration using harmony for batch correction

library(Seurat)
library(tidyverse)
library(harmony)

obj <- readRDS("")
output_dir <- ""

print("starting SCT")
obj <- obj %>%
  PercentageFeatureSet(., pattern = "^MT-", col.name = "percent.mt.int") %>% 
  SCTransform(.,
              vars.to.regress = c("percent.mt.int"),
              method="glmGamPoi",
              vst.flavor = "v2",
              variable.features.n = 3000,
              return.only.var.genes = F,
              verbose = T) %>%
  RunPCA(.,
         features = VariableFeatures(object = .),
         verbose = FALSE)

pdf(paste0(output_dir, "plots/","Elbow.pdf"), width = 8, height = 7)
print(ElbowPlot(obj, ndims = 50, reduction = "pca"))
dev.off()

dims = c(1:29) # select dims based on Elbow plot

obj <- obj %>%
  RunUMAP(., reduction = "pca",
          dims = dims
  ) %>%
  FindNeighbors(., reduction = "pca",
                dims = dims,
                distance.matrix = FALSE,
                k.param = 20,
                return.neighbor = F,
                compute.SNN = T,
                prune.SNN = 1/15,
                nn.method = "annoy",
                n.trees = 50,
                annoy.metric = "euclidean"
  ) %>%
  FindClusters(., resolution = 0.5)

pdff = function(out_file, ...,  width = 8, height = 7) {
  plots = list(...)
  print(paste0("saving ", length(plots), " plot(s)"))
  pdf(out_file, width, height)
  for (i in 1:length(plots)) {
    print(plots[[i]])
    print(paste0("saved plot ",i))
  }
  dev.off()
}

# plot to show marker gene expression
pdff(paste0(output_dir,"plots/","original.pdf"),
     DimPlot(obj, label = T) , DimPlot(obj, group.by = "orig.ident"),
     DimPlot(obj, group.by = "cell_type_v1"),
     FeaturePlot(obj, "percent.mt"),
     FeaturePlot(obj, c("nCount_RNA", "nCount_FB", "nCount_RNA_log", "nCount_FB_log")), 
     FeaturePlot(obj, c("nCount_RNA", "nCount_SCT", "nFeature_RNA", "nFeature_SCT")),
     FeaturePlot(obj, "CD3D"),
     FeaturePlot(obj, "CD4"),
     FeaturePlot(obj, "CD8A"),
     FeaturePlot(obj, "TRDC"),
     FeaturePlot(obj, "GZMK"),
     FeaturePlot(obj, "GZMB"),
     FeaturePlot(obj, "TCF7"),
     FeaturePlot(obj, "TRDV1"),
     FeaturePlot(obj, "TRDV2"),
     FeaturePlot(obj, "TRDV3"))

# run harmony to correct batch effect
obj <- RunHarmony(obj, group.by.vars = "orig.ident", theta = 0, assay.use = "SCT")

pdf(paste0(output_dir,"plots/","Elbow_harmony.pdf"), width = 8, height = 7)
ElbowPlot(obj, ndims = 40, reduction = "harmony")
dev.off()

dims = c(1:30)

obj <- obj %>%
  RunUMAP(., reduction = "harmony",
          #map.method = "uwot",
          dims = dims
  ) %>%
  FindNeighbors(., reduction = "harmony",
                dims = dims,
                distance.matrix = FALSE,
                k.param = 20,
                return.neighbor = F,
                compute.SNN = T,
                prune.SNN = 1/15,
                nn.method = "annoy",
                n.trees = 50,
                annoy.metric = "euclidean"
  ) %>%
  FindClusters(., resolution = 0.5)