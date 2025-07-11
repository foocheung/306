library(Seurat,lib="/hpcdata/scratch/cheungf")
library(docopt)
library(future)
source("dev.306_4_functions_cite.R")
library("dplyr")
library("matrixStats")
library('tidyverse')
library(tidyseurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(harmony)
library(dsb)
library(SingleR)
library(celldex)
library(RSpectra)

# Function to process each dataset
process_dataset <- function(config_file) {
  # Load configuration from YAML file
  config <- yaml::read_yaml(config_file)
  
  B1_US_data <- Read10X_h5(config$data_path)
  
  B1_US_SeuratObj <- CreateSeuratObject(counts = (B1_US_data$'Gene Expression'), assay = "RNA", min.feature = 20)
  B1_US_SeuratObj[["CITE"]] <- CreateAssayObject(counts = B1_US_data$'Antibody Capture'[1:24,])
  
  B1_US_SeuratObj$Lane <- gsub('.*-', "", x = colnames(B1_US_SeuratObj))
  
  B1_US_demuxbestList <- list()
  demuxbestList <- list()
  lanes <- config$lanes
  
  for (i in 1:length(lanes)) {
    demuxbestList[[i]] <- read.table(paste("DEMUXALOT", "/23_306_3_1_", lanes[[i]], "/assignments_refined.tsv", sep = ""), sep = "\t", header = TRUE, col.names = c("sBARCODE", "sCALL"))
  }
  
  names(demuxbestList) <- 1:length(lanes)
  
  for (i in 1:length(lanes)) {
    demuxbestList[[i]]$sNewBarcode <- paste(substr(demuxbestList[[i]]$sBARCODE, start = 1, stop = 17), i, sep = "")
    demuxbestList[[i]] <- demuxbestList[[i]] %>%
      mutate(sDROPLET.TYPE = ifelse(grepl("\\+", sCALL), "DBL", ifelse(grepl("AMB", sCALL), "AMB", "SNG")))
  }
  
  for (i in 1:length(lanes)) {
    length(which(colnames(B1_US_SeuratObj %>% dplyr::filter(Lane == i)) %in% demuxbestList[[i]]$sNewBarcode))
    print(setdiff(colnames(B1_US_SeuratObj %>% dplyr::filter(Lane == i)), demuxbestList[[i]]$sNewBarcode))
  }
  
  demuxbestdf <- plyr::ldply(demuxbestList, data.frame)
  rownames(demuxbestdf) <- demuxbestdf$sNewBarcode
  
  B1_US_SeuratObj <- AddMetaData(B1_US_SeuratObj, metadata = demuxbestdf[colnames(B1_US_SeuratObj),])
  
  mito.genes <- grep(pattern = "^MT-", x = rownames(B1_US_SeuratObj), value = TRUE)
  B1_US_SeuratObj <- AddMetaData(object = B1_US_SeuratObj, metadata = Matrix::colSums(B1_US_SeuratObj[mito.genes,])/Matrix::colSums(B1_US_SeuratObj), col.name = "percent.mito")
  B1_US_SeuratObj <- PercentageFeatureSet(B1_US_SeuratObj, "^RP[SL]", col.name = "percent_ribo")
  B1_US_SeuratObj <- PercentageFeatureSet(B1_US_SeuratObj, "^HB[^(P)]", col.name = "percent_hb")
  

ggsave(paste(config$output_prefix, "_pre-filtering.pdf", sep = ""), height = 10, plot = VlnPlot(B1_US_SeuratObj, c("percent.mito", "percent_ribo", "percent_hb", "nFeature_RNA", "nCount_RNA"), pt.size = 0, ncol = 1 ,group.by="Lane"))
 
  
  ## Save raw Seurat object
  saveRDS(B1_US_SeuratObj, paste(config$output_prefix, ".rds", sep = ""))
  
  B1_US_SeuratObj_neg <- subset(B1_US_SeuratObj, subset = sDROPLET.TYPE == "SNG" & nFeature_RNA < 500)
  
  B1_US_SeuratObj_pos <- B1_US_SeuratObj %>%
    dplyr::filter(nFeature_RNA > 1000 & nCount_RNA < 20000 & percent.mito < 0.30)
  
  DSB.ADT.list <- list()
  
  DefaultAssay(B1_US_SeuratObj_neg) <- "CITE"
  DefaultAssay(B1_US_SeuratObj_pos) <- "CITE"
  
  isotype.control.name.vec <- c("MouseIgG2a", "MouseIgG2b", "MouseIgG1")
  DSB.ADT.list <- DSBNormalizeProtein(cell_protein_matrix = as.matrix(B1_US_SeuratObj_pos[["CITE"]]@counts), empty_drop_matrix = as.matrix(B1_US_SeuratObj_neg[["CITE"]]@counts), denoise.counts = TRUE, isotype.control.name.vec = isotype.control.name.vec)
  
  B1_US_SeuratObj_pos <- SetAssayData(B1_US_SeuratObj_pos, assay = "CITE", new.data = DSB.ADT.list, slot = "data")
 

ggsave(paste(config$output_prefix, "_pos_post-filtering.pdf", sep = ""), height = 10, plot = VlnPlot(B1_US_SeuratObj_pos, c("percent.mito", "percent_ribo", "percent_hb", "nFeature_RNA", "nCount_RNA"), pt.size = 0, ncol = 1 ,group.by="Lane"))


  
  
  metadata <- B1_US_SeuratObj@meta.data
  
  sng_data <- metadata %>%
    group_by(Lane, sDROPLET.TYPE) %>%
    summarise(n = n())
  write.csv(sng_data, paste(config$output_prefix, "_sng_data_demuxalot.csv", sep = ""))
 


seur.list <- SplitObject(B1_US_SeuratObj_pos, split.by = "Lane")

# perform SCTransform normalization
if(length(seur.list) > 1){
seur.list <- lapply(X = seur.list, FUN = SCTransform)
seur <- merge(seur.list[[1]], seur.list[2:length(seur.list)])
}
else{
seur<-SCTransform(B1_US_SeuratObj_pos)
}
VariableFeatures(seur) <- rownames(seur[["RNA"]])
seur <- RunPCA(seur)
seur <- FindNeighbors(seur, dims = 1:30)
seur <- FindClusters(seur, resolution = 0.8, algorithm=3, verbose = FALSE)

seur$monaco_ann1 <- monaco_ann1(seur)
seur$monaco_ann2 <- monaco_ann2(seur)

seur <- RunUMAP(seur, reduction = 'pca', dims = 1:10, assay = 'RNA', 
               reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

ggsave(paste(config$output_prefix,"_UMAP.pdf" , sep=""),height =3, width=length(B1_US_SeuratObj_pos$Lane  %>% unique()) + 1 * 3, plot=scCustomize::DimPlot_scCustom(seur, split.by = "Lane",label = FALSE,num_columns=length(B1_US_SeuratObj_pos$Lane  %>% unique()) + 1, reduction = "rna.umap", group.by = "monaco_ann1", split_seurat = TRUE ))



  saveRDS(seur, paste(config$output_prefix, "_pos_monaco.rds", sep = "")) 
  saveRDS(B1_US_SeuratObj_pos, paste(config$output_prefix, "_pos.rds", sep = ""))
  saveRDS(B1_US_SeuratObj_neg, paste(config$output_prefix, "_neg.rds", sep = ""))
  saveRDS(B1_US_SeuratObj, paste(config$output_prefix, ".rds", sep = ""))
}

# Process each configuration file
process_dataset("config_4.yaml")
process_dataset("config_5.yaml")
process_dataset("config_1.yaml")
process_dataset("config_2.yaml")
process_dataset("config_3.yaml")
