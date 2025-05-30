# Note ------------------------------------------------------------------
# Tutorial for "Inference and analysis of cell-cell communication using CellChat" (https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html)


# Set environment: load the required libraries ---------------------------------
#install.packages("BiocManager")
#BiocManager::install("BiocNeighbors")
#install.packages("remotes")
#remotes::install_github("sqjin/CellChat")
library(CellChat)

library(dplyr)
library(Seurat)
#library(SeuratDisk)
library(SeuratData)
library(patchwork)
#dyn.load('/apps/hdf5/gnu/9.1/openmpi/4.0/1.12.0/lib/libhdf5_hl.so.200')
library(hdf5r)
library(ggplot2)
library(cowplot)
library(here)
library(Polychrome)
library(harmony)
library(dittoSeq)
library(RColorBrewer)
library(tidyr)
sample_color <- colorRampPalette(brewer.pal(12, "Paired"))(12)
library(R.utils)
library(qs)
library(stringr)
library("readxl")
library(hdf5r)
setwd("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/")
set.seed(100)



# Load data --------------------------------------------------------------------
GSE201586 <- qread("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Data/preprocessed.data/GSE201586.qs")
CellChatDB <- CellChatDB.human  # use the human database

aged_obj <- subset(GSE201586, subset = group == "Aged" & celltype != "Adipocyte")
young_obj <- subset(GSE201586, subset = group == "Young"& celltype != "Adipocyte")

# Run CCC for aged group -------------------------------------------------------
DefaultAssay(aged_obj) <- "RNA" # CellChat requires biologically meaningful expression values, i.e., log-normalized counts from the "RNA" or "SCT" assay
cellchat_aged <- createCellChat(object = aged_obj, group.by = "celltype")
cellchat_aged@DB <- CellChatDB
cellchat_aged <- subsetData(cellchat_aged)
cellchat_aged <- identifyOverExpressedGenes(cellchat_aged)
cellchat_aged <- identifyOverExpressedInteractions(cellchat_aged)
cellchat_aged <- computeCommunProb(cellchat_aged)
cellchat_aged <- filterCommunication(cellchat_aged)
cellchat_aged <- computeCommunProbPathway(cellchat_aged)
cellchat_aged <- aggregateNet(cellchat_aged)
qsave(cellchat_aged, file = "cellchat_aged.qs")


# Run CCC for young group ------------------------------------------------------
DefaultAssay(young_obj) <- "RNA" # CellChat requires biologically meaningful expression values, i.e., log-normalized counts from the "RNA" or "SCT" assay
cellchat_young <- createCellChat(object = young_obj, group.by = "celltype")
cellchat_young@DB <- CellChatDB
cellchat_young <- subsetData(cellchat_young)
cellchat_young <- identifyOverExpressedGenes(cellchat_young)
cellchat_young <- identifyOverExpressedInteractions(cellchat_young)
cellchat_young <- computeCommunProb(cellchat_young)
cellchat_young <- filterCommunication(cellchat_young)
cellchat_young <- computeCommunProbPathway(cellchat_young)
cellchat_young <- aggregateNet(cellchat_young)
qsave(cellchat_young, file = "cellchat_young.qs")


# Load cellchat object for each group ------------------------------------------
cellchat_aged<-qread("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/cellchat_aged.qs")
cellchat_young<-qread("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/cellchat_young.qs")


# CCC Visualization for Aged.group  -------------------------------------------
setwd("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/Aged.group")
groupSize <- as.numeric(table(cellchat_aged@idents)) # number of cells in each cell group

## 1. Calculate the aggregated ccc network for aged group
ptm <- Sys.time()
cellchat_aged <- aggregateNet(cellchat_aged)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs")) # 0.01440549

ptm <- Sys.time()
# Save to PDF (recommended)
pdf("cellchat_circleplot_aged.pdf", width = 8, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat_aged@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat_aged@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength")
dev.off()

## 2. CCC network: signaling sent from each cell group

###count
pdf("circleplots.count_each.celltype_aged.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), xpd = TRUE)

mat <- cellchat_aged@net$count  
groupSize <- as.numeric(table(cellchat_aged@idents))  # required for node sizing

for (i in 1:nrow(mat)) {
  sender <- rownames(mat)[i]
  interaction_count <- sum(mat[i, ])
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2,
                   #vertex.weight = groupSize,
                   weight.scale = TRUE,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i],
                   label.edge = TRUE)  # 👈 show count or weight between cell types
}

dev.off()


###strength
pdf("circleplots.strength_each.celltype_aged.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), xpd = TRUE)

mat <- cellchat_aged@net$weight  # or @net$count if you want raw counts
groupSize <- as.numeric(table(cellchat_aged@idents))  # required for node sizing

for (i in 1:nrow(mat)) {
  sender <- rownames(mat)[i]
  interaction_count <- sum(mat[i, ])
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2,
                   #vertex.weight = groupSize,
                   weight.scale = TRUE,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i],
                   label.edge = TRUE)  # 👈 show count or weight between cell types
}

dev.off()



# CCC Visualization for Young.group  -------------------------------------------
setwd("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/Young.group/")
groupSize <- as.numeric(table(cellchat_young@idents)) # number of cells in each cell group


## 1.Calculate the aggregated CCC network for young group ----------------------
ptm <- Sys.time()
cellchat_young <- aggregateNet(cellchat_young)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs")) # 0.01333904

ptm <- Sys.time()
# Save to PDF (recommended)
pdf("cellchat_circleplot_young.pdf", width = 8, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat_young@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat_young@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength")
dev.off()

## 2. CCC network: signaling sent from each cell group
###count
pdf("circleplots.count_each.celltype_young.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), xpd = TRUE)

mat <- cellchat_young@net$count  
groupSize <- as.numeric(table(cellchat_young@idents))  # required for node sizing

for (i in 1:nrow(mat)) {
  sender <- rownames(mat)[i]
  interaction_count <- sum(mat[i, ])
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2,
                   #vertex.weight = groupSize,
                   weight.scale = TRUE,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i],
                   label.edge = TRUE)  # 👈 show count or weight between cell types
}

dev.off()


###strength
pdf("circleplots.strength_each.celltype_young.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), xpd = TRUE)

mat <- cellchat_young@net$weight  # or @net$count if you want raw counts
groupSize <- as.numeric(table(cellchat_young@idents))  # required for node sizing

for (i in 1:nrow(mat)) {
  sender <- rownames(mat)[i]
  interaction_count <- sum(mat[i, ])
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2,
                   #vertex.weight = groupSize,
                   weight.scale = TRUE,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i],
                   label.edge = TRUE)  # 👈 show count or weight between cell types
}

dev.off()

