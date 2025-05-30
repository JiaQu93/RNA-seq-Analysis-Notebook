# Note ------------------------------------------------------------------
# Tutorial for "Comparison analysis of multiple datasets using CellChat" (https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html)

# Set environment: load the required libraries ---------------------------------
#install.packages("BiocManager")
#BiocManager::install("BiocNeighbors")
#install.packages("remotes")
#remotes::install_github("sqjin/CellChat")
library(CellChat)
library(Cairo)
library(showtext)
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
library(gridExtra)
library(grid)


setwd("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/Comparison.analysis/")
set.seed(100)


# Load cellchat object for each group and merge together -----------------------
cellchat_aged<-qread("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/cellchat_aged.qs")
cellchat_young<-qread("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/cellchat_young.qs")

object.list <- list(Young = cellchat_young, Aged = cellchat_aged)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 
#>  855 signaling genes.
#>  42940 cells. 
#> CellChat analysis of single cell RNA-seq data!
ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))


# Identify altered interactions and cell populations ---------------------------
setwd("/fs/ess/PAS2556/Bioinformatics_analysis/Aging.DRG_snRNA-seq/Analysis/CellChat/Comparison.analysis/altered.interactions")

## 1. Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))+
  theme(text = element_text(family = "Arial", size = 11))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")+
  theme(text = element_text(family = "Arial", size = 11))
compareInteractions <- gg1 + gg2
showtext_auto()  
ggsave("1_compareInteractions.pdf", compareInteractions, width = 8, height = 6)


## 2. Compare the number of interactions and interaction strength among different cell populations
###(A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
CairoPNG("2A_netVisual_diffInteraction.png", width = 12, height = 6, units = "in", res = 100)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, label.edge = F)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight", label.edge = F)
dev.off()


CairoPNG("2A_netVisual_diffInteraction_lable.png", width = 12, height = 6, units = "in", res = 100)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight", label.edge = T)
dev.off()



###(B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets

netVisual_heatmap_count <- netVisual_heatmap(cellchat)
CairoPNG("2B_netVisual_heatmap_count.png", width = 8, height = 6, units = "in", res = 100)
plot(netVisual_heatmap_count)
dev.off()

netVisual_heatmap_strength <- netVisual_heatmap(cellchat, measure = "weight")
CairoPNG("2B_netVisual_heatmap_strength.png", width = 8, height = 6, units = "in", res = 100)
plot(netVisual_heatmap_strength)
dev.off()

###(C) Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets
CairoPNG("2C_netVisual_circle.png", width = 12, height = 6, units = "in", res = 100)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

###(D) Circle plot showing the differential number of interactions or interaction strength among coarse cell types
#categorize the cell populations into 4 cell types
group.cellType <- c(
  "Satellite glia"              = "Glia",
  "Myelinating Schwann cell"    = "Myelinating Schwann cell",
  "Fibroblast"                  = "Stroma",
  "Endothelial"                 = "Stroma",
  "Pericyte/Smooth muscle"      = "Stroma",
  "Macrophage"                  = "Macrophage",
  "Lymphocyte"                  = "Lymphocyte",
  "DRG Neuron"                  = "DRG Neuron"
)

object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#show the number of interactions or interaction strength between any two cell types in each dataset.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#show the differential number or strength of interactions between any two cell types
CairoPNG("2D_diffInteraction.count_immune.neuron.png", width = 12, height = 6, units = "in", res = 100)
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
dev.off()

CairoPNG("2D_diffInteraction.strength_immune.neuron.png", width = 12, height = 6, units = "in", res = 100)
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
dev.off()



## (E) Circle plot showing the differential number or strength of interactions among specific cell types

celltypes_use <- c("Macrophage", "Lymphocyte", "DRG Neuron")
#extract the matrix
count_young <- cellchat@net$Young$count
count_aged  <- cellchat@net$Aged$count
weight_young <- cellchat@net$Young$weight
weight_aged  <- cellchat@net$Aged$weight
#Subset matrices to the 3 cell types
# Ensure ordering is consistent
celltypes_use <- intersect(celltypes_use, rownames(count_young))  # for safety
# Subset
count_young_sub  <- count_young[celltypes_use, celltypes_use]
count_aged_sub   <- count_aged[celltypes_use, celltypes_use]
weight_young_sub <- weight_young[celltypes_use, celltypes_use]
weight_aged_sub  <- weight_aged[celltypes_use, celltypes_use]

#Calculate the difference 
diff_count  <- count_aged_sub - count_young_sub
diff_weight <- weight_aged_sub - weight_young_sub

#Plot
CairoPNG("2E_diffInteraction.count_mac.lym.neuron.png", width = 12, height = 6, units = "in", res = 100)
groupSize <- rowSums(diff_count) + colSums(diff_count)  # size based on communication involvement
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(
  diff_count,
  #vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = TRUE,
  title.name = "Differential number of interactions among mac, lym and neuron"
)
dev.off()

CairoPNG("2E_diffInteraction.strength_mac.lym.neuron.png", width = 12, height = 6, units = "in", res = 100)
groupSize <- rowSums(diff_weight) + colSums(diff_weight)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(
  diff_weight,
  #vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = TRUE,
  title.name = "Differential strength of interactions among mac, lym and neuron"
)
dev.off()



## (F) Circle plot showing the differential number or strength of interactions sending from all sender cell types

#extract the matrix
count_young <- cellchat@net$Young$count
count_aged  <- cellchat@net$Aged$count
weight_young <- cellchat@net$Young$weight
weight_aged  <- cellchat@net$Aged$weight
#Calculate the difference 
diff_count  <- count_aged - count_young
diff_weight <- weight_aged - weight_young
#Calculate vertex weights (optional: use row + col sums)
groupSize <- rowSums(abs(diff_count)) + colSums(abs(diff_count))


### Plots for diff count 
pdf("2F_diffInteraction.count_EachSender.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), xpd = TRUE)
# Loop over each sender (row)
for (i in 1:nrow(diff_count)) {
  sender <- rownames(diff_count)[i]
  
  # Create matrix with only one sender row
  mat_single <- matrix(0, nrow = nrow(diff_count), ncol = ncol(diff_count),
                       dimnames = dimnames(diff_count))
  mat_single[i, ] <- diff_count[i, ]
  # Draw circle plot
  netVisual_circle(
    mat_single,
    #vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(abs(diff_count)),
    label.edge = TRUE,
    title.name = paste("Diff count from", sender)
  )
  # Optional: reset layout every 4 plots
  if (i %% 4 == 0) {
    par(mfrow = c(2, 2), xpd = TRUE)
  }
}
dev.off()


### Plots for diff strength 
groupSize <- rowSums(abs(diff_weight)) + colSums(abs(diff_weight))
pdf("2F_diffInteraction.strength_EachSender.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), xpd = TRUE)

#Loop over each sender cell type (row)
for (i in 1:nrow(diff_weight)) {
  sender <- rownames(diff_weight)[i]
  
  # Create matrix showing only outgoing edges from this sender
  mat_single <- matrix(0, nrow = nrow(diff_weight), ncol = ncol(diff_weight),
                       dimnames = dimnames(diff_weight))
  mat_single[i, ] <- diff_weight[i, ]
  # Plot circle
  netVisual_circle(
    mat_single,
    #vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(abs(diff_weight)),
    label.edge = TRUE,
    title.name = paste("Δ Strength from", sender)
  )
  
  # Optional: reset layout every 4 plots
  if (i %% 4 == 0) {
    par(mfrow = c(2, 2), xpd = TRUE)
  }
}

dev.off()


## (G) Circle plot showing the differential number or strength of interactions sending from DRG Neuron

#extract the matrix
count_young <- cellchat@net$Young$count
count_aged  <- cellchat@net$Aged$count
weight_young <- cellchat@net$Young$weight
weight_aged  <- cellchat@net$Aged$weight
#Calculate the difference 
diff_count  <- count_aged - count_young
diff_weight <- weight_aged - weight_young
#Calculate vertex weights (optional: use row + col sums)
groupSize <- rowSums(abs(diff_count)) + colSums(abs(diff_count))


### Plots for diff count of DRG Neuron
CairoPNG("2G_diffInteraction.count_DRGNeuron.png", width = 8, height = 6, units = "in", res = 100)
par(mfrow = c(1, 1), xpd = TRUE)
netVisual_diffInteraction(
  object = cellchat,
  label.edge = TRUE,
  weight.scale = TRUE,
  measure = "count",  # or "weight"
  sources.use = "DRG Neuron",
  targets.use = NULL,  # All targets
  title.name = "Diff count from DRG Neuron"
)
dev.off()


### Plots for diff count of DRG Neuron
CairoPNG("2G_diffInteraction.strength_DRGNeuron.png", width = 8, height = 6, units = "in", res = 100)
par(mfrow = c(1, 1), xpd = TRUE)
netVisual_diffInteraction(
  object = cellchat,
  label.edge = TRUE,
  weight.scale = TRUE,
  measure = "weight",
  sources.use = "DRG Neuron",
  targets.use = NULL,  # All targets
  title.name = "Diff count from DRG Neuron"
)
dev.off()




