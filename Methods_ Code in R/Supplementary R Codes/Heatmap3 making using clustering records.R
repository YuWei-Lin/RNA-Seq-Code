### GSE75688_Breast
CaseN <- "D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/Inputs/GSE75688_Breast_NoBulk_CDF.csv"
DataMem <- "D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/Inputs/GSE75688_BreastClusGene_Membership.csv"
Mel.stand <- read.csv(CaseN, header = TRUE, stringsAsFactors = F)
Mems <- read.csv(DataMem, header = T)
# Remove genes that contain zero more than 7/10
Mel.stand <- Mel.stand[apply(Mel.stand, 1, function(x) sum(is.na(x)) < (ncol(Mel.stand)*(0.7))), ]
Mel.stand <- as.matrix(Mel.stand)
Mel.stand <- t(Mel.stand)

# Colors Customization
library(cluster) #General color Idex
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe"
           , "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")
# Tumor Identity Colors labeling
labels <- strsplit(rownames(Mel.stand), "_")
labels <- unlist(lapply(labels, function(x) x[1]))
tumors_colors=rep(colors[1], nrow(Mel.stand))
for(w in 2:length(table(labels))){
  tumors_colors[labels==names(table(labels))[w]]=colors[w]
}
tumors_colors <- as.matrix(tumors_colors) # Legend Attached Form
colnames(tumors_colors)=c("Tumor_ID") # Legend Label Name 

library(RColorBrewer) 
grcol <- colorRampPalette(c("green","red"))(64) #heat maps color keys

# Duplicate input and convert NAs to "0"
AA <- Mel.stand
AA[is.na(AA)]=0 # "AA" ready for Kmeans
CC <- Mel.stand # "CC" a duplicate

### Get cluster colors and genes from previous clustering results
ColorRecord <- list.files("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/GeneColorsLists/")

p_time <- proc.time()
for(k in 2:12){
  K.Clustcolors <- read.csv(paste("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/GeneColorsLists/", ColorRecord[k], sep = ""), header = T, stringsAsFactors = F)
  BB <- AA[ , order(K.Clustcolors$Group)]
  main_title=paste("GSE75688_Breast_CCP_k=", k, sep = "")
  par(cex.main=0.5)
  K.Clustcolors <- K.Clustcolors[order(K.Clustcolors$Group), ]
  Cluster_Colors <- K.Clustcolors$Col

  plot(1:length(Cluster_Colors), col=Cluster_Colors, pch=16, cex=3)
  Cluster_Colors <- as.matrix(t(Cluster_Colors))
  rownames(Cluster_Colors)=c("Clusters")
  
  tiff(paste("D:/ALL_FINAL_RESULTS/GSE75688_Breast_DATA/K.Means Clustering/K_Heatmaps/", main_title, "NA00.tiff", sep=""), width=2200, height=1600, compression="lzw", res=300)
  heatmap.3(t(BB), na.rm = TRUE, scale="none", dendrogram="none", margins=c(6,12), RowSideColors=Cluster_Colors,
            Rowv=FALSE, Colv=FALSE, ColSideColors=tumors_colors, symbreaks=FALSE, key=TRUE, symkey=FALSE,
            density.info="none", trace="none", main=main_title, labCol=FALSE, labRow=FALSE, cexRow=1, col=grcol,
            ColSideColorsSize=1, RowSideColorsSize=1)
  
  legend("bottomleft",legend=c(paste("ConsenClus", 1:k, sep = "")), fill=unique(K.Clustcolors$Col), border=FALSE, bty="n", y.intersp = 1, cex=0.7)
  legend("topright",legend=c(names(table(labels))), fill=colors[1:12], border=FALSE, bty="n", y.intersp = 1, cex=0.7)
  dev.off()
}
t_time <- proc.time()-p_time
print(t_time)