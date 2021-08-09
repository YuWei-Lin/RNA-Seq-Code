### Below are process script after Running ConsensusClusterPlus Package.
### Location address and Cases Setting
CaseN <- "GSE75688_Breast"
CaseN <- "GSE72056_Melanoma"
CaseN <- "GSE70630_OG"
CaseN <- "Astrocytoma"
CaseN <- "GSE81383_Melanoma"
CaseN <- "GSE103322_HNSCC"

FileDes1 <- paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/ClusterGeneLists/", sep = "")
FileDes2 <- paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/GeneColorsLists/", sep = "")
FileRead <- paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/Inputs/", sep = "")

### Preserve all results from ConsensusClusterPlus Package
p_time <- proc.time()
for(k in 2:12){
  if(k < 10){
    k <- paste("0", k, sep = "")
  }
  main_title=paste(CaseN, "_CCP_k=", k, sep = "")
  for(y in 1:as.numeric(k)){
    k <- as.numeric(k)
    y <- as.numeric(y)
    CGset <- rownames(CC)[resultss[[k]][["consensusClass"]]==y]
    if(as.numeric(k) < 10){
      k <- as.numeric(k)
      k <- paste("0", k, sep = "")
    }
    if(as.numeric(y) < 10){
      y <- paste("0", y, sep = "")
    }
    write.table(CGset, file = paste(FileDes1, CaseN, "_ClusGenes_K", k, "C", y, ".csv", sep = ""), quote = F, row.names = F, sep = ",") 
  }
  k <- as.numeric(k)
  g1 <- resultss[[k]][[3]]
  dfg1 <- data.frame(Gene = names(g1), Group = as.numeric(g1), stringsAsFactors = FALSE)
  col1 <- unlist(resultss[[k]][[5]][1])
  dfg1$Col <- col1
  if(k < 10){
    k = paste("0",k,sep ="")
  }
  write.csv(dfg1, paste(FileDes2, main_title, ".csv", sep = ""), row.names = FALSE)
}
t_time <- proc.time()-p_time
print(t_time)

### Cluster Membership and Inheritance Matrix
ClusLi <- list.files(FileDes1, pattern = "ClusGenes_K")
Data <- read.table(paste(FileDes1, ClusLi[1], sep = ""), header = T, stringsAsFactors = FALSE)
BVC <- Data$x
for (y in 1:length(ClusLi)) {
  DA <- read.table(paste(FileDes1, ClusLi[y], sep = ""), header = T, stringsAsFactors = FALSE)
  ColRo <- Data$x%in%DA$x
  print(nrow(DA))
  print(sum(ColRo))
  BVC <- cbind(BVC, ColRo)
}
## Select one of the matching code to generate the colnames of BVC
#colnames(BVC) <- c("Gene", unlist(strsplit(unlist(strsplit(ClusLi, "_"))[(1:234)%%3==0], "\\."))[(1:156)%%2==1])
colnames(BVC) <- c("Gene", unlist(strsplit(unlist(strsplit(ClusLi, "_"))[(1:312)%%4==0], "\\."))[(1:156)%%2==1])

BVC <- BVC[-1, ]
BVC <- as.data.frame(BVC)
write.csv(BVC, file = paste(FileRead, CaseN, "_ClusGene_Membership.csv", sep = ""), row.names = F)
Clust.Num <- sapply(BVC[,2:ncol(BVC)], function(x) sum(as.logical(x)))
CluRatio <- as.data.frame(matrix(NA, nrow = 78, ncol = 78))

### Cluster Heatmaps (Run Heatmap.3 function first and get the sample color band ready!)
ColorRecord <- list.files(FileDes2)
labels <- names(table(as.numeric(meldat[1,])))
p_time <- proc.time()
for(k in 2:3){
  K.Clustcolors <- read.csv(paste(FileDes2, ColorRecord[k], sep = ""), header = T, stringsAsFactors = F)
  BB <- meldatOri[order(K.Clustcolors$Group), ]
  BB <- BB[4:nrow(BB), ]
  main_title=paste(CaseN, "_CCP_k=", k, sep = "")
  par(cex.main=0.5)
  K.Clustcolors <- K.Clustcolors[order(K.Clustcolors$Group), ]
  Cluster_Colors <- K.Clustcolors$Col
  
  plot(1:length(Cluster_Colors), col=Cluster_Colors, pch=16, cex=3)
  Cluster_Colors <- as.matrix(t(Cluster_Colors))
  rownames(Cluster_Colors)=c("Clusters")
  
  tiff(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/K_Heatmaps/", main_title, ".tiff", sep=""), width=2200, height=1600, compression="lzw", res=300)
  heatmap.3(BB, na.rm = TRUE, scale="none", dendrogram="none", margins=c(6,12), RowSideColors=Cluster_Colors,
            Rowv=FALSE, Colv=FALSE, ColSideColors=tumors_colors, symbreaks=FALSE, key=TRUE, symkey=FALSE,
            density.info="none", trace="none", main=main_title, labCol=FALSE, labRow=FALSE, cexRow=1, col=grcol,
            ColSideColorsSize=1, RowSideColorsSize=1)
  
  legend("bottomleft",legend=c(paste("ConsenClus", 1:k, sep = "")), fill=unique(K.Clustcolors$Col), border=FALSE, bty="n", y.intersp = 1, cex=0.7)
  legend("topright",legend=c(labels), fill=colors[1:length(labels)], border=FALSE, bty="n", y.intersp = 1, cex=0.7)
  dev.off()
}
t_time <- proc.time()-p_time
print(t_time)
### GOWhole Enrichment FDR Values Table
library(ggplot2)
library(cluster)
library(qvalue)
BVC <- read.csv(paste(FileRead, CaseN, "_ClusGene_Membership.csv", sep = ""), header = T, stringsAsFactors = F)
rownames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
colnames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
GOIDList <- read.csv("C:/Users/User/Desktop/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
ClusLi <- list.files(FileDes1, pattern = "ClusGenes_K")

# Pre-screening and generating of GOnames
N <- table(GOIDList$ProcessName)
UG <- length(unique(GOIDList$Gene)) # The number of Universal Genes = Toatal balls!

p_time <- proc.time()
for(k in 1:length(ClusLi)){
  Data <- read.table(paste(FileDes1, ClusLi[k], sep = ""), header = FALSE, stringsAsFactors = FALSE)
  DWs <- sum(Data$V1[-(1:2)]%in%unique(GOIDList$Gene)) # Actual Genes in comparison = Draws!
  pV <- NULL
  pV_Name <- NULL
  terN <- NULL
  Hits <- NULL
  for (i in 1:length(N)) {
    terN[i] <- as.numeric(N[i]) # The number of target Gene set
    Hits[i] <- length(intersect(GOIDList[which(GOIDList$ProcessName==(names(N)[i])), 1], Data$V1)) # How many genes are located in target Gene set
    pV_Name[i] <- names(N)[i]
    pV[i] <- phyper(Hits[i]-1, terN[i], UG-terN[i], DWs, lower.tail = FALSE)
  }
  fdr <- qvalue(pV)
  GeEnR <- cbind(pV_Name, terN, Hits, pV, fdr$qvalues)
  colnames(GeEnR) <- c("GO_Term", "Gene_Number", "Overlapping", "p-Value", "FDR")
  #GeEnR <- GeEnR[GeEnR[,5]<0.05, ]
  if(is.null(dim(GeEnR))==TRUE){
    write.table(GeEnR, file = paste(FileRead, strsplit(ClusLi[k], "\\.")[[1]][1], ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  else{
    GeEnR <- GeEnR[order(GeEnR[, 5]), ]
    write.table(GeEnR, file = paste(FileRead, strsplit(ClusLi[k], "\\.")[[1]][1], ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  print(k)
}
t_time <- proc.time()-p_time
print(t_time)

ClusLi <- list.files(FileRead, pattern = "ClusGenes_K")
# Create value interval check function
in_interval <- function(x, interval){ 
  for (o in 1:(length(interval)-1)){
    if((interval[o]+1<=x & x<=interval[o+1]) == TRUE){
      return(o)
    }
  }
}
# Type I Layout - Summary Table for individual K
p_time <- proc.time()
for(k in 1:12){
  if(k < 10){
    k = paste("0", k,sep ="")
  }
  ID <- grep(paste("K", k, sep = ""), ClusLi)
  RW <- NULL
  Data <- NULL
  Intval <- NULL
  for (y in 1:length(ID)) {
    Data <- read.csv(paste(FileRead, ClusLi[ID[y]], sep = ""), header = TRUE, stringsAsFactors = FALSE)
    RW <- rbind(RW, Data)
    Intval[y] <- nrow(RW)
  }
  Intval <- c(0, Intval)
  WIN <- data.frame(matrix(NA, nrow = length(unique(RW$GO_Term)), ncol = as.numeric(k)+1))
  GoTerms <- names(sort(table(RW$GO_Term), decreasing = T))
  WIN[ ,1] <- GoTerms
  for (i in 1:nrow(WIN)) {
    TAR <- NULL
    UU <- rep(NA, k)
    for (j in 1:nrow(RW[RW$GO_Term==GoTerms[i], ])) {
      TAR <- c(TAR, in_interval(as.numeric(rownames(RW[RW$GO_Term==GoTerms[i], ])[j]), Intval))
      UU[TAR[j]] <- RW[RW$GO_Term==GoTerms[i], ]$FDR[j]
    }
    WIN[i, 2:ncol(WIN)] <- UU
  }
  k <- as.numeric(k)
  if(k==1){
    colnames(WIN) <- c("Enrichment GOTerms", colnames(CluRatio)[1])
  }else{
    colnames(WIN) <- c("Enrichment GOTerms", colnames(CluRatio)[(sum(1:(k-1))+1):((sum(1:(k-1))+1)+k-1)])
  }
  if(k < 10){
    k = paste("0", k,sep ="")
  }
  write.table(WIN, file = paste(FileRead, CaseN, "_ClusGenes_K", k, "_Summary.csv", sep = ""), sep = ",", row.names = FALSE)
}
t_time <- proc.time()-p_time
print(t_time)
# Type II Layout - Concise GOTerm Table across K
Tracing <- list.files(FileRead, pattern = "Summary")
#Tracing <- Tracing[order(as.numeric(sub("K", "", unlist(strsplit(Tracing, "_"))[(1:60)%%5==4])))]
ALLTe <- NULL
ColN <- NULL
for(k in 1:12){
  Data <- read.csv(paste(FileRead, Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ALLTe <- c(ALLTe, Data$Enrichment.GOTerms)
  if(k==1){
    ColN <- c(ColN, colnames(CluRatio)[1])
  }else{
    ColN <- c(ColN, colnames(CluRatio)[(sum(1:(k-1))+1):((sum(1:(k-1))+1)+k-1)])
  }
}
WIN <- data.frame(matrix(NA, nrow = length(unique(ALLTe)), ncol = 79))
GoTerms <- names(sort(table(ALLTe), decreasing = T))
WIN[ ,1] <- GoTerms
FRE <- NULL
for(i in 1:length(GoTerms)){
  LM <- NULL
  for(k in 1:12){
    Data <- read.csv(paste(FileRead, Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
    if(names(sort(table(ALLTe), decreasing = T))[i]%in%Data$Enrichment.GOTerms == TRUE){
      LM <- c(LM, Data[ Data$Enrichment.GOTerms == names(sort(table(ALLTe), decreasing = T))[i], 2:(k+1)])
    }
    else{
      LM <- c(LM, rep(NA, length(2:(k+1))))
    }
  }
  FRE <- c(FRE, (78-sum(is.na(LM))))
  WIN[i, 2:ncol(WIN)] <- LM
  print(i)
}
WIN <- cbind(WIN[ , 1], FRE, WIN[ , -1])
colnames(WIN) <- c("Enrichment GOTerms", "Frequency", ColN)
WIN <- WIN[order(WIN$Frequency, decreasing = T), ]
write.table(WIN, file = paste(FileRead, CaseN, "_KClusters_GOWhole.csv", sep = ""), sep = ",", row.names = FALSE)