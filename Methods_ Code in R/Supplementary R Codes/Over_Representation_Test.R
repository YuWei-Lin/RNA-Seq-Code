
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
# Import Packages
library(ggplot2)
library(cluster)
library(qvalue)
# Import GO and processed ssc-RNAseq data
GOIDList <- read.csv("C:/Users/User/Desktop/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)

ClusLi <- list.files("D:/GSE75688_Breast_CCP_kClusters/BC_Version2", pattern = "ClusGenes_K")

# Pre-screening and generating of GOnames
N <- table(GOIDList$ProcessName)
UG <- length(unique(GOIDList$Gene)) # The number of Universal Genes = Toatal balls!

p_time <- proc.time()
for(k in 1:length(ClusLi)){
  Data <- read.table(paste("D:/GSE75688_Breast_CCP_kClusters/BC_Version2/", ClusLi[k], sep = ""), header = FALSE, stringsAsFactors = FALSE)
  DWs <- sum(Data$V1%in%unique(GOIDList$Gene)) # Actual Genes in comparison = Draws!
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
  #fdr2 <- p.adjust(pV, method = "BH")
  #fdr3 <- p.adjust(pV, method = "fdr")
  #fdr4 <- p.adjust(pV, method = "bonferroni")
  GeEnR <- cbind(pV_Name, terN, Hits, pV, fdr$qvalues)
  colnames(GeEnR) <- c("GO_Term", "Gene_Number", "Overlapping", "p-Value", "FDR")
  #GeEnR <- GeEnR[GeEnR[,5]<0.05, ]
  if(is.null(dim(GeEnR))==TRUE){
    write.table(GeEnR, file = paste(strsplit(ClusLi[k], "\\.")[[1]][1], ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  else{
    GeEnR <- GeEnR[order(GeEnR[, 5]), ]
    write.table(GeEnR, file = paste(strsplit(ClusLi[k], "\\.")[[1]][1], ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  print(k)
}
t_time <- proc.time()-p_time
print(t_time)
### For hallmark gene sets session
install.packages("GSA")
library(GSA)
Data <- GSA.read.gmt("C:/Users/User/Documents/h.all.v6.2.symbols.gmt")
FIN <- data.frame()
SIN <- data.frame()
for ( p in 1:50) {
  SIN <- cbind(Data$genesets[[p]], rep(Data$geneset.names[p], length(Data$genesets[[p]])))
  FIN <- rbind(FIN, SIN)
}
colnames(FIN) <- c("Gene", "ProcessName")
write.csv(FIN , "GSA_HALLMark_GENES_Category.csv",row.names = FALSE)
GOIDList <- read.csv("C:/Users/User/Documents/GSA_HALLMark_GENES_Category.csv", header = TRUE, stringsAsFactors = FALSE)
ClusLi <- list.files("D:/GSE103322_HNSCC_CCP_kClusters/", pattern = "ClusGenes_K")
# Pre-screening and generating of GOnames
N <- table(GOIDList$ProcessName)
UG <- length(unique(GOIDList$Gene)) # The number of Universal Genes = Toatal balls!
GeEnR <- NULL
p_time <- proc.time()
for(k in 1:length(ClusLi)){
  Data <- read.table(paste("D:/GSE103322_HNSCC_CCP_kClusters/", ClusLi[k], sep = ""), header = FALSE, stringsAsFactors = FALSE)
  DWs <- sum(Data$V1%in%unique(GOIDList$Gene)) # Actual Genes in comparison = Draws!
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
  colnames(GeEnR) <- c("HallMark_Term", "Gene_Number", "Overlapping", "p-Value", "FDR")
  GeEnR <- as.data.frame(GeEnR, stringsAsFactors = F)
  GeEnR <- GeEnR[as.numeric(GeEnR$FDR)<0.05, ]
  if(is.null(dim(GeEnR))==TRUE){
    write.table(GeEnR, file = paste(sub("ClusGenes", "HALLMARKs", strsplit(ClusLi[k], "\\.")[[1]][1]), ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  else{
    GeEnR <- GeEnR[order(GeEnR[, 5]), ]
    write.table(GeEnR, file = paste(sub("ClusGenes", "HALLMARKs", strsplit(ClusLi[k], "\\.")[[1]][1]), ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  print(k)
}
t_time <- proc.time()-p_time
print(t_time)
# Find a position of term in a list
match("GO_DNA_REPAIR", names(N))
length(ClusLi)
length(N)
for(i in 1:length(N)){
  if(length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand))) >= 50 & 
     length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand))) <= 500){
    AD <- Mel.stand[ ,colnames(Mel.stand)%in%GOIDList[which(GOIDList$GOID==(names(N)[i])),1]]
    GNlist <- append(GNlist, names(N)[i], after = length(GNlist))
    write.csv(AD, file = paste("C:/Users/tonyxp/Desktop/BreastCan_ZRMCDF/GSE75688_Mig_GOid/GO_", substr(names(N)[i],4,nchar(as.character(names(N)[i]))), sep="",".csv"))
  }
}
write.csv(GNlist , "GNlistGSE7568_ZRCDF.csv",row.names = FALSE)

phyper(49-1, 1977, 17106-1977, 379, lower.tail = F)
length(unique(GOIDList$Gene))
length(GOIDList$GOID[GOIDList$GOID=="GO:0006796"])