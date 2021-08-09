### Location address Setting
CaseN <- "GSE75688_Breast"
CaseN <- "GSE72056_Melanoma"
CaseN <- "GSE70630_OG"
CaseN <- "Astrocytoma"
CaseN <- "GSE81383_Melanoma"
CaseN <- "GSE103322_HNSCC"
### 2. Clusters at each K table (Use txt gene lists generated from ConsensusClusterPlus package)
ClusLi <- list.files(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/ClusterGeneLists", sep = ""), pattern = "ClusGenes_K")
for (i in 1:length(ClusLi)) {
  Data <- read.table(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/ClusterGeneLists/", ClusLi[i], sep = ""), header = FALSE, stringsAsFactors = FALSE)
  if(nchar(unlist(strsplit(ClusLi[i], "_"))[4]) == 2){
    A <- sub( "K", "K0", unlist(strsplit(ClusLi[i], "_"))[4])
  }else{
    A <- unlist(strsplit(ClusLi[i], "_"))[4]
  }
  if(nchar(unlist(strsplit(ClusLi[i], "_"))[5]) == 6){
    B <- sub( "C", "C0", unlist(strsplit(ClusLi[i], "_"))[5])
  }else{
    B <- unlist(strsplit(ClusLi[i], "_"))[5]
  }
  B <- sub("txt","csv", B)
  write.table( Data, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/ClusterGeneLists/", paste(CaseN, "_ClusGenes_", A, B, sep = ""), sep = ""), row.names = F)
}
# Go back to delete txt version of genelists, then reload following code
ClusLi <- list.files(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/ClusterGeneLists", sep = ""), pattern = "ClusGenes_K")
Data <- read.table(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/ClusterGeneLists/", ClusLi[1], sep = ""), header = T, stringsAsFactors = FALSE)
BVC <- Data$V1
for (y in 1:length(ClusLi)) {
  DA <- read.table(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/ClusterGeneLists/", ClusLi[y], sep = ""), header = T, stringsAsFactors = FALSE)
  ColRo <- Data$V1%in%DA$V1
  print(nrow(DA))
  print(sum(ColRo))
  BVC <- cbind(BVC, ColRo)
}
#colnames(BVC) <- c("Gene", unlist(strsplit(unlist(strsplit(ClusLi, "_"))[(1:234)%%3==0], "\\."))[(1:156)%%2==1])
colnames(BVC) <- c("Gene", unlist(strsplit(unlist(strsplit(ClusLi, "_"))[(1:312)%%4==0], "\\."))[(1:156)%%2==1])
BVC <- BVC[-1, ]
write.csv(BVC, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN,"ClusGene_Membership.csv", sep = ""), row.names = F)
Clust.Num <- sapply(BVC[,2:ncol(BVC)], function(x) sum(x))
### 3.  Overlap between clusters at adjacent levels of K.(Only overlaping number of size)
BVC <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN,"ClusGene_Membership.csv", sep = ""), header = T, stringsAsFactors = F)
CluRatio <- as.data.frame(matrix(NA, nrow = 78, ncol = 78))
rownames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
colnames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
for (k in 1:11) {
  Siz <- k+(k+1)
  if(k==1){
    for (m in 1:3) {
      for (n in 1:3) {
        TOP <- BVC$Gene[BVC[,m+1] == TRUE]
        BOT <- BVC$Gene[BVC[,n+1] == TRUE]
        CluRatio[m, n] <- length(intersect(TOP, BOT))
      }
    }
  }else{
    for (m in (sum(1:(k-1))+1):((sum(1:(k-1))+1)+Siz-1)) {
      for (n in (sum(1:(k-1))+1):((sum(1:(k-1))+1)+Siz-1)) {
        TOP <- BVC$Gene[BVC[,m+1] == TRUE]
        BOT <- BVC$Gene[BVC[,n+1] == TRUE]
        CluRatio[m, n] <- length(intersect(TOP, BOT))
      }
    }
  }
}
write.csv(CluRatio, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN,"ClusGene_Overlap(Size).csv", sep = ""), row.names = T)

### 4.  Overlap between clusters at adjacent levels of K.(overlaping sizes are divided by parent size)
BVC <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN,"ClusGene_Membership.csv", sep = ""), header = T, stringsAsFactors = F)
CluRatio <- as.data.frame(matrix(NA, nrow = 78, ncol = 78))
rownames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
colnames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
for (k in 1:11) {
  Siz <- k+(k+1)
  if(k==1){
    for (m in 1:3) {
      for (n in 1:3) {
        TOP <- BVC$Gene[BVC[,m+1] == TRUE]
        BOT <- BVC$Gene[BVC[,n+1] == TRUE]
        CluRatio[m, n] <- length(intersect(TOP, BOT))/length(TOP)
      }
    }
  }else{
    for (m in (sum(1:(k-1))+1):((sum(1:(k-1))+1)+Siz-1)) {
      for (n in (sum(1:(k-1))+1):((sum(1:(k-1))+1)+Siz-1)) {
        TOP <- BVC$Gene[BVC[,m+1] == TRUE]
        BOT <- BVC$Gene[BVC[,n+1] == TRUE]
        CluRatio[m, n] <- length(intersect(TOP, BOT))/length(TOP)
      }
    }
  }
}
CluRatio[lower.tri(CluRatio, diag = T)] <- NA
write.csv(CluRatio, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN,"ClusGene_Overlap(dvd Parent).csv", sep = ""), row.names = T)

### 5.  Overlap between clusters at adjacent levels of K.(overlaping sizes are divided by child size)
BVC <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/Inputs/", CaseN,"ClusGene_Membership.csv", sep = ""), header = T, stringsAsFactors = F)
CluRatio <- as.data.frame(matrix(NA, nrow = 78, ncol = 78))
rownames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
colnames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
for (k in 1:11) {
  Siz <- k+(k+1)
  if(k==1){
    for (m in 1:3) {
      for (n in 1:3) {
        TOP <- BVC$Gene[BVC[,m+1] == TRUE]
        BOT <- BVC$Gene[BVC[,n+1] == TRUE]
        CluRatio[m, n] <- length(intersect(TOP, BOT))/length(BOT)
      }
    }
  }else{
    for (m in (sum(1:(k-1))+1):((sum(1:(k-1))+1)+Siz-1)) {
      for (n in (sum(1:(k-1))+1):((sum(1:(k-1))+1)+Siz-1)) {
        TOP <- BVC$Gene[BVC[,m+1] == TRUE]
        BOT <- BVC$Gene[BVC[,n+1] == TRUE]
        CluRatio[m, n] <- length(intersect(TOP, BOT))/length(BOT)
      }
    }
  }
}
CluRatio[lower.tri(CluRatio, diag = T)] <- NA
write.csv(CluRatio, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN,"ClusGene_Overlap(dvd Child).csv", sep = ""), row.names = T)

### 6.1.1 Use Overlap(dvd Child) matrix as a start in this section to generate the cluster hierarchy matrix based on Child view

write.csv(CluRatio, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/TABLES/cluster hierarchy(Child).csv", sep = ""), row.names = T)
### 6.1.2 cluster hierarchy(Parent)


### 6.2.1 Child cluster hierarchy tree
### 6.2.2 Parent cluster hierarchy tree

