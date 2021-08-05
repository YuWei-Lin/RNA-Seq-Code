### Location address Setting
CaseN <- "GSE75688_Breast"
CaseN <- "GSE72056_Melanoma"
CaseN <- "GSE70630_OG"
CaseN <- "Astrocytoma"
CaseN <- "GSE81383_Melanoma"
CaseN <- "GSE103322_HNSCC"
CaseN <- "E_MTAB_6149_NSCLC"
CaseN <- "GSE81861_CRC"
CaseN <- "GSE76312_CML"

STAT <- "Normal"
STAT <- "Cancerous"
### Folder Routes
RouteNum1 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", STAT, "/K.Means/Filled with Zero/", sep = "")
RouteNum2 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", STAT, "/K.Means/Filled with Mean/", sep = "")

### Generate Genes and clusters membership table
ClusLi <- list.files(paste(RouteNum2, "ClusterGeneLists/", sep = ""), pattern = "ClusGenes_K")
Data <- read.table(paste(RouteNum2, "ClusterGeneLists/", ClusLi[1], sep = ""), header = T, stringsAsFactors = FALSE)
BVC <- Data$x
for (y in 1:length(ClusLi)) {
  DA <- read.table(paste(RouteNum2, "ClusterGeneLists/", ClusLi[y], sep = ""), header = T, stringsAsFactors = FALSE)
  ColRo <- Data$x%in%DA$x
  print(nrow(DA))
  print(sum(ColRo))
  BVC <- cbind(BVC, ColRo)
}
### Different Case has different hytpen labels
#colnames(BVC) <- c("Gene", unlist(strsplit(unlist(strsplit(ClusLi, "_"))[(1:234)%%3==0], "\\."))[(1:156)%%2==1])
colnames(BVC) <- c("Gene", unlist(strsplit(unlist(strsplit(ClusLi, "_"))[(1:312)%%4==0], "\\."))[(1:156)%%2==1])
#colnames(BVC) <- c("Gene", unlist(strsplit(unlist(strsplit(ClusLi, "_"))[(1:468)%%6==0], "\\."))[(1:156)%%2==1])
BVC <- BVC[-1, ]
BVC <- data.frame(BVC, stringsAsFactors = F)
write.csv(BVC, file = paste(RouteNum2, "Inputs/", CaseN, "_", STAT, "ClusGene_Membership.csv", sep = ""), row.names = F)
Clust.Num <- sapply(BVC[,2:ncol(BVC)], function(x) sum(as.logical(x)))
##########################################################################################

### CluRatio calculation and the types (dvd. by parental cluster gene size OR by child cluster gene size)
### Parent cluster hierarchy tree (Compute parent CluRatio beforehand as input)
# Convert upper triangle part of Parent CluRatio into Parent linear vector
## Overlap between clusters at adjacent levels of K.(overlaping sizes are divided by parent size)
BVC <- read.csv(paste(RouteNum2, "Inputs/", CaseN, "_", STAT, "ClusGene_Membership.csv", sep = ""), header = T, stringsAsFactors = F)
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
### Cluster Trees of parent and children view
# 18 colors for marking different clusters
thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
             "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f", "#000000", "#076f25", "#93cd7f", 
             "#4d0776", "#ffffff")

# Clusters Coordinates
X <- NULL
Y <- NULL
for (k in 1:12) {
  for (p in 1:k) {
    X <- append(X, 39/(k+1)*p, after = length(X))
    Y <- append(Y, 396-(k-1)*44, after = length(Y))
  }
}
# ClustLines: Clusters are within the same linage and  labeled  as same color
op <- list.files(paste(RouteNum2, "GeneColorsLists", sep = ""), pattern = "CPP_K")
Clust.Col <- NULL
for (i in 1:12) {
  dt <- read.csv(paste(RouteNum2, "GeneColorsLists/", op[i], sep = ""), header = T, stringsAsFactors = F)
  for (p in 1:(i+1)) {
    Clust.Col <- append(Clust.Col, unique(dt[ dt$Group == p, 3]), after = length(Clust.Col)) 
  }
}
ClustLine.D <- cbind(colnames(CluRatio), Clust.Col)
ClustLine.O <- rev(sort(table(Clust.Col)))
plot( 1:16, col = names(rev(sort(table(Clust.Col)))), pch = 16, cex = 5)
ClustLine.L <- NULL
for (t in 1:length(table(Clust.Col))) {
  ClustLine.L[[t]] <- ClustLine.D[ ClustLine.D[,2] == names(ClustLine.O)[t], 1]
}
ClustLine.L <- ClustLine.L[1:length(ClustLine.L)-1]
Legend.Colorder <- NULL
for (e in 1:length(ClustLine.L)) {
  names(ClustLine.L)[e] <- paste("ClustLine_", e, "(", length(ClustLine.L[[e]]), ")", sep = "")
  AA <- ClustLine.D[ClustLine.D[ ,1] == ClustLine.L[[e]][1], 2]
  Legend.Colorder <- c(Legend.Colorder, AA)
}

Global.ratio <- NULL
for (t in 1:11) {
  if(t == 1){
    GG <- CluRatio[t, (t+1):((t+1)+t)]
  }else{
    GG <- CluRatio[(sum(1:(t-1))+1):(sum(1:(t-1))+t) , (sum(1:(t-1))+t+1):(sum(1:(t-1))+t+1+t)]
  }
  Global.ratio <- append(Global.ratio, unlist(c(t(GG))), after = length(Global.ratio))
}
hist(Global.ratio)
tiff(paste(RouteNum2, "K_Maps/", CaseN,"_Whole_Overlapping_ParentRatio.tiff", sep = ""), width=1600, height=1600, compression="lzw", res=300)
hist(Global.ratio)
arrows(0.2, 250, 0.4, 350, lwd = 0.5)
arrows(0.4, 250, 0.6, 350, lwd = 2.5)
arrows(0.6, 250, 0.8, 350, lwd = 5.0)
arrows(0.8, 250, 1.0, 350, lwd = 8.0)
dev.off()

tiff(paste(RouteNum2, "K_Maps/", CaseN,"_General_Tree_ParentView.tiff", sep = ""), width=2400, height=3200, compression="lzw", res=300)
plot(NULL, type="n", xlab="", ylab="", xlim=c(0, 40), ylim=c(-80, 400), main = CaseN, axes = F)
#1 For generating cumulative numbers of knots and arrows ammount within certain K
CX <- NULL
MP <- NULL
for (k in 1:12) {
  CX <- c(CX, sum(1:k))
  MP <- c(MP, k*(k+1))
}
#2 Starting from first cluster to the "k-1" knots; Set k=1 as a special case to deal with.
for (b in 1:66) {
  rect(X[b]-1, Y[b]-5, X[b]+1, Y[b]+5, col = Clust.Col[b])
  k = which(c(T, b>CX)&(b<=CX))
  if(k==1){
    Kra <- Global.ratio[1:2]
  }
  else{
    Kra <- Global.ratio[(sum(MP[1:(k-1)])+1):sum(MP[1:k])]
    Kra <- split(Kra,as.numeric(gl(length(Kra),k+1,length(Kra)))) 
    # unlist(Kra[b-CX[k-1]]) serves as index to extract the coordinates set corresponding to that particular knot.
  }
  if(b==1){
    for (j in 1:2) {
      if(Kra[j]>=0.8&Kra[j]<=1.0){
        arrows(X[b], Y[b]-8, X[sum(1:k)+j], Y[sum(1:k)+j]+8, lwd = 8.0, length = 0.1)
      }
      if(Kra[j]>=0.6&Kra[j]<0.8){
        arrows(X[b], Y[b]-8, X[sum(1:k)+j], Y[sum(1:k)+j]+8, lwd = 5.0, length = 0.1)
      }
      if(Kra[j]>=0.4&Kra[j]<0.6){
        arrows(X[b], Y[b]-8, X[sum(1:k)+j], Y[sum(1:k)+j]+8, lwd = 2.5, length = 0.1)
      }
      if(Kra[j]>=0.2&Kra[j]<=0.4){
        arrows(X[b], Y[b]-8, X[sum(1:k)+j], Y[sum(1:k)+j]+8, lwd = 0.5, length = 0.1)
      }
    }
  }
  else{
    for (t in 1:(k+1)) {
      if(unlist(Kra[b-CX[k-1]])[t]>=0.8&unlist(Kra[b-CX[k-1]])[t]<=1.0){
        arrows(X[b], Y[b]-8, X[sum(1:k)+t], Y[sum(1:k)+t]+8, lwd = 8.0, length = 0.1)
      }
      if(unlist(Kra[b-CX[k-1]])[t]>=0.6&unlist(Kra[b-CX[k-1]])[t]<0.8){
        arrows(X[b], Y[b]-8, X[sum(1:k)+t], Y[sum(1:k)+t]+8, lwd = 5.0, length = 0.1)
      }
      if(unlist(Kra[b-CX[k-1]])[t]>=0.4&unlist(Kra[b-CX[k-1]])[t]<0.6){
        arrows(X[b], Y[b]-8, X[sum(1:k)+t], Y[sum(1:k)+t]+8, lwd = 2.5, length = 0.1)
      }
      if(unlist(Kra[b-CX[k-1]])[t]>=0.2&unlist(Kra[b-CX[k-1]])[t]<0.4){
        arrows(X[b], Y[b]-8, X[sum(1:k)+t], Y[sum(1:k)+t]+8, lwd = 0.5, length = 0.1)
      }
    }
  }
}
#3 Adding the last layer of K to close up the tree
for (b in 67:78) {
  rect(X[b]-1, Y[b]-5, X[b]+1, Y[b]+5, col = Clust.Col[b])
}
for (h in 1:12) {
  text( 0.1, as.numeric(rev(names(table(Y))))[h], labels = paste("k= ", h, sep = ""))
}
for (b in 1:78) {
  text(X[b], Y[b], labels = Clust.Num[b], cex = 0.5)
}
legend("topright",legend=names(ClustLine.L), fill=Legend.Colorder, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()
#################################################################################

## Overlap between clusters at adjacent levels of K.(overlaping sizes are divided by child size)
### Child cluster hierarchy tree (Compute Child CluRatio beforehand as input)
# Convert upper triangle part of Child CluRatio into Child linear vector
# ClustLines: Clusters are within the same linage and  labeled  as same color
BVC <- read.csv(paste(RouteNum2, "Inputs/", CaseN, "_", STAT, "ClusGene_Membership.csv", sep = ""), header = T, stringsAsFactors = F)
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
op <- list.files(paste(RouteNum2, "GeneColorsLists", sep = ""), pattern = "CPP_K")
Clust.Col <- NULL
for (i in 1:12) {
  dt <- read.csv(paste(RouteNum2, "GeneColorsLists/", op[i], sep = ""), header = T, stringsAsFactors = F)
  for (p in 1:(i+1)) {
    Clust.Col <- append(Clust.Col, unique(dt[ dt$Group == p, 3]), after = length(Clust.Col)) 
  }
}
ClustLine.D <- cbind(colnames(CluRatio), Clust.Col)
ClustLine.O <- rev(sort(table(Clust.Col)))
plot( 1:16, col = names(rev(sort(table(Clust.Col)))), pch = 16, cex = 5)
ClustLine.L <- NULL
for (t in 1:length(table(Clust.Col))) {
  ClustLine.L[[t]] <- ClustLine.D[ ClustLine.D[,2] == names(ClustLine.O)[t], 1]
}
ClustLine.L <- ClustLine.L[1:length(ClustLine.L)-1]
Legend.Colorder <- NULL
for (e in 1:length(ClustLine.L)) {
  names(ClustLine.L)[e] <- paste("ClustLine_", e, "(", length(ClustLine.L[[e]]), ")", sep = "")
  AA <- ClustLine.D[ClustLine.D[ ,1] == ClustLine.L[[e]][1], 2]
  Legend.Colorder <- c(Legend.Colorder, AA)
}

Global.ratio <- NULL
for (t in 1:11) {
  if(t == 1){
    GG <- CluRatio[t, (t+1):((t+1)+t)]
  }else{
    GG <- CluRatio[(sum(1:(t-1))+1):(sum(1:(t-1))+t) , (sum(1:(t-1))+t+1):(sum(1:(t-1))+t+1+t)]
  }
  Global.ratio <- append(Global.ratio, unlist(c(t(GG))), after = length(Global.ratio))
}
hist(Global.ratio)
tiff(paste(RouteNum2, "K_Maps/", CaseN,"_Whole_Overlapping_ChildRatio.tiff", sep = ""), width=1600, height=1600, compression="lzw", res=300)
hist(Global.ratio)
arrows(0.2, 250, 0.4, 350, lwd = 0.5)
arrows(0.4, 250, 0.6, 350, lwd = 2.5)
arrows(0.6, 250, 0.8, 350, lwd = 5.0)
arrows(0.8, 250, 1.0, 350, lwd = 8.0)
dev.off()

tiff(paste(RouteNum2, "K_Maps/", CaseN,"_General_Tree_ChildView.tiff", sep = ""), width=2400, height=3200, compression="lzw", res=300)
plot(NULL, type="n", xlab="", ylab="", xlim=c(0, 40), ylim=c(-80, 400), main = CaseN, axes = F)
#1 For generating cumulative numbers of knots and arrows ammount within certain K
CX <- NULL
MP <- NULL
for (k in 1:12) {
  CX <- c(CX, sum(1:k))
  MP <- c(MP, k*(k+1))
}
#2 Starting from first cluster to the "k-1" knots; Set k=1 as a special case to deal with.
for (b in 1:66) {
  rect(X[b]-1, Y[b]-5, X[b]+1, Y[b]+5, col = Clust.Col[b])
  k = which(c(T, b>CX)&(b<=CX))
  if(k==1){
    Kra <- Global.ratio[1:2]
  }
  else{
    Kra <- Global.ratio[(sum(MP[1:(k-1)])+1):sum(MP[1:k])]
    Kra <- split(Kra,as.numeric(gl(length(Kra),k+1,length(Kra)))) 
    # unlist(Kra[b-CX[k-1]]) serves as index to extract the coordinates set corresponding to that particular knot.
  }
  if(b==1){
    for (j in 1:2) {
      if(Kra[j]>=0.8&Kra[j]<=1.0){
        arrows(X[b], Y[b]-8, X[sum(1:k)+j], Y[sum(1:k)+j]+8, lwd = 8.0, length = 0.1)
      }
      if(Kra[j]>=0.6&Kra[j]<0.8){
        arrows(X[b], Y[b]-8, X[sum(1:k)+j], Y[sum(1:k)+j]+8, lwd = 5.0, length = 0.1)
      }
      if(Kra[j]>=0.4&Kra[j]<0.6){
        arrows(X[b], Y[b]-8, X[sum(1:k)+j], Y[sum(1:k)+j]+8, lwd = 2.5, length = 0.1)
      }
      if(Kra[j]>=0.2&Kra[j]<=0.4){
        arrows(X[b], Y[b]-8, X[sum(1:k)+j], Y[sum(1:k)+j]+8, lwd = 0.5, length = 0.1)
      }
    }
  }
  else{
    for (t in 1:(k+1)) {
      if(unlist(Kra[b-CX[k-1]])[t]>=0.8&unlist(Kra[b-CX[k-1]])[t]<=1.0){
        arrows(X[b], Y[b]-8, X[sum(1:k)+t], Y[sum(1:k)+t]+8, lwd = 8.0, length = 0.1)
      }
      if(unlist(Kra[b-CX[k-1]])[t]>=0.6&unlist(Kra[b-CX[k-1]])[t]<0.8){
        arrows(X[b], Y[b]-8, X[sum(1:k)+t], Y[sum(1:k)+t]+8, lwd = 5.0, length = 0.1)
      }
      if(unlist(Kra[b-CX[k-1]])[t]>=0.4&unlist(Kra[b-CX[k-1]])[t]<0.6){
        arrows(X[b], Y[b]-8, X[sum(1:k)+t], Y[sum(1:k)+t]+8, lwd = 2.5, length = 0.1)
      }
      if(unlist(Kra[b-CX[k-1]])[t]>=0.2&unlist(Kra[b-CX[k-1]])[t]<0.4){
        arrows(X[b], Y[b]-8, X[sum(1:k)+t], Y[sum(1:k)+t]+8, lwd = 0.5, length = 0.1)
      }
    }
  }
}
#3 Adding the last layer of K to close up the tree
for (b in 67:78) {
  rect(X[b]-1, Y[b]-5, X[b]+1, Y[b]+5, col = Clust.Col[b])
}
for (h in 1:12) {
  text( 0.1, as.numeric(rev(names(table(Y))))[h], labels = paste("k= ", h, sep = ""))
}
for (b in 1:78) {
  text(X[b], Y[b], labels = Clust.Num[b], cex = 0.5)
}
legend("topright",legend=names(ClustLine.L), fill=Legend.Colorder, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()
#################################################################################

### FDR score computation session
library(ggplot2)
library(cluster)
library(qvalue)
# Import GO and processed ssc-RNAseq data
GOIDList <- read.csv("C:/Users/User/Desktop/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
ClusLi <- list.files(paste(RouteNum2, "ClusterGeneLists/", sep = ""), pattern = "ClusGenes_K")
# Pre-screening and generating of GOnames
N <- table(GOIDList$ProcessName)
UG <- length(unique(GOIDList$Gene)) # The number of Universal Genes = Toatal balls!

p_time <- proc.time()
for(k in 1:length(ClusLi)){
  Data <- read.table(paste(RouteNum2, "ClusterGeneLists/", ClusLi[k], sep = ""), header = FALSE, stringsAsFactors = FALSE)
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
    write.table(GeEnR, file = paste(RouteNum2, "Inputs/", strsplit(ClusLi[k], "\\.")[[1]][1], ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  else{
    GeEnR <- GeEnR[order(GeEnR[, 5]), ]
    write.table(GeEnR, file = paste(RouteNum2, "Inputs/", strsplit(ClusLi[k], "\\.")[[1]][1], ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  print(k)
}
t_time <- proc.time()-p_time
print(t_time)
# Create value interval check function
in_interval <- function(x, interval){ 
  for (o in 1:(length(interval)-1)){
    if((interval[o]+1<=x & x<=interval[o+1]) == TRUE){
      return(o)
    }
  }
}
### Summary FDR Table for individual K
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
    Data <- read.csv(paste(RouteNum2, "Inputs/", ClusLi[ID[y]], sep = ""), header = TRUE, stringsAsFactors = FALSE)
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
  write.table(WIN, file = paste(RouteNum2, "Inputs/", CaseN, "_ClusGenes_K", k, "_Summary.csv", sep = ""), sep = ",", row.names = FALSE)
}
t_time <- proc.time()-p_time
print(t_time)
# Type II Layout - Concise GOTerm Table across K
Tracing <- list.files(paste(RouteNum2, "Inputs", sep = ""), pattern = "Summary")
#Tracing <- Tracing[order(as.numeric(sub("K", "", unlist(strsplit(Tracing, "_"))[(1:60)%%5==4])))]
ALLTe <- NULL
ColN <- NULL
for(k in 1:12){
  Data <- read.csv(paste(RouteNum2, "Inputs/", Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
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
    Data <- read.csv(paste(RouteNum2, "Inputs/", Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
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
write.table(WIN, file = paste(RouteNum2, "Inputs/", CaseN, "_KClusters_GOWhole.csv", sep = ""), sep = ",", row.names = FALSE)





