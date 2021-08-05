### Overepresentation test of K Clusters using GOterms
### Location address Setting
CaseN <- "GSE75688_Breast"
### 7. GO terms membership (Just provide the file: "GOID2MigDB.csv")
### 8. FDR scores of all clusters
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(ggplot2)
library(cluster)
library(qvalue)
# Import GO and processed ssc-RNAseq data
GOIDList <- read.csv("C:/Users/User/Desktop/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
ClusLi <- list.files(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/ClusterGeneLists", sep = ""), pattern = "ClusGenes_K")

# Pre-screening and generating of GOnames
N <- table(GOIDList$ProcessName)
UG <- length(unique(GOIDList$Gene)) # The number of Universal Genes = Toatal balls!

p_time <- proc.time()
for(k in 1:length(ClusLi)){
  Data <- read.table(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/ClusterGeneLists/", ClusLi[k], sep = ""), header = FALSE, stringsAsFactors = FALSE)
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
    write.table(GeEnR, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", strsplit(ClusLi[k], "\\.")[[1]][1], ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  else{
    GeEnR <- GeEnR[order(GeEnR[, 5]), ]
    write.table(GeEnR, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", strsplit(ClusLi[k], "\\.")[[1]][1], ".csv", sep = ""), sep = ",", row.names = FALSE)
  }
  print(k)
}
t_time <- proc.time()-p_time
print(t_time)

ClusLi <- list.files(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs", sep = ""), pattern = "ClusGenes_K")
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
    Data <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", ClusLi[ID[y]], sep = ""), header = TRUE, stringsAsFactors = FALSE)
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
  write.table(WIN, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN, "_ClusGenes_K", k, "_Summary.csv", sep = ""), sep = ",", row.names = FALSE)
}
t_time <- proc.time()-p_time
print(t_time)
# Type II Layout - Concise GOTerm Table across K
Tracing <- list.files(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs", sep = ""), pattern = "Summary")
#Tracing <- Tracing[order(as.numeric(sub("K", "", unlist(strsplit(Tracing, "_"))[(1:60)%%5==4])))]
ALLTe <- NULL
ColN <- NULL
for(k in 1:12){
  Data <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
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
    Data <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
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
write.table(WIN, file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN, "_KClusters_GOWhole.csv", sep = ""), sep = ",", row.names = FALSE)

### 9. Child Tree maps with FDR value texted on the figures.
DATA <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN, "_KClusters_GOWhole.csv", sep = ""), header = T, stringsAsFactors = F)
for (w in 1:5917) {
  tiff(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/TABLES/ChildView FDR Tree Maps/", DATA$Enrichment.GOTerms[w],".tiff", sep = ""), width=2400, height=3200, compression="lzw", res=300)
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
  CORR <- as.numeric(DATA[w, 3:ncol(DATA)])
  mtext(text = DATA$Enrichment.GOTerms[w], side = 3)
  #4 Use round(-log10(CORR), digits = 2) if there is a vast range of small FDR values
  for (b in 1:78) {
    if(round(-log10(CORR), digits = 2)[b]>= 1.3){
      text(X[b], Y[b], labels = round(-log10(CORR), digits = 2)[b], cex = 0.5)
    }
    if(round(-log10(CORR), digits = 2)[b]< 1.3){
      text(X[b], Y[b], labels = round(-log10(CORR), digits = 2)[b], cex = 0.2)
    }
  }
  legend("topright",legend=names(ClustLine.L), fill=thisPal[1:length(ClustLine.L)], border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
  dev.off()
  print(w)
}
### 10. stumps matrix based on Child view tree map
# Childview Stumps
BrP.Parent <- NULL
BrP.Child <- NULL
for (k in 2:11) {
  REL <- CluRatio[(sum(1:(k-1))+1):(sum(1:(k-1))+k) , (sum(1:(k-1))+k+1):(sum(1:(k-1))+k+1+k)]
  BrP.Child <- append(BrP.Child, rownames(REL)[apply(REL, 1, function(x) sum(x >= 0.4)>=2) == TRUE], after = length(BrP.Child))
  BrP.Parent <- append(BrP.Parent, colnames(REL)[apply(REL, 2, function(x) sum(x >= 0.4)>=2) == TRUE], after = length(BrP.Parent))
}
stumps.L <- NULL
for (i in 1:length(ClustLine.L)) {
  pos <- sort(c(which(ClustLine.L[[i]]%in%BrP.Parent), which(ClustLine.L[[i]]%in%BrP.Child)+1))
  pos <- c(1, pos)
  pat <- rep(seq_along(pos), times=diff(c(pos, length(ClustLine.L[[i]]) + 1)))
  stumps.L <- c(stumps.L, split(ClustLine.L[[i]], pat))
}
names(stumps.L) <- 1:length(stumps.L)
### Run a function that convert list of list into data frame witholding NAs for inequvalent length of sublists
as.data.frame.list <- function(x, row.names=NULL, optional=FALSE, ...) {
  if(!all(unlist(lapply(x, class)) %in% 
          c('raw','character','complex','numeric','integer','logical'))) {
    warning('All elements of the list must be a vector.')
    NextMethod(x, row.names=row.names, optional=optional, ...)
  }
  allequal <- all(unlist(lapply(x, length)) == length(x[[1]]))
  havenames <- all(unlist(lapply(x, FUN=function(x) !is.null(names(x)))))
  if(havenames) { #All the vectors in the list have names we can use
    colnames <- unique(unlist(lapply(x, names)))
    df <- data.frame(matrix(
      unlist(lapply(x, FUN=function(x) { x[colnames] })),
      nrow=length(x), byrow=TRUE))
    names(df) <- colnames
  } else if(allequal) { #No names, but are of the same length
    df <- data.frame(matrix(unlist(x), nrow=length(x), byrow=TRUE), ...)
    hasnames <- which(unlist(lapply(x, FUN=function(x) !is.null(names(x)))))
    if(length(hasnames) > 0) { #We'll use the first element that has names
      names(df) <- names(x[[ hasnames[1] ]])
    }
  } else { #No names and different lengths, we'll make our best guess here!
    warning(paste("The length of vectors are not the same and do not ",
                  "are not named, the results may not be correct.", sep=''))
    #Find the largest
    lsizes <- unlist(lapply(x, length))
    start <- which(lsizes == max(lsizes))[1]
    df <- x[[start]]
    for(i in (1:length(x))[-start]) {
      y <- x[[i]]
      if(length(y) < length(x[[start]])) {
        y <- c(y, rep(NA, length(x[[start]]) - length(y)))
      }
      if(i < start) {
        df <- rbind(y, df)
      } else {
        df <- rbind(df, y)
      }
    }
    df <- as.data.frame(df, row.names=1:length(x))
    names(df) <- paste('Col', 1:ncol(df), sep='')
  }
  if(missing(row.names)) {
    row.names(df) <- names(x)
  } else {
    row.names(df) <- row.names
  }
  return(df)
}
###
for (i in 1:length(ClustLine.L)) {
  names(ClustLine.L)[i] <- paste("ClustLine_", i, "(", length(ClustLine.L[[i]]),")", sep = "")
}
write.table(as.data.frame(ClustLine.L), file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/TABLES/", "ClustLine_composition.csv", sep = ""), sep = ",", row.names = T)
for (i in 1:length(stumps.L)) {
  names(stumps.L)[i] <- paste("stump_", i, "(", length(stumps.L[[i]]),")", sep = "")
}
write.table(as.data.frame(stumps.L), file = paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/TABLES/", "stumps_composition.csv", sep = ""), sep = ",", row.names = T)

### 11. stump tree map from Childview
# Generate Stumps Map with new defined colors
stumps.col <- Clust.Col
for (e in 1:length(stumps.L)) {
  pos <- which(colnames(CluRatio)%in%stumps.L[[e]])
  stumps.col[pos] <- col_vector[e]
}
stumps.O <- rev(sort(table(stumps.col)))
# Drawing cluster boxes with corresponding coordinates(stumps_Tree_ChildView)
tiff(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/TABLES/stumps_Tree_ChildView.tiff", sep=""), width=2400, height=3200, compression="lzw", res=300)
plot(NULL, type="n", xlab="", ylab="", xlim=c(0, 40), ylim=c(-80, 400), main = "Breast Cancer", axes = F)
CX <- NULL
MP <- NULL
for (k in 1:12) {
  CX <- c(CX, sum(1:k))
  MP <- c(MP, k*(k+1))
}
for (b in 1:66) {
  rect(X[b]-1, Y[b]-5, X[b]+1, Y[b]+5, col = stumps.col[b])
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
    }
  }
}
#3 Adding the last layer of K to close up the tree
for (b in 67:78) {
  rect(X[b]-1, Y[b]-5, X[b]+1, Y[b]+5, col = stumps.col[b])
}
for (h in 1:12) {
  text( 0.1, as.numeric(rev(names(table(Y))))[h], labels = paste("k= ", h, sep = ""))
}
legend("topright",legend=names(stumps.L), fill=col_vector[1:length(stumps.L)], border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

### 12. stump tree map from Childview for every GOterms
DATA <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering/Inputs/", CaseN, "_KClusters_GOWhole.csv", sep = ""), header = T, stringsAsFactors = F)
for (w in 1:5917) {
  tiff(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/TABLES/C.stumps Tree Maps with FDR/", DATA$Enrichment.GOTerms[w],".tiff", sep = ""), width=2400, height=3200, compression="lzw", res=300)
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
    rect(X[b]-1, Y[b]-5, X[b]+1, Y[b]+5, col = stumps.col[b])
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
  CORR <- as.numeric(DATA[w, 3:ncol(DATA)])
  mtext(text = DATA$Enrichment.GOTerms[w], side = 3)
  #4 Use round(-log10(CORR), digits = 2) if there is a vast range of small FDR values
  for (b in 1:78) {
    if(round(-log10(CORR), digits = 2)[b]>= 1.3){
      text(X[b], Y[b], labels = round(-log10(CORR), digits = 2)[b], cex = 0.5)
    }
    if(round(-log10(CORR), digits = 2)[b]< 1.3){
      text(X[b], Y[b], labels = round(-log10(CORR), digits = 2)[b], cex = 0.2)
    }
  }
  legend("topright",legend=names(stumps.L), fill=col_vector[1:length(stumps.L)], border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
  dev.off()
  print(w)
}