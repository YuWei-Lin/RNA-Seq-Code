### 18 colors for marking different clusters
thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
             "#bd18ea", #magenta
             "#2ef4ca", #aqua
             "#f4cced", #pink,
             "#f4cc03", #lightorange
             "#05188a", #navy,
             "#e5a25a", #light brown
             "#06f106", #bright green
             "#85848f", #med gray
             "#000000", #black
             "#076f25", #dark green
             "#93cd7f", #lime green
             "#4d0776", #dark purple
             "#ffffff"  #white
)
### Clusters Coordinates
X <- NULL
Y <- NULL
Clust.Names <- NULL
for (k in 1:12) {
  for (p in 1:k) {
    X <- append(X, 39/(k+1)*p, after = length(X))
    Y <- append(Y, 396-(k-1)*44, after = length(Y))
    Clust.Names <-  append(Clust.Names, paste("K", k, "C", p, sep = ""), after = length(Clust.Names))
  }
}
### Clusters' colors in tree order
op <- list.files(path = "D:/GSE75688_Breast_CCP_kClusters/BC_Version2/GSE75688_Breast_CPP", pattern = "CPP_K")
Clust.Col <- NULL
for (i in 1:12) {
  dt <- read.csv(paste("D:/GSE75688_Breast_CCP_kClusters/BC_Version2/GSE75688_Breast_CPP/", op[i], sep = ""), header = T, stringsAsFactors = F)
  for (k in 1:(i+1)) {
    Clust.Col <- append(Clust.Col, unique(dt[ dt$Group == k, 3]), after = length(Clust.Col)) 
  }
}
### Clusters changing matrix
Global.ratio <- NULL
for (k in 1:11) {
  dT <- read.csv(paste("D:/GSE75688_Breast_CCP_kClusters/BC_Version2/GSE75688_Breast_CPP/", op[k], sep = ""), header = T, stringsAsFactors = F)
  dB <- read.csv(paste("D:/GSE75688_Breast_CCP_kClusters/BC_Version2/GSE75688_Breast_CPP/", op[k+1], sep = ""), header = T, stringsAsFactors = F)
  CluRatio <- matrix(NA, nrow = k, ncol = (k+1))
  for (m in 1:k) {
    for (n in 1:(k+1)) {
      CluRatio[m, n] <- length(intersect(dB[dB$Group==n,1], dT[dT$Group==m,1]))/length(dT[dT$Group==m,1])
      Global.ratio <- append(Global.ratio, CluRatio[m, n], after = length(Global.ratio))
    }
  }
  CluRatio <- as.data.frame(CluRatio)
  rownames(CluRatio) <- Clust.Names[(sum(1:k)-k+1):sum(1:k)]
  colnames(CluRatio) <- Clust.Names[(sum(1:k)+1):sum(1:(k+1))]
  #write.table(CluRatio, file = paste("GSE75688_Breast_ClustRatio_K", k, "&K", k+1, ".csv", sep = ""), sep = ",", row.names = T)
}
hist(Global.ratio)
tiff(paste("C:/Users/User/Desktop/Clusters_Case Global Ratio.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
hist(Global.ratio)
arrows(0.2, 250, 0.4, 350, lwd = 0.5)
arrows(0.4, 250, 0.6, 350, lwd = 2.5)
arrows(0.6, 250, 0.8, 350, lwd = 5.0)
arrows(0.8, 250, 1.0, 350, lwd = 8.0)
dev.off()


### Drawing cluster boxes with corresponding coordinates
tiff(paste("C:/Users/User/Desktop/CCC.tiff", sep=""), width=2400, height=3200, compression="lzw", res=300)
plot(NULL, type="n", xlab="", ylab="", xlim=c(0, 40), ylim=c(-80, 400), main = "Breast Cancer", axes = F)
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
dev.off()
#4 Attach particular term and evaluate its cut-off
DATA <- read.csv(paste("C:/Users/User/Desktop/GSE75688_Breast_ClusGenes/BC_Summary/GSE75688_Breast_KClusters_GOWhole.csv", sep = ""), header = T, stringsAsFactors = F)
CORR <- as.numeric(DATA[DATA$Enrichment.GOTerms=="GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY", 3:ncol(DATA)])
mtext(text = "GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY", side = 1)
#5 Use round(-log10(CORR), digits = 2) if there is a vast range of small FDR values
for (b in 1:78) {
  if(round(-log10(CORR), digits = 2)[b]>= 1.3){
    text(X[b], Y[b], labels = round(-log10(CORR), digits = 2)[b], cex = 0.5)
  }
  if(round(-log10(CORR), digits = 2)[b]< 1.3){
    text(X[b], Y[b], labels = round(-log10(CORR), digits = 2)[b], cex = 0.2)
  }
}
dev.off()
### Find out ClustLines in the tree
op <- list.files(path = "D:/GSE75688_Breast_CCP_kClusters/BC_Version2/GSE75688_Breast_CPP", pattern = "_CPP_K")
Clust.Col <- NULL
for (i in 1:12) {
  dt <- read.csv(paste("D:/GSE75688_Breast_CCP_kClusters/BC_Version2/GSE75688_Breast_CPP/", op[i], sep = ""), header = T, stringsAsFactors = F)
  for (k in 1:(i+1)) {
    Clust.Col <- append(Clust.Col, unique(dt[ dt$Group == k, 3]), after = length(Clust.Col)) 
  }
}
ClustLine.D <- cbind(Clust.Names, Clust.Col)
ClustLine.O <- rev(sort(table(Clust.Col)))
plot( 1:16, col = names(rev(sort(table(Clust.Col)))), pch = 16, cex = 5)
ClustLine.L <- NULL
for (t in 1:length(table(Clust.Col))) {
  ClustLine.L[[t]] <- ClustLine.D[ ClustLine.D[,2] == names(ClustLine.O)[t], 1]
}

### Stumps
ChM <- list.files("D:/GSE75688_Breast_CCP_kClusters/BC_Version2/GSE75688_Breast_CPP", pattern = "&")
BrP.Parent <- NULL
BrP.Child <- NULL
for (k in 2:11) {
  REL <- read.csv(paste("D:/GSE75688_Breast_CCP_kClusters/BC_Version2/GSE75688_Breast_CPP/", ChM[k], sep = ""), header = T, stringsAsFactors = F)
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
lapply(stumps.L, function(x) write.table( data.frame(x), "test.csv"  , append= T, sep=','))
### Generate "A" form Data table: ClustLines vesus GO gene sets
DATA <- read.csv("C:/Users/User/Desktop/GSE75688_Breast_ClusGenes/BC_Summary/GSE75688_Breast_KClusters_GOWhole.csv", header = T, stringsAsFactors = F)
DATA[ ,3:ncol(DATA)] <- round(-log10(DATA[ ,3:ncol(DATA)]), digits = 2)
A_form <- as.data.frame(matrix(NA, nrow = nrow(DATA), ncol = length(stumps.L)), stringsAsFactors = FALSE)
rownames(A_form) <- DATA$Enrichment.GOTerms
for (i in 1:length(stumps.L)) {
  colnames(A_form)[i] <- paste("stump_", i, "(", length(stumps.L[[i]]),")", sep = "")
}
for (g in 1:nrow(DATA)) {
  for (s in 1:length(stumps.L)) {
    if(any(DATA[ g, colnames(DATA)%in%paste("FDR_", stumps.L[[s]], sep = "")] >= 1.13)){
      A_form[ g, s] <- paste("Y", "_(", sum(DATA[ g, colnames(DATA)%in%paste("FDR_", stumps.L[[s]], sep = "")]>=1.13),")", sep = "")
      if(DATA[ g, colnames(DATA)%in%paste("FDR_", stumps.L[[s]], sep = "")][1]>= 1.13){
        A_form[ g, s] <- sub("Y", "YF", A_form[ g, s])
      }
    }
    else{
      A_form[ g, s] <- "N"
    }
  }
}
FRE <- sapply(A_form, function(x) nrow(A_form)-table(x)[1])
A_form <- rbind(FRE, A_form)
rownames(A_form)[1] <- "Y Counts by Col."
FRE2 <-apply(A_form[2:nrow(A_form), 1:(ncol(A_form)-1)], 1, function(x) (ncol(A_form)-1)-table(x)[names(table(x))=="N"])
A_form <- cbind(A_form[ , 1:(ncol(A_form)-1)], c(NA, FRE2), A_form[ , ncol(A_form)])
colnames(A_form)[ncol(A_form)-1] <- "Y Counts by Row."
colnames(A_form)[ncol(A_form)] <- "stumps_23(1)"
AA_form <-A_form[!(A_form$`Y Counts by Row.` == 0), ]
AA_form[1, ] <- A_form[1, ]
rownames(AA_form)[1] <- "Y Counts by Col."
write.csv(AA_form, file = "GSE75688_Breast_KTREE_stumps.csv", row.names = T)

### Generate "B" form Data table: Ancestor clusters of the stumps vesus GO gene sets
#1 Import Ancestor-Stumps association dataframe generated manually outside of R (Using Excel)
AncSt <- read.csv("C:/Users/User/Desktop/GSE75688_Breast_KTREE_Ance&Stumps_Relation.csv", header = T, stringsAsFactors = F)
names(stumps.L) <- c(colnames(A_form)[1:(ncol(A_form)-2)], colnames(A_form)[ncol(A_form)])
dF <- NULL
for (n in 2:nrow(A_form)) {
  ID <- colnames(A_form)[grep("YF", A_form[n, ])]
  if(length(ID) != 0){
    NBB <- NULL
    for (m in 1:length(ID)) {
      BB <- NULL
      BB <- c(rownames(A_form)[n], unlist(stumps.L[names(stumps.L)==ID[m]])[1], AncSt[AncSt$X == paste(ID[m]), ])
      NBB <- rbind(NBB, BB)
    }
    dF <- rbind(dF, NBB)
    rownames(dF) <- 1:nrow(dF)
  }
}
dF <- as.data.frame(dF, stringsAsFactors = F)
#2 Get rid of the terms that both ancestors are "NA"
dF <- dF[!(is.na(dF$Ancestor_1)&is.na(dF$Ancestor_2)), ]
colnames(dF)[1:3] <- c("GOTerms", "Top_Cluster", "Stump")
CV <- NULL
for (p in 1:nrow(dF)) {
  FCV1 <- DATA[DATA$Enrichment.GOTerms == dF$GOTerms[p], sub("FDR_", "", colnames(DATA))%in%unlist(dF[p, c(2,4,5)])]
  if(ncol(FCV1)==3){
    FCV1 <- cbind(FCV1[ ,3], FCV1[ , 1:2])
  }
  if(ncol(FCV1)==2){
    FCV1 <- FCV1[ , rev(1:ncol(FCV1))]
    FCV1 <- cbind(FCV1, as.data.frame(matrix(NA, nrow = 1, ncol = 1)))
  }
  colnames(FCV1) <- c("Top_Cluster_FDR", "Ancestor_Val1", "Ancestor_Val2")
  CV <- rbind(CV, FCV1)
}
rownames(CV) <- 1:nrow(CV)
Difference_1 <- CV$Top_Cluster_FDR - CV$Ancestor_Val1
Difference_2 <- CV$Top_Cluster_FDR - CV$Ancestor_Val2
CV <- cbind(CV, Difference_1, Difference_2)
VX <- cbind(dF[ ,1:2], CV[ ,1], dF[ , 3:5], CV[ ,2:5])
VX <- VX[ order(VX$Difference_1, decreasing = T), ]
colnames(VX)[3] <- "Top_Cluster_FDR"
VX <- apply(VX,2,as.character)
write.table(VX, file = "GSE75688_Breast_KTREE_Significane_JUMP_Summary.csv", sep = ",", row.names = F)

### Generate Stumps Map with new defined colors
stumps.col <- Clust.Col
for (e in 1:length(stumps.L)) {
  pos <- which(Clust.Names%in%stumps.L[[e]])
  stumps.col[pos] <- thisPal[e]
}
stumps.O <- rev(sort(table(stumps.col)))
names(stumps.L) <- colnames(A_form)[-(ncol(A_form)-1)]
### Drawing cluster boxes with corresponding coordinates
tiff(paste("C:/Users/User/Desktop/CCCC.tiff", sep=""), width=2400, height=3200, compression="lzw", res=300)
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
legend("topright",legend=names(stumps.L), fill=thisPal[1:length(stumps.L)], border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

DATA <- read.csv(paste("C:/Users/User/Desktop/GSE75688_Breast_ClusGenes/BC_Summary/GSE75688_Breast_KClusters_GOWhole.csv", sep = ""), header = T, stringsAsFactors = F)
CORR <- as.numeric(DATA[DATA$Enrichment.GOTerms=="GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY", 3:ncol(DATA)])
mtext(text = "GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY", side = 3)
#5 Use round(-log10(CORR), digits = 2) if there is a vast range of small FDR values
for (b in 1:78) {
  if(round(-log10(CORR), digits = 2)[b]>= 1.3){
    text(X[b], Y[b], labels = round(-log10(CORR), digits = 2)[b], cex = 0.5)
  }
  if(round(-log10(CORR), digits = 2)[b]< 1.3){
    text(X[b], Y[b], labels = round(-log10(CORR), digits = 2)[b], cex = 0.2)
  }
}
dev.off()
