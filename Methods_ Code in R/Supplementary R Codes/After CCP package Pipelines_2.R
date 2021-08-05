### Location address Setting
CaseN <- "GSE75688_Breast"
CaseN <- "GSE72056_Melanoma"
CaseN <- "GSE70630_OG"
CaseN <- "Astrocytoma"
CaseN <- "GSE81383_Melanoma"
CaseN <- "GSE103322_HNSCC"
# 18 colors for marking different clusters
thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
             "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f", "#000000", "#076f25", "#93cd7f", 
             "#4d0776", "#ffffff")
# N distinctive colors for labeling stumps
library(RColorBrewer)
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

# Clusters Coordinates
X <- NULL
Y <- NULL
for (k in 1:12) {
  for (p in 1:k) {
    X <- append(X, 39/(k+1)*p, after = length(X))
    Y <- append(Y, 396-(k-1)*44, after = length(Y))
  }
}
# Find out ClustLines in the tree (One type)
op <- list.files(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/GeneColorsLists", sep = ""), pattern = "CPP_K")
Clust.Col <- NULL
for (i in 1:12) {
  dt <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/GeneColorsLists/", op[i], sep = ""), header = T, stringsAsFactors = F)
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
Clust.Num <- sapply(BVC[,2:ncol(BVC)], function(x) sum(x))

#################################################################################
### ### 6.2.1 Child cluster hierarchy tree (Use CluRatio dvd by Child as input)
# Convert upper triangle part of Child CluRatio into Child linear vector
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
tiff(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/K_Maps/", CaseN,"_Whole_Overlapping_ChildRatio.tiff", sep = ""), width=1600, height=1600, compression="lzw", res=300)
hist(Global.ratio)
arrows(0.2, 250, 0.4, 350, lwd = 0.5)
arrows(0.4, 250, 0.6, 350, lwd = 2.5)
arrows(0.6, 250, 0.8, 350, lwd = 5.0)
arrows(0.8, 250, 1.0, 350, lwd = 8.0)
dev.off()

tiff(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/K_Maps/", CaseN,"_General_Tree_ChildView.tiff", sep = ""), width=2400, height=3200, compression="lzw", res=300)
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
##########################################################################################
### ### 6.2.2 Parent cluster hierarchy tree (Use CluRatio dvd by Parent as input)
# Convert upper triangle part of Parent CluRatio into Parent linear vector
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
tiff(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/K_Maps/", CaseN,"_Whole_Overlapping_ParentRatio.tiff", sep = ""), width=1600, height=1600, compression="lzw", res=300)
hist(Global.ratio)
arrows(0.2, 250, 0.4, 350, lwd = 0.5)
arrows(0.4, 250, 0.6, 350, lwd = 2.5)
arrows(0.6, 250, 0.8, 350, lwd = 5.0)
arrows(0.8, 250, 1.0, 350, lwd = 8.0)
dev.off()

tiff(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/K_Maps/", CaseN,"_General_Tree_ParentView.tiff", sep = ""), width=2400, height=3200, compression="lzw", res=300)
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