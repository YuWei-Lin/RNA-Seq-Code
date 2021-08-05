MEL <- read.table(file.choose(), header = TRUE, stringsAsFactors = FALSE) 

write.table(MEL, file = "GSE72056_Melanoma_RAW.csv", sep = ",", row.names = FALSE)

BC <- read.table(file.choose(), header = TRUE, stringsAsFactors = FALSE) 

write.table(BC, file = "GSE75688_BreastCan_RAW.csv", sep = ",", row.names = FALSE)

MigDB <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)

GG <- NULL
for(i in unique(MigDB$ProcessName)){
  AB <- matrix(0, 1, 17106)
  colnames(AB) <- unique(MigDB$Gene)
  rownames(AB) <- i
  AB[ ,colnames(AB)%in%MigDB[MigDB$ProcessName == i, 1]] <- 1
  GG <- rbind(GG, AB)
}

write.table(GG, file = "MigDB_Binary_Table.csv", sep = ",")
GG <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
PV05 <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
PV95 <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
X <- PV05[PV05$Total >= 4, 1]
Y <- PV95[PV95$Total >= 4, 1]
Z <- GG[X, ]
V <- GG[Y, ]

#Remove Zero columns

Z <- Z[ , !(apply(Z, 2, function(x) all(x==0)))]
V <- V[ , !(apply(V, 2, function(x) all(x==0)))]

#Color customization - see help for brewer.pal for more information
library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(9, "Oranges"))(100)
image(1:11,1,as.matrix(1:11), col = brewer.pal(11, "Paired"))
#Color customization - see help for brewer.pal for more information
library(RColorBrewer) 
grcol <- colorRampPalette(c("lightyellow","green","red"))(100)
#Column labels - patients
labels <- colnames(Z)
#Customized heatmap
library(gplots) 
library(rafalib)
Z <- as.matrix(Z)
V <- as.matrix(V)
tiff(paste("C:/Users/tonyxp/Desktop/PVal.05_terms_hp.tiff", sep=""), width=3000, height=3000, compression="lzw", res=300)
heatmap.2(Z, labCol=labels, col=grcol, hclustfun = hclust, dendrogram = "both", trace = "none", main = "PVal.05", cexRow = 0.3)
dev.off()
heatmap.2(M, col=grcol, hclustfun = hclust, dendrogram = "column", trace = "none")
###trace = "none", margins = c(7,5), scale = "row"
# Install "NMF"
install.packages('NMF')
library(NMF)

p_time <- proc.time()
MSE <- NULL
K.SSB <- NULL
K.SSW <- NULL
for(j in 2:3){
  res <- nmf(V, j)
  w<- basis(res)
  h<- coef(res)
  #MSE[j-1] <- sum(apply((Z - (w %*% h))**2, 1, sum))
  print(j)
  
  #Alter col sequence according to "h" matrix
  COE <- as.data.frame(apply(h, 2, function(x) max(x)))
  COE$GROUP <- max.col(t(h))
  colnames(COE) <- c("Max_V", "K_Group")
  YY <- COE[order(COE$K_Group), ]
  hh <- h[ ,match(rownames(YY), colnames(h))]
  tiff(paste("C:/Users/tonyxp/Desktop/NMF_FILES/NMF_H/PV95_hhcoef_K", j,".tiff", sep=""), width=6000, height=6000, compression="lzw", res=300)
  heatmap.2(hh, labCol=colnames(hh), col=grcol, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", main = "hhcoef", cexRow = 0.6, margins = c(5,15))
  dev.off()
  #Alter row sequence according to "w" matrix
  IND <- as.data.frame(apply(w, 1, function(x) max(x)))
  IND$GROUP <- max.col(w)
  colnames(IND) <- c("Max_V", "K_Group")
  SS <- IND[order(IND$K_Group), ]
  ww <- w[match(rownames(SS), rownames(w)), ]
  
  tiff(paste("C:/Users/tonyxp/Desktop/NMF_FILES/NMF_W/PV95_wwbasis_K", j,".tiff", sep=""), width=6000, height=6000, compression="lzw", res=300)
  heatmap.2(ww, labCol=colnames(ww), col=grcol, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", main = "wwbasis", cexRow = 0.6, margins = c(5,15))
  dev.off()
  
  #Original table heatmaps in wbasis order
  New_V <- V[match(rownames(SS), rownames(V)), ]
  
  tiff(paste("C:/Users/tonyxp/Desktop/NMF_FILES/NMF_X/PV95_Xori_K", j,".tiff", sep=""), width=6000, height=6000, compression="lzw", res=300)
  heatmap.2(New_V, labCol=colnames(New_Z), col=grcol, Rowv = FALSE, dendrogram = "column", trace = "none", main = "New_V", cexRow = 0.6, margins = c(5,15))
  dev.off()
  
  #Compute ANOVA
  Inner <- NULL
  Inter <- NULL
  for(n in 1:j){
    GM <- mean(IND$Max_V)
    Inner[n] <- nrow(w[IND$K_Group == n, ])*((mean(w[IND$K_Group == n, n])-GM)**2)
    SSB <- sum(Inner)
    Inter[n] <- sum((w[IND$K_Group == n, n]-mean(w[IND$K_Group == n, n]))**2)
    SSW <- sum(Inter)
  }
  K.SSB <- append(K.SSB, SSB, after = length(K.SSB))
  K.SSW <- append(K.SSW, SSW, after = length(K.SSW))
}
t_time <- proc.time()-p_time
print(t_time)

# Check Values
table(apply(w, 1, function(x) all(x < 2.220446e-15)))

tiff(paste("C:/Users/tonyxp/Desktop/PVal.05_MSE.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
plot(2:30, MSE, main = "PVal.05_MSE", xlab = "k", pch = 16, cex = 0.4)
dev.off()

tiff(paste("C:/Users/tonyxp/Desktop/PVal.95_Wbasis.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
basismap(res, main = "PVal.95_Wbasis Components")
dev.off()

tiff(paste("C:/Users/tonyxp/Desktop/PVal.95_Hcoef.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
coefmap(res, main = "PVal.95_Hcoef Mix")
dev.off()



