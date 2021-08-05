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

# Install "NNLM"
install.packages("NNLM")
library(NNLM)

MSE <- NULL
for(k in 2:5){
  p_time <- proc.time()
  NNLM_res <- nnmf(V, k)
  t_time <- proc.time()-p_time
  print(t_time)
  MSE[k] <- sum(apply((V - (NNLM_res$W %*% NNLM_res$H))**2, 1, sum))
}

tiff(paste("C:/Users/tonyxp/Desktop/PVal.05_MSE.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
plot(2:30, MSE, main = "PVal.05_MSE", xlab = "k", pch = 16, cex = 0.4)
dev.off()

# Install "NMF"
install.packages('NMF')
library(NMF)

#Variable of iterative subset GOterms
A <- V[rownames(INDV[INDV$K_Group == 1, ]), ]
B <- V[rownames(INDV[INDV$K_Group == 2, ]), ]

C <- A[rownames(INDA[INDA$K_Group == 1, ]), ]
D <- A[rownames(INDA[INDA$K_Group == 2, ]), ]
E <- B[rownames(INDB[INDB$K_Group == 1, ]), ]
G <- B[rownames(INDB[INDB$K_Group == 2, ]), ]

H <- C[rownames(INDC[INDC$K_Group == 1, ]), ]
I <- C[rownames(INDC[INDC$K_Group == 2, ]), ]
J <- D[rownames(INDD[INDD$K_Group == 1, ]), ]
K <- D[rownames(INDD[INDD$K_Group == 2, ]), ]
L <- E[rownames(INDE[INDE$K_Group == 1, ]), ]
M <- E[rownames(INDE[INDE$K_Group == 2, ]), ]
N <- G[rownames(INDG[INDG$K_Group == 1, ]), ]
O <- G[rownames(INDG[INDG$K_Group == 2, ]), ]

A <- A[ , !(apply(A, 2, function(x) all(x==0)))]
B <- B[ , !(apply(B, 2, function(x) all(x==0)))]

C <- C[ , !(apply(C, 2, function(x) all(x==0)))]
D <- D[ , !(apply(D, 2, function(x) all(x==0)))]
E <- E[ , !(apply(E, 2, function(x) all(x==0)))]
G <- G[ , !(apply(G, 2, function(x) all(x==0)))]

H <- H[ , !(apply(H, 2, function(x) all(x==0)))]
I <- I[ , !(apply(I, 2, function(x) all(x==0)))]
J <- J[ , !(apply(J, 2, function(x) all(x==0)))]
K <- K[ , !(apply(K, 2, function(x) all(x==0)))]
L <- L[ , !(apply(L, 2, function(x) all(x==0)))]
M <- M[ , !(apply(M, 2, function(x) all(x==0)))]
N <- N[ , !(apply(N, 2, function(x) all(x==0)))]
O <- O[ , !(apply(O, 2, function(x) all(x==0)))]

p_time <- proc.time()
res <- nmf(V, 2)
w_V<- basis(res)
h_V<- coef(res)
t_time <- proc.time()-p_time
print(t_time)

p_time <- proc.time()
res <- nnmf(V, 16)
nw_V<- res$W
nh_V<- res$H
t_time <- proc.time()-p_time
print(t_time)

#Alter row sequence according to "w" matrix
p_time <- proc.time()
INDVV <- as.data.frame(apply(nw_V, 1, function(x) max(x)))
INDVV$GROUP <- max.col(nw_V)
colnames(INDVV) <- c("Max_V", "K_Group")
SS <- INDVV[order(INDVV$K_Group), ]
ww <- nw_V[match(rownames(SS), rownames(nw_V)), ]
tiff(paste("C:/Users/tonyxp/Desktop/NMF_FILES/V_Split_NMF_W/PV95_W_V_nnmf.tiff", sep=""), width=6000, height=6000, compression="lzw", res=300)
heatmap.2(ww, labRow = rownames(ww), col=grcol, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", main = "W_nnmf", cexRow = 0.6, margins = c(5,25))
dev.off()
t_time <- proc.time()-p_time
print(t_time)

#Alter col sequence according to "h" matrix
p_time <- proc.time()
COE <- as.data.frame(apply(h_V, 2, function(x) max(x)))
COE$GROUP <- max.col(t(h_V))
colnames(COE) <- c("Max_V", "K_Group")
YY <- COE[order(COE$K_Group), ]
hh <- h_V[ ,match(rownames(YY), colnames(h_V))]
tiff(paste("C:/Users/tonyxp/Desktop/NMF_FILES/V_Split_NMF_H/PV95_H_V_nmf.tiff", sep=""), width=3000, height=3000, compression="lzw", res=300)
heatmap.2(hh, labCol=colnames(hh), col=grcol, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", main = "H_V", cexRow = 1.2)
dev.off()
t_time <- proc.time()-p_time
print(t_time)

#Original table heatmaps in wbasis order
p_time <- proc.time()
New_V <- V[match(rownames(SS), rownames(V)), ]
tiff(paste("C:/Users/tonyxp/Desktop/NMF_FILES/NMF_X/PV95_Xori_K", j,".tiff", sep=""), width=6000, height=6000, compression="lzw", res=300)
heatmap.2(New_V, labCol=colnames(New_V), col=grcol, Rowv = FALSE, dendrogram = "column", trace = "none", main = "New_V", cexRow = 0.6, margins = c(5,15))
dev.off()
t_time <- proc.time()-p_time
print(t_time)

#Writeout by hand
SS <- IND16[order(IND16$K_Group), ]
SS$TERMS <- rownames(SS)
write.table(SS, file = "IND16.csv", sep = ",", row.names = FALSE)
w_V$TERMS <- rownames(w_V)
write.table(nw_V, file = "w_V_nnmf.csv", sep = ",", row.names = TRUE)
h_V$TERMS <- rownames(h_V)
write.table(nh_V, file = "h_V_nnmf.csv", sep = ",", row.names = TRUE)

intersect(rownames(INDVV[ INDVV$K_Group == 1, ]), rownames(INDV[ INDV$K_Group == 1, ]))
###
SSS <- c(INDV, INDA, INDB, INDC, INDD, INDE, INDG, INDH, INDI, INDJ, INDK, INDL, INDM, INDN, INDO)
DV <- NULL
for(i in 1:1){
  write.table(i, file = "IND16.csv", sep = ",", row.names = FALSE)
}
###
length(rownames(INDV[INDV$K_Group == 1, ]))
length(rownames(INDV[INDV$K_Group == 2, ]))
length(rownames(INDVV[INDVV$K_Group == 1, ]))
length(rownames(INDVV[INDVV$K_Group == 2, ]))
length(intersect(rownames(INDV[INDV$K_Group == 1, ]), rownames(INDVV[INDVV$K_Group == 1, ])))
length(intersect(rownames(INDV[INDV$K_Group == 2, ]), rownames(INDVV[INDVV$K_Group == 2, ])))
length(intersect(rownames(INDV[INDV$K_Group == 1, ]), rownames(INDVV[INDVV$K_Group == 2, ])))
length(intersect(rownames(INDV[INDV$K_Group == 2, ]), rownames(INDVV[INDVV$K_Group == 1, ])))