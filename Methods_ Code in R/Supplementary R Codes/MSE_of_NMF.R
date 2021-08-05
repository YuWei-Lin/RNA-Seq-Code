# Original NMF in R package
p_time <- proc.time()
res174 <- nmf(V, 174)
w_174<- basis(res174)
h_174<- coef(res174)
t_time <- proc.time()-p_time
print(t_time)

write.table(V, file = "pV.95_Binary_GOTerms.csv", sep = ",", row.names = TRUE)
write.table(Z, file = "pV.05_Binary_GOTerms.csv", sep = ",", row.names = TRUE)
# Install "NNLM"
install.packages("NNLM")
library(NNLM)

MSE <- NULL
for(k in 2:174){
  
  NNLM_res <- nnmf(V, k)

  MSE[k-1] <- sum(apply((V - (NNLM_res$W %*% NNLM_res$H))**2, 1, sum))
  print(k)
}

tiff(paste("C:/Users/tonyxp/Desktop/PVal.95_NNLM_MSE.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
plot(2:174, MSE, main = "PVal.95_NNLM_MSE", xlab = "k", pch = 16, cex = 0.4)
dev.off()

#Alter row sequence according to "w" matrix
p_time <- proc.time()
IND8 <- as.data.frame(apply(w_8, 1, function(x) max(x)))
IND8$GROUP <- max.col(w_8)
colnames(IND8) <- c("Max_V", "K_Group")
SS <- IND8[order(IND8$K_Group), ]
ww <- w_8[match(rownames(SS), rownames(w_8)), ]
tiff(paste("C:/Users/tonyxp/Desktop/NMF_FILES/V_Split_NMF_W/PV95_W_8.tiff", sep=""), width=6000, height=6000, compression="lzw", res=300)
heatmap.2(ww, labRow = rownames(ww), col=grcol, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", main = "W_8", cexRow = 0.6, margins = c(5,30))
dev.off()
t_time <- proc.time()-p_time
print(t_time)

#Alter col sequence according to "h" matrix
p_time <- proc.time()
COE <- as.data.frame(apply(h_8, 2, function(x) max(x)))
COE$GROUP <- max.col(t(h_8))
colnames(COE) <- c("Max_V", "K_Group")
YY <- COE[order(COE$K_Group), ]
hh <- h_8[ ,match(rownames(YY), colnames(h_8))]
tiff(paste("C:/Users/tonyxp/Desktop/NMF_FILES/V_Split_NMF_H/PV95_H_8.tiff", sep=""), width=3000, height=3000, compression="lzw", res=300)
heatmap.2(hh, labCol=colnames(hh), col=grcol, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", main = "H_8", cexRow = 1.2)
dev.off()
t_time <- proc.time()-p_time
print(t_time)

#Writeout by hand

SS <- IND8[order(IND8$K_Group), ]
SS$TERMS <- rownames(SS)
write.table(SS, file = "IND8.csv", sep = ",", row.names = FALSE)

for(p in 1:5){
  XX <- NULL
  if(p == p){
    XX[p] <- "NaN"
  }
}

