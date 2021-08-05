# Import Packages
library(cluster)
library(Rtsne)
library(gplots)

# Import GO and processed ssc-RNAseq data
GOIDName <- read.csv("C:/Users/tonyxp/Desktop/CDF_PureCancer_SC_ReadyFiles/GSE75688_Breast/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
GOIDList <- read.csv("C:/Users/tonyxp/Desktop/CDF_PureCancer_SC_ReadyFiles/GSE75688_Breast/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
Mel.stand <- read.csv("C:/Users/tonyxp/Desktop/CDF_PureCancer_SC_ReadyFiles/GSE75688_Breast/CDF_GSE75688_BreastCan_ZRM30.csv", header = TRUE)


Mel.stand <- t(Mel.stand)
# Pre-screening and generating of GOnames
N <- table(GOIDList$GOID)
GNlist <- NULL
for(i in 1:length(N)){
  
  if(length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand))) >= 50 & 
     length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand))) <= 500){
    AD <- Mel.stand[ ,colnames(Mel.stand)%in%GOIDList[which(GOIDList$GOID==(names(N)[i])),1]]
    GNlist <- append(GNlist, names(N)[i], after = length(GNlist))
    write.csv(AD, file = paste("C:/Users/tonyxp/Desktop/BreastCan_ZRMCDF/GSE75688_Mig_GOid/GO_", substr(names(N)[i],4,nchar(as.character(names(N)[i]))), sep="",".csv"))
  }
}
write.csv(GNlist , "GNlistGSE7568_ZRCDF.csv",row.names = FALSE)
#Union the genes as a pool after GO screening
gepoo <- list.files("C:/Users/tonyxp/Desktop/BreastCan_ZRMCDF/GSE75688_Mig_GOid")
gepo <- gsub("_", ":", gepoo)
gepo <- substr(gepo, 4, 10)
tar <- NULL
for (i in 1:length(gepo)){
  tar <- rbind(tar, GOIDList[substr(GOIDList$GOID, 4, 10) == gepo[i], ])
}
tartest <- tar[!duplicated(tar$Gene), ]
write.table(tartest, file = "Genepool_of_GSE75688_Mig.csv", sep = ",", row.names = FALSE)

# Screening and generate GOnames, GOtsne, and GOsil 
tartest <- read.table("C:/Users/tonyxp/Desktop/R00006/Genepool of GO screening of GSE75688 FIX.csv", header = TRUE, sep = ",")
Mel.stand2 <- Mel.stand[ , colnames(Mel.stand)%in%(tartest$Gene)]

Mel.stand2[is.nan(Mel.stand2)] <- 0

N <- table(GOIDList$GOID)
GNlist <- NULL
GNVa <- NULL
GNum <- NULL
Indsil <- NULL
GNAvgss <- NULL
p_time <- proc.time()
for(i in 1:length(N)){
  
  if(length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand2))) >= 50 & 
     length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand2))) <= 500){
    AD <- Mel.stand2[ ,colnames(Mel.stand2)%in%GOIDList[which(GOIDList$GOID==(names(N)[i])),1]]
    AD <- AD[!apply(AD, 1,function(x) all(abs(x-x[1])<0.00000000001)), ]
    tsne <- Rtsne(AD, dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
    cormat <- cor(t(AD), method = "pearson", use = "pairwise.complete.obs")
    cortrans <- 2*(1-cormat) 
    sil = silhouette (gps, cortrans)
    tt <- sil[ ,"sil_width"]
    
    GNsubsil <- NULL
    GNlist <- append(GNlist, names(N)[i], after = length(GNlist))
    GNVa <- append(GNVa, length(tt[tt<=0])/length(tt), after = length(GNVa))
    for(k in 1:length(table(gps))){
      KK <- sil[ ,"cluster"] == k
      GNsubsil <- append(GNsubsil, sum(sil[KK,"sil_width"]<=0)/sum(KK), after = length(GNsubsil))
    }
    Indsil <- rbind(Indsil, GNsubsil)
    GNAvgss <- append(GNAvgss, mean(GNsubsil), after = length(GNAvgss))
    GNum <- append(GNum, ncol(AD), after = length(GNum))
    
    tiff(paste("C:/Users/tonyxp/Desktop/BreastCan_ZRMCDF/BreastCan_ZRMCDF_tsne/", substr(names(N)[i],4,nchar(as.character(names(N)[i])))," Autotsne.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(tsne$Y, main= paste(names(N)[i], format(length(tt[tt<=0])/length(tt), digits = 2, format = T), unique(GOIDName[GOIDName$GOID==names(N)[i], ][3])), cex.main = 0.8, col= cc, pch = 16, xlim = c(-35,35), ylim = c(-35,35), cex = 0.4)
    dev.off()
    
    tiff(paste("C:/Users/tonyxp/Desktop/BreastCan_ZRMCDF/BreastCan_ZRMCDF_sil/", substr(names(N)[i],4,nchar(as.character(names(N)[i])))," Autosil.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(sil, main = paste(names(N)[i], format(length(tt[tt<=0])/length(tt), digits = 2, format = T), unique(GOIDName[GOIDName$GOID==names(N)[i], ][3])), size = 0.8, col= colors[1:13], border=NA)
    dev.off()
  }
}
t_time <- proc.time()-p_time
print(t_time)

# Randomly runs for computing p-value
p_time <- proc.time()
GNVaAcu <- NULL
GNVaAcuInd <- NULL
for(i in 1:length(GNum)){
  GNVaR <- NULL
  GNAvgssR <- NULL
  for(j in 1:1000){
    ADR <- Mel.stand2[ ,sample(ncol(Mel.stand2), GNum[i], replace = FALSE)]
    cormatR <- cor(t(ADR), method = "pearson")
    cortransR <- 2*(1-cormatR) 
    silR = silhouette (gps, cortransR)
    ttR <- silR[ ,"sil_width"]
    GNVaR <- append(GNVaR, length(ttR[ttR<=0])/length(ttR), after = length(GNVaR))
    GNsubsilR <- NULL
    for(k in 1:length(table(gps))){
      KKR <- silR[ ,"cluster"] == k
      GNsubsilR <- append(GNsubsilR, sum(silR[KKR,"sil_width"]<=0)/sum(KKR), after = length(GNsubsilR))
    }
    GNAvgssR <- append(GNAvgssR, mean(GNsubsilR), after = length(GNAvgssR))  
  }
  GNVaAcu <- append(GNVaAcu, 1-(sum(as.numeric(GNVa[i] > GNVaR))/1000), after = length(GNVaAcu))
  GNVaAcuInd <- append(GNVaAcuInd, 1-(sum(as.numeric(GNAvgss[i] > GNAvgssR))/1000), after = length(GNVaAcuInd))
}
t_time <- proc.time()-p_time
print(t_time)
E1 <- cbind(GNlist, GNum, Indsil)
E2 <- cbind(GNAvgss, GNVa, GNVaAcuInd, GNVaAcu)
E <- cbind(E1, E2)
E <- as.data.frame(E, stringsAsfactor = FALSE)
E <- E[rev(order(E$GNVaAcu)), ]
NN <- GOIDName[GOIDName$GOID%in%E$GNlist, ]
E <- cbind(NN[match(E$GNlist, NN$GOID), 2], E)
colnames(E) <- c("GO terms", "GO-ID", "Gene Numbers", paste("BC0", 1:9, sep = ""), "BC10", "BC11", "BC03LN","BC07LN", "Clust_NSV", "Whole_NSV", "Clust_p-Val", "Whole_p-Val")
write.table(E, file = "GSE75688_BC_PatientID FIX Go-Scoring-of Clust&Whole_p-val.csv", sep = ",", row.names = FALSE)

