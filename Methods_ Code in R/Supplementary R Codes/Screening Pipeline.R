# Import Packages
library(cluster)
library(Rtsne)
library(gplots)

# Import GO and processed ssc-RNAseq data
GOIDName <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
GOIDList <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
Mel.stand <- read.csv("C:/Users/tonyxp/Desktop/R00006/R0006 TNBC standardized table.CSV", header = TRUE)

# Pre-screening and generating of GOnames
N <- table(GOIDList$GOID)
GNlist <- NULL
for(i in 1:length(N)){
  
  if(length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand))) >= 50 & 
     length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand))) <= 500){
    AD <- Mel.stand[ ,colnames(Mel.stand)%in%GOIDList[which(GOIDList$GOID==(names(N)[i])),1]]
    GNlist <- append(GNlist, names(N)[i], after = length(GNlist))
    write.csv(AD, file = paste("C:/Users/tonyxp/Desktop/ALL GO Runs/GSE72056 FIX AuGOid/GO_", substr(names(N)[i],4,nchar(as.character(names(N)[i]))), sep="",".csv"))
  }
}

#Union the genes as a pool after GO screening
gepoo <- list.files("C:/Users/tonyxp/Desktop/ALL GO Runs/GSE72056 FIX AuGOid")
gepo <- gsub("_", ":", gepoo)
gepo <- substr(gepo, 4, 10)
tar <- NULL
for (i in 1:length(gepo)){
  tar <- rbind(tar, GOIDList[substr(GOIDList$GOID, 4, 10) == gepo[i], ])
}
tartest <- tar[!duplicated(tar$Primary.symbol), ]
write.table(tartest, file = "Genepool of GO screening of GSE72056 FIX.csv", sep = ",", row.names = FALSE)

# Screening and generate GOnames, GOtsne, and GOsil 
tartest <- read.table("C:/Users/tonyxp/Desktop/R00006/Genepool of GO screening of GSE75688 FIX.csv", header = TRUE, sep = ",")
Mel.stand2 <- Mel.stand[ , colnames(Mel.stand)%in%(tartest$Primary.symbol)]
N <- table(GOIDList$GOID)
GNlist <- NULL
GNVa <- NULL
GNum <- NULL
p_time <- proc.time()
for(i in 1:length(N)){
  
  if(length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand2))) >= 50 & 
     length(intersect(GOIDList[which(GOIDList$GOID==(names(N)[i])),1],colnames(Mel.stand2))) <= 500){
    AD <- Mel.stand2[ ,colnames(Mel.stand2)%in%GOIDList[which(GOIDList$GOID==(names(N)[i])),1]]
    tsne <- Rtsne(AD, dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
    cormat <- cor(t(AD), method = "pearson")
    cortrans <- 2*(1-cormat) 
    sil = silhouette (gps, cortrans)
    tt <- sil[ ,"sil_width"]
    
    GNlist <- append(GNlist, names(N)[i], after = length(GNlist))
    GNVa <- append(GNVa, length(tt[tt<=0])/length(tt), after = length(GNVa))
    GNum <- append(GNum, ncol(AD), after = length(GNum))
    
    tiff(paste("C:/Users/tonyxp/Desktop/ALL GO Runs/GSE75688(BC_PatientID) FIX AuGOtsne/", substr(names(N)[i],4,nchar(as.character(names(N)[i])))," Autotsne.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(tsne$Y, main= paste(names(N)[i], format(length(tt[tt<=0])/length(tt), digits = 2, format = T), GOIDName[GOIDName$GOID==names(N)[i], ][2]), col= cc, pch = 16, xlim = c(-35,35), ylim = c(-35,35), cex = 0.4)
    dev.off()
    
    tiff(paste("C:/Users/tonyxp/Desktop/ALL GO Runs/GSE75688(BC_PatientID) FIX AuGOsil/", substr(names(N)[i],4,nchar(as.character(names(N)[i])))," Autosil.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(sil, main = paste(names(N)[i], format(length(tt[tt<=0])/length(tt), digits = 2, format = T), GOIDName[GOIDName$GOID==names(N)[i], ][2]), col= colors[1:11], border=NA)
    dev.off()
  }
}
t_time <- proc.time()-p_time
print(t_time)
p_time <- proc.time()
GNVaAcu <- NULL
for(i in 1:length(GNum)){
  GNVaR <- NULL
  for(j in 1:1000){
    ADR <- Mel.stand2[ ,sample(ncol(Mel.stand2), GNum[i], replace = FALSE)]
    cormatR <- cor(t(ADR), method = "pearson")
    cortransR <- 2*(1-cormatR) 
    silR = silhouette (gps, cortransR)
    ttR <- silR[ ,"sil_width"]
    GNVaR <- append(GNVaR, length(ttR[ttR<=0])/length(ttR), after = length(GNVaR))
  }
  GNVaAcu <- append(GNVaAcu, 1-(sum(as.numeric(GNVa[i] > GNVaR))/1000), after = length(GNVaAcu))
}
t_time <- proc.time()-p_time
print(t_time)
E <- cbind(GNlist, GNum, GNVa, GNVaAcu)
E <- as.data.frame(E, stringsAsfactor = FALSE)
E <- E[rev(order(E$GNVaAcu)), ]
NN <- GOIDName[GOIDName$GOID%in%E$GNlist, ]
E <- cbind(NN[match(E$GNlist, NN$GOID), 2], E)
colnames(E) <- c("GO terms", "GO-ID", "Gene Numbers", "NSV Percentage", "p-Value")
write.table(E, file = "GSE75688_BC_PatientID FIX Go-Scoring-of p-value.csv", sep = ",", row.names = FALSE)

