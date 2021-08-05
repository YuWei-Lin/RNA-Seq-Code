### General Pipe needs a input case name
CaseN <- "Astrocytoma"

### Import Packages
library(cluster)
library(Rtsne)

### Sample labeling and colors
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe"
           , "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")
colorss = c("#ffe119", "#f58231", "#42d4f4", "#469990", "#800000", "#aaffc3", "#e6beff", "#bfef45", "#e6194B")

### File Directories
Filepath1 <- paste("/NA2/tonyxp/", CaseN, "_DATA/Cancerous/", sep = "")

### Data sets and MSigDB import
GOIDList <- read.csv("/NA2/tonyxp/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
TARFILE2 <- list.files(Filepath1, pattern = "AnnoRows")
ADD <- read.csv(paste(Filepath1, TARFILE2, sep = ""), header = T, stringsAsFactors = F)
TXT <- ADD$Sample[2]
ADD <- ADD[ ,-1]
ADD <- ADD[ , order(ADD[1, ], decreasing = F)]
TARFILE <- list.files(Filepath1, pattern = "RAW_CDF_DEL")
Mel.stand <- read.csv(paste(Filepath1, TARFILE, sep = ""), header = T, stringsAsFactors = F)
Mel.stand <- Mel.stand[!duplicated(Mel.stand$Gene), ]
rownames(Mel.stand) <- Mel.stand$Gene
Mel.stand <- Mel.stand[ ,-1]
Mel.stand <- Mel.stand[apply(Mel.stand, 1, function(x) sum(is.na(x)) < (ncol(Mel.stand)*(0.7))), ]
Mel.stand <- Mel.stand[ , order(ADD[1, ], decreasing = F)]


### Patient Identity Colors labeling
titCP <- names(table(as.numeric(ADD[1, ])))
Patient_colors=rep(colors[1], ncol(ADD))
for(i in 2:length(titCP)){
  Patient_colors[ADD[1, ]==titCP[i]]=colors[i]
}
### Grouping based on Patient Identity
gps <- rep(titCP[1], ncol(ADD))
for(i in 2:length(titCP)){
  gps[ADD[1, ]==titCP[i]]=titCP[i]
}
gps <- as.numeric(gps)

#######################################################################
# If you want to use cell types as your groups for Normal cells, do following code part.
#######################################################################
### Cell type Colors labeling
# titCT <- names(table(as.numeric(ADD[2, ])))
# Cell_type <- rep(colorss[1], ncol(ADD))
# for(i in 2:length(titCT)){ Cell_type[ADD[2, ]==titCT[i]]=colorss[i] }
### Grouping based on Cell types
# gps <- rep(titCT[1], ncol(ADD))
# for(i in 2:length(titCT)){ gps[ADD[2, ]==titCT[i]]=titCT[i] }
# gps <- as.numeric(gps)


#######################################################################
# Union genes from all GOterms that pass the criteria ">=50 & <=500"
#######################################################################
### Pre-screening and generating of GOnames
N <- table(GOIDList$ProcessName)
GNlist <- NULL
tartest <- NULL
for(i in 1:length(N)){
  if(length(intersect(GOIDList[which(GOIDList$ProcessName==(names(N)[i])),1],rownames(Mel.stand))) >= 50 & 
     length(intersect(GOIDList[which(GOIDList$ProcessName==(names(N)[i])),1],rownames(Mel.stand))) <= 500){
    GNlist <-  append(GNlist, names(N)[i], after = length(GNlist))
    tartest <- append(tartest, GOIDList[which(GOIDList$ProcessName==(names(N)[i])),1], after = length(tartest))
  }
}
tartest <- unique(tartest)
Mel.stand2 <- Mel.stand[rownames(Mel.stand)%in%tartest, ]

#######################################################################
# Determine the NA values; either with "Mean" OR "Zero"
#######################################################################
### Convert NAs to "0"
# Mel.stand2[is.na(Mel.stand2)] = 0
# Filepath2 <- paste("/NA2/tonyxp/", CaseN, "_DATA/Cancerous/GOSiltsne/Filled with Zero/", sep = "")

### Convert NAs to "Mean"
for (i in 1:nrow(Mel.stand2)) { Mel.stand2[ i, is.na(Mel.stand2[i, ])] <- mean(na.omit(as.numeric(Mel.stand2[i, ]))) }
Filepath2 <- paste("/NA2/tonyxp/", CaseN, "_DATA/Cancerous/GOSiltsne/Filled with Mean/", sep = "")

### GOSiltsne Pipe and Negative Sil Value, NSV

GNVa <- NULL
GNum <- NULL
Indsil <- NULL
GNAvgss <- NULL
p_time <- proc.time()
for(i in 1:1){
  AD <- Mel.stand2[rownames(Mel.stand2)%in%GOIDList[which(GOIDList$ProcessName==(GNlist[i])), 1], ]
  AD <- AD[!apply(AD, 1,function(x) all(abs(x-x[1])<0.00000000001)), ]
  tsne <- Rtsne(t(AD), dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
  cormat <- cor(AD, method = "pearson", use = "pairwise.complete.obs")
  cortrans <- 2*(1-cormat) 
  sil = silhouette (gps, cortrans)
  tt <- sil[ ,"sil_width"]
  
  GNsubsil <- NULL
 
  GNVa <- append(GNVa, length(tt[tt<=0])/length(tt), after = length(GNVa))
  for(k in as.numeric(titCP)){
    KK <- sil[ ,"cluster"] == k
    GNsubsil <- append(GNsubsil, sum(sil[KK,"sil_width"]<=0)/sum(KK), after = length(GNsubsil))
  }
  Indsil <- rbind(Indsil, GNsubsil)
  GNAvgss <- append(GNAvgss, mean(GNsubsil), after = length(GNAvgss))
  GNum <- append(GNum, nrow(AD), after = length(GNum))
  
  tiff(paste(Filepath2, "tsneplots/", GNlist[i],".tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
  plot(tsne$Y, main= paste(format(length(tt[tt<=0])/length(tt), digits = 2, format = T), GNlist[i]), cex.main = 0.8, col= Patient_colors, pch = 16, cex = 0.4)
  par(xpd=T)
  legend("topright",legend=c(titCP), fill=c(colors[1:length(titCP)]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
  dev.off()
  
  tiff(paste(Filepath2, "Silplots/", GNlist[i],".tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
  plot(sil, main = paste(format(length(tt[tt<=0])/length(tt), digits = 2, format = T), GNlist[i]), col = Patient_colors, border=NA, cex = 0.1, cex.main = 0.1)
  dev.off()
  
  print(i)
}
t_time <- proc.time()-p_time
print(t_time)

### Randomly runs for computing p-value
p_time <- proc.time()
GNVaAcu <- NULL
GNVaAcuInd <- NULL
for(i in 1:length(GNum)){
  GNVaR <- NULL
  GNAvgssR <- NULL
  for(j in 1:500){
    ADR <- Mel.stand2[sample(nrow(Mel.stand2), GNum[i], replace = FALSE), ]
    cormatR <- cor(ADR, method = "pearson")
    cortransR <- 2*(1-cormatR) 
    silR = silhouette(gps, cortransR)
    ttR <- silR[ ,"sil_width"]
    GNVaR <- append(GNVaR, length(ttR[ttR<=0])/length(ttR), after = length(GNVaR))
    GNsubsilR <- NULL
    for(k in as.numeric(titCP)){
      KKR <- silR[ ,"cluster"] == k
      GNsubsilR <- append(GNsubsilR, sum(silR[KKR,"sil_width"]<=0)/sum(KKR), after = length(GNsubsilR))
    }
    GNAvgssR <- append(GNAvgssR, mean(GNsubsilR), after = length(GNAvgssR))  
  }
  GNVaAcu <- append(GNVaAcu, 1-(sum(as.numeric(GNVa[i] > GNVaR))/500), after = length(GNVaAcu))
  GNVaAcuInd <- append(GNVaAcuInd, 1-(sum(as.numeric(GNAvgss[i] > GNAvgssR))/500), after = length(GNVaAcuInd))
  print(i)
}
t_time <- proc.time()-p_time
print(t_time)

### Resulting table E
E <- cbind(GNlist[1:100], GNum, Indsil, GNAvgss, GNVa, GNVaAcuInd, GNVaAcu)
E <- as.data.frame(E, stringsAsfactor = FALSE)
E <- E[rev(order(E$GNVaAcu)), ]
for (m in 1:length(titCP)) {
  if(as.numeric(titCP[m])<10){
    titCP[m] <- paste("0", titCP[m], sep = "")
  }
}
colnames(E) <- c("GO terms", "Gene Numbers", paste(strsplit(CaseN, "_")[[1]][length(strsplit(CaseN, "_")[[1]])], "_", titCP, sep = ""), "Clust_NSV", "Whole_NSV", "Clust_p-Val", "Whole_p-Val")
write.table(E, file = paste(Filepath2, CaseN, "_GOSilScore_YY.csv", sep = ""), sep = ",", row.names = FALSE)

 
  