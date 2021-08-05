### Pipeline for SC cancer datasets tsne plot comparisons
CaseN <- c("GSE75688_Breast", "GSE72056_Melanoma", "GSE70630_OG", "Astrocytoma", "GSE81383_Melanoma", "GSE103322_HNSCC", "E_MTAB_6149_NSCLC", "GSE81861_CRC")

### Import Packages
library(cluster)
library(Rtsne)

### Sample labeling and colors
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe"
           , "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")

colorss = c("#ffe119", "#f58231", "#42d4f4", "#469990", "#800000", "#aaffc3", "#e6beff", "#bfef45", "#e6194B")

### For cases that have both normal and cancer cells mixed tables 
CaseN <- CaseN[N]
Filepath <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", sep = "")
Filepath1 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", CaseN, "_RAW_CDF.csv", sep = "")
Filepath2 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/Normal/", CaseN, "Nor_RAW_CDF.csv", sep = "")
Filepath3 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/Cancerous/", CaseN, "Can_RAW_CDF.csv", sep = "")
AD <- read.csv(Filepath1, header = T, stringsAsFactors = F)
AD <- AD[!duplicated(AD[ , 1]), ]
rownames(AD) <- AD$X
AD <- AD[ ,-1]
AD <- AD[!apply(AD,1,function(x) all(is.na(x))), ]

# Patient Identity Colors labeling
Patient_colors=rep(colors[1], ncol(AD))
for(i in 2:length(table(as.numeric(AD[1, ])))){
  Patient_colors[AD[1, ]==names(table(as.numeric(AD[1, ])))[i]]=colors[i]
}
# Cell type Colors labeling
Cell_type <- rep(colorss[1], ncol(AD))
for(i in 2:length(table(as.numeric(AD[6, ])))){
  Cell_type[AD[6, ]==names(table(as.numeric(AD[6, ])))[i]]=colorss[i]
}
# For GSE103322_HNSCC: Cancer&unknown=0 ; Fibroblast=1 ; T cell=2 ; Endothelial=3 ; B cell=4 ; Mast=5 ; Macrophage=6 ; Dendritic=7 ; myocyte=8
REPLAS <- c("Zero", "Mean")

for (g in 1:length(REPLAS)) {
  if(g==1){
    AA <- AD[7:nrow(AD), ] 
    AA[is.na(AA)] = 0 # Convert NAs to "0" and "AA" ready for tsne
    tsne <- Rtsne(t(AA), dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
    tiff(paste(Filepath, CaseN, "_ZeroCell.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(tsne$Y, main= paste(CaseN, "_Celltype_ZeroCDF", sep=""), cex.main = 0.8, col= Cell_type, pch = 16, cex = 0.4)
    par(xpd=T)
    mtext("Cancer=0;Fibro=1;T=2;Endo=3;B=4;Mast=5;Macro=6;Dendri=7;myo=8", side = 3)
    legend("topright",legend=c(names(table(as.numeric(AD[6,])))), fill=c(colorss[1:length(table(as.numeric(AD[6, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
    dev.off()
    tiff(paste(Filepath, CaseN, "_ZeroPat.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(tsne$Y, main= paste(CaseN, "_PatientID_ZeroCDF", sep=""), cex.main = 0.8, col= Patient_colors, pch = 16, cex = 0.4)
    par(xpd=T)
    legend("topright",legend=c(names(table(as.numeric(AD[1,])))), fill=c(colors[1:length(table(as.numeric(AD[1, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
    dev.off()
  }else{
    AA <- AD[7:nrow(AD), ]
    for (i in 1:nrow(AA)) {
      AA[ i, is.na(AA[i, ])] <- mean(na.omit(as.numeric(AA[i, ])))
    } # Convert NAs to "0" for "Mean version"
    tsne <- Rtsne(t(AA), dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
    tiff(paste(Filepath, CaseN, "_MeanCell.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(tsne$Y, main= paste(CaseN, "_Celltype_MeanCDF", sep=""), cex.main = 0.8, col= Cell_type, pch = 16, cex = 0.4)
    par(xpd=T)
    mtext("Cancer=0;Fibro=1;T=2;Endo=3;B=4;Mast=5;Macro=6;Dendri=7;myo=8", side = 3)
    legend("topright",legend=c(names(table(as.numeric(AD[6,])))), fill=c(colorss[1:length(table(as.numeric(AD[6, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
    dev.off()
    tiff(paste(Filepath, CaseN, "_MeanPat.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(tsne$Y, main= paste(CaseN, "_PatientID_MeanCDF", sep=""), cex.main = 0.8, col= Patient_colors, pch = 16, cex = 0.4)
    par(xpd=T)
    legend("topright",legend=c(names(table(as.numeric(AD[1,])))), fill=c(colors[1:length(table(as.numeric(AD[1, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
    dev.off()
  }
}