### Pipeline for merging the number of remaining genes after various cutoffs of NA levels in all 9 cases
CaseN <- c("Astrocytoma","GSE70630_OG","GSE72056_Melanoma","GSE75688_Breast","GSE81861_CRC","GSE103322_HNSCC","GSE81383_Melanoma","E_MTAB_6149_NSCLC","GSE76312_CML")
### Import Packages
library(cluster)
library(Rtsne)
### Sample labeling and colors
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe"
           , "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")

colorss = c("#ffe119", "#f58231", "#42d4f4", "#469990", "#800000", "#aaffc3", "#e6beff", "#bfef45", "#e6194B")

for (R in 1:length(CaseN)) {
  if(R<=6){
    CaseN <- CaseN[R]
    STAT <- c("Normal", "Cancerous")
    Filepath1 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", sep = "")
    Filepath2 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", STAT[1], "/", sep = "")
    Filepath3 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", STAT[2], "/", sep = "")
    for (WWW in c(Filepath1, Filepath2, Filepath3)) {
      TARFILE <- list.files(WWW, pattern = "RAW_CDF_DEL")
      SEE <- read.csv(paste(WWW, TARFILE, sep = ""), header = T, stringsAsFactors = F)
      SEE <- SEE[apply(SEE, 1, function(x) sum(is.na(x)) < (ncol(SEE)*(0.7))), ]
      TARFILE2 <- list.files(WWW, pattern = "AnnoRows")
      AD <- read.csv(paste(WWW, TARFILE2, sep = ""), header = T, stringsAsFactors = F)
      TXT <- AD$Sample[2]
      AD <- AD[ ,-1]
      # Patient Identity Colors labeling
      Patient_colors=rep(colors[1], ncol(AD))
      for(i in 2:length(table(as.numeric(AD[1, ])))){
        Patient_colors[AD[1, ]==names(table(as.numeric(AD[1, ])))[i]]=colors[i]
      }
      # Cell type Colors labeling
      Cell_type <- rep(colorss[1], ncol(AD))
      for(i in 2:length(table(as.numeric(AD[2, ])))){
        Cell_type[AD[2, ]==names(table(as.numeric(AD[2, ])))[i]]=colorss[i]
      }
      
      REPLAS <- c("Zero", "Mean")
      for (g in 1:length(REPLAS)) {
        if(g==1){
          AA <- SEE[ ,-1]
          AA[is.na(AA)] = 0 # Convert NAs to "0" and "AA" ready for tsne
          tsne <- Rtsne(t(AA), dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
          tiff(paste(WWW, CaseN, "_ZeroCell70.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
          plot(tsne$Y, main= paste(CaseN, "_Celltype_ZeroCDF", sep=""), cex.main = 0.8, col= Cell_type, pch = 16, cex = 0.4)
          par(xpd=T)
          mtext(TXT, side = 3)
          legend("topright",legend=c(names(table(as.numeric(AD[2,])))), fill=c(colorss[1:length(table(as.numeric(AD[2, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
          dev.off()
          tiff(paste(WWW, CaseN, "_ZeroPat70.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
          plot(tsne$Y, main= paste(CaseN, "_PatientID_ZeroCDF", sep=""), cex.main = 0.8, col= Patient_colors, pch = 16, cex = 0.4)
          par(xpd=T)
          legend("topright",legend=c(names(table(as.numeric(AD[1,])))), fill=c(colors[1:length(table(as.numeric(AD[1, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
          dev.off()
        }else{
          AA <- SEE[ ,-1]
          for (i in 1:nrow(AA)) {
            AA[ i, is.na(AA[i, ])] <- mean(na.omit(as.numeric(AA[i, ])))
          } # Convert NAs to "0" for "Mean version"
          tsne <- Rtsne(t(AA), dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
          tiff(paste(WWW, CaseN, "_MeanCell70.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
          plot(tsne$Y, main= paste(CaseN, "_Celltype_MeanCDF", sep=""), cex.main = 0.8, col= Cell_type, pch = 16, cex = 0.4)
          par(xpd=T)
          mtext(TXT, side = 3)
          legend("topright",legend=c(names(table(as.numeric(AD[2,])))), fill=c(colorss[1:length(table(as.numeric(AD[2, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
          dev.off()
          tiff(paste(WWW, CaseN, "_MeanPat70.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
          plot(tsne$Y, main= paste(CaseN, "_PatientID_MeanCDF", sep=""), cex.main = 0.8, col= Patient_colors, pch = 16, cex = 0.4)
          par(xpd=T)
          legend("topright",legend=c(names(table(as.numeric(AD[1,])))), fill=c(colors[1:length(table(as.numeric(AD[1, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
          dev.off()
        }
      }
    }
  }else{
    CaseN <- CaseN[R]
    STAT <- c("Normal", "Cancerous")
    Filepath3 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", STAT[2], "/", sep = "")
    WWW <- Filepath3
    TARFILE <- list.files(WWW, pattern = "RAW_CDF_DEL")
    SEE <- read.csv(paste(WWW, TARFILE, sep = ""), header = T, stringsAsFactors = F)
    SEE <- SEE[apply(SEE, 1, function(x) sum(is.na(x)) < (ncol(SEE)*(0.7))), ]
    TARFILE2 <- list.files(WWW, pattern = "AnnoRows")
    AD <- read.csv(paste(WWW, TARFILE2, sep = ""), header = T, stringsAsFactors = F)
    TXT <- AD$Sample[2]
    AD <- AD[ ,-1]
    
    # Patient Identity Colors labeling
    Patient_colors=rep(colors[1], ncol(AD))
    for(i in 2:length(table(as.character(AD[1, ])))){
      Patient_colors[AD[1, ]==names(table(as.character(AD[1, ])))[i]]=colors[i]
    }
    
    # Cell type Colors labeling
    Cell_type <- rep(colorss[1], ncol(AD))
    for(i in 2:length(table(as.numeric(AD[2, ])))){
      Cell_type[AD[2, ]==names(table(as.numeric(AD[2, ])))[i]]=colorss[i]
    }
    REPLAS <- c("Zero", "Mean")
    for (g in 1:length(REPLAS)){
      if(g==1){
        AA <- SEE[ ,-1]
        AA[is.na(AA)] = 0 # Convert NAs to "0" and "AA" ready for tsne
        tsne <- Rtsne(t(AA), dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
        tiff(paste(WWW, CaseN, "_ZeroCell70.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
        plot(tsne$Y, main= paste(CaseN, "_Celltype_ZeroCDF", sep=""), cex.main = 0.8, col= Cell_type, pch = 16, cex = 0.4)
        par(xpd=T)
        mtext(TXT, side = 3)
        legend("topright",legend=c(names(table(as.numeric(AD[2,])))), fill=c(colorss[1:length(table(as.numeric(AD[2, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
        dev.off()
        tiff(paste(WWW, CaseN, "_ZeroPat70.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
        plot(tsne$Y, main= paste(CaseN, "_PatientID_ZeroCDF", sep=""), cex.main = 0.8, col= Patient_colors, pch = 16, cex = 0.4)
        par(xpd=T)
        legend("topright",legend=c(names(table(as.character(AD[1,])))), fill=c(colors[1:length(table(as.character(AD[1, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
        dev.off()
      }else{
        AA <- SEE[ ,-1]
        for (i in 1:nrow(AA)) {
          AA[ i, is.na(AA[i, ])] <- mean(na.omit(as.numeric(AA[i, ])))
        } # Convert NAs to "0" for "Mean version"
        tsne <- Rtsne(t(AA), dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
        tiff(paste(WWW, CaseN, "_MeanCell70.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
        plot(tsne$Y, main= paste(CaseN, "_Celltype_MeanCDF", sep=""), cex.main = 0.8, col= Cell_type, pch = 16, cex = 0.4)
        par(xpd=T)
        mtext(TXT, side = 3)
        legend("topright",legend=c(names(table(as.numeric(AD[2,])))), fill=c(colorss[1:length(table(as.numeric(AD[2, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
        dev.off()
        tiff(paste(WWW, CaseN, "_MeanPat70.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
        plot(tsne$Y, main= paste(CaseN, "_PatientID_MeanCDF", sep=""), cex.main = 0.8, col= Patient_colors, pch = 16, cex = 0.4)
        par(xpd=T)
        legend("topright",legend=c(names(table(as.character(AD[1,])))), fill=c(colors[1:length(table(as.character(AD[1, ])))]), border=FALSE, bty="n", y.intersp = 1, cex=0.7)  
        dev.off()
      }
    }
  }
  print(CaseN)
  CaseN <- c("Astrocytoma","GSE70630_OG","GSE72056_Melanoma","GSE75688_Breast","GSE81861_CRC","GSE103322_HNSCC","GSE81383_Melanoma","E_MTAB_6149_NSCLC","GSE76312_CML")
}