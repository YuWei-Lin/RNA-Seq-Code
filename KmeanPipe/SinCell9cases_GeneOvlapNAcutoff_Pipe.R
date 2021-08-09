### Pipeline for merging the number of remaining genes after various cutoffs of NA levels in all 9 cases
CaseN <- c("Astrocytoma","GSE70630_OG","GSE72056_Melanoma","GSE75688_Breast","GSE81861_CRC","GSE103322_HNSCC","GSE81383_Melanoma","E_MTAB_6149_NSCLC","GSE76312_CML")

for (R in 2:length(CaseN)) {
  if(R<=6){
    CaseN <- CaseN[R]
    STAT <- c("Normal", "Cancerous")
    Filepath1 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", sep = "")
    Filepath2 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", STAT[1], "/", sep = "")
    Filepath3 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", STAT[2], "/", sep = "")
    for (WWW in c(Filepath1, Filepath2, Filepath3)) {
      TARFILE <- list.files(WWW, pattern = "RAW_CDF_DEL")
      SEE <- read.csv(paste(WWW, TARFILE, sep = ""), header = T, stringsAsFactors = F)
      if(WWW==Filepath1){
        GeneSetF1 <- NULL
        GeneNumF1 <- NULL
        for (CUT in 1:10) {
          SEE1 <- SEE[apply(SEE, 1, function(x) sum(is.na(x)) < (ncol(SEE)*(0.1*CUT))), ]
          GeneSetF1[[CUT]] <- SEE1$Gene
          GeneNumF1[CUT] <- length(SEE1$Gene)
        }
        print(WWW)
      }
      if(WWW==Filepath2){
        GeneSetF2 <- NULL
        GeneNumF2 <- NULL
        for (CUT in 1:10) {
          SEE2 <- SEE[apply(SEE, 1, function(x) sum(is.na(x)) < (ncol(SEE)*(0.1*CUT))), ]
          GeneSetF2[[CUT]] <- SEE2$Gene
          GeneNumF2[CUT] <- length(SEE2$Gene)
        }
        print(WWW)
      }
      if(WWW==Filepath3){
        GeneSetF3 <- NULL
        GeneNumF3 <- NULL
        for (CUT in 1:10) {
          SEE3 <- SEE[apply(SEE, 1, function(x) sum(is.na(x)) < (ncol(SEE)*(0.1*CUT))), ]
          GeneSetF3[[CUT]] <- SEE3$Gene
          GeneNumF3[CUT] <- length(SEE3$Gene)
        }
        print(WWW)
      }
    }
    df <- matrix(nrow = 10, ncol = 3)
    for (i in 1:10) {
      for (j in 1:3) {
        if(j==1){
          df[i,j] <- length(intersect(GeneSetF1[[i]], GeneSetF2[[i]]))
        }
        if(j==2){
          df[i,j] <- length(intersect(GeneSetF1[[i]], GeneSetF3[[i]]))
        }
        if(j==3){
          df[i,j] <- length(intersect(GeneSetF2[[i]], GeneSetF3[[i]]))
        }
      }
    }
    df <- cbind((1:10)*10, GeneNumF1, GeneNumF2, GeneNumF3, df)  
    colnames(df) <- c("Remove%NA","MixGenesNum","NorGenesNum","CanGenesNum","MixNorIntSecNum","MixCanIntSecNum","NorCanIntSecNum")
    write.csv(df, file = paste(Filepath1, CaseN,"_GeneOvpinNAcutoff.csv", sep = ""), row.names = F)
  }else{
    CaseN <- CaseN[R]
    STAT <- c("Normal", "Cancerous")
    Filepath1 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", sep = "")
    Filepath3 <- paste("D:/SC Cases Completed/", CaseN, "_DATA/", STAT[2], "/", sep = "")
    TARFILE <- list.files(WWW, pattern = "RAW_CDF_DEL")
    SEE <- read.csv(paste(WWW, TARFILE, sep = ""), header = T, stringsAsFactors = F)
    WWW <- Filepath3
    if(WWW==Filepath3){
      GeneSetF4 <- NULL
      GeneNumF4 <- NULL
      for (CUT in 1:10) {
        SEE4 <- SEE[apply(SEE, 1, function(x) sum(is.na(x)) < (ncol(SEE)*(0.1*CUT))), ]
        GeneSetF4[[CUT]] <- SEE4$Gene
        GeneNumF4[CUT] <- length(SEE4$Gene)
      }
      print(WWW)
    }
    df <- cbind((1:10)*10, GeneNumF4)  
    colnames(df) <- c("Remove%NA","CanGenesNum")
    write.csv(df, file = paste(Filepath1, CaseN,"_GeneOvpinNAcutoff.csv", sep = ""), row.names = F)
  }
  print(CaseN)
  CaseN <- c("Astrocytoma","GSE70630_OG","GSE72056_Melanoma","GSE75688_Breast","GSE81861_CRC","GSE103322_HNSCC","GSE81383_Melanoma","E_MTAB_6149_NSCLC","GSE76312_CML")
}

### Pipeline for evaluating how many genes remain after using cutoff of NA considering gene sparsity for each patient in all 9 cases







### Tumor Identity Colors labeling
tumors_colors=rep(colors[1], 295)
for(i in 2:length(table(as.numeric(SEE[1,])))){
  tumors_colors[SEE[1, ]==names(table(as.numeric(SEE[1,])))[i]]=colors[i]
}
tumors_colors <- as.matrix(tumors_colors) # Legend Attached Form

colorss = c("#ffe119", "#f58231", "#42d4f4", "#469990", "#800000", "#aaffc3")
Cell_type <- rep(colorss[1], 295)
for(i in 2:length(table(as.numeric(SEE[2,])))){
  Cell_type[SEE[2, ]==names(table(as.numeric(SEE[2,])))[i]]=colorss[i]
}
Cell_type <- as.matrix(Cell_type) # Legend Attached Form
clab <- cbind(Cell_type, tumors_colors)
colnames(clab)=c("Cell_Type", "Tumor_ID") # Legend Label Name 
Anno <- rownames(SEE)[3]
### Duplicate input and convert NAs to "0" for "Zero version"
AA <- SEE[3:nrow(SEE), ]
AA[is.na(AA)]=0 # "AA" ready for Kmeans
CC <- SEE[3:nrow(SEE), ]
### Duplicate input and convert NAs to "0" for "Mean version"
AA <- SEE[3:nrow(SEE), ]
for (i in 1:nrow(AA)) {
  AA[ i, is.na(AA[i, ])] <- mean(na.omit(as.numeric(AA[i, ])))
}
CC <- SEE[3:nrow(SEE), ]