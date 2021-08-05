# Import Significant GO Terms 
GOIDList <- read.csv("C:/Users/User/Desktop/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)

ClusLi <- list.files("C:/Users/User/Desktop/GSE75688_Breast_ClusGenes", pattern = "ClusGenes_K")
ClusLi <- ClusLi[order(as.numeric(sub("K", "", unlist(strsplit(ClusLi, "_"))[(1:390)%%5==4])), as.numeric(sub( "csv","", sub("C", "", unlist(strsplit(ClusLi, "_"))[(1:390)%%5==0]))))]
# Create value interval check function
in_interval <- function(x, interval){ 
  for (o in 1:(length(interval)-1)){
    if((interval[o]+1<=x & x<=interval[o+1]) == TRUE){
      return(o)
    }
  }
}
# Type I Layout - Summary Table for individual K
p_time <- proc.time()
for(k in 1:12){
  ID <- grep(paste("K", k, sep = ""), ClusLi)
  if(k==1){
    ID <- grep(paste("K", k, "_C1", sep = ""), ClusLi)
  }
  RW <- NULL
  Data <- NULL
  Intval <- NULL
  for (y in 1:length(ID)) {
    Data <- read.csv(paste("C:/Users/User/Desktop/GSE75688_Breast_ClusGenes/", ClusLi[ID[y]], sep = ""), header = TRUE, stringsAsFactors = FALSE)
    RW <- rbind(RW, Data)
    Intval[y] <- nrow(RW)
  }
  Intval <- c(0, Intval)
  WIN <- data.frame(matrix(NA, nrow = length(unique(RW$GO_Term)), ncol = k+1))
  GoTerms <- names(sort(table(RW$GO_Term), decreasing = T))
  WIN[ ,1] <- GoTerms
  for (i in 1:nrow(WIN)) {
    TAR <- NULL
    UU <- rep(NA, k)
    for (j in 1:nrow(RW[RW$GO_Term==GoTerms[i], ])) {
      TAR <- c(TAR, in_interval(as.numeric(rownames(RW[RW$GO_Term==GoTerms[i], ])[j]), Intval))
      UU[TAR[j]] <- RW[RW$GO_Term==GoTerms[i], ]$FDR[j]
    }
    WIN[i, 2:ncol(WIN)] <- UU
  }
  colnames(WIN) <- c("Enrichment GOTerms", paste("FDR_K", k, "C", 1:k, sep = ""))
  write.table(WIN, file = paste("GSE75688_Breast_ClusGenes_K", k, "_Summary.csv", sep = ""), sep = ",", row.names = FALSE)
}
t_time <- proc.time()-p_time
print(t_time)
# Type II Layout - Concise GOTerm Table across K
Tracing <- list.files("C:/Users/User/Desktop/GSE75688_Breast_ClusGenes/BC_Summary")
Tracing <- Tracing[order(as.numeric(sub("K", "", unlist(strsplit(Tracing, "_"))[(1:60)%%5==4])))]

ALLTe <- NULL
ColN <- NULL
for(k in 1:12){
  Data <- read.csv(paste("C:/Users/User/Desktop/GSE75688_Breast_ClusGenes/BC_Summary/", Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ALLTe <- c(ALLTe, Data$Enrichment.GOTerms)
  ColN <- c(ColN, paste("FDR_K", k, "C", 1:k, sep = ""))
}
WIN <- data.frame(matrix(NA, nrow = length(unique(ALLTe)), ncol = 79))
GoTerms <- names(sort(table(ALLTe), decreasing = T))
WIN[ ,1] <- GoTerms
FRE <- NULL
for(i in 1:length(GoTerms)){
  LM <- NULL
  for(k in 1:12){
    Data <- read.csv(paste("C:/Users/User/Desktop/GSE75688_Breast_ClusGenes/BC_Summary/", Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
    if(names(sort(table(ALLTe), decreasing = T))[i]%in%Data$Enrichment.GOTerms == TRUE){
      LM <- c(LM, Data[ Data$Enrichment.GOTerms == names(sort(table(ALLTe), decreasing = T))[i], 2:(k+1)])
    }
    else{
      LM <- c(LM, rep(NA, length(2:(k+1))))
    }
  }
  FRE <- c(FRE, (78-sum(is.na(LM))))
  WIN[i, 2:ncol(WIN)] <- LM
  print(i)
}

WIN <- cbind(WIN[ , 1], FRE, WIN[ , -1])
colnames(WIN) <- c("Enrichment GOTerms", "Frequency", ColN)
WIN <- WIN[order(WIN$Frequency, decreasing = T), ]
write.table(WIN, file = paste("GSE75688_Breast_KClusters_GOWhole.csv", sep = ""), sep = ",", row.names = FALSE)
# Type III Layout - GOTerm Tracing Table across K
### Run rbind.na & cbind.na function files if neccessary
RPmatrix <- data.frame()
for (m in 1:nrow(WIN)) {
  ST <- NULL
  for (n in 3:80) {
    if(is.na(WIN[ m, n]) == FALSE){
      ST <- c(ST, sub("FDR_", "", colnames(WIN)[n]))
    }
  }
  RPmatrix <- rbind.na(RPmatrix, ST)
  print(m)
}
RPmatrix <- cbind(WIN[ ,1:2], RPmatrix)
write.table(RPmatrix, file = paste("GSE75688_Breast_KClusters_Tracing.csv", sep = ""), sep = ",", row.names = FALSE)
