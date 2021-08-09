# Import Significant GO Terms 
GOIDList <- read.csv("C:/Users/User/Desktop/HALL_MARKs/GSA_HALLMark_GENES_Category.csv", header = TRUE, stringsAsFactors = FALSE)
ClusLi <- list.files("C:/Users/User/Desktop/HALL_MARKs", pattern = "Breast")
ClusLi <- ClusLi[order(as.numeric(sub("K", "", unlist(strsplit(ClusLi, "_"))[(1:385)%%5==4])), as.numeric(sub( "csv","", sub("C", "", unlist(strsplit(ClusLi, "_"))[(1:385)%%5==0]))))]
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
for(k in 2:12){
  ID <- grep(paste("K", k, sep = ""), ClusLi)
  RW <- NULL
  Data <- NULL
  Intval <- NULL
  for (y in 1:length(ID)) {
    Data <- read.csv(paste("C:/Users/User/Desktop/HALL_MARKs/", ClusLi[ID[y]], sep = ""), header = TRUE, stringsAsFactors = FALSE)
    RW <- rbind(RW, Data)
    Intval[y] <- nrow(RW)
  }
  Intval <- c(0, Intval)
  WIN <- data.frame(matrix(NA, nrow = length(unique(RW$HallMark_Term)), ncol = k+1))
  GoTerms <- names(sort(table(RW$HallMark_Term), decreasing = T))
  WIN[ ,1] <- GoTerms
  for (i in 1:nrow(WIN)) {
    TAR <- NULL
    UU <- rep(NA, k)
    for (j in 1:nrow(RW[RW$HallMark_Term==GoTerms[i], ])) {
      TAR <- c(TAR, in_interval(as.numeric(rownames(RW[RW$HallMark_Term==GoTerms[i], ])[j]), Intval))
      UU[TAR[j]] <- RW[RW$HallMark_Term==GoTerms[i], ]$FDR[j]
    }
    WIN[i, 2:ncol(WIN)] <- UU
  }
  colnames(WIN) <- c("Enrichment HMTerms", paste("FDR_K", k, "C", 1:k, sep = ""))
  write.table(WIN, file = paste("GSE75688_Breast_HALLMARKs_K", k, "_Summary.csv", sep = ""), sep = ",", row.names = FALSE)
}
t_time <- proc.time()-p_time
print(t_time)
# Type II Layout - Concise GOTerm Table across K
Tracing <- list.files("C:/Users/User/Desktop/GSE75688_Breast_HALLMARKs_Summary")
Tracing <- Tracing[order(as.numeric(sub("K", "", unlist(strsplit(Tracing, "_"))[(1:55)%%5==4])))]

ALLTe <- NULL
ColN <- NULL
for(k in 1:11){
  Data <- read.csv(paste("C:/Users/User/Desktop/GSE75688_Breast_HALLMARKs_Summary/", Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
  ALLTe <- c(ALLTe, Data$Enrichment.HMTerms)
  ColN <- c(ColN, paste("FDR_K", k+1, "C", 1:(k+1), sep = ""))
}
WIN <- data.frame(matrix(NA, nrow = length(unique(ALLTe)), ncol = 78))
GoTerms <- names(sort(table(ALLTe), decreasing = T))
WIN[ ,1] <- GoTerms
FRE <- NULL
for(i in 1:length(GoTerms)){
  LM <- NULL
  for(k in 1:11){
    Data <- read.csv(paste("C:/Users/User/Desktop/GSE75688_Breast_HALLMARKs_Summary/", Tracing[k], sep = ""), header = TRUE, stringsAsFactors = FALSE)
    if(names(sort(table(ALLTe), decreasing = T))[i]%in%Data$Enrichment.HMTerms == TRUE){
      LM <- c(LM, Data[ Data$Enrichment.HMTerms == names(sort(table(ALLTe), decreasing = T))[i], 2:(k+2)])
    }
    else{
      LM <- c(LM, rep(NA, length(2:(k+2))))
    }
  }
  FRE <- c(FRE, (77-sum(is.na(LM))))
  WIN[i, 2:ncol(WIN)] <- LM
  print(i)
}

WIN <- cbind(WIN[ , 1], FRE, WIN[ , -1])
colnames(WIN) <- c("Enrichment HMTerms", "Frequency", ColN)
WIN <- WIN[order(WIN$Frequency, decreasing = T), ]
write.table(WIN, file = paste("GSE75688_Breast_KClus_HALLMARKs.csv", sep = ""), sep = ",", row.names = FALSE)
# Type III Layout - GOTerm Tracing Table across K
### Run rbind.na & cbind.na function files if neccessary
RPmatrix <- data.frame()
for (m in 1:nrow(WIN)) {
  ST <- NULL
  for (n in 3:79) {
    if(is.na(WIN[ m, n]) == FALSE){
      ST <- c(ST, sub("FDR_", "", colnames(WIN)[n]))
    }
  }
  RPmatrix <- rbind.na(RPmatrix, ST)
  print(m)
}
RPmatrix <- cbind(WIN[ ,1:2], RPmatrix)
write.table(RPmatrix, file = paste("GSE75688_Breast_KClus_HALLMARKs_Tracing.csv", sep = ""), sep = ",", row.names = FALSE)
