CBS <- read.csv("C:/Users/tonyxp/Desktop/geneset_CBSofBC/GSE75688_BC_PatientID FIX Go-Scoring-of Clust&Whole_p-val.csv", header = TRUE) 

CB <- scale(CBS[ , 4:14])
CB <- cbind(CB, rowSums(CB))
CB <- cbind(CBS[ , 1:3], CB)
CB <- CB[order(CB$V12), ]

bina <- CBS[ , 4:14]
binana <- bina <= 0.2
binana[binana == "TRUE"] <- 1
binana[binana == "FALSE"] <- 0
binana <- cbind(binana, rowSums(binana))
binana <- cbind(CBS[ , 1:3], binana)
binana <- binana[rev(order(binana$V12)), ]

CadID <- CB[CB$V12 <= -13, 2]
CadNum <- CB[CB$V12 <= -13, 3]

CadID <- gsub(":", "_", CadID)


COL.OMIT <- function(u){
  any(colSums(binana[u, , drop = FALSE]) == 0)
}

g=function(z){
  AA <- binana[ ,colSums(binana[z, , drop = FALSE]) == 0, drop = FALSE]
  CC <- AA[(rowSums(AA) == max(rowSums(AA))), , drop=FALSE]
  rownames(CC[!duplicated(CC), ,drop = FALSE])

  #CC <- AA[(rowSums(AA) == max(rowSums(AA))), , drop=FALSE]
  #PaT[!duplicated(PaT), ,drop = FALSE]
  #which(rowSums(AA)==max(rowSums(AA)))
}

recursive.Gf <- function(x){
  
  if(any(unlist(lapply(x, function(z) COL.OMIT(z)))) == FALSE){
    
    return(x)
  }
  else{
    temp=lapply(x, function(z) g(z))
    xx=list()
    nn=1
    for(i in 1:length(x)){
      for(j in 1:length(temp[[i]])){
        xx[[nn]]=c(x[[i]],temp[[i]][j])
        nn=nn+1
      }
    }
    #print(xx)
    #print(length(x))
    return(recursive.Gf(xx))
  }
}

YY=recursive.Gf(PosL)

CoT <- NULL
Lyer <- NULL
CuO <- NULL
for(i in 2){
  binana <- bina <= (i*0.05)
  binana[binana == "TRUE"] <- 1
  binana[binana == "FALSE"] <- 0
  Q <- sum(colSums(binana) == 722)

  binana <- binana[ ,!apply(binana,2,function(x) all(x==x[1]))]
  rownames(binana) <- 1:nrow(binana)
  PaT <- binana[(rowSums(binana) == max(rowSums(binana))), , drop=FALSE]
  PaT <- PaT[!duplicated(PaT), ,drop = FALSE]
  PosL <- list()
  for(k in 1:length(rownames(PaT))){
    PosL[[k]] <- rownames(PaT)[k]
  }
  YY <- recursive.Gf(PosL)
  CoT <- append(CoT, (ncol(PaT)+Q), after = length(CoT))
  Lyer <- append(Lyer, length(YY[[1]]), after = length(Lyer))
  CuO <- append(CuO, i*0.05, after = length(CuO))
  print(i)
}

write.csv(CB, file = paste("C:/Users/tonyxp/Desktop/CBSofBC_P_zscore_rank.csv", sep=""), row.names = FALSE)
write.csv(binana, file = paste("C:/Users/tonyxp/Desktop/CBSofBC_P_binary_rank.csv", sep=""), row.names = FALSE)
