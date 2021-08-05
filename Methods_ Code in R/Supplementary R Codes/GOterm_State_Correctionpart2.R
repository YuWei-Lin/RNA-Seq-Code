### Location address Setting
CaseN <- "GSE75688_Breast"
CaseN <- "GSE72056_Melanoma"
CaseN <- "GSE70630_OG"
CaseN <- "Astrocytoma"
CaseN <- "GSE81383_Melanoma"
CaseN <- "GSE103322_HNSCC"
Filepath <- paste("D:/ALL_FINAL_RESULTS/Manual_Check_States/", CaseN, "/FromMatlabOut/", CaseN, sep = "")
### Read Data
GOterms <- read.csv(paste("D:/ALL_FINAL_RESULTS/Manual_Check_States/", CaseN, "/FromMatlabOut/ALL_GOterms.csv", sep = ""), header = F, stringsAsFactors = F)
GOdex.M <- read.csv(paste(Filepath, "_GOindex_M.csv", sep = ""), header = F, stringsAsFactors = F)
GOstate.M <- read.csv(paste(Filepath, "_state_M.csv", sep = ""), header = F, stringsAsFactors = F)
GOFDRs.M <- read.csv(paste(Filepath, "_clusFDRs_M.csv", sep = ""), header = F, stringsAsFactors = F)

GOterms <- GOterms$V1
GOdex.M <- as.numeric(GOdex.M)
GOFDRs.M <- GOFDRs.M[GOdex.M, ]
### Process states with mean value replacement
BVC <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/Inputs/", CaseN,"_ClusGene_Membership.csv", sep = ""), header = T, stringsAsFactors = F)
for (i in 1:length(GOdex.M)) {
  GOname <- GOterms[as.numeric(GOdex.M[i])]
  Chefile <- rbind(colnames(BVC)[2:ncol(BVC)], 1:78, round(-log10(GOFDRs.M[i, ]), digits = 2), GOstate.M[i, ], GOstate.M[i, ])
  Chefile <- t(Chefile)
  colnames(Chefile) <- c("Cluster", "ClustNumber", "neglog10FDRs", "Predicted State", "Modified State")
  write.table(GOname, file = paste("D:/ALL_FINAL_RESULTS/Manual_Check_States/", CaseN, "/", CaseN, "_UsingMean/", GOname, "_UsingMean_stateChange.csv", sep = ""), sep = ",", row.names = F, col.names = F)
  write.table(Chefile, file = paste("D:/ALL_FINAL_RESULTS/Manual_Check_States/", CaseN, "/", CaseN, "_UsingMean/", GOname, "_UsingMean_stateChange.csv", sep = ""), sep = ",", row.names = F, col.names = T, append = T)  
}

### Enumerate Routes of the childview tree
### 5.  Overlap between clusters at adjacent levels of K.(overlaping sizes are divided by child size)
BVC <- read.csv(paste("D:/ALL_FINAL_RESULTS/", CaseN, "_DATA/K.Means Clustering_UsingMean/Inputs/", CaseN,"_ClusGene_Membership.csv", sep = ""), header = T, stringsAsFactors = F)
CluRatio <- as.data.frame(matrix(NA, nrow = 78, ncol = 78))
rownames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
colnames(CluRatio) <- colnames(BVC)[2:ncol(BVC)]
for (k in 1:11) {
  Siz <- k+(k+1)
  if(k==1){
    for (m in 1:3) {
      for (n in 1:3) {
        TOP <- BVC$Gene[BVC[,m+1] == TRUE]
        BOT <- BVC$Gene[BVC[,n+1] == TRUE]
        CluRatio[m, n] <- length(intersect(TOP, BOT))/length(BOT)
      }
    }
  }else{
    for (m in (sum(1:(k-1))+1):((sum(1:(k-1))+1)+Siz-1)) {
      for (n in (sum(1:(k-1))+1):((sum(1:(k-1))+1)+Siz-1)) {
        TOP <- BVC$Gene[BVC[,m+1] == TRUE]
        BOT <- BVC$Gene[BVC[,n+1] == TRUE]
        CluRatio[m, n] <- length(intersect(TOP, BOT))/length(BOT)
      }
    }
  }
}
CluRatio[lower.tri(CluRatio, diag = T)] <- NA

### Find out all possible route in the tree
CluRatio[is.na(CluRatio)] <- 0
Curpath <- list(1)
for (k in 2:12) {
  for (y in 1:length(Curpath)) {
    for (n in 1:k) {
      if(CluRatio[Curpath[[y]][length(Curpath[[y]])], (sum(1:(k-1))+n)]>=0.2){
        Curpath[[length(Curpath)+1]] <- c(Curpath[[y]], (sum(1:(k-1))+n))
      }
    }
  }
}
for (i in 1:(length(Curpath)-1)) {
  for (j in (i+1):length(Curpath)) {
    if(all(Curpath[[i]]%in%Curpath[[j]])==T){
      Curpath[[i]] <- 0
      break
    }
  }
}
TF <- lapply(Curpath, function(x) length(x)>1)
Curpath <- Curpath[TF==T]

### "state 1" linage Check in every route and reasonably elongate the length of "state 1" linage under specifying conditions.
### Gap could be filled up or it could really shows the indepence between "state 1" linages
destFolder <- paste("D:/ALL_FINAL_RESULTS/Manual_Check_States/", CaseN, "/", CaseN, "_UsingMean", sep = "")
Allstateterms <- list.files(destFolder)
for (f in 1:length(Allstateterms)) {
  DAmat <- read.csv(paste(destFolder, "/", Allstateterms[f], sep = ""), stringsAsFactors = F, header = F)
  DAmat <- DAmat[3:nrow(DAmat), ]
  rownames(DAmat) <- DAmat$V2
  colnames(DAmat) <- c("Cluster", "ClustNumber", "neglog10FDRs", "Predicted State", "Modified State")
  DAmat[,2:ncol(DAmat)] <- apply(DAmat[,2:ncol(DAmat)], 2, as.numeric)
  QP <- 0.7
  if(all(DAmat$`Predicted State`==0)==T){
    next
  }else{
    Sum1Rank <- list()
    for (E in 1:length(Curpath)) {
      CAD <- DAmat$`Predicted State`[Curpath[[E]]]
      names(CAD) <- Curpath[[E]]
      CAD <- CAD[CAD!=0]
      if(length(CAD)==0){
        CAD <- NA
        Sum1Rank <- append(Sum1Rank, CAD, after = length(Sum1Rank))
      }else{
        Sum1Rank <- append(Sum1Rank, list(CAD), after = length(Sum1Rank))
      }
      names(Sum1Rank)[length(Sum1Rank)] <- E
    }
    Sum1Rank <- Sum1Rank[is.na(Sum1Rank)==F]
    RMindex <- NULL
    for (i in 1:length(unique(Sum1Rank))) {
      for (j in 1:length(unique(Sum1Rank))) {
        if(all(names(unique(Sum1Rank)[[i]])%in%names(unique(Sum1Rank)[[j]])==T)&(i!=j)){
          RMindex <- c(RMindex, i)
        }
      }
    }
    RMindex <- unique(RMindex)
    if(is.null(RMindex)==T){
      PAT <- unique(Sum1Rank)
    }else{
      PAT <- unique(Sum1Rank)[-(RMindex)]
    }
    Refine.Path <- function(x) {
      if(all(as.numeric(names(PAT[[P]]))%in%x == T)) {
        New.Curpath <- c(New.Curpath, x)
      }
    }
    Now.Curpath <- list()
    for (P in 1:length(PAT)) {
      New.Curpath <- NULL
      New.Curpath <- lapply(Curpath, Refine.Path)
      New.Curpath <- Filter(Negate(is.null), setNames(New.Curpath,seq_along(New.Curpath)))
      Now.Curpath <- append(Now.Curpath, list(New.Curpath), after=length(Now.Curpath))
    }
    Path.NEDex <- NULL
    for (H in 1:length(Now.Curpath)) {
      Path.NEDex <- c(Path.NEDex, as.numeric(names(Now.Curpath[[H]])))
    }
    Route.Sel <- Curpath[Path.NEDex]
    for (A in QP) {
      for (R in 1:length(Route.Sel)) {
        if(any(DAmat$`Predicted State`[Route.Sel[[R]]] != 0) == T){
          exist1path <- DAmat$`Predicted State`[Route.Sel[[R]]]
          if(sum(exist1path) > 1){
            temptable <- rbind(cumsum(rle(exist1path)$lengths), rle(exist1path)$lengths, rle(exist1path)$values)
            ClusLineOnly1 <- rep(0, temptable[1, ncol(temptable)])
            for (P in which(temptable[3,]==1)) {
              ClusLineOnly1 <- rep(0, temptable[1, ncol(temptable)])
              ClusLineOnly1[(temptable[ ,P][1]-temptable[ ,P][2]+1):temptable[ ,P][1]]<-1
              if(sum(ClusLineOnly1) == 1){
                next
              }
              if(length(rle(ClusLineOnly1)$values)==2&rle(ClusLineOnly1)$values[1]==1){
                Pos <- rle(ClusLineOnly1)$lengths[1]
                idx <- Route.Sel[[R]][which(ClusLineOnly1==1)]
                idx <- idx[1:(length(idx)-1)]
                thres.elong <- mean(DAmat[ idx, 3])*A
                while(DAmat[Route.Sel[[R]][Pos], 3] <= thres.elong) {
                  DAmat[Route.Sel[[R]][Pos], 5] <- 0
                  ClusLineOnly1 <- DAmat$`Modified State`[Route.Sel[[R]]]
                  idx <- Route.Sel[[R]][which(ClusLineOnly1==1)]
                  if(length(idx)==0){
                    break()
                  }
                  idx <- idx[1:(length(idx)-1)]
                  thres.elong <- mean(DAmat[ idx, 3])*A
                  Pos <- Pos-1
                  if(Pos == 0){
                    break
                  }
                }
              }else if(length(rle(ClusLineOnly1)$values)==2&rle(ClusLineOnly1)$values[2]==1){
                Pos <- sum(rle(ClusLineOnly1)$lengths[1:(length(rle(ClusLineOnly1)$lengths)-1)])+1
                idx <- Route.Sel[[R]][which(ClusLineOnly1==1)]
                idx <- idx[2:length(idx)]
                thres.elong <- mean(DAmat[ idx, 3])*A
                while(DAmat[Route.Sel[[R]][Pos], 3] <= thres.elong) {
                  DAmat[Route.Sel[[R]][Pos], 5] <- 0
                  ClusLineOnly1 <- DAmat$`Modified State`[Route.Sel[[R]]]
                  idx <- Route.Sel[[R]][which(ClusLineOnly1==1)]
                  if(length(idx)==0){
                    break()
                  }
                  idx <- idx[2:length(idx)]
                  thres.elong <- mean(DAmat[ idx, 3])*A
                  Pos <- Pos+1
                  if(Pos > length(Route.Sel[[R]])){
                    break
                  }
                }
              }else{
                if(all(exist1path==1)==F){
                  PosFront <- rle(ClusLineOnly1)$lengths[1]+1
                  idxFront <- Route.Sel[[R]][which(ClusLineOnly1==1)]
                  if(length(idxFront)!=1){
                    idxFront <- idxFront[2:length(idxFront)]
                  }
                  thres.elong.F <- mean(DAmat[ idxFront, 3])*A
                  PosEnd <- sum(rle(ClusLineOnly1)$lengths[1:2])
                  idxEnd <- Route.Sel[[R]][which(ClusLineOnly1==1)]
                  idxEnd <- idxEnd[1:(length(idxEnd)-1)]
                  thres.elong.E <- mean(DAmat[ idxEnd, 3])*A
                  while(DAmat[Route.Sel[[R]][PosFront], 3] <= thres.elong.F) {
                    DAmat[Route.Sel[[R]][PosFront], 5] <- 0
                    ClusLineOnly1 <- DAmat$`Modified State`[Route.Sel[[R]]]
                    idxFront <- Route.Sel[[R]][which(ClusLineOnly1==1)]
                    if(length(idxFront)==0){
                      break()
                    }
                    if(length(idxFront)!=1){
                      idxFront <- idxFront[2:length(idxFront)]
                    }
                    thres.elong.F <- mean(DAmat[ idxFront, 3])*A
                    PosFront <- PosFront+1
                    if(PosFront > length(Route.Sel[[R]])){
                      break
                    }
                  }
                  while(DAmat[Route.Sel[[R]][PosEnd], 3] <= thres.elong.E) {
                    DAmat[Route.Sel[[R]][PosEnd], 5] <- 0
                    ClusLineOnly1 <- DAmat$`Modified State`[Route.Sel[[R]]]
                    idxEnd <- Route.Sel[[R]][which(ClusLineOnly1==1)]
                    if(length(idxEnd)==0){
                      break()
                    }
                    idxEnd <- idxEnd[1:(length(idxEnd)-1)]
                    thres.elong.E <- mean(DAmat[ idxEnd, 3])*A
                    PosEnd <- PosEnd-1
                    if(PosEnd == 0){
                      break
                    }
                  }
                }
              }
            }
          }
          exist1path <- DAmat$`Modified State`[Route.Sel[[R]]]
          if(all(exist1path==1)==T){
            next
          }else{
            temptable <- rbind(cumsum(rle(exist1path)$lengths), rle(exist1path)$lengths, rle(exist1path)$values)
            
            for (Q in which(temptable[3,]==1)) {
              ClusLineOnly1 <- rep(0, temptable[1, ncol(temptable)])
              ClusLineOnly1[(temptable[ ,Q][1]-temptable[ ,Q][2]+1):temptable[ ,Q][1]]<-1
              thres.elong <- mean(DAmat[ Route.Sel[[R]][which(ClusLineOnly1==1)], 3])*A
              
              if(length(rle(ClusLineOnly1)$values)==2&rle(ClusLineOnly1)$values[1]==1){
                Pos <- rle(ClusLineOnly1)$lengths[1]+1
                while(DAmat[Route.Sel[[R]][Pos], 3] >= thres.elong) {
                  DAmat[Route.Sel[[R]][Pos], 5] <- 1
                  ClusLineOnly1 <- DAmat$`Modified State`[Route.Sel[[R]]]
                  thres.elong <- mean(DAmat[ Route.Sel[[R]][which(ClusLineOnly1==1)], 3])*A
                  Pos <- Pos+1
                  if(Pos > length(Route.Sel[[R]])){
                    break
                  }
                }
              }else if(length(rle(ClusLineOnly1)$values)==2&rle(ClusLineOnly1)$values[2]==1){
                Pos <- sum(rle(ClusLineOnly1)$lengths[1:(length(rle(ClusLineOnly1)$lengths)-1)])
                while(DAmat[Route.Sel[[R]][Pos], 3] >= thres.elong) {
                  DAmat[Route.Sel[[R]][Pos], 5] <- 1
                  ClusLineOnly1 <- DAmat$`Modified State`[Route.Sel[[R]]]
                  thres.elong <- mean(DAmat[ Route.Sel[[R]][which(ClusLineOnly1==1)], 3])*A
                  Pos <- Pos-1
                  if(Pos == 0){
                    break
                  }
                }
              }else{
                PosFront <- rle(ClusLineOnly1)$lengths[1]
                PosEnd <- sum(rle(ClusLineOnly1)$lengths[1:2])+1
                while(DAmat[Route.Sel[[R]][PosFront], 3] >= thres.elong) {
                  DAmat[Route.Sel[[R]][PosFront], 5] <- 1
                  ClusLineOnly1 <- DAmat$`Modified State`[Route.Sel[[R]]]
                  thres.elong <- mean(DAmat[ Route.Sel[[R]][which(ClusLineOnly1==1)], 3])*A
                  PosFront <- PosFront-1
                  if(PosFront == 0){
                    break
                  }
                }
                while(DAmat[Route.Sel[[R]][PosEnd], 3] >= thres.elong) {
                  DAmat[Route.Sel[[R]][PosEnd], 5] <- 1
                  ClusLineOnly1 <- DAmat$`Modified State`[Route.Sel[[R]]]
                  thres.elong <- mean(DAmat[ Route.Sel[[R]][which(ClusLineOnly1==1)], 3])*A
                  PosEnd <- PosEnd+1
                  if(PosEnd > length(Route.Sel[[R]])){
                    break
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  write.table(Allstateterms[f], file=paste(destFolder, "/", Allstateterms[f], sep = ""), sep = ",", row.names = F, col.names = F)
  write.table(DAmat, file=paste(destFolder, "/", Allstateterms[f], sep = ""), sep = ",", row.names = F, col.names = T, append = T)
}


### Collect all FDR scores to see the distribution (Random 50 GOterms files)
deFolder <- "C:/Users/User/Desktop/StateTestRcoding"
All <- list.files(deFolder)
AllFDRs <- NULL
for (y in 1:length(All)) {
  FILES <- read.csv(paste(deFolder, "/", All[y], sep = ""), stringsAsFactors = F, header = F)
  FILES <- FILES[3:nrow(FILES), ]
  rownames(FILES) <- FILES$V2
  colnames(FILES) <- c("Cluster", "ClustNumber", "neglog10FDRs", "Predicted State", "Modified State")
  FILES[,2:ncol(FILES)] <- apply(FILES[,2:ncol(FILES)], 2, as.numeric)
  AllFDRs <- c(AllFDRs, FILES[which(FILES[ ,4]==1),3])
}
cad <- NULL
Myfunction <- function(x) {
  if(19%in%x == T) {
    cad <- c(cad, x)
  }
}
cad <- lapply(Route.Sel, Myfunction)
cad <- Filter(Negate(is.null), setNames(cad,seq_along(cad)))



