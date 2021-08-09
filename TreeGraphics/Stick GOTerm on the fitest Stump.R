### A lot of  reference is needed in R code-"Making Cluster Tree.R"; Like stump.L, AA_form...
### Import and modify needed tables
T_form <- AA_form[ , -23]
colnames(T_form)[23] <- "stump_23(1)"
JKK <- read.csv("C:/Users/User/Desktop/stumps relations.csv", stringsAsFactors = F, header = T)
### Select the most significant GOterms
TAR.G <- NULL
for (i in 2:nrow(T_form)) {
  stuNames <- which(T_form[i, ] != "N")
  STUV <- NULL # this is the average FDR values regarding each stumps in a GOterm; We can deem it as overall scale of FDR within a GOterm 
  for (y in 1:length(stuNames)) {
    stuDA <- DATA[ DATA$Enrichment.GOTerms == rownames(T_form[i, ]), colnames(DATA)%in%paste("FDR_", stumps.L[[stuNames[y]]], sep = "")]
    STUV <- append( STUV, sum(stuDA)/length(stumps.L[[stuNames[y]]]), after = length(STUV))
  }
  if(length(STUV)==1&any(STUV >= 10)){
    TAR.G <- append(TAR.G, rownames(T_form[i, ]))
  }
  if(length(STUV)>=2){
    if(any(STUV >= 10)&any(diff(combn(STUV, 2)) >= 2)){
      TAR.G <- append(TAR.G, rownames(T_form[i, ]))
    }
  }
}
#WQ <- cbind(colnames(T_form)[stuNames], STUV)

### Determine the fitest stumps for every significant GOterm
LOC <- NULL
for (i in 1:length(TAR.G)) {
  stuNames <- which(T_form[rownames(T_form)==TAR.G[i], ] != "N")
  STUV <- NULL # this is the average FDR values regarding each stumps in a GOterm; We can deem it as overall scale of FDR within a GOterm 
  for (y in 1:length(stuNames)) {
    stuDA <- DATA[ DATA$Enrichment.GOTerms == TAR.G[i], colnames(DATA)%in%paste("FDR_", stumps.L[[stuNames[y]]], sep = "")]
    STUV <- append( STUV, sum(stuDA)/length(stumps.L[[stuNames[y]]]), after = length(STUV))
  }
  WQ <- as.data.frame(cbind(colnames(T_form)[stuNames], STUV), stringsAsFactors = F)
  WQ[ ,2] <- as.numeric(WQ$STUV)
  WQ <- WQ[order(WQ[,2], decreasing = T), ]
  if(is.na(var(STUV))==T){
    LOC <- append(LOC, WQ$V1, after = length(LOC))
  }
  if(length(STUV)==2){
    LOC <- append(LOC, WQ[WQ$STUV==max(STUV), ]$V1, after = length(LOC))
  }
  else{
    if(var(sort(STUV, decreasing = T)[1:2])/var(STUV) >= 0.2){
      LOC <- append(LOC, WQ[WQ$STUV==max(STUV), ]$V1, after = length(LOC))
    }
    else{
      S1 <- WQ[1:3, ]$V1
      if(all(S1%in%JKK$X == F)){
        LOC <- append(LOC, WQ[WQ$STUV==max(STUV), ]$V1, after = length(LOC))
      }
      else{
        S1 <- JKK$X[min(which(JKK$X%in%S1))]
        S2 <- JKK[JKK$X==S1,2]
        S3 <- JKK[JKK$X==S1,3]
        if((is.na(S3)==T)&(S2%in%WQ$V1 == F)){
          LOC <- append(LOC, S1, after = length(LOC))
        }
        if((is.na(S3)==T)&(S2%in%WQ$V1 == T)){
          if((WQ[WQ$V1==S2, 2]-WQ[WQ$V1==S1, 2]) >= 2){
            LOC <- append(LOC, S2, after = length(LOC))
          }
          else{
            LOC <- append(LOC, S1, after = length(LOC))
          }
        }
      }
    }
  }
}
