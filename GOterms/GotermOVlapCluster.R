GOIDList <- read.csv("C:/Users/User/Desktop/GOID2GeneMSigDB.csv", header = TRUE, stringsAsFactors = FALSE)
Files <- list.files("C:/Users/User/Desktop/GSE75688 Breast_GotermOVlapCluster")
CoLink <- myDATA$GO_Term
for (i in 2:78) {
  myDATA <- read.csv(paste("C:/Users/User/Desktop/GSE75688 Breast_GotermOVlapCluster/", Files[i], sep = ""), stringsAsFactors = F, header = T)
  myDATA <- myDATA[order(myDATA$GO_Term), ]
  CoLink <- cbind(CoLink, myDATA$Overlapping)
  Clust <- strsplit(Files[i], "_")[[1]][4]
  colnames(CoLink)[i+1] <- gsub(".csv", "", Clust)
}
colnames(CoLink)[1:2] <- c("GOterms", "K01C01")
CoLinkk <- cbind(myDATA[,1:2], CoLink[,2:ncol(CoLink)])
write.table(CoLinkk, file = "C:/Users/User/Desktop/GSE75688 Breast_GotermOVlapCluster/GSE75688 Breast_GotermOVlapCluster.csv", sep = ",", col.names = T, row.names = F)
