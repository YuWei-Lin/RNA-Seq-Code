SC_N <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
SC_O <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
SC_N <- SC_N[order(SC_N$Clust_NSV), ]
SC_O <- SC_O[order(SC_O$Clust_NSV), ]

ID <- intersect(SC_N$GO.ID, SC_O$GO.ID)
N <- SC_N[SC_N$GO.ID%in%ID, colnames(SC_N) == "Clust_NSV"]
O <- SC_O[SC_O$GO.ID%in%ID, colnames(SC_O) == "Clust_NSV"]

plot(N, O, main = "Clust_NSV_of_BC",col = c("red", "green"), pch = 16)

MigDB <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
A <- MigDB[MigDB$ProcessName == "GO_SISTER_CHROMATID_SEGREGATION", ]
B <- MigDB[MigDB$ProcessName == "GO_CHROMOSOME_SEGREGATION", ]
C <- MigDB[MigDB$ProcessName == "GO_DNA_REPLICATION", ]

AB <- intersect(A$Gene, B$Gene)
ABC <- intersect(AB, C$Gene)

GOIDName <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
GOIDList <- read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)

NN <- GOIDName[match(GOIDList$GOID, GOIDName$GOID), 2]
#There are 36 missing matches. 
GOIDList$TERMS <- NN
XX <- GOIDList[!(is.na(NN)), ]
write.table(XX, file = "GOIDNAME_Combined.csv", sep = ",", row.names = FALSE)

KK <- SC_O$GO_Name
YY <- SC_N$GOName
MA <- NULL
for(i in 1:nrow(SC_O)){
  for(j in 1:nrow(SC_N)){
    OV <- length(intersect(XX[XX$TERMS == KK[i], 1], MigDB[MigDB$ProcessName == YY[j], 1]))/length(XX[XX$TERMS == KK[i], 1])
    if(OV >= 0.9){
      MA <- rbind(MA, c(i, KK[i], YY[j], j))
    }
  }
  print(i)
}

colnames(MA) <- c("Old_Ranking", "Old_Name", "New_Name", "New_Ranking")
MA <- as.data.frame(MA, stringsAsFactors = FALSE)
write.table(MA, file = "GOName_Translation.05_Full.csv", sep = ",", row.names = FALSE)

PF <- MA[rownames(unique(MA[ ,1:2])), ]
write.table(PF, file = "GOName_Translation.05_Iconic.csv", sep = ",", row.names = FALSE)

length(MigDB[MigDB$ProcessName == "GO_REGULATION_OF_CELL_MATRIX_ADHESION", 1])
length(MigDB[MigDB$ProcessName == "GO_CELL_MATRIX_ADHESION", 1])
length(intersect(MigDB[MigDB$ProcessName == "GO_CELL_MATRIX_ADHESION", 1], MigDB[MigDB$ProcessName == "GO_REGULATION_OF_CELL_MATRIX_ADHESION", 1]))
