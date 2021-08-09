#Customized heatmap
library(gplots) 
library(rafalib)
# Create labels
labels <- strsplit(rownames(RD), "_")
labels <- unlist(lapply(labels, function(x) x[1]))
#labels <- substr(rownames(RD), 1, 4)

#Color customization - see help for brewer.pal for more information
library(RColorBrewer) 
grcol <- colorRampPalette(c("green","lightyellow","red"))(100)

# Molecular Subtypes labeling
library(Rtsne)
library(cluster)
colors = c(hsv(h=0, s=1, v=1),hsv(h=1/3, s=1, v=1),hsv(h=0.75, s=1, v=1),hsv(h=0.5, s=1, v=1),hsv(h=1/6, s=1, v=1),hsv(h=, s=0, v=0),
           hsv(h=, s=0, v=3/4),hsv(h=2/3, s=1, v=1),hsv(h=1/12, s=1, v=1),hsv(h=5/6, s=1, v=1),hsv(h=7/12, s=1, v=1),hsv(h=11/12, s=1/2, v=1), hsv(h=1/6, s=1, v=1/2))
cc=rep(colors[1], 549)
cc[substr(rownames( RD ),1,4)=="BC02"]=colors[2]
cc[substr(rownames( RD ),1,4)=="BC03"]=colors[3]
cc[substr(rownames( RD ),1,4)=="BC04"]=colors[4]
cc[substr(rownames( RD ),1,4)=="BC05"]=colors[5]
cc[substr(rownames( RD ),1,4)=="BC06"]=colors[6]
cc[substr(rownames( RD ),1,4)=="BC07"]=colors[7]
cc[substr(rownames( RD ),1,4)=="BC08"]=colors[8]
cc[substr(rownames( RD ),1,4)=="BC09"]=colors[9]
cc[substr(rownames( RD ),1,4)=="BC10"]=colors[10]
cc[substr(rownames( RD ),1,4)=="BC11"]=colors[11]
cc[substr(rownames( RD ),1,6)=="BC03LN"]=colors[3]
cc[substr(rownames( RD ),1,6)=="BC07LN"]=colors[7]

gps <- rep(1, 549)
gps[substr(rownames( RD ),1,4)=="BC02"]=2
gps[substr(rownames( RD ),1,4)=="BC03"]=3
gps[substr(rownames( RD ),1,4)=="BC04"]=4
gps[substr(rownames( RD ),1,4)=="BC05"]=5
gps[substr(rownames( RD ),1,4)=="BC06"]=6
gps[substr(rownames( RD ),1,4)=="BC07"]=7
gps[substr(rownames( RD ),1,4)=="BC08"]=8
gps[substr(rownames( RD ),1,4)=="BC09"]=9
gps[substr(rownames( RD ),1,4)=="BC10"]=10
gps[substr(rownames( RD ),1,4)=="BC11"]=11
gps[substr(rownames( RD ),1,6)=="BC03LN"]=3
gps[substr(rownames( RD ),1,6)=="BC07LN"]=7

# Targeted conbined GOID files
p_time <- proc.time()
for(i in 1:length(CadID)-1){
  RD_1 <- read.csv(paste("C:/Users/tonyxp/Desktop/ALL GO Runs/GSE75688(BC) FIX AuGOid/", CadID[i], ".csv", sep = ""), header = TRUE)
  row.names(RD_1) <- RD_1$X
  RD_1 <- RD_1[ ,-1]
  for(j in (i+1):length(CadID)){
    RD_2 <- read.csv(paste("C:/Users/tonyxp/Desktop/ALL GO Runs/GSE75688(BC) FIX AuGOid/", CadID[j], ".csv", sep = ""), header = TRUE)
    RD <- cbind(RD_1, RD_2[ ,-1])
    RD <- as.matrix(RD)
    tsne <- Rtsne(RD, dims = 2, perplexity=30, max_iter = 5000, check_duplicates = FALSE)
    cormat <- cor(t(RD), method = "pearson")
    cortrans <- 2*(1-cormat) 
    sil = silhouette (gps, cortrans)
    tt <- sil[ ,"sil_width"]
    
    tiff(paste("C:/Users/tonyxp/Desktop/geneset_CBSofBC/CBSofBC_P_tsne/", CadID[i],"+", CadID[j]," CBSofBC_P.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(tsne$Y, main= paste(CadID[i],"+", CadID[j], sep = ""), col= cc, pch = 16, cex = 0.4)
    dev.off()
    
    tiff(paste("C:/Users/tonyxp/Desktop/geneset_CBSofBC/CBSofBC_P_sil/", CadID[i],"+", CadID[j]," CBSofBC_P.tiff", sep=""), width=1600, height=1600, compression="lzw", res=300)
    plot(sil, main = paste(CadID[i],"+", CadID[j], sep = ""), col= colors[1:11], border=NA)
    dev.off()
    
    tiff(paste("C:/Users/tonyxp/Desktop/geneset_CBSofBC/CBSofBC_P_heatmap/", CadID[i],"+", CadID[j]," CBSofBC_P.tiff", sep=""), width=10000, height=5000, compression="lzw", res=300)
    heatmap.2(t(RD), labCol=labels, ColSideColors=cc, col=grcol,trace = "none", cexCol = 0.8, hclustfun = hclust, dendrogram = "column", cexRow = 0.25)
    dev.off()
    
    print(paste(i, "AND",j))
  }
}
t_time <- proc.time()-p_time
print(t_time)
