
theori <- t(data.frame(rep(0, 24762)))
rownames(theori) <- paste("The Origin")
Mel.standOri <- rbind(theori, Mel.stand)
BC0102oridist <- as.matrix(dist(Mel.standOri, method = "euclidean"))[1, ]

B01ori <- BC0102oridist[substr(names(BC0102oridist),1,4)=="BC01"]
B02ori <- BC0102oridist[substr(names(BC0102oridist),1,4)=="BC02"]

tiff(paste("C:/Users/tonyxp/Desktop/R00006 BC01BC02/","BC0102 the origin distance.tiff", sep=""), width=10000, height=5000, compression="lzw", res=300)
plot(density(BC0102oridist), main = "BC0102 the origin distance", col = "green")
lines(density(B01ori, adjust = 5), col = "blue")
lines(density(B02ori, adjust = 5), col = "red")
dev.off()

Mel.stand <- Mel.stand[,!apply(Mel.stand,2,function(x) all(x==x[1]))]
PCA1 <- prcomp(scale(Mel.stand))
plot(prcomp(Mel.stand))
summary(PCA1)
biplot(PCA1)
plot(PCA1$x[,1:2])

