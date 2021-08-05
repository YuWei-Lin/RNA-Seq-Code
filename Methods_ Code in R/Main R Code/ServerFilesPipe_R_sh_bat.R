### Pipeline for merging the number of remaining genes after various cutoffs of NA levels in all 9 cases
CaseN <- c("Astrocytoma","GSE70630_OG","GSE103322_HNSCC","E_MTAB_6149_NSCLC")
CaseN <- "Astrocytoma"

Figpath=paste("D:/ISS_Cluster/", CaseN, "_SeedM/", sep = "")

### Generate .R files
RUN <- round(length(GNlist)/100)+1
for(y in 1:RUN){
  temp=readLines("D:/ISS_Cluster/NegSilhoutte Score Pipe Case Template.R")
  temp[93]=gsub("1:1", paste((y-1)*100+1, ":", min(y*100, length(GNlist)), sep=""),temp[93])
  temp[157]=gsub("1:100", paste((y-1)*100+1, ":", min(y*100, length(GNlist)), sep=""),temp[157])
  if(y < 10){ y <- paste("0", y, sep = "") }
  temp[166]=gsub("YY", y,temp[166])
  write(temp, file=paste(Figpath, CaseN,"_GOSil_M", y, ".r", sep = ""))
}
# GOSil_M stands for using "Mean"; GOSil_Z stands for using "Zero"

### Generate .sh files

for(y in 1:RUN){
  temp=readLines("D:/ISS_Cluster/SCCase_sh_draft.txt")
  temp[4]=gsub("SCCase", CaseN, temp[4])
  if(y < 10){ y <- paste("0", y, sep = "") }
  temp[4]=gsub("XX", paste("M", y, sep = ""), temp[4])
  write(temp, file=paste(Figpath, CaseN,"_GOSil_M", y, ".sh", sep = ""))
}
# GOSil_M stands for using "Mean"; GOSil_Z stands for using "Zero"

### Generate .bat files

lish <- list.files(Figpath, pattern = ".sh")
x <- NULL
for(e in 1:length(lish)){
  x <- append(x, paste("qsub -q hmque -l nodes=node39  /NA2/tonyxp/", CaseN, "_SeedM/", lish[e], sep = ""), after = length(x))
  write(x, file=paste(Figpath, CaseN, "_GOSilshM.sh", sep = "")) 
}
### Convert this window's ".sh" into Unix's real ".sh" by using linux commend "dos2unix"
### Then use "sh" commend to execute it.

#qsub -l nodes=1 -q hmque /NA2/tonyxp/GSE75688_Breast_Seed/GSE75688_BC_pipeZRM.sh
#llii=list.files(path = "D:/GSE75688_Breast_Seed")
#sp <- lapply(llii, function(x) strsplit(x, "\\.")[[1]]



