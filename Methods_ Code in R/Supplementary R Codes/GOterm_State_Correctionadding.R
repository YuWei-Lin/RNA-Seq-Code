
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
for (i in 1:(length(unique(Sum1Rank))-1)) {
  for (j in (i+1):length(unique(Sum1Rank))) {
    if(all(names(unique(Sum1Rank)[[i]])%in%names(unique(Sum1Rank)[[j]])==T)){
      RMindex <- c(RMindex, i)
    }
  }
}
PAT <- unique(Sum1Rank)[-c(2,3)]
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
