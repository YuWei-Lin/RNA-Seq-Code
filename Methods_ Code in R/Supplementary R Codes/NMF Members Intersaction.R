O1 <- NULL
O2 <- NULL
N1 <- NULL
N2 <- NULL
IneA <- NULL
IneB <- NULL
for(i in 1:10){
  p_time <- proc.time()
  res <- nmf(V, 2)
  w_V<- basis(res)
  h_V<- coef(res)
  t_time <- proc.time()-p_time
  print(t_time)
  
  INDVV <- as.data.frame(apply(w_V, 1, function(x) max(x)))
  INDVV$GROUP <- max.col(w_V)
  colnames(INDVV) <- c("Max_V", "K_Group")
  
  O1[i] <- length(rownames(INDV[INDV$K_Group == 1, ]))
  O2[i] <- length(rownames(INDV[INDV$K_Group == 2, ]))
  N1[i] <- length(rownames(INDVV[INDVV$K_Group == 1, ]))
  N2[i] <- length(rownames(INDVV[INDVV$K_Group == 2, ]))
  
  if(abs(O1[i] - N1[i]) < abs(O1[i] - N2[i])){
    IneA[i] <- length(intersect(rownames(INDV[INDV$K_Group == 1, ]), rownames(INDVV[INDVV$K_Group == 1, ])))/mean(c(O1[i], N1[i]))
    IneB[i] <- length(intersect(rownames(INDV[INDV$K_Group == 2, ]), rownames(INDVV[INDVV$K_Group == 2, ])))/mean(c(O2[i], N2[i]))
  } else if(abs(O1[i] - N1[i]) == abs(O1[i] - N2[i])) {
    IneA[i] <- "NaN"
    IneB[i] <- "NaN"
  } else {
    IneA[i] <- length(intersect(rownames(INDV[INDV$K_Group == 1, ]), rownames(INDVV[INDVV$K_Group == 2, ])))/mean(c(O1[i], N2[i]))
    IneB[i] <- length(intersect(rownames(INDV[INDV$K_Group == 2, ]), rownames(INDVV[INDVV$K_Group == 1, ])))/mean(c(O2[i], N1[i]))
  }
  INDV <- INDVV 
}


sum(!(rownames(INDVV[INDVV$K_Group == 1, ])%in%rownames(INDV[INDV$K_Group == 1, ])))
sum(!(rownames(INDV[INDV$K_Group == 1, ])%in%rownames(INDVV[INDVV$K_Group == 1, ])))