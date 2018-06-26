rare_func=function(otu, from=1, to=11, by=1, times=100){
  rf = data.frame()
  
  for (j in 1:times){
  
    r = seq(from=from, to = to, by = by)
  
    for (i in r){
  
      temp=otu[rowSums(otu) > i, ]
  
      temp=rrarefy(temp, i)
  
      temp=data.frame(t(data.frame(estimateR(temp))))
  
      temp$depth = rep(i)
      temp$Sample = row.names(temp)
      rf = rbind(rf, temp)
    }
  }
  
  rf = rf[,which(names(rf) %in% c("S.obs","Sample","depth"))]
  rf = aggregate(S.obs ~ (Sample + depth), data = rf, mean)
  rf$S.obs=round(rf$S.obs, 0)
  names(rf)[3]="observed"
  return(rf)
}
