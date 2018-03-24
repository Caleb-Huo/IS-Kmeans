GetWCSS <- function(x, Cs, ws=NULL){
  wcss.perfeature <- numeric(ncol(x))
  for(k in unique(Cs)){
    whichers <- (Cs==k)
    if(sum(whichers)>1) wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  tss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)
  bcss.perfeature <- tss.perfeature-wcss.perfeature
  r <- bcss.perfeature
  
  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature, r=r))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature, r=r))
}
