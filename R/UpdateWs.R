UpdateWs <- function(x, Cs, l1bound){
  wcss.perfeature <- GetWCSS(x, Cs)$wcss.perfeature
  tss.perfeature <- GetWCSS(x, rep(1, nrow(x)))$wcss.perfeature
  lam <- BinarySearch(-wcss.perfeature+tss.perfeature, l1bound)
  ws.unscaled <- soft(-wcss.perfeature+tss.perfeature,lam)
  return(ws.unscaled/l2n(ws.unscaled))
}


