UpdateCs <- function(x, K, ws, Cs){
  x <- x[,ws!=0]
  if(sum(ws!=0)==1){ 
	  z <- x * sqrt(ws[ws!=0])
  } else {
	  z <- sweep(x, 2, sqrt(ws[ws!=0]), "*")  	
  }
  nrowz <- nrow(z)
  km <- kmeans(z, centers=K, nstart=10)
  return(km$cluster)
}

