UpdateWsADMM <- function(d, Cs, ws, currentY=NULL, groupInfo){


  r <- GetWCSS(d, Cs)$r

  J <- groupInfo$J
  L <- groupInfo$L
  G <- groupInfo$G
  groupLevel <- groupInfo$groupLevel 
  genePos <- groupInfo$genePos - 1
  coef <- groupInfo$coef
  
  
  stopifnot(L == length(groupLevel))
  stopifnot(L == length(genePos))
  stopifnot(L == length(coef))
  stopifnot(1:L == order(groupLevel))
  stopifnot(max(groupLevel) == G)
  stopifnot(max(genePos) == J-1)
    
  if(is.null(currentY)) currentY<-numeric(L)
  
  x <- numeric(L)
  z <- ws
	  
  ADMMobj <- .C('ADMM_updatew_R', 
    x = as.double(x),
	currentY = as.double(currentY),
	z = as.double(z),
	r = as.double(r),
	objective = as.double(0),
	groupLevel = as.integer(groupLevel),
	genePos = as.integer(genePos),
	coef = as.double(coef),
	J = as.integer(J),
	G = as.integer(G),
	L = as.integer(L)
  )
  
  ADMMobj$x <- NULL
  ADMMobj$r <- NULL
  
  ADMMobj$groupLevel <- NULL
  ADMMobj$genePos <- NULL
  ADMMobj$coef <- NULL
  ADMMobj$J <- NULL
  ADMMobj$G <- NULL
  ADMMobj$L <- NULL
  
  return(ADMMobj)
}

