updateISKmeans <- function(d, K, groupInfo, Cs, ws, silent=FALSE, maxiter=20){
  J <- ncol(d)
  ws.old <- rnorm(J)
  nonTrivialFlag = 1
  niter <- 0
  currentY <- NULL
	
  while((sum(abs(ws-ws.old))/sum(abs(ws.old)))>1e-4 && niter<maxiter){
    if(!silent) cat('Iteration',niter, ':\n', fill=FALSE)
    niter <- niter+1
    ws.old <- ws
	if(sum(ws!=0)<1){
		nonTrivialFlag<-0
		wsPre <- ws
		objective <- 0
		obj0 <- 0
		break
	}
	if(!silent) cat('Updating CS...\n', fill=FALSE)
    if(niter>1) Cs <- UpdateCs(d, K, ws, Cs) # if niter=1, no need to update!!
 	if(!silent) cat('Updating WS...\n', fill=FALSE)
	if(is.null(groupInfo)){
		nonTrivialFlag<-0			
		wcss=GetWCSS(d, Cs)
		ws <- wcss$r/sqrt(sum(wcss$r^2))
		wsPre <- ws
		objective <- - sum(ws * wcss$r)
		obj0 <-  - sum(ws * wcss$r)			
		print(objective)
	} else {
		ADMMobject <- UpdateWsADMM(d, Cs, ws, currentY=currentY, groupInfo)
		ws <- ADMMobject$z
		print(sum(ws!=0))
		currentY <- ADMMobject$currentY	
		print(ADMMobject$objective)
				
	}
  }
  if(nonTrivialFlag){
	  ##ws[ws<sum(ws)/ncol(d)] <- 0
	  Cs <- UpdateCs(d, K, ws, Cs)
	  wcss=GetWCSS(d, Cs)
	  wsPre <- ws
	  objective = ADMMobject$objective
	  obj0 <-  - sum(ws * wcss$r)
	  ## original implementation
	  ## BIC <- (n - 1) * sum(ws * wcss$r) - log(n) * sum(ws)	  	
  }
  res <- list(ws=ws, Cs=Cs, obj0 = obj0, objective=objective, gamma=gamma, alpha=alpha)	  
  return(res)
}