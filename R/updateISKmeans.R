updateISKmeans <- function(d, K, groupInfo, Css, wss, silent=silent, maxiter=maxiter){
  J <- ncol(d)
  wss.old <- rnorm(J)
  nonTrivialFlag = 1
  niter <- 0
  currentY <- NULL
cat('initilizaing results using alpha = 1, step 3\n')
	
  while((sum(abs(wss-wss.old))/sum(abs(wss.old)))>1e-4 && niter<maxiter){
    if(!silent) cat('Iteration',niter, ':\n', fill=FALSE)
    niter <- niter+1
    wss.old <- wss
	if(sum(wss!=0)<1){
		nonTrivialFlag<-0
		wssPre <- wss
		objective <- 0
		obj0 <- 0
		break
	}
	cat('initilizaing results using alpha = 1, step 4\n')
	
	if(!silent) cat('Updating CS...\n', fill=FALSE)
    if(niter>1) Css <- UpdateCss(d, K, wss, Css) # if niter=1, no need to update!!
 	if(!silent) cat('Updating WS...\n', fill=FALSE)
	if(is.null(groupInfo)){
		nonTrivialFlag<-0			
		wcss=GetWCSS(d, Css)
		wss <- wcss$r/sqrt(sum(wcss$r^2))
		wssPre <- wss
		objective <- - sum(wss * wcss$r)
		obj0 <-  - sum(wss * wcss$r)			
		print(objective)
	} else {
		ADMMobject <- UpdateWsADMM(d, Css, wss, currentY=currentY, groupInfo)
		wss <- ADMMobject$z
		currentY <- ADMMobject$currentY	
		print(ADMMobject$objective)
				
	}
  }
  if(nonTrivialFlag){
	  ##wss[wss<sum(wss)/ncol(d)] <- 0
	  Css <- UpdateCss(d, K, wss, Css)
	  wcss=GetWCSS(d, Css)
	  wssPre <- wss
	  objective = ADMMobject$objective
	  obj0 <-  - sum(wss * wcss$r)
	  ## original implementation
	  ## BIC <- (n - 1) * sum(wss * wcss$r) - log(n) * sum(wss)	  	
  }
  res <- list(wss=wss, Css=Css, obj0 = obj0, objective=objective, gamma=groupInfo$gamma,alpha=groupInfo$alpha)	  
  return(res)
}