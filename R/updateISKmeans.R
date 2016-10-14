updateISKmeans <- function(d, K, groupInfo, Css, wss, silent=silent, maxiter=maxiter){
  J <- ncol(d)
  wss.old <- rnorm(J)
  nonTrivialFlag = 1
  niter <- 0
  currentY <- NULL
cat('initilizaing results using alpha = 1, step 3\n')
	
  while((sum(abs(wss-wss.old))/sum(abs(wss.old)))>1e-4 && niter<maxiter){
  	cat('initilizaing results using alpha = 1, step 4\n')
    if(!silent) cat('Iteration',niter, ':\n', fill=FALSE)
		cat('initilizaing results using alpha = 1, step 4\n')
    niter <- niter+1
	cat('initilizaing results using alpha = 1, step 5\n')
    wss.old <- wss
	cat('initilizaing results using alpha = 1, step 6\n')
	if(sum(wss!=0)<1){
		cat('initilizaing results using alpha = 1, step 7\n')
		nonTrivialFlag<-0
		cat('initilizaing results using alpha = 1, step 8\n')
		wssPre <- wss
		cat('initilizaing results using alpha = 1, step 9\n')
		objective <- 0
		cat('initilizaing results using alpha = 1, step 10\n')
		obj0 <- 0
		cat('initilizaing results using alpha = 1, step 11\n')
		break
	}
	cat('initilizaing results using alpha = 1, step 12\n')
	
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