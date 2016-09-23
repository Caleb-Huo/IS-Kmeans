## iRes <- ISKmeans(d, K=3, gamma=seq(0,1,0.1), alpha=0.5, group=group, nstart=20, sparseStart=TRUE, silent=FALSE, maxiter=6)


ISKmeans <-
function(d, K=NULL, gamma=NULL, alpha=0.5, group=NULL, nstart=20, wsPre=NULL ,sparseStart=TRUE ,silent=FALSE, maxiter=6){
  # The criterion is : minimize_{w, C} sum_j w_j (R_j) + gamma_1*\sum_group penalty + gamma_2*||w||_1 s.t. ||w||_2=1, w_j>=0
  # x is the data, nxp
  # K is the number of clusters desired
  # gamma is tuning parameter controlling number of selected features. Larger gamma will decrease total number of genes.
  # alpha is ratio between individual penalty and group penalty. alpha=1 means there only exist individual penalty.
  # rho is initial value for ADMM, should be omitted later
  # group is the overlapping group information. Data structure of module should be a list. Each element of the list represents a group, 
  ## which contains feature index information.
  ## 
  
  ## check the input variables are complete.
  # wbounds is a vector of L1 constraints on w, of the form  sum(abs(w))<=wbounds[i]
  if(is.null(K)) stop("Must provide either K or centers.")
  	  	  	  
  ## obtain basic information
  J <- ncol(d)
  G0 <- length(group)
  
  ## initialization, ws, Y
  if(!is.null(wsPre)){
	## pre-specify ws  
  	iniws <- wsPre
	Cs <- UpdateCs(d, K, iniws, Cs=NULL)
  } else if(sparseStart){
  	## sparse start
  	skm <- sparcl::KMeansSparseCluster(d,K=K,wbounds=10,nstart=nstart)[[1]]
	Cs <- skm$Cs
	iniws <- skm$ws
  }else{
    ## uniform start
    Cs <- kmeans(d, centers=K, nstart=nstart)$cluster	
	iniws <- rep(1/sqrt(J),J)
  }
  
  
  ## iteratively update CS, WS
  
  gamma <- sort(gamma,decreasing=TRUE)
  out <- replicate(length(gamma),list())
  for(i in 1:length(gamma)){
	  agamma <- gamma[i]
	  if(agamma <= 0){
	  	agamma <- 1e-7
	  }

	  ws.old <- rnorm(J)
	  ws <- iniws
	  niter <- 0
	  currentY <- NULL

      groupInfo <- prepareGroup(group, J, G0, agamma, alpha)
	  while((sum(abs(ws-ws.old))/sum(abs(ws.old)))>1e-4 && niter<maxiter){
	    if(!silent) cat('Iteration',niter, ':\n', fill=FALSE)
	    niter <- niter+1
	    ws.old <- ws
		  cat('Updating CS...\n', fill=FALSE)
	    if(niter>1) Cs <- UpdateCs(d, K, ws, Cs) # if niter=1, no need to update!!
	 	  cat('Updating WS...\n', fill=FALSE)
		ADMMobject <- UpdateWsADMM(d, Cs, ws, currentY=currentY, groupInfo)
		ws <- ADMMobject$z
		currentY <- ADMMobject$currentY
	  }
  
	  ##ws[ws<sum(ws)/ncol(d)] <- 0
	  Cs <- UpdateCs(d, K, ws, Cs)
	  wcss=GetWCSS(d, Cs)
	  wsPre <- ws
	  n <- nrow(d)
	  BIC <- (n - 1) * sum(ws * wcss$r) - log(n) * sum(ws)
  
	  out[[i]] <- list(ws=ws, Cs=Cs, objective=ADMMobject$objective, BIC=BIC, gamma=agamma,alpha=alpha)  
  
  }
  
  if(length(gamma)==1) out <- out[[i]]	            
  class(out) <- "ISKmeans"
  return(out)
}

