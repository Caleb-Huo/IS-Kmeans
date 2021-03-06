##' Integrative Sparse KMeans
##'
##'  Integrative Sparse KMeans, integrating multiple omics dataset, using prior group information.
##' @title IS-Kmeans
##' @param d combined data matrix n*J, where n is number of subjects, J=J1+J2+... and J1 is number of features in omics dataset 1 and J2 is number of features in omics dataset 2...
##' @param K number of clusters
##' @param gamma Penalty on total number of features. Larger gamma will yeild small number of selected features.
##' @param alpha balance between group sparsity and individual sparsity. alpha=1 yeilds no group sparsity. alpha=0 yeilds no individual penalty.
##' @param group Prior group information. Potentially these group can contain overlap features. group is a list and each element of the list is feature index.
##' @param nstart Number of initials for Kmeans for sparse Kmeans
##' @param wsPre Initial feature weight.
##' @param penaltyInfo only for the purpose of gap statitics. Here we will fix the penalty design to perform gap statistics. The input should be a list of groupInfo. See groupInfo for details.
##' @param sparseStart Use Sparse Kmeans to do initialization.
##' @param silent Output progress.
##' @param maxiter Maximum numbre of iteration between ws and Cs.
##' @return m lists, m is length of gamma parameters. Each list is consisting of
##' (ws=ws, Cs=Cs, objective=ADMMobject$objective, BIC=BIC, gamma=agamma,alpha=alpha
##' \item{ws}{weight for each feature. Zero weight means the feature is not selected.}
##' \item{Cs}{Cluster Assignment}
##' \item{objective}{objective value}
##' \item{obj0}{sum of weighted separation ability. This term is for the purpose of gap statistics.}
##' \item{groupInfo}{a list containing group design, alpha, gamma}
##' @export
##' @useDynLib ISKmeans
##' @author Caleb
##' @examples
##' set.seed(123)
##'
##' # Generate two random omics datasets
##' mu <- c(-3,1,3)
##' Simu1_mRNA <- rbind(cbind(matrix(rnorm(40*5, mu[1], 0.1),40,5),
##'                           matrix(rnorm(40*5, mu[2], 0.1),40,5),
##'                           matrix(rnorm(40*5, mu[3], 0.1),40,5)),
##'                     matrix(rnorm(10*15,0,0.1),10,15))
##'
##' mu <- c(1,3,-3)
##' Simu1_methyl <- rbind(cbind(matrix(rnorm(40*5, mu[1], 0.1),40,5),
##'                             matrix(rnorm(40*5, mu[2], 0.1),40,5),
##'                             matrix(rnorm(40*5, mu[3], 0.1),40,5)),
##'                       matrix(rnorm(10*15,0,0.1),10,15))
##'
##' ## feature modules across two datasets
##' group <- list(c(1:10,51:60), c(11:20,61:70), c(21:30,71:80), c(31:40,81:90))
##'
##' DList <- rbind(Simu1_mRNA, Simu1_methyl)
##' dim(DList)
##'
##' ## colSums(module)
##' K=3
##' nstart=20
##' silent=FALSE
##' maxiter=6
##' centers=NULL
##' error=1e-4
##' sparseStart = TRUE
##' gamma = 0.2
##' alpha = 0.5
##' d <- t(DList)
##'
##' iRes <- ISKmeans(d, K=3, gamma=0.5, alpha=0.5, group=group)
##'
ISKmeans <-
function(d, K=NULL, gamma=NULL, alpha=0.5, group=NULL, nstart=20, wsPre=NULL ,penaltyInfo=NULL ,sparseStart=TRUE ,silent=FALSE, maxiter=20){
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
  if(!is.null(penaltyInfo)){
	  if(!(length(gamma) == length(penaltyInfo))){
		  stop('gamma and penaltyInfo must have the same length.')
	  }
  }
  if(is.null(K)) stop("Must provide either K or centers.")

  ## obtain basic information
  J <- ncol(d)
  G0 <- length(group)
  
  ## convert BCSS to BCSS/TSS
  tss.x <- apply(scale(d, center = TRUE, scale = FALSE)^2, 2, sum)
  d <- sweep(d, 2, sqrt(1/tss.x), "*")
  d[,tss.x==0] <- 0
  
  ## initialization, ws, Y
  if(!is.null(wsPre)){
	## pre-specify ws
  	iniws <- wsPre
	Cs <- UpdateCs(d, K, iniws, Cs=NULL)
  } else if(sparseStart){
  	## sparse start
  	skm <- KMeansSparseCluster(d,K=K,wbounds=10,nstart=nstart,silent=silent)[[1]]
	Cs <- skm$Cs
	iniws <- skm$ws
  }else{
    ## uniform start
    Cs <- kmeans(d, centers=K, nstart=nstart)$cluster
	iniws <- rep(1/sqrt(J),J)
  }
  wsPre <- iniws

  ## iteratively update CS, WS
  out <- replicate(length(gamma),list())
  for(i in 1:length(gamma)){
	agamma <- gamma[i]
	if(is.null(penaltyInfo)){
		cat('initilizaing results using alpha = 1\n')
		groupInfoIni <- prepareGroup(group, J, G0, agamma, 1, wsPre)
	   	ADMMobjectIni <- updateISKmeans(d, K, groupInfoIni, Cs, wsPre)
		cat('initilizaing groups\n')
	    groupInfo <- prepareGroup(group, J, G0, agamma, alpha, ADMMobjectIni$ws)		
		ADMMobject <- updateISKmeans(d, K, groupInfo, ADMMobjectIni$Cs, ADMMobjectIni$ws)
	} else {
		cat('using defined groups\n')	
		groupInfo <- penaltyInfo[[i]]
		ADMMobject <- updateISKmeans(d, K, groupInfo, Cs, wsPre)
	}
	out[[i]] <- ADMMobject
  }

  if(length(gamma)==1) out <- out[[i]]
  class(out) <- "ISKmeans"
  return(out)
}

