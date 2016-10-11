## prepare group information   
## here module should be G-lists. Each element g contain feature indexes whose domain is P.
## L is the expanded length of non-zero element in J by G0 matrix.
## groupLevel (L): 1,1,1,1,2,2,2,3,3,3,... increasing, same number indicate same group.
## genePos (L): position.
## coef (L): coef for the expanded features.
## z (J): is the feature weight.
## x (J): primal variable.
## y (J): dual variable.
## ws: feature weight of previous iteration

prepareGroup <- function(group, J, G0, gamma, alpha, ws){  
  ## take care of trivial class
  if(gamma==0){
	  return(NULL)
  }
  
  groupFeatureCounts <- numeric(J)
  for(g in 1:G0){
  	groupFeatureCounts[group[[g]]] <- groupFeatureCounts[group[[g]]] + 1
  }
  
  curPos <- 1
  preCoef <- gamma * (1 - alpha)
  
  if(alpha==0){
	  J0logic <- groupFeatureCounts==0
	  J0 = sum(J0logic)
  	
	  L <- sum(groupFeatureCounts) + J0
	  groupLevel <- numeric(L)
	  genePos <- numeric(L)
	  coef <- numeric(L)
	
	  for(g in 1:G0){
		  agroup <- group[[g]]
		  aws <- ws[agroup]
		  alen <- length(agroup)
		  endPos <- curPos + alen - 1
		  groupLevel[curPos:endPos] <- g
		  genePos[curPos:endPos] <- agroup
		  a_inv_groupFeatureCounts <- 1/groupFeatureCounts[agroup]
		  agroupPenalty <- min(sum(a_inv_groupFeatureCounts[aws!=0]),1)
		  cat(agroupPenalty)
		  cat(' ')
		  coef[curPos:endPos] <- preCoef*sqrt(a_inv_groupFeatureCounts)*sqrt(agroupPenalty)
		  curPos <- curPos + alen
	  }  	
	  cat('\n')
	
	  endPos <- curPos + J0 - 1
	  groupLevel[curPos:endPos] <- (G0+1):(G0+J0)
	  genePos[curPos:endPos] <- (1:J)[J0logic]
	  coef[curPos:endPos] <- gamma  	
	
  } else if(alpha==1){
	  J0 = J

	  L <- J
	  groupLevel <- numeric(J)
	  genePos <- numeric(J)
	  coef <- numeric(J)  	
	  G0=0
  	
	  endPos <- curPos + J - 1
	  groupLevel[curPos:endPos] <- (G0+1):(G0+J)
	  genePos[curPos:endPos] <- 1:J
	  coef[curPos:endPos] <- alpha * gamma
  } else {
	  J0=J
	  L <- sum(groupFeatureCounts) + J
	  groupLevel <- numeric(L)
	  genePos <- numeric(L)
	  coef <- numeric(L)

	  for(g in 1:G0){
		  agroup <- group[[g]]
		  alen <- length(agroup)
		  endPos <- curPos + alen - 1
		  groupLevel[curPos:endPos] <- g
		  genePos[curPos:endPos] <- agroup
		  a_inv_groupFeatureCounts <- 1/groupFeatureCounts[agroup]
		  coef[curPos:endPos] <- preCoef*sqrt(a_inv_groupFeatureCounts)*sqrt(sum(a_inv_groupFeatureCounts))
		  curPos <- curPos + alen
	  }  	
 
	  endPos <- curPos + J - 1
	  groupLevel[curPos:endPos] <- (G0+1):(G0+J)
	  genePos[curPos:endPos] <- 1:J
	  coef[curPos:endPos] <- alpha * gamma
	  coef[curPos:endPos][groupFeatureCounts==0] <- gamma  	 
  }  

  groupInfo <- list(groupLevel=groupLevel, genePos=genePos, coef=coef, L=L, G=G0+J0 ,J=J)
  return(groupInfo)
}
