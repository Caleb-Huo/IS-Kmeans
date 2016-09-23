## prepare group information   
## here module should be G-lists. Each element g contain feature indexes whose domain is P.
## L is the expanded length of non-zero element in J by G0 matrix.
## groupLevel (L): 1,1,1,1,2,2,2,3,3,3,... increasing, same number indicate same group.
## genePos (L): position.
## coef (L): coef for the expanded features.
## z (J): is the feature weight.
## x (J): primal variable.
## y (J): dual variable.

prepareGroup <- function(group, J, G0, gamma, alpha){  
  groupFeatureCounts <- numeric(J)
  for(g in 1:G0){
  	groupFeatureCounts[group[[g]]] <- groupFeatureCounts[group[[g]]] + 1	
  }

  L <- sum(groupFeatureCounts) + J
  groupLevel <- numeric(L)
  genePos <- numeric(L)
  coef <- numeric(L)
  
  curPos <- 1
  preCoef <- gamma * (1 - alpha)
  for(g in 1:G0){
	  agroup <- group[[g]]
	  alen <- length(agroup)
	  endPos <- curPos + alen - 1
	  groupLevel[curPos:endPos] <- g
	  genePos[curPos:endPos] <- agroup
	  coef[curPos:endPos] <- preCoef/groupFeatureCounts[agroup]
	  curPos <- curPos + alen
  }

  ## count all individual features: including both l1 and l2 norm penalty
  endPos <- curPos + J - 1
  groupLevel[curPos:endPos] <- (G0+1):(G0+J)
  genePos[curPos:endPos] <- 1:J
  coef[curPos:endPos] <- alpha * gamma
  coef[curPos:endPos][groupFeatureCounts==0] <- gamma

  groupInfo <- list(groupLevel=groupLevel, genePos=genePos, coef=coef, L=L, G=G0+J ,J=J)
  return(groupInfo)
}
