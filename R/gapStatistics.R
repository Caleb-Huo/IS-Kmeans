##' gap statistics to tune gamma
##'
##' gap statistics to tune gamma, the tuning parameter to control number of features.
##' @title gap statistics
##' @param d combined data matrix n*J, where n is number of subjects, J=J1+J2+... and J1 is number of features in omics dataset 1 and J2 is number of features in omics dataset 2...
##' @param K number of clusters
##' @param B number of permutations.
##' @param gamma Penalty on total number of features. Larger gamma will yeild small number of selected features.
##' @param alpha balance between group sparsity and individual sparsity. alpha=1 yeilds no group sparsity. alpha=0 yeilds no individual penalty.
##' @param group Prior group information. Potentially these group can contain overlap features. group is a list and each element of the list is feature index.
##' @param seed random seed
##' @param silence do not print progress in silence = TRUE.
##' @return a table. Each row represents a gamma. 
##' \item{gamma}{input gammas}
##' \item{score}{objective score: see obj0 in ISKmeans function}
##' \item{E.score}{mean value of permutated score}
##' \item{se.score}{standard error of permutated score}
##' \item{gapStat}{gap statistics, score - E.score}
##' \item{numF}{number of selected features}
##' @export
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
##' d <- t(DList)
##'
##' gapResult <- gapStatistics(d,K=3,B=10,group=group)
##' print(gapResult)
gapStatistics <-
function(d,K=3,B=10,gamma=NULL,alpha=1, group=NULL,seed=15213,silence=FALSE){
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0)
     stop("'B' has to be a positive integer")
  if (is.null(gamma))
	  gamma <- seq(0,0.8, 0.05)
  if (min(gamma) < 0)
      stop("gamma should be greater than or equal to 0")

  ## get true objective score
  if (!silence)
  	cat("calculating true score...\n")
  set.seed(seed)
  trueRes <- ISKmeans(d,K=K,gamma=gamma,alpha=alpha,group=group,silent=TRUE)

  numF <- sapply(trueRes,function(x) sum(x$ws!=0))
  nGamma <- sapply(trueRes,function(x) x$gamma)
  score <- sapply(trueRes,function(x) x$obj0)

  if (!silence)
     cat("calculating permutated score, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", sep = "")
   E.score.full <- NULL
  for(b in 1:B){
	  cat(".", if (b%%50 == 0)
	         paste(b, "\n"))

	  set.seed(15213+b)
      ad = t(apply(d,1,function(x) sample(x)))
	  permRes <- ISKmeans(ad,K=K,gamma=gamma,alpha=alpha,group=group,silent=TRUE)
	  scoreperm <- sapply(permRes,function(x) x$obj0)
	  E.score.full <- rbind(E.score.full, scoreperm)
  }

  if (!silence && (B%%50 != 0))
          cat("", B, "\n")

  E.score <- colMeans(E.score.full)
  se.score <- apply(E.score.full,2,sd)/nrow(E.score.full)

  gapStat <- score - E.score

  result <- data.frame(gamma, score, E.score, se.score, gapStat,numF=numF)
  return(result)
}
