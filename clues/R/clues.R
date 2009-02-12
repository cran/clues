# v0.2.4 created on Feb. 4, 2009 by Weiliang Qiu
#  (1) added code to output final results
#  
# y - data matrix. Rows are observations; Columns are variables
# n0 - initial guess of the number of clusters
clues<-function(y, 
                n0=5, 
                alpha=0.05, 
                eps=1.0e-4, 
                itmax=500, 
                K2.vec=n0, 
                strengthMethod="sil", 
                strengthIni=-3, 
                disMethod="Euclidean",
                plotFlag=FALSE,
                plot.dim=c(1,2), 
                quiet=FALSE)
                
{
  strengthMethod=match.arg(strengthMethod, c("sil", "CH"))
  disMethod=match.arg(disMethod, c("Euclidean", "1-corr"))

  # y-- data matrix; rows are observations; columns are variables
  if(!is.matrix(y))
  { y<-matrix(y, ncol=1) }
  if(strengthMethod == "CH")
  {
    res<-clues.CH(y, n0, alpha, eps, itmax, K2.vec, 
                  strengthIni, disMethod, plotFlag, plot.dim, quiet)
    if(!quiet)
    {
      cat("********* output results ********\n")
      cat("CH = ", res$CH, "\n")
      cat("Estimated number of clusters = ", res$g, "\n")
      cat("*********************************\n")
    }

  } else if (strengthMethod == "sil") {
    res<-clues.sil(y, n0, alpha, eps, itmax, K2.vec, 
                   strengthIni, disMethod, plotFlag, plot.dim, quiet)
    if(!quiet)
    {
      cat("********* output results ********\n")
      cat("avg Silhouette = ", res$avg.s, "\n")
      cat("Estimated number of clusters = ", res$g, "\n")
      cat("*********************************\n")
    }

  } else {
    stop("'strengthMethod' does not match 'CH' or 'sil'!\n")
  }

  invisible(res)
}



