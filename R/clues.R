# v0.3.1 created on Feb. 15, 2009 by Weiliang Qiu
#  (1) if number of clusters = number of data points,
#      then we set number of clusters=1, 
#      and CH index (or silhouette index)=NULL
#
# v0.2.4 created on Feb. 4, 2009 by Weiliang Qiu
#  (1) added code to output final results
#  
# y - data matrix. Rows are observations; Columns are variables
# n0 - initial guess of the number of clusters
clues<-function(y, 
                n0=5, 
                alpha=0.05, 
                eps=1.0e-4, 
                itmax=20, 
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

  nr<-nrow(y)

  if(strengthMethod == "CH")
  {
    res<-clues.CH(y, n0, alpha, eps, itmax, K2.vec, 
                  strengthIni, disMethod, plotFlag, plot.dim, quiet)
  } else if (strengthMethod == "sil") {
    res<-clues.sil(y, n0, alpha, eps, itmax, K2.vec, 
                   strengthIni, disMethod, plotFlag, plot.dim, quiet)
  } else {
    stop("'strengthMethod' does not match 'CH' or 'sil'!\n")
  }

  if(nr==res$g)
  { res$g=1
    res$mem<-rep(1, nr)
    res$size<-nr
    res$neighbor<-NULL
    if(strengthMethod=="CH")
    { 
      res$CH=NULL
    } else {
      res$avg.s=NULL
      res$s=rep(NULL, nr)
    }
  }

  if(!quiet)
  {
    cat("********* output results ********\n")
    if(strengthMethod=="CH")
    {
      cat("CH = ", res$CH, "\n")
    } else {
      cat("avg Silhouette = ", res$avg.s, "\n")
    }
    cat("Estimated number of clusters = ", res$g, "\n")
    cat("*********************************\n")
  }



  invisible(res)
}



