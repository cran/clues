\name{clues}
\alias{clues}
\title{Clustering Method Based on Local Shrinking}
\description{
Automatically estimate the number of clusters for a given data set and 
get a partition.
}
\usage{
clues(y, n0 = 5, alpha = 0.05, eps = 1.0e-4, itmax = 20, 
      K2.vec = n0, strengthMethod = "sil", strengthIni = -3, 
      disMethod ="Euclidean", plotFlag = FALSE, 
      plot.dim = c(1, 2), quiet = FALSE)
}
\arguments{
  \item{y}{
  data matrix which is a R matrix object (for dimension > 1) or vector 
  object (for dimension=1) with rows being observations and columns being variables.
  }
  \item{n0}{a guess for the number of clusters.}
  \item{alpha}{speed factor.}
  \item{eps}{a small positive number. A value is regarded as zero if it is
less than \code{eps}.}
  \item{itmax}{maximum number of iterations allowed.}
  \item{K2.vec}{range for the number of nearest neighbors for the second pass of the iteration.}
  \item{strengthMethod}{
  specifies the prefered measure of the strength of the clusters 
(i.e., compactness of the clusters). Two available methods are \dQuote{sil} 
(silhouette index) and \dQuote{CH} (CH index).
  }
  \item{strengthIni}{
  initial value for the lower bound of the measure of the strength for the clusters.  Any negative values will do.
  }
  \item{disMethod}{specification of the dissimilarity measure. The available measures are \dQuote{Euclidean} and \dQuote{1-corr}.}
  \item{plotFlag}{logical. Indicates if a scatter plot of clusters should be output.}
  \item{plot.dim}{
  specifies the two dimensions to be plot. 
  }
  \item{quiet}{logical. Indicates if intermediate results should be output.}
}
\value{
  \item{K}{number of nearest neighbors can be used to get final clustering.}
  \item{size}{vector of the number of data points for clusters.}
  \item{mem}{vector of the cluster membership of data points. 
    The cluster member ship takes values: \eqn{1}, \eqn{2}, \eqn{\ldots}, 
    \eqn{g}, where \eqn{g} is the estimated number of clusters.}
  \item{g}{an estimate of the number of clusters.}
  \item{CH}{CH index value for the final partition if \code{strengthMethod} is 
\dQuote{CH}.}
  \item{avg.s}{average of the Silhoutte index value for the final partition if \code{strengthMethod} is \dQuote{sil}.}
  \item{s}{vector of Silhoutte indices for data points if \code{strengthMethod} is \dQuote{sil}.}
  \item{neighbor}{nearest neighbor clusters for data points if \code{strengthMethod} is \code{"sil"}.}
  \item{K.vec}{number of nearest neighbors used for each iteration.}
  \item{g.vec}{number of clusters obtained in each iteration.}
  \item{myupdate}{logical. Indicates if the partition obtained in the first pass is the same as that obtained in the second pass.}
  \item{y.old1}{data used for shrinking and clustering.}
  \item{y.old2}{data returned after shrinking and clustering.}
}
\references{
Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.
}
\note{ 
Occasionally, the number of clusters estimated by \code{clues} will be 
equal to the number of data points (that is, each data point forms a cluster). 
In this case, the estimated number of clusters was set to be equal to one. 
And the CH index or silhouette index will be set to be equal to NULL 
since CH index and silhouette index are not defined when the number of clusters is equal to one.
}
\examples{
  # ruspini data
  data(Ruspini)
  # data matrix
  ruspini <- Ruspini$ruspini
    
  res <- clues(ruspini)
}
\keyword{cluster}

