\name{Get.Silhouette}
\alias{Get.Silhouette}
\title{Compute Silhouette Index}
\description{Compute Silhouette index for a given partition of a data set.}
\usage{Get.Silhouette(y, mem, disMethod = "Euclidean")}
\arguments{
  \item{y}{
  data matrix which is a R matrix object (for dimension > 1) or vector 
  object (for dimension=1) with rows be observations and columns be variables.  
  }
  \item{mem}{vector of the cluster membership of data points.
    The cluster member ship takes values: \eqn{1}, \eqn{2}, \eqn{\ldots}, 
    \eqn{g}, where \eqn{g} is the estimated number of clusters.}
  \item{disMethod}{specification of the dissimilarity measure. The available measures are \dQuote{Euclidean} and \dQuote{1-corr}.}
}
\value{A list of 3 elements:

  \item{avg.s}{average Sihouette index.}
  \item{s}{vector of Sihouette indices for data points.}
  \item{neighbor}{a vector, the $i$-th element of which indicates 
which cluster is the nearest neighbor cluster of the $i$-th data point.}
}
\references{ 
Kaufman, L., Rousseeuw, P.J., (1990). 
\emph{Finding Groups in Data: An Introduction to Cluster Analysis}. 
Wiley, New York.

Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.
}
\examples{
# ruspini data
data(Ruspini)
# data matrix
ruspini <- Ruspini$ruspini
mem <- Ruspini$ruspini.mem
 
tt <- Get.Silhouette(ruspini, mem)
tt$avg.s
}
\keyword{cluster}
