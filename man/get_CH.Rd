\name{get_CH}
\alias{get_CH}
\title{Compute CH Index}
\description{Compute CH index for a given partition of a data set.}
\usage{get_CH(y, mem, disMethod = "Euclidean")}
\arguments{
  \item{y}{
  data matrix which is an R matrix object (for dimension > 1) or vector 
  object (for dimension=1) with rows be observations and columns be variables.  
  }
  \item{mem}{vector of the cluster membership of data points.
    The cluster membership takes values: \eqn{1}, \eqn{2}, \eqn{\ldots}, 
    \eqn{g}, where \eqn{g} is the estimated number of clusters.}
  \item{disMethod}{specification of the dissimilarity measure. The available measures are \dQuote{Euclidean} and \dQuote{1-corr}.}
}
\value{The value of the CH index.}
\references{ 
Calinski, R.B., Harabasz, J., (1974). 
A dendrite method for cluster analysis.
\emph{Communications in Statistics}, Vol. 3, pages 1-27.

Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.
}
\examples{
    data(Maronna)
    # data matrix
    maronna <- Maronna$maronna

    # cluster membership
    maronna.mem <- Maronna$maronna.mem
     
    get_CH(maronna, maronna.mem)
}
\keyword{cluster}
