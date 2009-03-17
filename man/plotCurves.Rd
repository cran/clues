\name{plotCurves}
\alias{plotCurves}
\title{Plot Trajectories for Each Cluster}
\description{Plot trajectories for each cluster.}
\usage{plotCurves(y, mem, myxlab = "variable", myylab = "observation")}
\arguments{
  \item{y}{
  data matrix which is a R matrix object (for dimension > 1) or vector 
  object (for dimension=1) with rows being observations and columns being variables.
  }
  \item{mem}{vector of the cluster membership for data points.
    The cluster member ship takes values: \eqn{1}, \eqn{2}, \eqn{\ldots}, 
    \eqn{g}, where \eqn{g} is the estimated number of clusters.}
  \item{myxlab}{a title for the x axis.}
  \item{myylab}{a title for the y-axis.}
}
\references{
Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.
}
\examples{
data(Curve)

# data matrix
curve <- Curve$curve

res <- clues(curve, disMethod = "1-corr")
plotCurves(curve, res$mem)
}
\keyword{cluster}

