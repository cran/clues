\name{plotAvgCurves}
\alias{plotAvgCurves}
\title{Average Trajectory Plot for Each Cluster}
\description{Average trajectory plot for each cluster.}
\usage{
  plotAvgCurves(y, mem, xlab = NULL, ylab = NULL,
      xlim = NULL, ylim = NULL, las = NULL, lwd = NULL, ...)
}
\arguments{
  \item{y}{
  data matrix which is an R matrix object (for dimension > 1) or vector 
  object (for dimension=1) with rows being observations and columns being variables.
  }
  \item{mem}{vector of the cluster membership of data points.
    The cluster membership takes values: \eqn{1}, \eqn{2}, \eqn{\ldots}, 
    \eqn{g}, where \eqn{g} is the estimated number of clusters.}
  \item{xlab}{a title for the x axis. If \code{xlab = NULL}, then
    \code{xlab} will be set to be \dQuote{variable}.}
  \item{ylab}{a title for the y-axis. If \code{ylab = NULL}, then
    \code{ylab} will be set to be \dQuote{average observation}.}
  \item{xlim}{range of x-axis. If \code{xlim = NULL}, then \code{xlim} will
    be set to be \code{(0, nc + 1)}, where \code{nc} is the number of
    columns of the matrix \code{y}.}
  \item{ylim}{range of y-axis. If \code{ylim = NULL}, then \code{ylim} will
    be set to be the range of the matrix \code{y}.} 
  \item{las}{numeric in \eqn{0,1,2,3}; the style of axis labels (see \code{\link{par}}. 
0 means that the labels always parallel to the axis [default]; 
1 means that the labels always horizontal; 
2 means that the labels always perpendicular to the axis, 
3 means that the lables always vertical. 
If \code{las = NULL}, then \code{las} will be set to be 2.
}
  \item{lwd}{the line width, a positive number, defaulting to 1 (see
\code{\link{par}}. If \code{lwd = NULL}, the \code{lwd} will be set to be 3.}
  \item{...}{graphical parameters (see \code{\link{par}}).}
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
    
    # 'true' cluster membership
    curve.mem <- Curve$curve.mem

    # 'true' number of clusters
    nClust <- length(unique(curve.mem))
 
    # plot average trajectories
    plotAvgCurves(curve, curve.mem)
}
\keyword{cluster}
