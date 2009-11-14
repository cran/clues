\name{plotClusters}
\alias{plotClusters}
\title{Scatter Plots of Clusters in Specified Dimensions}
\description{Scatter plots of clusters in specified dimensions. If more than two dimensions are specified, pair-wise scatter plots will be plotted.}
\usage{
  plotClusters(y, mem, plot.dim = NULL,
      xlab = NULL, ylab = NULL,
      xlim = NULL, ylim = NULL, cex = NULL, ...)

}
\arguments{
  \item{y}{data matrix with rows being a collection of observations and columns being a list of variables.}
  \item{mem}{vector of the cluster membership of data points.
    The cluster membership takes values: \eqn{1}, \eqn{2}, \eqn{\ldots}, 
    \eqn{g}, where \eqn{g} is the estimated number of clusters.}
  \item{plot.dim}{
  specifies the dimensions to be plot. If \code{plot.dim} is \code{NULL},
then pair-wise scatter plots of all dimensions will be plotted.
  }
  \item{xlab}{a title for the x axis. If \code{xlab = NULL}, then
    \code{xlab} will be set to be \dQuote{}.}
  \item{ylab}{a title for the y-axis. If \code{ylab = NULL}, then
    \code{ylab} will be set to be \dQuote{}.}
  \item{xlim}{range of x-axis. If \code{xlim = NULL}, then \code{xlim} will
    be set to be the range of the matrix \code{y} in the dimensions specified by
    \code{plot.dim}.}
  \item{ylim}{range of y-axis. If \code{ylim = NULL}, then \code{ylim} will
    be set to be the range of the matrix \code{y} in the dimensions 
    specified by \code{plot.dim}.} 
  \item{cex}{A numerical value giving the amount by which plotting text
          and symbols should be magnified relative to the default. (see \code{\link{par}}). If \code{cex = NULL}, then \code{cex} will be set to be 2.}
  \item{...}{graphical parameters (see \code{\link{par}}).}
}
\references{
Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.
}
\examples{
    # Maronna data set
    data(Maronna)

    # data matrix
    maronna <- Maronna$maronna

    # cluster membership
    maronna.mem <- Maronna$maronna.mem
    
    plotClusters(maronna, maronna.mem)
}
\keyword{cluster}

