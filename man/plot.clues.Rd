\name{plot.clues}
\alias{plot.clues}
\title{Draw Scatter Plots or Plot of Average Trajectory per Cluster, or Both}
\description{Draw scatter plots, or plot of average trajectory per cluster,
or both. A simple menu will prompt for the user to choose what will be plotted.}
\usage{
\method{plot}{clues}(
    x, ask = TRUE, plot.dim = NULL, 
    xlab = NULL, ylab = NULL,
    xlim = NULL, ylim = NULL, cex = NULL,
    las = NULL, lwd = NULL, 
    xlab.avg.curve = "variable",
    ylab.avg.curve = "average observation", ...)
}
\arguments{
  \item{x}{an object returend by the function \code{clues}.}
  \item{ask}{logical; if true,
    \code{plot.clues} operates in interactive mode, via \code{\link{menu}}.}
  \item{plot.dim}{
  specifies the dimensions to be plot in pair-wise scatter plots. 
The number of dimensions specified
can be more than 2. If \code{plot.dim} is \code{NULL},
then pair-wise scatter plots of all dimensions will be plotted.}
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
  \item{las}{numeric in \eqn{0,1,2,3}; the style of axis labels (see \code{\link{par}}. 
0 means that the labels always parallel to the axis [default]; 
1 means that the labels always horizontal; 
2 means that the labels always perpendicular to the axis, 
3 means that the lables always vertical. 
If \code{las = NULL}, then \code{las} will be set to be 2.
}
  \item{lwd}{the line width, a positive number, defaulting to 1 (see
\code{\link{par}}. If \code{lwd = NULL}, the \code{lwd} will be set to be 3.}
  \item{xlab.avg.curve}{label for x-axis in the plot of average
trajectory per cluster.}
  \item{ylab.avg.curve}{label for y-axis in the plot of average
trajectory per cluster.}
  \item{...}{graphical parameters (see \code{\link{par}}).}
}
\examples{
    # Maronna data set
    data(Maronna)
    # data matrix
    maronna <- Maronna$maronna
    # apply clues to maronna data set
    res <- clues(maronna)
    # plot
    \dontrun{ plot(res) }
}
\keyword{cluster}

