\name{plotClusters}
\alias{plotClusters}
\title{Plot Scatter Plot of Clusters in a Specified Two Dimensional Space}
\description{Plot scatter plot of clusters in a specified two dimensional space.}
\usage{plotClusters(y, mem, plot.dim = c(1, 2))}
\arguments{
  \item{y}{data matrix with rows being a collection of observations and columns being a list of variables.}
  \item{mem}{vector of the cluster membership of data points.}
  \item{plot.dim}{
  specifies the two dimensions to be plot. 
  }
}
\references{
Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.
}
\examples{
  # ruspini data
  data(Ruspini)
  ruspini <- Ruspini$ruspini
  ruspini.mem <- Ruspini$ruspini.mem
  
  plotClusters(ruspini, ruspini.mem, plot.dim = c(1, 2))
}
\keyword{cluster}

