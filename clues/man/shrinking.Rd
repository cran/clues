\name{shrinking}
\alias{shrinking}
\title{Data Sharpening Based on K-nearest Neighbors}
\description{Data sharpening based on K-nearest neighbors.}
\usage{shrinking(y, K, disMethod = "Euclidean", eps = 1e-04, itmax = 500)}
\arguments{
  \item{y}{data matrix with rows being the observations and columns being variables.}
  \item{K}{number of nearest neighbors.}
  \item{disMethod}{specification of the dissimilarity measure. The available measures are \dQuote{Euclidean} and \dQuote{1-corr}.}
  \item{eps}{a small positive number. A value is regarded as zero if it is less than \code{eps}.}
  \item{itmax}{maximum number of iterations allowed.}
}
\details{Within each iteration, each data point is replaced by the vector of the coordinate-wise medians of its \code{K} nearest neighbors. Data points will
move toward the locally most dense data point by this shrinking process.
}
\value{Sharpened data set.}
\seealso{
\code{\link{clustering}}
}
\examples{
  # ruspini data
  data(Ruspini)
  # data matrix
  ruspini <- Ruspini$ruspini
  
  tt <- shrinking(ruspini, K = 25)
  tt2 <- clustering(tt)
  plotClusters(ruspini, tt2$mem)
}
\keyword{cluster}
