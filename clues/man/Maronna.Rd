\name{Maronna}
\alias{Maronna}
\alias{maronna}
\alias{maronna.mem}
\title{The Maronna Data Set}
\description{The Maronna data set.}
\usage{data(Maronna)}
\value{
An object contains a 200 by 2 data matrix and a 200 by 1 cluster membership vector.
There are 4 equally distributed clusters in a two-dimensional space.
The first cluster consists of the first 50 data points, 
the second cluster corresponds to the next 50 data points and so on.
}
\references{
Maronna, R., Jacovkis, P.M, 1974. Multivariate clustering procedures with variable metrics.
Biometrics 30, 499-505.
}
\examples{
  data(Maronna)
  # data matrix
  maronna <- Maronna$maronna
  # cluster membership
  maronna.mem <- Maronna$maronna.mem
}
\keyword{cluster}

