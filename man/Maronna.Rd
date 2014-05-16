\name{Maronna}
\alias{Maronna}
\alias{maronna}
\alias{maronna.mem}
\title{The Maronna Data Set}
\description{The Maronna data set.}
\usage{data(Maronna)}
\value{
A list contains a 200 by 2 data matrix and a 200 by 1 cluster membership vector.
There are 4 equally distributed clusters in a two-dimensional space.
The first cluster consists of the first 50 data points, 
the second cluster corresponds to the next 50 data points and so on.
}
\note{
We generated this data set by using the models and parameter settings
described in Maronna and Jacovkis (1974).
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
    # 'true' number of clusters
    nClust <- length(unique(maronna.mem))

    # scatter plots
    plotClusters(maronna, maronna.mem)

}
\keyword{cluster}

