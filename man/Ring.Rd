\name{Ring}
\alias{Ring}
\alias{ring}
\alias{ring.mem}
\title{The Broken-Ring Data Set}
\description{The broken-ring data set studied in Wang et al. (2007).}
\usage{data(Ring)}
\value{
A list contains a 800 by 2 data matrix and a 800 by 1 cluster membership vector. 
There are 5 clusters containing 200, 157, 151, 151 and 141 data points, respectively, 
in a two-dimensional space.
}
\references{
Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.
}
\examples{
    data(Ring)
 
    # data matrix
    ring <- Ring$ring
 
    # 'true' cluster membership
    ring.mem <- Ring$ring.mem

    # 'true' number of clusters
    nClust <- length(unique(ring.mem))

    # scatter plots
    plotClusters(ring, ring.mem)

}
\keyword{cluster}

