\name{adjustedRand}
\alias{adjustedRand}
\title{Calculate Agreement Indices Between Two Partitions for a Data Set}
\description{
Calculate the five agreement indices: Rand index, 
Hubert and Arabie's adjusted Rand index, 
Morey and Agresti's adjusted Rand index, 
Fowlkes and Mallows's index,
and Jaccard index,
which measure the agreement between any two partitions for a data set. 
}
\usage{adjustedRand(cl1, cl2, randMethod = c("Rand", "HA", "MA", "FM", "Jaccard"))}
\arguments{
  \item{cl1}{partition 1 of a data set.}
  \item{cl2}{partition 2 of a data set. \code{cl2} must have the same length as \code{cl1}, but could have different number of clusters.}
  \item{randMethod}{
    specifies the prefered external index to meaure the agreement between the
    two partitions \code{cl1} and \code{cl2}. Available indices are:
    \dQuote{Rand}, 
    \dQuote{HA}  (Hubert and Arabie's adjusted Rand index),
    \dQuote{MA}  (Morey and Agresti's adjusted Rand index),
    \dQuote{FM} (Fowlkes and Mallows's index),
    \dQuote{Jaccard} (Jaccard index). 
    By default, all 5 indices will be output.
  }
}
\value{Returns a vector of the index values.}
\references{
Milligan, G.W. and Cooper, M.C. (1986) 
A study of the comparability of external criteria for hierarchical
cluster analysis. 
\emph{Multivariate Behavioral Research \bold{21}, 441--458.}

Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.

}
\examples{
    cl1 <- c(1, 1, 1, 2, 2, 2, 2, 2)
    cl2 <- c(1, 2, 1, 2, 1, 2, 1, 3)
    adjustedRand(cl1, cl2)
    
    # perfect agreement
    cl1 <- c(1, 1, 2, 2)
    cl2 <- cl1
    adjustedRand(cl1, cl2)
    
}
\keyword{cluster}

