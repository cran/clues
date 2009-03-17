\name{compClust}
\alias{compClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare different partitions for a data set}
\description{
Compare different partitions for a data set based on 
agreement indices, average sihouette index and CH index.
}
\usage{
compClust(y, memMat, disMethod = "Euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ 
  data matrix which is a R matrix object (for dimension > 1) or vector 
  object (for dimension=1) with rows being observations and columns being variables.
  }
  \item{memMat}{cluster membership matrix. Each column corresponds to a partition
  of the matrix \code{y}. The numbers of clusters for different partitions can
  be different. The cluster membership of a \eqn{g}-cluster data set should
    take values: \eqn{1}, \eqn{2}, \eqn{\ldots}, \eqn{g}.}
  \item{disMethod}{specification of the dissimilarity measure. 
  The available measures are \dQuote{Euclidean} and \dQuote{1-corr}.
  }
}
\value{
  \item{avg.s}{a vector of average sihouette indices for the different partitions
    in \code{memMat}.}
  \item{CH}{a vector of CH indices for the different partitions in \code{memMat}.}
  \item{Rand}{a matrix of Rand indices measuring the pair-wise agreement among
   the different partitions in \code{memMat}.}
  \item{HA}{a matrix of Hubert and Arabie's adjusted Rand indices measuring 
   the pair-wise agreement among the different partitions in \code{memMat}.}
  \item{MA}{a matrix of Morey and Agresti's adjusted Rand indices measuring 
   the pair-wise agreement among the different partitions in \code{memMat}.}
  \item{FM}{a matrix of Fowlkes and Mallows's indices measuring 
   the pair-wise agreement among the different partitions in \code{memMat}.}
  \item{Jaccard}{a matrix of Jaccard indices measuring 
   the pair-wise agreement among the different partitions in \code{memMat}.}
}
\references{ 
Calinski, R.B., Harabasz, J., (1974). 
A dendrite method for cluster analysis.
\emph{Communications in Statistics}, Vol. 3, pages 1-27.


Kaufman, L., Rousseeuw, P.J., (1990). 
\emph{Finding Groups in Data: An Introduction to Cluster Analysis}. 
Wiley, New York.


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

  # ruspini data
  data(Ruspini)
  # data matrix
  ruspini <- Ruspini$ruspini
  ruspini.mem <- Ruspini$ruspini.mem
    
  res.CH <- clues(ruspini, strengthMethod = "CH", quiet = TRUE)
  res.sil <- clues(ruspini, strengthMethod = "sil", quiet = TRUE)
  res.km.HW <- kmeans(ruspini, 4, algorithm = "Hartigan-Wong")
  res.km.L <- kmeans(ruspini, 4, algorithm = "Lloyd")
  res.km.F <- kmeans(ruspini, 4, algorithm = "Forgy")
  res.km.M <- kmeans(ruspini, 4, algorithm = "MacQueen")

  memMat <- cbind(ruspini.mem, res.CH$mem, res.sil$mem, res.km.HW$cluster,
    res.km.L$cluster, res.km.F$cluster, res.km.M$cluster) 

  colnames(memMat) <- c("true", "clues.CH", "clues.sil", "km.HW", "km.L", "km.F",
    "km.M")

  res <- compClust(ruspini, memMat)
  round(res$avg.s, 1)
  round(res$CH, 1)
  round(res$Rand, 1)
  round(res$HA, 1)
  round(res$MA, 1)
  round(res$FM, 1)
  round(res$Jaccard, 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
