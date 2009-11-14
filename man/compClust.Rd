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
  data matrix which is an R matrix object (for dimension > 1) or vector 
  object (for dimension = 1) with rows being observations and columns being variables.
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

    # Maronna data set
    data(Maronna)
    # data matrix
    maronna <- Maronna$maronna
    # cluster membership
    maronna.mem <- Maronna$maronna.mem    

    # partition by clues and kmeans
    res_CH <- clues(maronna, strengthMethod = "CH", quiet = TRUE)
    res_sil <- clues(maronna, strengthMethod = "sil", quiet = TRUE)
    res_km_HW <- kmeans(maronna, 4, algorithm = "Hartigan-Wong")
    res_km_L <- kmeans(maronna, 4, algorithm = "Lloyd")
    res_km_F <- kmeans(maronna, 4, algorithm = "Forgy")
    res_km_M <- kmeans(maronna, 4, algorithm = "MacQueen")
 
    memMat <- cbind(maronna.mem, res_CH$mem, res_sil$mem, 
        res_km_HW$cluster, res_km_L$cluster, 
        res_km_F$cluster, res_km_M$cluster) 
 
    colnames(memMat) <- c("true", "clues_CH", "clues_sil", 
        "km_HW", "km_L", "km_F", "km_M")
 
    res <- compClust(maronna, memMat)

    print(sapply(res, function(x) {round(x,1)}))

}
\keyword{ cluster }
