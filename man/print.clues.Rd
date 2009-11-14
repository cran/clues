\name{print.clues}
\alias{print.clues}
\title{Print Method for CLUES Objects}
\description{
  Prints the number of data points, number of variables, number of clusters,
  cluster sizes, strength method, strength value, dissimilarity measurement,
  components of the object retruend by \code{clues}.

  This is a method for the function \code{\link{print}()} for objects
  inheriting from class \code{\link{clues}}.
}
\usage{
\method{print}{clues}(x, \dots)
}
\arguments{
  \item{x}{a clues object.}
  \item{\dots}{potential further arguments (require by generic).}
}
\seealso{
  \code{\link{clues}}, \code{\link{summary.clues}}.
}
\examples{
    # Maronna data set
    data(Maronna)
    # data matrix
    maronna <- Maronna$maronna
    # apply clues to maronna data set
    res <- clues(maronna)

    print(res)
}
\keyword{cluster}
\keyword{print}
