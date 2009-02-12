\name{Iris}
\alias{Iris}
\alias{iris}
\alias{iris.mem}
\title{The Iris Data Set}
\description{The Iris data set.}
\usage{data(Iris)}
\value{
An object contains a 150 by 4 data matrix and a 150 by 1 cluster membership vector.
There are 3 balanced clusters in a 4 dimensional space spanned by the 4 columns.
}
\references{
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) 
\emph{The New S Language}. 
Wadsworth & Brooks/Cole. (has 'iris3' as 'iris'.)
}
\examples{
  data(Iris)
  # data matrix
  iris <- Iris$iris
  # cluster membership
  iris.mem <- Iris$iris.mem
}
\keyword{cluster}

