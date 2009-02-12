\name{Vowel}
\alias{Vowel}
\alias{vowel}
\alias{vowel.mem}
\title{The Vowel Data Set}
\description{The vowel data set.}
\usage{data(Vowel)}
\value{
An object consists of a 528 by 2 data matrix and a 528 by 1 cluster membership vector. 
There are 11 clusters each containing 48 data points in a two-dimensional space.
}
\details{
In the original vowel data set, the 11 different words, with each one representing a vowel sound,
correspond to 11 different clusters. There are 528 training observations and 462 testing observations.  
The training observations are employed to asset the performance of our algorithm and refered as the Vowel data set.
We projected the 10 dimensional data into a 2-dimensional space for examination purpose 
according to the method described in Hastie, Tibshirani and Friedman(pages 84, 92, 2001).
}
\references{
Hastie, T., Tibshirani, R., Friedman, J. (2001)
\emph{The elements of statistical learning.}
Springer-Verlag.
}
\examples{
  data(Vowel)

  # data matrix
  vowel <- Vowel$vowel

  # cluster membership
  vowel.mem <- Vowel$vowel.mem
}
\keyword{cluster}

