\name{Vowel}
\alias{Vowel}
\alias{vowel}
\alias{vowel.mem}
\title{The Vowel Data Set}
\description{The vowel data set studied in Hastie, Tibshirani and
Friedman (2001).}
\usage{data(Vowel)}
\value{
A list consists of a 528 by 2 data matrix and a 528 by 1 cluster membership vector. 
There are 11 clusters each containing 48 data points in a two-dimensional space.
}
\details{
In the original vowel data set, the 11 different words, with each one representing a vowel sound,
correspond to 11 different clusters. There are 528 training observations and 462 testing observations.  
The training observations are employed to asset the performance of 
the \code{clues} algorithm in Wang et al. (2007)
and refered as the Vowel data set.
The 10 dimensional data was projected into a 2-dimensional space for examination purpose 
according to the method described in Hastie, Tibshirani and Friedman(pages 84, 92, 2001).
}
\references{
Hastie, T., Tibshirani, R., Friedman, J. (2001)
\emph{The elements of statistical learning.}
Springer-Verlag.

Wang, S., Qiu, W., and Zamar, R. H. (2007).
CLUES: A non-parametric clustering method based on local shrinking.
\emph{Computational Statistics & Data Analysis}, Vol. 52, issue 1, 
pages 286-298.
}
\examples{
    data(Vowel)
 
    # data matrix
    vowel <- Vowel$vowel
 
    # 'true' cluster membership
    vowel.mem <- Vowel$vowel.mem

    # scatter plots
    plotClusters(vowel, vowel.mem)

}
\keyword{cluster}

