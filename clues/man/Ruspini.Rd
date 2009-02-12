\name{Ruspini}
\alias{Ruspini}
\alias{ruspini}
\alias{ruspini.mem}
\title{The Ruspini Data Set}
\description{The Ruspini data set.}
\usage{data(Ruspini)}
\value{
An object consists of a 75 by 2 data matrix and a 75 by 1 cluster membership vector.  
There are 4 clusters in the data set.}
\references{
Anja Struyf, Mia Hubert & Peter J. Rousseeuw (1996)
Clustering in an Object-Oriented Environment. 
\emph{Journal of Statistical Software},
*1*. <URL: http://www.stat.ucla.edu/journals/jss/>
}
\examples{
  data(Ruspini)
  # data matrix
  ruspini <- Ruspini$ruspini
  # cluster membership
  ruspini.mem <- Ruspini$ruspini.mem
}
\keyword{cluster}

