\name{clustering}
\alias{clustering}
\title{Data Clustering (After Data Shrinking)}
\description{Data clustering (after data shrinking).}
\usage{clustering(y, disMethod = "Euclidean")}
\arguments{
  \item{y}{
  data matrix which is an R matrix object (for dimension > 1) or vector 
  object (for dimension=1) with rows be observations and columns be variables.}
  \item{disMethod}{specification of the dissimilarity measure. The available measures are \dQuote{Euclidean} and \dQuote{1-corr}.}
  }
\details{
We first store the first observation (data point) in \code{point[1]}. 
We then get the nearest neighbor of \code{point[1]}. Store it in 
\code{point[2]}. Store the dissimilarity between \code{point[1]} and
\code{point[2]} to \code{db[1]}. We next remove \code{point[1]}. 
We then find the nearest neighbor of \code{point[2]}. 
Store it in \code{point[3]}. Store the dissimilarity between \code{point[2]}
and \code{point[3]} to \code{db[2]}.  We then remove \code{point[2]} 
and find the nearest neighbor of \code{point[3]}. We repeat this procudure
until we find \code{point[n]} and \code{db[n-1]} where \code{n} is the 
total number of data points.

Next, we calculate the interquartile range (IQR) of the vector \code{db}.
We then check which elements of \code{db} are larger than \code{avg+1.5IQR}
where \code{avg} is the average of the vector \code{db}. The mininum value of
these outlier dissimilarities will be stored in \code{omin}. 
An estimate of the number of clusters is \code{g} where \code{g-1} is the number
of the outlier dissimilarities.
The position of an outlier dissimilarity 
indicates the end of a cluster and the start of a new cluster.

To get a reasonable clustering result, 
data sharpening (shrinking) is recommended before data clustering.
}
\value{
  \item{mem}{vector of the cluster membership of data points.
    The cluster membership takes values: \eqn{1}, \eqn{2}, \eqn{\ldots}, 
    \eqn{g}, where \eqn{g} is the estimated number of clusters.}
  \item{size}{vector of the number of data points for clusters.}
  \item{g}{an estimate of the number of clusters.}
  \item{db}{vector of dissimilarities between sorted consecutive data points 
    (c.f. details).}
  \item{point}{vector of sorted consecutive data points (c.f. details).}
  \item{omin}{The minimum value of the outlier dissimilarities (c.f. details).}
}
\references{ 
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

    tt <- shrinking(maronna, K = 50, itmax = 20)
    tt2 <- clustering(tt)

    # Plot of disimilarities between the sorted consecutive data points
    #     versus the sorted consecutive data points
    # This plot can be used to assess the estimated number of clusters
    db <- tt2$db
    point <- tt2$point
    n <- length(point)
    plot(1:(n - 1), db, type = "l",
        xlab = "sorted consecutive data points", 
        ylab = "disimilarities between the sorted consecutive data points", 
        xlim = c(0, n), axes = FALSE)
    box()
    axis(side = 2)
    axis(side = 1, at = c(0, 1:(n - 1)), labels = point)

}
\keyword{cluster}
