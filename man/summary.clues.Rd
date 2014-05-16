\name{summary.clues}
\alias{summary.clues}
\alias{print.summary.clues}
\title{Summary Method for CLUES Objects}
\description{Summarize a \code{\link{clues}} object and return an object
  of class \code{summary.clues}.
  There's a \code{\link{print}} method for the latter.
}
\usage{
\method{summary}{clues}(object, \dots)
\method{print}{summary.clues}(x, \dots)
}
\arguments{
  \item{x, object}{a \code{\link{clues}} object.}
  \item{\dots}{potential further arguments (require by generic).}
}
\seealso{\code{\link{clues}}, \code{\link{print.clues}}.
}
\examples{
    # Maronna data set
    data(Maronna)
    # data matrix
    maronna <- Maronna$maronna
    # apply clues to maronna data set
    res <- clues(maronna)

    summary(res)
}
\keyword{cluster}

