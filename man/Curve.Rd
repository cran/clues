\name{Curve}
\alias{Curve}
\alias{curve}
\alias{curve.mem}
\docType{data}
\title{Curve Data Set}
\description{
  A toy example used to illustrate curve clustering.
}
\usage{data(Curve)}
\format{
  A list contains a 300 by 10 data matrix (\code{curve}) and a 300 by 1 cluster membership vector (\code{curve.mem}). There are 3 clusters, each containing 100 data points, respectively, in a 10-dimensional space.

}
\details{
cluster one is generated from the model:
\deqn{y_{ik} = sin(2 * PI * x_k) + e_{ik}, x_k \sim N(0, 1),  e_{ik} \sim N(0, 0.1), i = 1, \ldots, 100, k = 1, \ldots, 10.}{%
yik = sin(2 * PI * xk) + eik, xk ~ N(0, 1),  eik ~ N(0, 0.1), i = 1, \ldots, 100, k = 1, \ldots, 10.
}

cluster two is generated from the model:
\deqn{y_{ik} = cos(2 * PI * x_k) + e_{ik}, x_k \sim N(0, 1),  e_{ik} \sim N(0, 0.1), i = 1, \ldots, 100, k = 1, \ldots, 10.}{%
yik = cos(2 * PI * xk) + eik, xk ~ N(0, 1),  eik ~ N(0, 0.1), i = 1, ..., 100, k = 1, ..., 10.
}

cluster three is generated from the model:
\deqn{y_{ik} = e_{ik}, e_{ik} \sim N(0, 1), i = 1,\ldots, 100, k = 1, \ldots, 10.}{%
yik = eik, eik ~ N(0, 1), i = 1, ..., 100, k = 1, ..., 10.
}
}
\examples{
    data(Curve)
 
    # data matrix
    curve <- Curve$curve
 
    # 'true' cluster membership
    curve.mem <- Curve$curve.mem
    # 'true' number of clusters
    nClust <- length(unique(curve.mem))
 
    # plot average trajectories
    plotAvgCurves(curve, curve.mem)

}
\keyword{datasets}
