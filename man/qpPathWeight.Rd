\name{qpPathWeight}
\alias{qpPathWeight}
\alias{qpPathWeight,character-method}

\title{
Calculation of path weights
}
\description{
Calculates the path weight for a path of an undirected graph.
}
\usage{
\S4method{qpPathWeight}{character}(path, sigma, Q=integer(0), normalized=TRUE,
                                   R.code.only=TRUE)
}
\arguments{
  \item{path}{character vector of consecutive vertex names defining a path in an              undirected graph.}
  \item{sigma}{covariance matrix.}
  \item{Q}{indexes or names of the variables in \code{sigma} that should be
           used for conditioning.}
  \item{normalized}{logical; TRUE (default) when the calculated path weight should be normalized so that weights are comparable between paths with different endpoints, and false otherwise.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used                      (not yet available); if TRUE then only R code is executed
                     (default).}
}
\details{
Calculation of path weights. This implementation is still under development and will give only correct results with either population covariance matrices or sample covariance matrices estimated from data with n >> p. Consult Rovearto and Castelo (2016) for further details.

}
\value{
The calculated path weight for the given path.
}
\references{
Roverato, A. and Castelo, R. The networked partial correlation and its application to the analysis of genetic interactions. \emph{J. R. Stat. Soc. Ser. C-Appl. Stat., accepted}.
}
\author{R. Castelo and A. Roverato}
\examples{
## coming soon
}
\keyword{models}
\keyword{multivariate}
