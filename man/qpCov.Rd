\name{qpCov}
\alias{qpCov}

\title{
Calculation of the sample covariance matrix
}
\description{
Calculates the sample covariance matrix, just as the function \code{cov()}
but returning a \code{\link{dspMatrix-class}} object which efficiently
stores such a dense symmetric matrix.
}
\usage{
qpCov(X)
}
\arguments{
  \item{X}{data set from where to calculate the sample covariance matrix.
           As the \code{cov()} function, it assumes the columns correspond
           to random variables and the rows to multivariate observations.}
}
\details{
The calculations made by this function are the same as the ones made for
a single pair of variables by the function \code{\link{cor.test}} but for
all the pairs of variables in the data set.
}
\value{
A sample covariance matrix stored as a \code{\link{dspMatrix-class}} object.
See the \code{Matrix} package for full details on this object class.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpPCC}}
}
\examples{
require(graph)
require(mvtnorm)

nVar <- 50 ## number of variables
nObs <- 10 ## number of observations to simulate

set.seed(123)

g <- randomEGraph(as.character(1:nVar), p=0.15)

Sigma <- qpG2Sigma(g, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

S <- qpCov(X)

## estimate Pearson correlation coefficients by scaling the sample covariance matrix
R <- cov2cor(as(S, "matrix"))

## get the corresponding boolean adjacency matrix
A <- as(g, "matrix") == 1

## Pearson correlation coefficients of the present edges
summary(abs(R[upper.tri(R) & A]))

## Pearson correlation coefficients of the missing edges
summary(abs(R[upper.tri(R) & !A]))

}
\keyword{models}
\keyword{multivariate}
