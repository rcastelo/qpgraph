\name{qpPCC}
\alias{qpPCC}
\alias{qpPCC,ExpressionSet-method}
\alias{qpPCC,data.frame-method}
\alias{qpPCC,matrix-method}

\title{
Estimation of Pearson correlation coefficients
}
\description{
Estimates Pearson correlation coefficients (PCCs) and their corresponding
P-values between all pairs of variables from an input data set.
}
\usage{
\S4method{qpPCC}{ExpressionSet}(X)
\S4method{qpPCC}{data.frame}(X, long.dim.are.variables=TRUE)
\S4method{qpPCC}{matrix}(X, long.dim.are.variables=TRUE)
}
\arguments{
  \item{X}{data set from where to estimate the Pearson correlation
       coefficients. It can be an ExpressionSet object, a data frame or a
       matrix.}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when \code{X} is a data frame or a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
}
\details{
The calculations made by this function are the same as the ones made for
a single pair of variables by the function \code{\link{cor.test}} but for
all the pairs of variables in the data set, with the exception of the treatment
of missing values, since only complete observations across all variables in
\code{X} are used.
}
\value{
A list with two matrices, one with the estimates of the PCCs and the other with
their P-values.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpPAC}}
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

pcc.estimates <- qpPCC(X)

## get the corresponding boolean adjacency matrix
A <- as(g, "matrix") == 1

## Pearson correlation coefficients of the present edges
summary(abs(pcc.estimates$R[upper.tri(pcc.estimates$R) & A]))

## Pearson correlation coefficients of the missing edges
summary(abs(pcc.estimates$R[upper.tri(pcc.estimates$R) & !A]))

}
\keyword{models}
\keyword{multivariate}
