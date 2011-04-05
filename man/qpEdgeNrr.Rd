\name{qpEdgeNrr}
\alias{qpEdgeNrr}
\alias{qpEdgeNrr,ExpressionSet-method}
\alias{qpEdgeNrr,data.frame-method}
\alias{qpEdgeNrr,matrix-method}

\title{
Non-rejection rate estimation for a pair of variables
}
\description{
Estimates the non-rejection rate for one pair of variables.
}
\usage{
\S4method{qpEdgeNrr}{ExpressionSet}(X, i=1, j=2, q=1, I=NULL, restrict.Q=NULL, nTests=100,
                                    alpha=0.05, exact.test=TRUE, R.code.only=FALSE)
\S4method{qpEdgeNrr}{data.frame}(X, i=1, j=2, q=1, I=NULL, restrict.Q=NULL, nTests=100,
                                 alpha=0.05, long.dim.are.variables=TRUE, exact.test=TRUE,
                                 R.code.only=FALSE)
\S4method{qpEdgeNrr}{matrix}(X, i=1, j=2, q=1, I=NULL, restrict.Q=NULL, n=NULL, nTests=100,
                             alpha=0.05, long.dim.are.variables=TRUE, exact.test=TRUE,
                             R.code.only=FALSE)
}
\arguments{
  \item{X}{data set from where the non-rejection rate should be estimated. It
       can be either an \code{ExpressionSet} object, a data frame, or a matrix.
       If it is a matrix and the matrix is squared then this function assumes the
       matrix is the sample covariance matrix of the data and the sample size
       parameter \code{n} should be provided.}
  \item{i}{index or name of one of the two variables in \code{X} to test.}
  \item{j}{index or name of the other variable in \code{X} to test.}
  \item{q}{order of the conditioning subsets employed in the calculation.}
  \item{I}{indexes or names of the variables in \code{X} that are discrete.}
  \item{restrict.Q}{indexes or names of the variables in \code{X} that
       restrict the sample space of conditioning subsets Q.}
  \item{n}{number of observations in the data set. Only necessary when the
       sample covariance matrix is provided through the \code{X} parameter.}
  \item{nTests}{number of tests to perform for each pair for variables.}
  \item{alpha}{significance level of each test.}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when data are in a data frame or in a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
  \item{exact.test}{logical; if \code{FALSE} an asymptotic conditional independence
       test is employed with mixed (i.e., continuous and discrete) data;
       if \code{TRUE} (default) then an exact conditional independence test with
       mixed data is employed.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
The estimation of the non-rejection rate for a pair of variables is calculated
as the fraction of tests that accept the null hypothesis of independence given a
set of randomly sampled q-order conditionals.

Note that the possible values of \code{q} should be in the range 1 to
\code{min(p,n-3)}, where \code{p} is the number of variables and \code{n}
the number of observations. The computational cost increases linearly with
\code{q}.
}
\value{
An estimate of the non-rejection rate for the particular given pair of
variables.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
  \code{\link{qpAvgNrr}}
  \code{\link{qpHist}}
  \code{\link{qpGraphDensity}}
  \code{\link{qpClique}}
}
\examples{
require(mvtnorm)

nObs <- 100 ## number of observations to simulate

## the following adjacency matrix describes an undirected graph
## where vertex 3 is conditionally independent of 4 given 1 AND 2
A <- matrix(c(FALSE,  TRUE,  TRUE,  TRUE,
              TRUE,  FALSE,  TRUE,  TRUE,
              TRUE,   TRUE, FALSE, FALSE,
              TRUE,   TRUE, FALSE, FALSE), nrow=4, ncol=4, byrow=TRUE)
Sigma <- qpG2Sigma(A, rho=0.5)

X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

qpEdgeNrr(X, i=3, j=4, q=1, long.dim.are.variables=FALSE)

qpEdgeNrr(X, i=3, j=4, q=2, long.dim.are.variables=FALSE)
}
\keyword{models}
\keyword{multivariate}
