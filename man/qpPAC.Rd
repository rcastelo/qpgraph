\name{qpPAC}
\alias{qpPAC}
\alias{qpPAC,ExpressionSet-method}
\alias{qpPAC,data.frame-method}
\alias{qpPAC,matrix-method}

\title{
Estimation of partial correlation coefficients
}
\description{
Estimates partial correlation coefficients (PACs) for a Gaussian graphical
model with undirected graph G and their corresponding P-values for the
hypothesis of zero partial correlations.
}
\usage{
\S4method{qpPAC}{ExpressionSet}(X, g, return.K=FALSE, tol=0.001,
                                matrix.completion=c("HTF", "IPF"), verbose=TRUE,
                                R.code.only=FALSE)
\S4method{qpPAC}{data.frame}(X, g, return.K=FALSE, long.dim.are.variables=TRUE,
                             tol=0.001, matrix.completion=c("HTF", "IPF"),
                             verbose=TRUE, R.code.only=FALSE)
\S4method{qpPAC}{matrix}(X, g, return.K=FALSE, long.dim.are.variables=TRUE,
                         tol=0.001, matrix.completion=c("HTF", "IPF"),
                         verbose=TRUE, R.code.only=FALSE)
}
\arguments{
  \item{X}{data set from where to estimate the partial correlation
       coefficients. It can be an ExpressionSet object, a data frame or a
       matrix.}
  \item{g}{either a \code{graphNEL} object or an adjacency matrix of the given
       undirected graph.}
  \item{return.K}{logical; if TRUE this function also returns the concentration
       matrix \code{K}; if FALSE it does not return it (default).}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when \code{X} is a data frame or a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
  \item{tol}{maximum tolerance in the application of the IPF algorithm.}
  \item{matrix.completion}{algorithm to employ in the matrix completion operations
                           employed to construct a positive definite matrix with the
                           zero pattern specified in \code{g}}
  \item{verbose}{show progress on the calculations.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
In the context of maximum likelihood estimation (MLE) of PACs it is a necessary
condition for the existence of MLEs that the sample size \code{n} is larger
than the clique number \code{w(G)} of the graph \code{G}.

The PAC estimation is done by first obtaining a MLE of the covariance matrix
using the \code{\link{qpIPF}} function and the P-values are calculated based on
the estimation of the standard errors (see Roverato and Whittaker, 1996).
}
\value{
A list with two matrices, one with the estimates of the PACs and the other with
their P-values.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
networks from microarray data with qp-graphs. \emph{J. Comp. Biol.},
16(2):213-227, 2009.

Roverato, A. and Whittaker, J. Standard errors for the parameters of graphical
Gaussian models. \emph{Stat. Comput.}, 6:297-302, 1996.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpGraph}}
  \code{\link{qpCliqueNumber}}
  \code{\link{qpClique}}
  \code{\link{qpGetCliques}}
  \code{\link{qpIPF}}
}
\examples{
require(mvtnorm)

nVar <- 50  ## number of variables
maxCon <- 5 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

A <- qpRndGraph(p=nVar, d=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

nrr.estimates <- qpNrr(X, verbose=FALSE)

g <- qpGraph(nrr.estimates, threshold=0.5)

pac.estimates <- qpPAC(X, g=g, verbose=FALSE)

## distribution absolute values of the estimated
## partial correlation coefficients of the present edges
summary(abs(pac.estimates$R[upper.tri(pac.estimates$R) & A]))

## distribution absolute values of the estimated
## partial correlation coefficients of the missing edges
summary(abs(pac.estimates$R[upper.tri(pac.estimates$R) & !A]))
}
\keyword{models}
\keyword{multivariate}
