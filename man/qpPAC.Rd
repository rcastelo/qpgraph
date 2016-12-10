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
model with undirected graph G and their corresponding p-values for the
null hypothesis of zero-partial correlation.
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
  \item{g}{either a \code{qpGraph} object, or a \code{graphNEL}, \code{graphAM}
           or \code{graphBAM} object, or an adjacency matrix of an undirected graph.}
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
than the clique number \code{w(G)} of the graph \code{G}. If the sample size
\code{n} is larger than the maximum boundary of the input graph \code{bd(G)},
then the default matrix completion algorithm HTF by Hastie, Tibshirani and
Friedman (2009) can be used (see the function \code{\link{qpHTF}()} for details),
which has the avantage that is faster than IPF (see the function
\code{\link{qpIPF}()} for details).

The PAC estimation is done by first obtaining a MLE of the covariance matrix
using the \code{\link{qpIPF}} function and the p-values are calculated based on
the estimation of the standard errors (see Roverato and Whittaker, 1996) and
performing Wald tests based on the asymptotic chi-squared distribution.
}
\value{
A list with two matrices, one with the estimates of the PACs and the other with
their p-values. If \code{return.K=TRUE} then the MLE of the inverse covariance is
also returned as part of the list.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
networks from microarray data with qp-graphs. \emph{J. Comp. Biol.},
16(2):213-227, 2009.

Hastie, T., Tibshirani, R. and Friedman, J.H. \emph{The Elements of Statistical Learning},
Springer, 2009.

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

qpg <- qpGraph(nrr.estimates, epsilon=0.5)
qpg$g

pac.estimates <- qpPAC(X, g=qpg, verbose=FALSE)

## distribution absolute values of the estimated
## partial correlation coefficients of the present edges
summary(abs(pac.estimates$R[upper.tri(pac.estimates$R) & A]))

## distribution absolute values of the estimated
## partial correlation coefficients of the missing edges
summary(abs(pac.estimates$R[upper.tri(pac.estimates$R) & !A]))
}
\keyword{models}
\keyword{multivariate}
