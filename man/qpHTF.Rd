\name{qpHTF}
\alias{qpHTF}

\title{
Hastie Tibshirani Friedman algorithm
}
\description{
Performs maximum likelihood estimation of a covariance matrix
given the independence constraints from an input undirected graph.
}
\usage{
qpHTF(S, g, tol = 0.001, verbose = FALSE, R.code.only = FALSE)
}
\arguments{
  \item{S}{input matrix, in the context of this package, the sample covariance
       matrix.}
  \item{g}{input undirected graph.}
  \item{tol}{tolerance under which the iterative algorithm stops.}
  \item{verbose}{show progress on calculations.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
This is an alternative to the Iterative Proportional Fitting (IPF) algorithm
(see, Whittaker, 1990, pp. 182-185 and \code{\link{qpIPF}}) which also
adjusts the input matrix to the independence constraints in the input undirected
graph. However, differently to the IPF, it works by going through each of the
vertices fitting the marginal distribution over the corresponding vertex boundary.
It stops when the adjusted matrix at the current iteration differs from the matrix
at the previous iteration in less or equal than a given tolerance value. This
algorithm is described by Hastie, Tibshirani and Friedman (2009, pg. 634), hence we
name it here HTF, and it has the advantage over the IPF that it does not require the
list of maximal cliques of the graph which may be exponentially large. In
contrast, it requires that the maximum boundary size of the graph is below the
number of samples where the input sample covariance matrix \code{S} was estimated.
For the purpose of exploring qp-graphs that meet such a requirement, one can use
the function \code{\link{qpBoundary}}.
}
\value{
The input matrix adjusted to the constraints imposed by the input undirected graph,
i.e., a maximum likelihood estimate of the sample covariance matrix that
includes the independence constraints encoded in the undirected graph.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Hastie, T., Tibshirani, R. and Friedman, J.H. \emph{The Elements of Statistical Learning}, Springer, 2009.

Whittaker, J. \emph{Graphical Models in Applied Multivariate Statistics.}
Wiley, 1990.

}
\note{
Thanks to Giovanni Marchetti for bringing us our attention to this algorithm and
sharing an early version of its implementation on the R package \code{ggm}.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpBoundary}}
  \code{\link{qpIPF}}
  \code{\link{qpPAC}}
}
\examples{
require(graph)
require(mvtnorm)

nVar <- 50  ## number of variables
nObs <- 100 ## number of observations to simulate

set.seed(123)

g <- randomEGraph(as.character(1:nVar), p=0.15)

Sigma <- qpG2Sigma(g, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

## MLE of the sample covariance matrix
S <- cov(X)

## more efficient MLE of the sample covariance matrix using HTF
S_htf <- qpHTF(S, g)

## get the adjacency matrix and put the diagonal to one
A <- as(g, "matrix")
diag(A) <- 1

## entries in S and S_htf for present edges in g should coincide
max(abs(S_htf[A==1] - S[A==1]))

## entries in the inverse of S_htf for missing edges in g should be zero
max(solve(S_htf)[A==0])
}
\keyword{models}
\keyword{multivariate}
