\name{qpIPF}
\alias{qpIPF}

\title{
Iterative proportional fitting algorithm
}
\description{
Performs maximum likelihood estimation of a covariance matrix
given the independence constraints from an input list of (maximal)
cliques.
}
\usage{
qpIPF(vv, clqlst, tol = 0.001, verbose = FALSE, R.code.only = FALSE)
}
\arguments{
  \item{vv}{input matrix, in the context of this package, the sample covariance
       matrix.}
  \item{clqlst}{list of maximal cliques obtained from an undirected graph
       by using the function \code{\link{qpGetCliques}}.}
  \item{tol}{tolerance under which the iterative algorithm stops.}
  \item{verbose}{show progress on calculations.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
The Iterative proportional fitting algorithm (see, Whittaker, 1990, pp. 182-185)
adjusts the input matrix to the independence constraints in the undirected graph
from where the input list of cliques belongs to, by going through each of the
cliques fitting the marginal distribution over the clique for the fixed
conditional distribution of the clique. It stops when the adjusted matrix at
the current iteration differs from the matrix at the previous iteration in less
or equal than a given tolerance value.
}
\value{
The input matrix adjusted to the constraints imposed by the list of cliques,
i.e., a maximum likelihood estimate of the sample covariance matrix that
includes the independence constraints encoded in the undirected graph formed
by the given list of cliques.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Tur, I., Roverato, A. and Castelo, R. Mapping eQTL networks with mixed graphical Markov models.
\emph{Genetics}, 198(4):1377-1393, 2014.

Whittaker, J. \emph{Graphical models in applied multivariate statistics.}
Wiley, 1990.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpGetCliques}}
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

## more efficient MLE of the sample covariance matrix using IPF
clqs <- qpGetCliques(g, verbose=FALSE)
S_ipf <- qpIPF(S, clqs)

## get the adjacency matrix and put the diagonal to one
A <- as(g, "matrix")
diag(A) <- 1

## entries in S and S_ipf for present edges in g should coincide
max(abs(S_ipf[A==1] - S[A==1]))

## entries in the inverse of S_ipf for missing edges in g should be zero
max(solve(S_ipf)[A==0])
}
\keyword{models}
\keyword{multivariate}
