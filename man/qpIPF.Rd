\name{qpIPF}
\alias{qpIPF}

\title{
Iterative proportional fitting algorithm
}
\description{
Performs maximum likelihood estimation of a sample covariance matrix
given the independence constraints from in input list of (maximal)
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

Whittaker, J. \emph{Graphical models in applied multivariate statistics.}
Wiley, 1990.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpGetCliques}}
  \code{\link{qpPAC}}
}
\examples{
nVar <- 50 # number of variables
maxCon <- 5  # maximum connectivity per variable
nObs <- 10 # number of observations to simulate

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)
Sigma <- qpDscale(solve(K)) # true covariance matrix

X <- qpSampleMvnorm(K, nObs)

# scaled sample covariance matrix
S <- qpDscale(cov(X))

# more efficient scaled sample covariance matrix
clqs <- qpGetCliques(I, verbose=FALSE)
S2 <- qpIPF(S, clqs)
S2 <- qpDscale(S2)

# mean squared error of S
mean((abs(Sigma-S)^2)[upper.tri(Sigma)])

# mean squared error of S2
mean((abs(Sigma-S2)^2)[upper.tri(Sigma)])

}
\keyword{models}
\keyword{multivariate}
