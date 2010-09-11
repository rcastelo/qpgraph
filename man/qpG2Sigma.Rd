\name{qpG2Sigma}
\alias{qpG2Sigma}

\title{
Random covariance matrix
}
\description{
Builds a positive definite matrix from an undirected graph G that can
be used as a covariance matrix for a Gaussian graphical model with
graph G. The inverse of the resulting matrix contains zeroes at the missing
edges of the given undirected graph G.
}
\usage{
qpG2Sigma(g, rho=0, verbose=FALSE, R.code.only=FALSE)
}
\arguments{
  \item{g}{undirected graph specified either as a \code{graphNEL} object or as
           an adjacency matrix.}
  \item{rho}{real number between -1/(n.var-1) and 1.}
  \item{verbose}{show progress on the calculations.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       in the internal call to the IPF algorithm (default); if TRUE then only R code is executed.}
}
\details{
The random covariance matrix is built by first generating a random matrix
with the function \code{\link{qpRndWishart}} from a Wishart distribution
whose expected value is a matrix with unit diagonal and constant off-diagonal
entries equal to \code{rho}.
}
\value{
A random positive definite matrix that can be used as a covariance matrix
for a Gaussian graphical model with graph \code{G}.
}
\references{
Castelo, R. and Roverato, A. Utilities for large Gaussian graphical model
inference and simulation with the R package qpgraph, submitted.
}
\author{A. Roverato}
\seealso{
  \code{\link{qpGetCliques}}
  \code{\link{qpIPF}}
  \code{\link{qpRndWishart}}
  \code{\link[mvtnorm]{rmvnorm}}
}
\examples{
require(graph)

n.var <- 5 # number of variables
set.seed(123)
g <- randomEGraph(as.character(1:n.var), p=0.15)

Sigma <- qpG2Sigma(g, rho=0.5)

round(solve(Sigma), digits=2)

as(g, "matrix")
}
\keyword{models}
\keyword{multivariate}
