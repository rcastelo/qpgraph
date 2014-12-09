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
qpG2Sigma(g, rho=0, matrix.completion=c("HTF", "IPF"), tol=0.001,
          verbose=FALSE, R.code.only=FALSE)
}
\arguments{
  \item{g}{undirected graph specified either as a \code{graphNEL} object or as
           an adjacency matrix.}
  \item{rho}{real number between -1/(n.var-1) and 1 corresponding to the mean
             marginal correlation}
  \item{matrix.completion}{algorithm to employ in the matrix completion operations
                           employed to construct a positive definite matrix with the
                           zero pattern specified in \code{g}}
  \item{tol}{tolerance under which the matrix completion algorithm stops.}
  \item{verbose}{show progress on the calculations.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       in the internal call to the HTF, or IPF, algorithm (default); if TRUE then
       only R code is executed.}
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
Tur, I., Roverato, A. and Castelo, R. Mapping eQTL networks with mixed graphical Markov models.
\emph{Genetics}, 198(4):1377-1393, 2014.
}
\author{A. Roverato}
\seealso{
  \code{\link{qpRndGraph}}
  \code{\link{qpGetCliques}}
  \code{\link{qpIPF}}
  \code{\link{qpRndWishart}}
  \code{\link[mvtnorm]{rmvnorm}}
}
\examples{
set.seed(123)
G <- qpRndGraph(p=5, d=2)

Sigma <- qpG2Sigma(G, rho=0.5)

round(solve(Sigma), digits=2)

as(G, "matrix")
}
\keyword{models}
\keyword{multivariate}
