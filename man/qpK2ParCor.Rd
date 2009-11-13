\name{qpK2ParCor}
\alias{qpK2ParCor}

\title{
Partial correlation coefficients
}
\description{
Obtains partial correlation coefficients from a given concentration matrix.
}
\usage{
qpK2ParCor(K)
}
\arguments{
  \item{K}{positive definite matrix, typically a concentration matrix.}
}
\details{
This function applies \code{\link[stats]{cov2cor}} to the given concentration
matrix and then changes the sign of the off-diagonal entries in order
to obtain a partial correlation matrix.
}
\value{
A partial correlation matrix.
}
\references{
Lauritzen, S.L. \emph{Graphical models}. Oxford University Press, 1996.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpG2Sigma}}
}
\examples{
require(graph)

n.var <- 5 # number of variables
set.seed(123)
g <- randomEGraph(as.character(1:n.var), p=0.15)

Sigma <- qpG2Sigma(g, rho=0.5)
K <- solve(Sigma)

round(qpK2ParCor(K), digits=2)

as(g, "matrix")
}
\keyword{models}
\keyword{multivariate}
