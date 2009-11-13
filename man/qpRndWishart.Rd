\name{qpRndWishart}
\alias{qpRndWishart}

\title{
Random Wishart distribution
}
\description{
Random generation for the (\code{n.var * n.var}) Wishart distribution (see
Press, 1972) with matrix parameter \code{A=diag(delta)\%*\%P\%*\%diag(delta)} and
degrees of freedom \code{df}.
}
\usage{
qpRndWishart(delta=1, P=0, df=NULL, n.var=NULL)
}
\arguments{
  \item{delta}{a numeric vector of \code{n.var} positive values. If a scalar
               is provided then this is extended to form a vector.}
  \item{P}{a (\code{n.var * n.var}) positive definite matrix with unit diagonal.
           If a scalar is provided then this number is used as constant
           off-diagonal entry for P.}
  \item{df}{degrees of freedom.}
  \item{n.var}{dimension of the Wishart matrix. It is required only when both
               delata and P are scalar.}
}
\details{
The degrees of freedom are \code{df > n.var-1} and the expected value of the
distribution is equal to \code{df * A}. The random generator is based on the
algorithm of Odell and Feiveson (1966).
}
\value{
A list of two \code{n.var * n.var} matrices \code{rW} and \code{meanW} where
\code{rW} is a random value from the Wishart and \code{meanW} is the expected
value of the distribution.
}
\references{
Odell, P.L. and Feiveson, A.G. A numerical procedure to generate a sample
covariance matrix. \emph{J. Am. Statist. Assoc.} 61, 199-203, 1966.

Press, S.J. \emph{Applied Multivariate Analysis: Using Bayesian and Frequentist
Methods of Inference}. New York: Holt, Rinehalt and Winston, 1972.
}
\author{A. Roverato}
\seealso{
  \code{\link{qpG2Sigma}}
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
