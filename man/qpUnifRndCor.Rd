\name{qpUnifRndCor}
\alias{qpUnifRndCor}

\title{
Uniformly random correlation values
}
\description{
Builds a matrix of uniformly random correlation values between -1 and +1
for all pairs of variables that follow from the number of variables given
as input argument.
}
\usage{
qpUnifRndCor(n.var, var.names=1:n.var)
}
\arguments{
  \item{n.var}{number of variables.}
  \item{var.names}{names of the variables to use as row and column names in the resulting matrix.}
}
\details{
This function simply generates random correlation values with no independence
pattern associated to them. For generating a random covariance matrix that reflects
inksuch a pattern use the function \code{\link{qpG2Sigma}}.
}
\value{
A symmetric matrix of uniformly random correlation values between -1 and +1.
}
\references{
}
\author{R. Castelo}
\seealso{
  \code{\link{qpG2Sigma}}
}
\examples{
rndcor <- qpUnifRndCor(100)
summary(rndcor[upper.tri(rndcor)])
}
\keyword{models}
\keyword{multivariate}
