\name{qpUnifRndAssociation}
\alias{qpUnifRndAssociation}

\title{
Uniformly random association values
}
\description{
Builds a matrix of uniformly random association values between -1 and +1
for all pairs of variables that follow from the number of variables given
as input argument.
}
\usage{
qpUnifRndAssociation(n.var, var.names=as.character(1:n.var))
}
\arguments{
  \item{n.var}{number of variables.}
  \item{var.names}{names of the variables to use as row and column names in the resulting matrix.}
}
\details{
This function simply generates uniformly random association values with no independence
pattern associated to them. For generating a random covariance matrix that reflects
such a pattern use the function \code{\link{qpG2Sigma}}.
}
\value{
A symmetric matrix of uniformly random association values between -1 and +1.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpG2Sigma}}
}
\examples{
rndassociation <- qpUnifRndAssociation(100)
summary(rndassociation[upper.tri(rndassociation)])
}
\keyword{models}
\keyword{multivariate}
