\name{qpPCC}
\alias{qpPCC}
\alias{qpPCC,ExpressionSet-method}
\alias{qpPCC,data.frame-method}
\alias{qpPCC,matrix-method}

\title{
Estimation of Pearson correlation coefficients
}
\description{
Estimates Pearson correlation coefficients (PCCs) and their corresponding
P-values between all pairs of variables from an input data set.
}
\usage{
\S4method{qpPCC}{ExpressionSet}(data, long.dim.are.variables=TRUE)
\S4method{qpPCC}{data.frame}(data, long.dim.are.variables=TRUE)
\S4method{qpPCC}{matrix}(data, long.dim.are.variables=TRUE)
}
\arguments{
  \item{data}{data set from where to estimate the Pearson correlation
       coefficients. It can be an ExpressionSet object, a data frame or a
       matrix.}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when data are in a data frame or in a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
}
\details{
The calculations made by this function are the same as the ones made for
a single pair of variables by the function \code{\link{cor.test}} but for
all the pairs of variables in the data set.
}
\value{
A list with two matrices, one with the estimates of the PCCs and the other with
their P-values.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpPAC}}
}
\examples{
nVar <- 50 # number of variables
maxCon <- 5 # maximum connectivity per variable
nObs <- 30 # number of observations to simulate

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

X <- qpSampleMvnorm(K, nObs)

pcc.estimates <- qpPCC(X)

# Pearson correlation coefficients of the present edges
summary(abs(pcc.estimates$R[upper.tri(pcc.estimates$R) & I]))

# Pearson correlation coefficients of the missing edges
summary(abs(pcc.estimates$R[upper.tri(pcc.estimates$R) & !I]))

}
\keyword{models}
\keyword{multivariate}
