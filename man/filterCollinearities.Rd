\name{filterCollinearities}
\alias{filterCollinearities}

\title{
Filter collinearities
}
\description{
Filters out variables or features that lead to collinearities in the
input data.
}
\usage{
filterCollinearities(X, soft.filter=FALSE, long.dim.are.variables=TRUE)
}
\arguments{
  \item{X}{data set where collinearities are identified.}
  \item{soft.filter}{logical; if FALSE (default) then the intput object
       \code{X} is returned without the variables or features that lead
       to collinearities, i.e., a hard-filtered version of \code{X}; if
       TRUE, then a logical mask is returned with as many positions as
       variables or features and those that lead to collinearities are
       set to \code{TRUE}, while the rest are set to \code{FALSE}.}
  \item{long.dim.are.variables}{logical; if TRUE (default) it is assumed
       that when \code{X} is a \code{data.frame} or a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the \code{data.frame} or matrix.}
}
\details{
The input object \code{X} can be either a \code{matrix} object, a
\code{data.frame} object or any other class of object that can be
handled by the function \code{qpPCC()}, which is internally called,
such as an \code{ExpressionSet} object.
}
\value{
The input object \code{X} without the variables or features that
lead to collinearities when \code{soft.filter=FALSE}, its default
value. Otherwise, when \code{soft.filter=TRUE} then a logical mask
is returned.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpPCC}}
}
\examples{
## build an undirected GMM model with
## average correlation 0.99 on the present edges
set.seed(1234)
gmm <- rUGgmm(dRegularGraphParam(), rho=0.99)
gmm

## sample n=100 observations from this GMM
X <- rmvnorm(100, gmm)
dim(X)
head(X)

## notice some variables lead to collinearities (r > 0.99)
cor(X)

## mask those variables
mask <- filterCollinearities(X, long.dim.are.variables=FALSE,
                             soft.filter=TRUE)
mask
head(X[, !mask])
}
\keyword{models}
\keyword{multivariate}
