\name{qpPathWeight}
\alias{qpPathWeight}
\alias{qpPathWeight,matrix-method}

\title{
Calculation of path weights
}
\description{
Calculates the path weight for a path of an undirected graph.
}
\usage{
\S4method{qpPathWeight}{matrix}(X, path, Q=integer(0), normalized=TRUE,
                                R.code.only=TRUE)
}
\arguments{
  \item{X}{covariance matrix.}
  \item{path}{character vector of consecutive vertex names defining a path in an
              undirected graph.}
  \item{Q}{indexes or names of the variables in \code{sigma} that should be
           used for conditioning.}
  \item{normalized}{logical; TRUE (default) when the calculated path weight should be
        normalized so that weights are comparable between paths with different endpoints,
        and false otherwise.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
        (not yet available); if TRUE then only R code is executed (default).} } 

\details{
Calculation of path weights. This implementation is still under development and will give only correct results with either population covariance matrices or sample covariance matrices estimated from data with n >> p. Consult (Roverato and Castelo, accepted) for further details.
}

\value{
The calculated path weight for the given path.
}
\references{
Roverato, A. and Castelo, R. The networked partial correlation and its application to the analysis of genetic interactions. \emph{J. R. Stat. Soc. Ser. C-Appl. Stat., accepted}.
}
\author{R. Castelo and A. Roverato}
\examples{
## example in Figure 1 from (Castelo and Roverato, accepted)

## undirected graph on 9 vertices
edg <- matrix(c(1, 4,
                2, 4,
                3, 4,
                4, 5,
                5, 6,
                5, 7,
                8, 9),
            ncol=2, byrow=TRUE)

## create a corresponding synthetic precision matrix with
## partial correlation values set to -0.4 for all present edges
K <- matrix(0, nrow=9, ncol=9, dimnames=list(1:9, 1:9))
K[edg] <- -0.4
K <- K + t(K)
diag(K) <- 1

## calculate the corresponding covariance matrix
S <- solve(K)

## calculate networked partial correlations for all present
## edges
npc <- sapply(1:nrow(edg), function(i) qpPathWeight(S, edg[i, ]))

## note that while all partial correlations are zero for missing
## edges and are equal to -0.4 for present edges, the corresponding
## networked partial ## correlations are also zero for missing edges
## but may be different between them for present edges, depending on
## the connections between the vertices
cbind(edg, npc)
}
\keyword{models}
\keyword{multivariate}
