\name{qpRndGraph}
\alias{qpRndGraph}

\title{
Undirected random d-regular graphs
}
\description{
Samples an undirected d-regular graph approximately uniformly at random.
}
\usage{
qpRndGraph(p=6, d=2)
}
\arguments{
  \item{p}{number of vertices.}
  \item{d}{degree of every vertex.}
}
\details{
This function implements the algorithm from Steger and Wormald (1999) for
sampling undirected d-regular graphs from a probability distribution of
all d-regular graphs on p vertices which is approximately uniform. More
concretely, for all vertex degree values d that grow as a small power of
p, all d-regular graphs on p vertices will have in the limit the same
probability as p grows large. Steger and Wormald (1999, pg. 396) believe
that for d >> sqrt(p) the resulting probability distribution will no
longer be approximately uniform.

This function is provided in order to generate a random undirected graph
as input to the function \code{\link{qpG2Sigma}} which samples a random
covariance matrix whose inverse (aka, precision matrix) has zeroes on those
cells corresponding to the missing edges in the input graph. d-regular
graphs are useful for working with synthetic graphical models for two
reasons: one is that d-regular graph density is a linear function of d and
the other is that the minimum connectivity degree of two disconnected vertices
is an upper bound of their outer connectivity (see Castelo and Roverato,
2006, pg. 2646).

}
\value{
The adjacency matrix of the resulting graph.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Steger, A. and Wormald, N.C. Generating random regular graphs quickly,
\emph{Combinatorics, Probab. and Comput.}, 8:377-396.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpG2Sigma}}
}
\examples{
nVar <- 50  ## number of vertices
maxCon <- 5 ## maximum connectivity per vertex

set.seed(123)

A <- qpRndGraph(p=nVar, d=maxCon)

summary(apply(A, 1, sum))
}
\keyword{models}
\keyword{multivariate}
