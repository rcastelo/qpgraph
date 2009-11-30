\name{qpCliqueNumber}
\alias{qpCliqueNumber}

\title{
Clique number
}
\description{
Calculates the size of the largest maximal clique (the so-called clique number
or maximum clique size) in a given undirected graph.
}
\usage{
qpCliqueNumber(g, exact.calculation=TRUE, return.vertices=FALSE,
               approx.iter=100, verbose=TRUE, R.code.only)
}
\arguments{
  \item{g}{either a \code{graphNEL} object or an adjacency matrix of the given
       undirected graph.}
  \item{exact.calculation}{logical; if TRUE then the exact clique number is
       calculated; if FALSE then a lower bound is given instead.}
  \item{return.vertices}{logical; if TRUE a set of vertices forming a maximal
       clique of maximum size is returned; if FALSE only the maximum clique
       size is returned.}
  \item{approx.iter}{number of iterations to be employed in the calculation of
       the lower bound (i.e., only applies when \code{exact.calculation=FALSE}.}
  \item{verbose}{show progress on calculations.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
The calculation of the clique number of an undirected graph is one of the basic
NP-complete problems (Karp, 1972) which means that its computational cost is
bounded by an exponential running time (Pardalos and Xue, 1994). The current
implementation uses C code from the GNU GPL Cliquer library by Niskanen and
Ostergard (2003) based on the, probably the fastest to date, algorithm by
Ostergard (2002).

The lower bound on the maximum clique size is calculated by ranking the
vertices by their connectivity degree, put the first vertex in a set and
go through the rest of the ranking adding those vertices to the set that
form a clique with the vertices currently within the set. Once the entire
ranking has been examined a large clique should have been built and eventually
one of the largests ones. This process is repeated a number of times
(\code{approx.iter}) each of which the ranking is altered with increasing
levels of randomness acyclically (altering 1 to $p$ vertices and again). Larger
values of \code{approx.iter} should provide tighter lower bounds although it has
been proven that no polynomial time algorithm can approximate the maximum
clique size within a factor of $n^\epsilon$ ($\epsilon > 0$), unless $P=NP$
(Feige et al, 1991; Pardalos and Xue, 1994).
}
\value{
a lower bound of the size of the largest maximal clique in the given graph, also
known as its clique number.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Feige, U., Goldwasser, S., Lov\'asz, L., Safra, S. and Szegedy, M.
Approximating the maximum clique is almost NP-Complete. \emph{Proc. 32nd IEEE
Symp. on Foundations of Computer Science}, 2-12, 1991.

Karp, R.M. Reducibility among combinatorial problems. \emph{Complexity of computer
computations}, 43:85-103, 1972.

Niskanen, S. Ostergard, P. Cliquer User's Guide, Version 1.0.
Communications Laboratory, Helsinki University of Technology, Espoo, Finland,
Tech. Rep. T48, 2003. (\url{http://users.tkk.fi/~pat/cliquer.html})

Ostergard, P. A fast algorithm for the maximum clique problem.
Discrete Appl. Math. 120:197-207, 2002.

Pardalos, P.M. and Xue, J. The maximum clique problem.
\emph{J. Global Optim.}, 4:301-328, 1994.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpClique}}
}
\examples{
require(graph)

nVar <- 50

set.seed(123)

g1 <- randomEGraph(V=as.character(1:nVar), p=0.3)
qpCliqueNumber(g1, verbose=FALSE)

g2 <- randomEGraph(V=as.character(1:nVar), p=0.7)
qpCliqueNumber(g2, verbose=FALSE)
}
\keyword{models}
\keyword{multivariate}
