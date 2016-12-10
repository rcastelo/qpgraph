\name{qpGraph-class}
\docType{class}
\alias{qpGraph-class}
\alias{qpGraph}
\alias{show,qpGraph-method}
\alias{qpGraph,dspMatrix-method}
\alias{qpGraph,matrix-method}
\alias{names,qpGraph-method}
\alias{$,qpGraph-method}

\title{The qp-graph}

\description{
  The \code{"qpGraph"} class is the class to store and manipulate
  q-order (partial) correlation graphs, or qp-graphs for short. See
  Castelo and Roverato (2006, 2009) for a mathematical and statistical
  definition of a qp-graph.

  In earlier versions 1.x of the \code{\link{qpgraph}} package there
  was a function called \code{qpGraph()} to obtain a qp-graph from a
  matrix of non-rejection rates. This function, as it was written,
  has been deprecated and replaced by this class and corresponding
  constructor methods of the same name. The main difference with respect
  to earlier 1.x versions is that the argument \code{threshold} is now
  called \code{epsilon}, the argument \code{return.type} has been
  removed and the current version returns an object of this class
  \code{qpGraph} described in this manual page.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{qpGraph(nrrMatrix, ...)}
  corresponding to constructor methods that take as input a matrix of
  non-rejection rates, calculated with \code{\link{qpNrr}}.
}
\section{Slots}{
  \describe{
    \item{\code{p}:}{number of vertices, in one-to-one correspondence with random variables.}
    \item{\code{q}:}{order of the qp-graph, always smaller than \code{p-2}.}
    \item{\code{n}:}{when the qp-graph has been estimated from data, this is the number of observations in the data set, which must be larger than \code{q+2}.}
    \item{\code{epsilon}:}{maximum cutoff value of the non-rejection rate met by the edges that are present in the qp-graph.}
    \item{\code{g}:}{undirected graph structure of the qp-graph stored as a \code{\link[graph]{graphBAM-class}} object.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{qpGraph(nrrMatrix, ...)}}{Constructor method where \code{nrrMatrix} is a matrix of non-rejection rates.}
    \item{\code{show(object)}}{Method to display some bits of information about
      the qp-graph stored in the input argument \code{object}.}
  }
}
\arguments{
  \item{epsilon}{maximum cutoff value met by the edges present in the
                   qp-graph.}
  \item{topPairs}{number of edges from the top of the ranking, defined by the
                  non-rejection rates in \code{nrrMatrix}, to use to form the
                  resulting qp-graph. This parameter is incompatible with a
                  value different from \code{NULL} in \code{epsilon}.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{q}{q-order employed to derive the input matrix of non-rejection rates
           \code{nrrMatrix}.}
  \item{n}{when the input matrix of non-rejection rates \code{nrrMatrix} has
           been estimated from data, this is the number of observations in the
           data set.}
}
\references{
  Castelo, R. and Roverato, A. A robust procedure for
  Gaussian graphical model search from microarray data with p larger than n,
  \emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

  Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
  networks from microarray data with qp-graphs. \emph{J. Comp. Biol.},
  16(2):213-227, 2009.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpNrr}}
}
\keyword{classes}
