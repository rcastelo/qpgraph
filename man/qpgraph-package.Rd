\name{qpgraph-package}
\alias{qpgraph-package}
\alias{qpgraph}
\docType{package}
\title{
The q-order partial correlation graph learning software, qpgraph.
}
\description{
q-order partial correlation graphs, or qp-graphs for short, are
undirected Gaussian graphical Markov models built from q-order
partial correlations. They are useful for learning undirected
graphical Gaussian Markov models from data sets where the number of
random variables p exceeds the available sample size n as, for
instance, in the case of microarray data where they can be employed
to reverse engineer a molecular regulatory network.
}
\details{
  \tabular{ll}{
  Package: \tab qp\cr
  Version: \tab 1.3.4\cr
  Date: \tab 30-11-2009\cr
  biocViews: \tab Microarray, GeneExpression, Transcription, Pathways, Bioinformatics, GraphsAndNetworks\cr
  Suggests: \tab mvtnorm, graph, Rgraphviz, annotate, genefilter, Category (>= 2.9.7), org.EcK12.eg.db (>= 2.2.6), GOstats\cr
  License: \tab GPL (>= 2)\cr
  URL: \tab \url{http://functionalgenomics.upf.edu/qpgraph}\cr
  }
}

\section{Functions}{
  \itemize{
    \item \code{\link{qpNrr}} estimates non-rejection rates for every pair
          of variables.
    \item \code{\link{qpAvgNrr}} estimates average non-rejection rates for
          every pair of variables.
    \item \code{\link{qpEdgeNrr}} estimate the non-rejection rate of one
          pair of variables.
    \item \code{\link{qpCItest}} performs a conditional independence test
          between two variables given a conditioning set.
    \item \code{\link{qpHist}} plots the distribution of non-rejection rates.
    \item \code{\link{qpGraph}} obtains a qp-graph from a matrix of
          non-rejection rates.
    \item \code{\link{qpAnyGraph}} obtains an undirected graph from a matrix of
          pairwise measurements.
    \item \code{\link{qpGraphDensity}} calculates and plots the graph density
          as function of the non-rejection rate.
    \item \code{\link{qpCliqueNumber}} calculates the size of the largest
          maximal clique (the so-called clique number or maximum clique size) in
          a given undirected graph.
    \item \code{\link{qpClique}} calculates and plots the size of the largest
          maximal clique (the so-called clique number or maximum clique size)
          as function of the non-rejection rate.
    \item \code{\link{qpGetCliques}} finds the set of (maximal) cliques of
          a given undirected graph.
    \item \code{\link{qpRndWishart}} random generation for the Wishart
          distribution.
    \item \code{\link{qpG2Sigma}} builds a random covariance matrix from an
          undrected graph. The inverse of the resulting matrix contains zeroes
          at the missing edges of the given undirected graph.
    \item \code{\link{qpK2ParCor}} obtains the partial correlation coefficients
          from a given concentration matrix.
    \item \code{\link{qpIPF}} performs maximum likelihood estimation of a
          sample covariance matrix given the independence constraints from
          an input list of (maximal) cliques.
    \item \code{\link{qpPAC}} estimates partial correlation coefficients and
          corresponding P-values for each edge in a given undirected graph,
          from an input data set.
    \item \code{\link{qpPCC}} estimates pairwise Pearson correlation coefficients
          and their corresponding P-values between all pairs of variables from an
          input data set.
    \item \code{\link{qpRndGraph}} builds a random undirected graph with a
          bounded maximum connectivity degree on every vertex.
    \item \code{\link{qpPrecisionRecall}} calculates the precision-recall curve
          for a given measure of association between all pairs of variables in a
           matrix.
    \item \code{\link{qpPRscoreThreshold}} calculates the score threshold at a
          given precision or recall level from a given precision-recall curve.
    \item \code{\link{qpImportNrr}} imports non-rejection rates.
    \item \code{\link{qpFunctionalCoherence}} estimates functional coherence of
          a given transcriptional regulatory network using Gene Ontology
          annotations.
  }

This package provides an implementation of the procedures described in (Castelo
and Roverato, 2006, 2009). An example of its use for reverse-engineering of
transcriptional regulatory networks from microarray data is available in the
vignette \code{qpTxRegNet}. This package is a contribution to the Bioconductor
(Gentleman et al., 2004) and gR (Lauritzen, 2002) projects.
}

\author{
R. Castelo and A. Roverato

Maintainer: R. Castelo <robert.castelo@upf.edu>
}

\references{
Castelo, R. and Roverato, A. A robust procedure for Gaussian graphical
model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
networks from microarray data with qp-graphs. \emph{J. Comput. Biol.},
16(2):213-227, 2009.

Gentleman, R.C., Carey, V.J., Bates, D.M., Bolstad, B., Dettling, M.,
Dudoit, S., Ellis, B., Gautier, L., Ge, Y., Gentry, J., Hornik, K. Hothorn, T.,
Huber, W., Iacus, S., Irizarry, R., Leisch, F., Li, C., Maechler, M. Rosinni,
A.J., Sawitzki, G., Smith, C., Smyth, G., Tierney, L., Yang, T.Y.H. and
Zhang, J. Bioconductor: open software development for computational biology
and bioinformatics. \emph{Genome Biol.}, 5:R80, 2004.

Lauritzen, S.L. (2002). gRaphical Models in R. \emph{R News}, 3(2)39.

}
\keyword{package}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
