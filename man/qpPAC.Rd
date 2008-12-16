\name{qpPAC}
\alias{qpPAC}
\alias{qpPAC,ExpressionSet-method}
\alias{qpPAC,data.frame-method}
\alias{qpPAC,matrix-method}

\title{
Estimation of partial correlation coefficients
}
\description{
Estimates partial correlation coefficients (PACs) and their corresponding
P-values in a given undirected graph, from an input data set.
}
\usage{
\S4method{qpPAC}{ExpressionSet}(data, g, return.K=FALSE,
                                long.dim.are.variables=TRUE, verbose=TRUE,
                                R.code.only=FALSE)
\S4method{qpPAC}{data.frame}(data, g, return.K=FALSE,
                             long.dim.are.variables=TRUE, verbose=TRUE,
                             R.code.only=FALSE)
\S4method{qpPAC}{matrix}(data, g, return.K=FALSE,
                         long.dim.are.variables=TRUE, verbose=TRUE,
                         R.code.only=FALSE)
}
\arguments{
  \item{data}{data set from where to estimate the partial correlation
       coefficients. It can be an ExpressionSet object, a data frame or a
       matrix.}
  \item{g}{either a \code{graphNEL} object or an incidence matrix of the given
       undirected graph.}
  \item{return.K}{logical; if TRUE this function also returns the concentration
       matrix \code{K}; if FALSE it does not return it (default).}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when data are in a data frame or in a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
  \item{verbose}{show progress on the calculations.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
The estimation of PACs requires that the sample size \code{n} is strictly larger
than the number of variables \code{p}. In the context of microarray data and
regulatory networks, genes play the role of variables and thus normally
\code{p >> n}. For this reason, we can estimate PACs from the edges of a
regulatory network represented by an undirected graph \code{G} if and only if the
maximum clique size of the graph, noted \code{w(G)}, is strictly smaller than the
sample size \code{n} (number of experiments in the microarray data context).

In the context of this package, the undirected graph should correspond to a
qp-graph (see function \code{\link{qpGraph}}) we have selected by thresholding
on the (average) non-rejection rate calculated from this same data set using
the functions \code{\link{qpNrr}} or \code{\link{qpAvgNrr}}. If the resulting
graph is sparse enough we may have a chance to meet the requirement of
\code{w(G) < n} and the function \code{\link{qpClique}} can be useful to
investigate this. In the context of transcriptional regulatory networks we may
consider to remove edges between non-transcription factor genes which will
substantially increase the sparseness of the network.

The PAC estimation is done by first obtaining a maximum likelihood estimate of
the sample covariance matrix of the input data set using the \code{\{link{qpIPF}}
function and the P-values are calculated based on the estimation of the standard
errors of the edges following the procedure by Roverato and Whittaker (1996).
}
\value{
A list with two matrices, one with the estimates of the PACs and the other with
their P-values.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
networks from microarray data with qp-graphs.
\emph{J. Comp. Biol., accepted}, 2008.

Roverato, A. and Whittaker, J. Standard errors for the parameters of graphical
Gaussian models. \emph{Stat. Comput.}, 6:297-302, 1996.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpGraph}}
  \code{\link{qpCliqueNumber}}
  \code{\link{qpClique}}
  \code{\link{qpGetCliques}}
  \code{\link{qpIPF}}
}
\examples{
nVar <- 50 # number of variables
maxCon <- 5  # maximum connectivity per variable
nObs <- 30 # number of observations to simulate

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

X <- qpSampleMvnorm(K, nObs)

nrr.estimates <- qpNrr(X, verbose=FALSE)

g <- qpGraph(nrr.estimates, 0.5)

pac.estimates <- qpPAC(X, g=g, verbose=FALSE)

# estimated partial correlation coefficients of the present edges
summary(abs(pac.estimates$R[upper.tri(pac.estimates$R) & I]))

# estimated partial correlation coefficients of the missing edges
summary(abs(pac.estimates$R[upper.tri(pac.estimates$R) & !I]))

}
\keyword{models}
\keyword{multivariate}
