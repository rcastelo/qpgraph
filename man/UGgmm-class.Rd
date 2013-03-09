\name{UGgmm-class}
\docType{class}
\alias{UGgmm-class}
\alias{UGgmm,missing-method}
\alias{UGgmm,matrix-method}
\alias{UGgmm,graphBAM-method}
\alias{UGgmm}
\alias{rUGgmm,graphParam,missing-method}
\alias{rUGgmm,missing,graphParam-method}
\alias{rUGgmm,numeric,graphParam-method}
\alias{rUGgmm,integer,graphParam-method}
\alias{rUGgmm,matrix,missing-method}
\alias{rUGgmm,missing,matrix-method}
\alias{rUGgmm,numeric,matrix-method}
\alias{rUGgmm,integer,matrix-method}
\alias{rUGgmm,graphBAM,missing-method}
\alias{rUGgmm,missing,graphBAM-method}
\alias{rUGgmm,numeric,graphBAM-method}
\alias{rUGgmm,integer,graphBAM-method}
\alias{rUGgmm}
\alias{show,UGgmm-method}
\alias{show,UGgmmSummary-method}
\alias{dim,UGgmm-method}
\alias{dimnames,UGgmm-method}
\alias{names,UGgmm-method}
\alias{$,UGgmm-method}
\alias{summary,UGgmm-method}
\alias{plot,UGgmm,ANY-method}
\alias{rmvnorm,integer,numeric-method}
\alias{rmvnorm,numeric,UGgmm-method}
\alias{rmvnorm,numeric,numeric-method}
\alias{rmvnorm}

\title{Undirected Gaussian graphical Markov model}
\description{
  The \code{"UGgmm"} class is the class of undirected Gaussian graphical
  Markov models defined within the \code{\link{qpgraph}} package to store
  simulate and manipulate this type of graphical Markov models (GMMs).

  An undirected Gaussian GMM is a family of multivariate normal distributions
  sharing a set of conditional independences encoded by means of an undirected
  graph. Further details can be found in the book of Lauritzen (1996).
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{UGgmm(g, ...)} corresponding
  to constructor methods or \code{rUGgmm(n, g, ...)} corresponding to random
  simulation methods.
}
\section{Slots}{
  \describe{
    \item{\code{p}:}{Object of class \code{"integer"} storing the dimension of the
      undirected Gaussian GMM corresponding to the number of random variables.}
    \item{\code{g}:}{Object of class \code{\link[graph]{graphBAM-class}} storing
      the associated undirected labeled graph.}
    \item{\code{mean}:}{Object of class \code{"numeric"} storing the mean vector.}
    \item{\code{sigma}:}{Object of class \code{\link[Matrix]{dspMatrix-class}}
      storing the covariance matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{UGgmm(g)}}{Constructor method where \code{g} can be either an
      adjacency matrix or a \code{\link[graph]{graphBAM-class}} object.}
    \item{\code{rUGgmm(n, g)}}{Constructor simulation method that allows one to
      simulate undirected Gaussian GMMs where \code{n} is the number of GMMs to
      simulate and \code{g} can be either a \linkS4class{graphParam} object,
      an adjacency matrix or a \code{\link[graph]{graphBAM-class}} object.}
    \item{\code{names(x)}}{Accessor method to obtain the names of the
      elements in the object \code{x} that can be retrieved with the \code{$}
      accessor operator.}
    \item{\code{$}}{Accessor operator to retrieve elements of the object
      in an analogous way to a \code{list}.}
    \item{\code{dim(x)}}{Dimension of the undirected Gaussian GMM corresponding
      to the total number of random variables.}
    \item{\code{dimnames(x)}}{Names of the random variables in the undirected
      Gaussian GMM.}
    \item{\code{show(object)}}{Method to display some bits of information about
      the input undirected Gaussian GMM specified in \code{object}.}
    \item{\code{summary(object)}}{Method to display a sumamry of the main features
      of the input undirected Gaussian GMM specified in \code{object}.}
    \item{\code{plot(x, ...)}}{Method to plot the undirected graph associated to the
      the input undirected Gaussian GMM specified in \code{x}. It uses the plotting
      capabilities from the \code{Rgraphviz} library to which further arguments
      specified in \code{...} are further passed.}
  }
}
\references{
  Lauritzen, S.L. \emph{Graphical models}. Oxford University Press, 1996.
}
\author{R. Castelo}
\seealso{
  \linkS4class{HMgmm}
}
\keyword{classes}
