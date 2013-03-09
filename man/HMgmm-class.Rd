\name{HMgmm-class}
\docType{class}
\alias{HMgmm-class}
\alias{HMgmm,missing-method}
\alias{HMgmm,matrix-method}
\alias{HMgmm,graphBAM-method}
\alias{HMgmm}
\alias{rHMgmm,markedGraphParam,missing-method}
\alias{rHMgmm,missing,markedGraphParam-method}
\alias{rHMgmm,numeric,markedGraphParam-method}
\alias{rHMgmm,integer,markedGraphParam-method}
\alias{rHMgmm,matrix,missing-method}
\alias{rHMgmm,missing,matrix-method}
\alias{rHMgmm,numeric,matrix-method}
\alias{rHMgmm,integer,matrix-method}
\alias{rHMgmm,graphBAM,missing-method}
\alias{rHMgmm,missing,graphBAM-method}
\alias{rHMgmm,numeric,graphBAM-method}
\alias{rHMgmm,integer,graphBAM-method}
\alias{rHMgmm}
\alias{show,HMgmm-method}
\alias{show,HMgmmSummary-method}
\alias{dim,HMgmm-method}
\alias{dimnames,HMgmm-method}
\alias{names,HMgmm-method}
\alias{$,HMgmm-method}
\alias{summary,HMgmm-method}
\alias{plot,HMgmm,ANY-method}
\alias{rcmvnorm,ANY,HMgmm-method}
\alias{rcmvnorm}


\title{Homogeneous mixed graphical Markov model}
\description{
  The \code{"HMgmm"} class is the class of homogeneous mixed graphical
  Markov models defined within the \code{\link{qpgraph}} package to store
  simulate and manipulate this type of graphical Markov models (GMMs).

  An homogeneous mixed GMM is a family of multivariate conditional Gaussian
  distributions on mixed discrete and continuous variables sharing a set of
  conditional independences encoded by means of a marked graph. Further
  details can be found in the book of Lauritzen (1996).
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{HMgmm(g, ...)} corresponding
  to constructor methods or \code{rHMgmm(n, g, ...)} corresponding to random
  simulation methods.
}
\section{Slots}{
  \describe{
    \item{\code{pI}:}{Object of class \code{"integer"} storing the number of
      discrete random variables.}
    \item{\code{pY}:}{Object of class \code{"integer"} storing the number of
      continuous random variables.}
    \item{\code{g}:}{Object of class \code{\link[graph]{graphBAM-class}} storing
      the associated marked graph.}
    \item{\code{vtype}:}{Object of class \code{"factor"} storing the type (discrete
      or continuous) of each random variable.}
    \item{\code{dLevels}:}{Object of class \code{"integer"} storing the number of
      levels of each discrete random variable.}
    \item{\code{a}:}{Object of class \code{"numeric"} storing the vector of additive
      linear effects on continuous variables connected to discrete ones.}
    \item{\code{rho}:}{Object of class \code{"numeric"} storing the value of the
      marginal correlation between two continuous random variables.}
    \item{\code{sigma}:}{Object of class \code{\link[Matrix]{dspMatrix-class}}
      storing the covariance matrix.}
    \item{\code{mean}:}{Object of class \code{"numeric"} storing the mean vector.}
    \item{\code{eta2}:}{Object of class \code{"numeric"} storing for each continuous
      variable connected to a discrete one, the fraction of variance of the
      continuous variable explained by the discrete one.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{HMgmm(g)}}{Constructor method where \code{g} can be either an
      adjacency matrix or a \code{\link[graph]{graphBAM-class}} object.}
    \item{\code{rHMgmm(n, g)}}{Constructor simulation method that allows one to
      simulate homogeneous mixed GMMs where \code{n} is the number of GMMs to
      simulate and \code{g} can be either a \linkS4class{markedGraphParam} object,
      an adjacency matrix or a \code{\link[graph]{graphBAM-class}} object.}
    \item{\code{names(x)}}{Accessor method to obtain the names of the
      elements in the object \code{x} that can be retrieved with the \code{$}
      accessor operator.}
    \item{\code{$}}{Accessor operator to retrieve elements of the object
      in an analogous way to a \code{list}.}
    \item{\code{dim(x)}}{Dimension of the homogeneous mixed GMM corresponding
      to the number of discrete and continuous random variables.}
    \item{\code{dimnames(x)}}{Names of the discrete and continuous random
      variables in the homogeneous mixed GMM.}
    \item{\code{show(object)}}{Method to display some bits of information about
      the input homogeneous mixed GMM specified in \code{object}.}
    \item{\code{summary(object)}}{Method to display a sumamry of the main features
      of the input homogeneous mixed GMM specified in \code{object}.}
    \item{\code{plot(x, ...)}}{Method to plot the undirected graph associated to the
      the input homogeneous mixed GMM specified in \code{x}. It uses the plotting
      capabilities from the \code{Rgraphviz} library to which further arguments
      specified in \code{...} are further passed.}
  }
}
\references{
  Lauritzen, S.L. \emph{Graphical models}. Oxford University Press, 1996.
}
\author{R. Castelo}
\seealso{
  \linkS4class{UGgmm}
}
\keyword{classes}
