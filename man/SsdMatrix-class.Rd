\name{SsdMatrix-class}
\docType{class}
\alias{SsdMatrix-class}
\alias{show,SsdMatrix-method}
\alias{determinant,SsdMatrix-method}
\alias{dim,SsdMatrix-method}
\alias{dimnames,SsdMatrix-method}

\title{Sum of squares and deviations Matrices}
\description{
  The \code{"SsdMatrix"} class is the class of symmetric, dense matrices
  in packed storage (just as a \code{\link[Matrix]{dspMatrix-class}}, i.e.,
  only the upper triangle is stored) defined within the \code{\link{qpgraph}}
  package to store corrected, or uncorrected, matrices of the sum of squares
  and deviations (SSD) of pairs of random variables. A corrected SSD matrix
  corresponds to a sample covariance matrix.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("SsdMatrix", ...)}
  or by using \code{\link{qpCov}()} which estimates a sample covariance
  matrix from data returning an object of this class.
}
\section{Slots}{
  \describe{
    \item{\code{ssd}:}{Object of class \code{\link[Matrix]{dspMatrix-class}}
      storing the SSD matrix.}
    \item{\code{n}:}{Object of class \code{"numeric"} storing the sample
      size employed to estimate the SSD matrix stored in the slot \code{ssd}.
      This is specially relevant when the SSD matrix was estimated from data
      with missing values by using complete observations only, which is the
      default mode of operation of \code{\link{qpCov}()}.}
  }
}
\section{Extends}{
  \code{"SsdMatrix"} extends class \code{"dspMatrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "SsdMatrix")}}
    \item{dimnames}{\code{signature(x = "SsdMatrix")}}
    \item{show}{\code{signature(object = "SsdMatrix")}}
    \item{determinant}{\code{signature(object = "SsdMatrix")}}
  }
}
\keyword{classes}
